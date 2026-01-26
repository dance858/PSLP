/*
 * Copyright 2025 Daniel Cederberg
 *
 * This file is part of the PSLP project (LP Presolver).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "SimpleReductions.h"
#include "Activity.h"
#include "Bounds.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "Locks.h"
#include "Matrix.h"
#include "Numerics.h"
#include "Postsolver.h"
#include "Problem.h"
#include "RowColViews.h"
#include "State.h"

void delete_fixed_cols_from_problem(Problem *prob)
{
    const Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT;
    const Bound *bounds = constraints->bounds;
    const RowTag *row_tags = constraints->row_tags;
    double *lhs = constraints->lhs;
    double *rhs = constraints->rhs;
    Activity *activities = constraints->state->activities;
    const iVec *fixed_cols_to_delete =
        prob->constraints->state->fixed_cols_to_delete;

    int i, j, col, len, row;
    int *rows;
    double *vals;

    for (i = 0; i < fixed_cols_to_delete->len; ++i)
    {
        col = fixed_cols_to_delete->data[i];
        assert(HAS_TAG(constraints->col_tags[col], C_TAG_FIXED));

        // if a variable has been fixed to INF it should not have been pushed to
        // fixed_cols_to_delete, so it shouldn't appear here
        assert(!IS_ABS_INF(bounds[col].ub));
        assert(bounds[col].lb == bounds[col].ub);

        // If the variable has been fixed to 0 or INF we must not do anything.
        // It is not necessary to update n_inf_min and n_inf_max (this has
        // already been done by boundChange).
        double ub = bounds[col].ub;
        if (ub == 0)
        {
            continue;
        }

        vals = AT->x + AT->p[col].start;
        rows = AT->i + AT->p[col].start;
        len = AT->p[col].end - AT->p[col].start;

        // ------------------------------------------------------------------------
        //           update left-and-right hand sides of constraints
        // ------------------------------------------------------------------------
        for (j = 0; j < len; ++j)
        {
            row = rows[j];

            // When a singleton row is used to fix a variable, the singleton row
            // will be inactive when trying to remove the contribution of the
            // fixed column. We want to skip it.
            if (HAS_TAG(row_tags[row], R_TAG_INACTIVE))
            {
                continue;
            }

            double contribution = vals[j] * ub;

            if (!HAS_TAG(row_tags[row], R_TAG_LHS_INF))
            {
                lhs[row] -= contribution;
            }

            if (!HAS_TAG(row_tags[row], R_TAG_RHS_INF))
            {
                rhs[row] -= contribution;
            }

            // not necessary to update n_inf_min and n_in_max since a fixed
            // variable has lb = ub = finite value so it doesn't contribute with
            // infinite bounds to the activity
            if (activities[row].n_inf_max == 0)
            {
                activities[row].max -= contribution;
            }

            if (activities[row].n_inf_min == 0)
            {
                activities[row].min -= contribution;
            }

            // Can a ranged row become an inequality due to numerics?
            assert(HAS_TAG(row_tags[row], R_TAG_LHS_INF) ||
                   HAS_TAG(row_tags[row], R_TAG_RHS_INF) ||
                   HAS_TAG(row_tags[row], R_TAG_EQ) ||
                   !IS_EQUAL_FEAS_TOL(lhs[row], rhs[row]));
        }

        // fix variable in the objective
        fix_var_in_obj(prob->obj, col, ub);
    }
}

static inline PresolveStatus remove_stonrow(Problem *prob, int row)
{
    // assume row 'row' is a singleton with variable k
    Constraints *constrs = prob->constraints;
    double lhs = constrs->lhs[row];
    double rhs = constrs->rhs[row];
    RowTag *row_tag = constrs->row_tags + row;
    RowRange *row_range = constrs->A->p + row;
    int k = constrs->A->i[row_range->start];
    double aik = constrs->A->x[row_range->start];
    ColTag col_tag = constrs->col_tags[k];
    int *row_size = constrs->state->row_sizes + row;
    PostsolveInfo *postsolve_info = constrs->state->postsolve_info;

    assert(!HAS_TAG(*row_tag, R_TAG_INACTIVE) &&
           !HAS_TAG(col_tag, C_TAG_SUBSTITUTED));

    // if the same variable appears in two singleton rows and the first one
    // that is processed is an equality row, then the variable has been fixed
    // and we should do nothing
    if (HAS_TAG(col_tag, C_TAG_FIXED))
    {
        return UNCHANGED;
    }

    assert(!HAS_TAG(col_tag, C_TAG_INACTIVE));

    // ----------------------------------------------------------------------
    //        if the row is an equality, we fix the variable
    // ----------------------------------------------------------------------
    if (HAS_TAG(*row_tag, R_TAG_EQ))
    {
        assert(lhs == rhs);

        // printf("fixing col %d to %f from row %d \n", k, rhs / aik, row);
        if (fix_col(constrs, k, rhs / aik, prob->obj->c[k]) == INFEASIBLE)
        {
            return INFEASIBLE;
        }

        // We manually mark the row inactive, since it is not necessary to
        // append it to the list of rows to delete. (It is not necessary
        // because the column has been appended to fixed_cols_to_delete.)
        *row_size = SIZE_INACTIVE_ROW;
        RESET_TAG(*row_tag, R_TAG_INACTIVE);
        row_range->end = row_range->start;
        constrs->A->nnz--;

        // dual postsolve:
        const Matrix *AT = constrs->AT;
        const int *rows = AT->i + AT->p[k].start;
        const double *vals = AT->x + AT->p[k].start;
        int len = AT->p[k].end - AT->p[k].start;
        save_retrieval_added_rows(postsolve_info, row, rows, vals, len, aik);
        save_retrieval_deleted_row(postsolve_info, row, prob->obj->c[k] / aik);
    }
    // ----------------------------------------------------------------------
    //  if the row is an inequality, we update the bounds (if they are tighter
    //  than the current bounds)
    // ----------------------------------------------------------------------
    else
    {
        bool rhs_inf = HAS_TAG(*row_tag, R_TAG_RHS_INF);
        bool lhs_inf = HAS_TAG(*row_tag, R_TAG_LHS_INF);

        if (aik > 0)
        {
            if ((!rhs_inf && (update_ub(constrs, k, rhs / aik,
                                        row HUGE_BOUND_IS_OK) == INFEASIBLE)) ||
                (!lhs_inf && (update_lb(constrs, k, lhs / aik,
                                        row HUGE_BOUND_IS_OK) == INFEASIBLE)))
            {
                return INFEASIBLE;
            }
        }
        else
        {
            if ((!rhs_inf && update_lb(constrs, k, rhs / aik,
                                       row HUGE_BOUND_IS_OK) == INFEASIBLE) ||
                (!lhs_inf && update_ub(constrs, k, lhs / aik,
                                       row HUGE_BOUND_IS_OK) == INFEASIBLE))
            {
                return INFEASIBLE;
            }
        }

        // must call set_row_to_inactive since the column is not appended to
        // fixed_cols_to_delete.
        set_row_to_inactive(row, row_tag, constrs->state->rows_to_delete,
                            postsolve_info, 0.0);
    }

    return REDUCED;
}

PresolveStatus remove_ston_rows(Problem *prob)
{
    iVec *ston_rows = prob->constraints->state->ston_rows;
    const int *ston_rows_data = ston_rows->data;
    int len = (int) ston_rows->len;
    RowTag *row_tags = prob->constraints->row_tags;
    DEBUG(verify_no_duplicates_sort(ston_rows));

    if (len == 0)
    {
        return UNCHANGED;
    }

    // first we process singleton equality rows, then singleton inequalities
    // (this ordering simplifies the dual postsolve when the same variable
    // appears in an equality and an inequality row)
    for (int i = 0; i < len; ++i)
    {
        int row = ston_rows_data[i];

        if (!HAS_TAG(row_tags[row], R_TAG_EQ) ||
            HAS_TAG(row_tags[row], R_TAG_INACTIVE))
        {
            continue;
        }

        assert(!HAS_TAG(row_tags[row], R_TAG_INACTIVE));

        if (INFEASIBLE == remove_stonrow(prob, row))
        {
            return INFEASIBLE;
        }
    }

    for (int i = 0; i < len; ++i)
    {
        int row = ston_rows_data[i];

        // eliminated equality rows have been marked as inactive
        // (or one might remain if the variable in it has been fixed
        // by another equality row)
        if (HAS_TAG(row_tags[row], (R_TAG_INACTIVE | R_TAG_EQ)))
        {
            continue;
        }

        assert(!HAS_TAG(row_tags[row], R_TAG_EQ));

        if (INFEASIBLE == remove_stonrow(prob, row))
        {
            return INFEASIBLE;
        }
    }

    // we only updated A->nnz inside remove_stonrow so we need to
    // update AT->nnz for consistency
    prob->constraints->AT->nnz = prob->constraints->A->nnz;

    iVec_clear_no_resize(ston_rows);

    // singleton inequality rows have been marked as inactive so we
    // must remove them
    delete_inactive_rows(prob->constraints);

    delete_fixed_cols_from_problem(prob);

    // variables might have been fixed (by ston equality rows) so we
    // must remove them
    delete_inactive_cols_from_A_and_AT(prob->constraints);

    DEBUG(verify_problem_up_to_date(prob->constraints));
    return REDUCED;
}

// If the returned row tag is R_TAG_INFEAS, then the constraint is feasible.
// If HAS_TAG(return_tag, R_TAG_RHS_INF), then the rhs is finite and redundant.
// If HAS_TAG(return_tag, R_TAG_LHS_INF), then the lhs is finite and redundant.
static inline RowTag check_activity(const Activity *act, double lhs, double rhs,
                                    RowTag row_tag)
{
    assert(!HAS_TAG(row_tag, R_TAG_EQ));
    RowTag return_tag = R_TAG_NONE;

    // -------------------------------------------------------------------------
    // Compare the activity with the rhs of the constraint. We check for both
    // infeasibility and redundancy.
    // -------------------------------------------------------------------------
    if (!HAS_TAG(row_tag, R_TAG_RHS_INF))
    {
        if (act->n_inf_min == 0 && act->min >= rhs + FEAS_TOL)
        {
            return R_TAG_INFEAS;
        }
        else if (act->n_inf_max == 0 && act->max <= rhs + FEAS_TOL)
        {
            UPDATE_TAG(return_tag, R_TAG_RHS_INF);
        }
    }

    // -------------------------------------------------------------------------
    // Compare the activity with the lhs of the constraint. We check for both
    // infeasibility and redundancy.
    // -------------------------------------------------------------------------
    if (!HAS_TAG(row_tag, R_TAG_LHS_INF))
    {
        if (act->n_inf_max == 0 && act->max <= lhs - FEAS_TOL)
        {
            return R_TAG_INFEAS;
        }
        else if (act->n_inf_min == 0 && act->min >= lhs - FEAS_TOL)
        {
            UPDATE_TAG(return_tag, R_TAG_LHS_INF);
        }
    }

#ifndef NDEBUG
    // can it happen that an equality constraint becomes an inequality or
    // redundant?
    if (HAS_TAG(row_tag, R_TAG_EQ))
    {
        assert(!HAS_TAG(return_tag, R_TAG_LHS_INF) &&
               !HAS_TAG(return_tag, R_TAG_RHS_INF));
    }
#endif

    return return_tag;
}

PresolveStatus check_activities(Problem *prob)
{
    const iVec *acts_to_check = prob->constraints->state->updated_activities;

    // verify that the same row hasn't been appended twice and that the
    // activities are correct
    DEBUG(verify_no_duplicates_sort(acts_to_check));
    DEBUG(verify_activities(prob->constraints));

    const Activity *activities = prob->constraints->state->activities;
    const int *row_sizes = prob->constraints->state->row_sizes;
    Lock *locks = prob->constraints->state->col_locks;
    iVec *rows_to_delete = prob->constraints->state->rows_to_delete;
    RowTag *row_tags = prob->constraints->row_tags;
    const double *lhs = prob->constraints->lhs;
    const double *rhs = prob->constraints->rhs;
    const Matrix *A = prob->constraints->A;

    int i, j, row;

    for (i = 0; i < acts_to_check->len; ++i)
    {
        row = acts_to_check->data[i];

        // Skip inactive rows, empty rows, and singleton rows. We must include
        // the check for an inactive row, since the row size of an inactive
        // row may not yet have been updated when this function is called.
        // IDEA: it is a design choice to have skip equality rows here, but
        // we should try to remove this restriction and see if it matters.
        if (HAS_TAG(row_tags[row], (R_TAG_INACTIVE | R_TAG_EQ)) ||
            row_sizes[row] <= 1)
        {
            continue;
        }

        assert(activities[row].n_inf_min == 0 || activities[row].n_inf_max == 0);

        RowTag status_tag =
            check_activity(activities + row, lhs[row], rhs[row], row_tags[row]);

        // ----------------------------------------------------------------
        //              if the constraint is redundant
        // ----------------------------------------------------------------
        if (HAS_TAG((status_tag | row_tags[row]), R_TAG_LHS_INF) &&
            HAS_TAG((status_tag | row_tags[row]), R_TAG_RHS_INF))
        {
            set_row_to_inactive(row, row_tags + row, rows_to_delete,
                                prob->constraints->state->postsolve_info, 0.0);
        }
        // ----------------------------------------------------------------
        //              if lhs is finite but redundant
        // (it can't happen that both HAS_TAG(status_tag, R_TAG_LHS_INF)
        //  and HAS_TAG(row_tags[row], R_TAG_LHS_INF) are true because of
        //  the design of check_activity)
        // ----------------------------------------------------------------
        else if (HAS_TAG(status_tag, R_TAG_LHS_INF))
        {
            assert(!HAS_TAG(row_tags[row], R_TAG_LHS_INF));
            RESET_TAG(row_tags[row], R_TAG_LHS_INF);

            // update locks
            for (j = A->p[row].start; j < A->p[row].end; ++j)
            {
                assert(A->x[j] != 0);
                if (A->x[j] > 0)
                {
                    locks[A->i[j]].down--;
                }
                else
                {
                    locks[A->i[j]].up--;
                }
            }

            prob->constraints->lhs[row] = -INF;
        }
        // ----------------------------------------------------------------
        //             if rhs is finite but redundant
        // (it can't happen that both HAS_TAG(status_tag, R_TAG_RHS_INF)
        //  and HAS_TAG(row_tags[row], R_TAG_RHS_INF) are true because of
        //  the design of check_activity)
        // ----------------------------------------------------------------
        else if (HAS_TAG(status_tag, R_TAG_RHS_INF))
        {
            assert(!HAS_TAG(row_tags[row], R_TAG_RHS_INF));
            RESET_TAG(row_tags[row], R_TAG_RHS_INF);

            // update locks
            for (j = A->p[row].start; j < A->p[row].end; ++j)
            {
                assert(A->x[j] != 0);
                if (A->x[j] > 0)
                {
                    locks[A->i[j]].up--;
                }
                else
                {
                    locks[A->i[j]].down--;
                }
            }

            prob->constraints->rhs[row] = INF;
        }
        // ----------------------------------------------------------------
        //             if the constraint is infeasible
        // ----------------------------------------------------------------
        else if (HAS_TAG(status_tag, R_TAG_INFEAS))
        {
            return INFEASIBLE;
        }
    }

    delete_inactive_rows(prob->constraints);

    // we should NOT clear updated activities here since we need them in
    // primal propagation
    DEBUG(verify_problem_up_to_date(prob->constraints));
    return UNCHANGED;
}

PresolveStatus remove_empty_rows(const Constraints *constraints)
{
    iVec *empty_rows = constraints->state->empty_rows;
    const int *empty_rows_data = empty_rows->data;
    int len = (int) empty_rows->len;

    if (len == 0)
    {
        return UNCHANGED;
    }

    int *row_sizes = constraints->state->row_sizes;
    RowTag *rtags = constraints->row_tags;
    const double *lhs = constraints->lhs;
    const double *rhs = constraints->rhs;
    PostsolveInfo *postsolve_info = constraints->state->postsolve_info;

    for (int i = 0; i < len; ++i)
    {
        int row = empty_rows_data[i];
        assert(row_sizes[row] == 0);

        bool infeasible =
            (!HAS_TAG(rtags[row], R_TAG_LHS_INF) && lhs[row] >= FEAS_TOL) ||
            (!HAS_TAG(rtags[row], R_TAG_RHS_INF) && rhs[row] <= -FEAS_TOL);

        if (infeasible)
        {
            return INFEASIBLE;
        }

        UPDATE_TAG(rtags[row], R_TAG_INACTIVE);
        row_sizes[row] = SIZE_INACTIVE_ROW;
        save_retrieval_deleted_row(postsolve_info, row, 0.0);
    }

    iVec_clear_no_resize(empty_rows);
    return UNCHANGED;
}

PresolveStatus remove_empty_cols(Problem *prob)
{
    iVec *empty_cols = prob->constraints->state->empty_cols;
    const int *empty_cols_data = empty_cols->data;
    int len = empty_cols->len;

    double val;
    int i, k;

    if (len == 0)
    {
        return UNCHANGED;
    }

    PostsolveInfo *postsolve_info = prob->constraints->state->postsolve_info;
    int *col_sizes = prob->constraints->state->col_sizes;
    ColTag *col_tags = prob->constraints->col_tags;
    Bound *bounds = prob->constraints->bounds;
    double *c = prob->obj->c;
    double *offset = &(prob->obj->offset);

    for (i = 0; i < len; ++i)
    {
        k = empty_cols_data[i];

        assert(col_sizes[k] == 0 && !HAS_TAG(col_tags[k], C_TAG_INACTIVE));

        // --------------------------------------------------------------------
        // if the objective coefficient is 0 we set the variable to 0 if
        // feasible, otherwise we set it equal to one of the bounds
        // --------------------------------------------------------------------
        if (c[k] == 0)
        {
            if (!HAS_TAG(col_tags[k], C_TAG_LB_INF) && bounds[k].lb > 0)
            {
                val = bounds[k].lb;
            }
            else if (!HAS_TAG(col_tags[k], C_TAG_UB_INF) && bounds[k].ub < 0)
            {
                val = bounds[k].ub;
            }
            else
            {
                val = 0;
            }
        }
        // ------------------------------------------------------------------
        //                  fix variable to upper bound
        // ------------------------------------------------------------------
        else if (c[k] < 0)
        {
            if (HAS_TAG(col_tags[k], C_TAG_UB_INF))
            {
                return UNBNDORINFEAS;
            }

            val = bounds[k].ub;
            bounds[k].lb = bounds[k].ub;
        }
        // ------------------------------------------------------------------
        //                   fix variable to lower bound
        // ------------------------------------------------------------------
        else
        {
            assert(c[k] > 0);
            if (HAS_TAG(col_tags[k], C_TAG_LB_INF))
            {
                return UNBNDORINFEAS;
            }

            val = bounds[k].lb;
            bounds[k].ub = bounds[k].lb;
        }

        // add offset to objective
        *offset += c[k] * val;

        // mark column as fixed
        UPDATE_TAG(col_tags[k], C_TAG_FIXED);
        col_sizes[k] = SIZE_INACTIVE_COL;

        // store postsolve information (for empty cols we want zk = ck)
        // passing NULL ptrs are safe since the length is 0
        save_retrieval_fixed_col(postsolve_info, k, val, prob->obj->c[k], NULL, NULL,
                                 0);
    }

    iVec_clear_no_resize(empty_cols);
    return UNCHANGED;
}

void clean_small_coeff_A(Matrix *A, const Bound *bounds, const RowTag *row_tags,
                         const ColTag *col_tags, double *rhs, double *lhs)
{
    int ii, row, col, len, n_rows, *cols;
    double Aik_abs, cum_changes, diff, *vals;
    bool has_finite_bounds, set_coeff_to_zero;
    n_rows = A->m;

    for (row = 0; row < n_rows; ++row)
    {
        cols = A->i + A->p[row].start;
        vals = A->x + A->p[row].start;
        len = A->p[row].end - A->p[row].start;
        cum_changes = 0.0;

        for (ii = 0; ii < len; ++ii)
        {
            set_coeff_to_zero = false;
            col = cols[ii];
            Aik_abs = ABS(vals[ii]);

            has_finite_bounds =
                !HAS_TAG(col_tags[col], (C_TAG_LB_INF | C_TAG_UB_INF));

            if (has_finite_bounds)
            {
                assert(!IS_ABS_INF(bounds[col].lb) && !IS_ABS_INF(bounds[col].ub));
                diff = bounds[col].ub - bounds[col].lb;
                assert(diff >= 0);

                // we take care of fixed variables later
                if (diff == 0)
                {
                    continue;
                }

                if (Aik_abs < CLEAN1 && Aik_abs * diff * len < 1e-2 * FEAS_TOL)
                {
                    set_coeff_to_zero = true;
                }
                else if (cum_changes + Aik_abs * diff < 1e-1 * FEAS_TOL)
                {
                    cum_changes += Aik_abs * diff;
                    set_coeff_to_zero = true;
                }
            }

            if (set_coeff_to_zero)
            {
                if (!HAS_TAG(row_tags[row], R_TAG_LHS_INF))
                {
                    assert(bounds[col].lb != -INF);
                    lhs[row] -= vals[ii] * bounds[col].lb;
                }

                if (!HAS_TAG(row_tags[row], R_TAG_RHS_INF))
                {
                    assert(bounds[col].lb != -INF);
                    rhs[row] -= vals[ii] * bounds[col].lb;
                }
            }

            if (Aik_abs < CLEAN2 || set_coeff_to_zero)
            {
                RowView row_view =
                    new_rowview(vals, cols, &len, A->p + row, NULL, NULL, NULL, -1);
                remove_coeff(&row_view, col);
                ii--;
                A->nnz--;
            }
        }
    }
}

PresolveStatus remove_variables_with_close_bounds(Problem *prob)
{
    Constraints *constraints = prob->constraints;
    const double *c = prob->obj->c;
    const Bound *bounds = constraints->bounds;
    const ColTag *col_tags = constraints->col_tags;
    int n_cols = constraints->n;
    const int *col_sizes = constraints->state->col_sizes;

    for (int ii = 0; ii < n_cols; ++ii)
    {
        if (col_sizes[ii] < 0 ||
            HAS_TAG(col_tags[ii], (C_TAG_LB_INF | C_TAG_UB_INF)))
        {
            continue;
        }

        assert(!HAS_TAG(col_tags[ii], C_TAG_INACTIVE));

        if (IS_EQUAL_FEAS_TOL(bounds[ii].lb, bounds[ii].ub))
        {
            // no need to check return value since bounds are equal
            fix_col(constraints, ii, bounds[ii].lb, c[ii]);
        }
        else if (bounds[ii].lb > bounds[ii].ub + FEAS_TOL)
        {
            return INFEASIBLE;
        }
    }

    // remove fixed columns
    delete_fixed_cols_from_problem(prob);
    delete_inactive_cols_from_A_and_AT(constraints);
    return UNCHANGED;
}
