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

#include "StonCols.h"
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

/* Helper for handling the case where a column-singleton variable in an
   equality row is only implied free from above. */
static void handle_impl_free_from_above_eq(RowView *row, int k, double Aik,
                                           double lb, double ub, ColTag *col_tag,
                                           Activity *act, Lock *locks, int *A_nnz,
                                           const Bound *bounds,
                                           PostsolveInfo *postsolve_info, double ck)
{
    /* dual postsolve */
    save_retrieval_eq_to_ineq(postsolve_info, row->i, ck / Aik);
    save_retrieval_sub_col(postsolve_info, k, row->cols, row->vals, *row->len,
                           *row->rhs, row->i, 0.0);

    assert(!IS_ABS_INF(lb) && !HAS_TAG(*col_tag, C_TAG_LB_INF));

    /* remove coefficient from A */
    remove_coeff(row, k);
    (*A_nnz)--;

    if (Aik > 0)
    {
        /* modify constraint */
        *row->rhs -= lb * Aik;
        *row->lhs = -INF;
        RESET_TAG(*row->tag, R_TAG_LHS_INF);

        /* update locks */
        update_locks_eq_to_ineq(locks, row->vals, row->cols, *row->len, true);

        /* update activity */
        if (act->n_inf_min == 0)
        {
            act->min -= lb * Aik;
        }

        if (HAS_TAG(*col_tag, C_TAG_UB_INF))
        {
            act->n_inf_max--;

            if (act->n_inf_max == 0)
            {
                act->max =
                    compute_max_act_no_tags(row->vals, row->cols, *row->len, bounds);
                assert(!IS_ABS_INF(act->max));
            }
        }
        else
        {
            assert(!IS_ABS_INF(ub));
            if (act->n_inf_max == 0)
            {
                act->max -= ub * Aik;
            }
        }
    }
    else
    {
        /* modify constraint and update locks */
        *row->lhs -= lb * Aik;
        *row->rhs = INF;
        RESET_TAG(*row->tag, R_TAG_RHS_INF);
        update_locks_eq_to_ineq(locks, row->vals, row->cols, *row->len, false);

        /* update activity */
        if (act->n_inf_max == 0)
        {
            act->max -= lb * Aik;
        }

        if (HAS_TAG(*col_tag, C_TAG_UB_INF))
        {
            act->n_inf_min--;

            if (act->n_inf_min == 0)
            {
                act->min =
                    compute_min_act_no_tags(row->vals, row->cols, *row->len, bounds);
                assert(!IS_ABS_INF(act->min));
            }
        }
        else
        {
            assert(!IS_ABS_INF(ub));
            if (act->n_inf_min == 0)
            {
                act->min -= ub * Aik;
            }
        }
    }
}

/* Helper for handling the case where a column-singleton variable in an
   equality row is only implied free from below. */
static void handle_impl_free_from_below_eq(RowView *row, int k, double Aik,
                                           double lb, double ub, ColTag *col_tag,
                                           Activity *act, Lock *locks, int *A_nnz,
                                           const Bound *bounds,
                                           PostsolveInfo *postsolve_info, double ck)
{
    save_retrieval_eq_to_ineq(postsolve_info, row->i, ck / Aik);
    save_retrieval_sub_col(postsolve_info, k, row->cols, row->vals, *row->len,
                           *row->rhs, row->i, 0.0);
    assert(!IS_ABS_INF(ub) && !HAS_TAG(*col_tag, C_TAG_UB_INF));

    /* remove coefficient from A */
    remove_coeff(row, k);
    (*A_nnz)--;

    if (Aik > 0)
    {
        /* modify constraint */
        *row->lhs -= ub * Aik;
        *row->rhs = INF;
        RESET_TAG(*row->tag, R_TAG_RHS_INF);

        /* update locks */
        update_locks_eq_to_ineq(locks, row->vals, row->cols, *row->len, false);

        /* update activity */
        if (act->n_inf_max == 0)
        {
            act->max -= ub * Aik;
        }

        if (HAS_TAG(*col_tag, C_TAG_LB_INF))
        {
            act->n_inf_min--;

            if (act->n_inf_min == 0)
            {
                act->min =
                    compute_min_act_no_tags(row->vals, row->cols, *row->len, bounds);
                assert(!IS_ABS_INF(act->min));
            }
        }
        else
        {
            assert(!IS_ABS_INF(lb));
            if (act->n_inf_min == 0)
            {
                act->min -= lb * Aik;
            }
        }
    }
    else
    {
        /* modify constraint and update locks*/
        *row->rhs -= ub * Aik;
        *row->lhs = -INF;
        RESET_TAG(*row->tag, R_TAG_LHS_INF);
        update_locks_eq_to_ineq(locks, row->vals, row->cols, *row->len, true);

        /* update activity */
        if (act->n_inf_min == 0)
        {
            act->min -= ub * Aik;
        }
        if (HAS_TAG(*col_tag, C_TAG_LB_INF))
        {
            act->n_inf_max--;

            if (act->n_inf_max == 0)
            {
                act->max =
                    compute_max_act_no_tags(row->vals, row->cols, *row->len, bounds);
                assert(!IS_ABS_INF(act->max));
            }
        }
        else
        {
            assert(!IS_ABS_INF(lb));
            if (act->n_inf_max == 0)
            {
                act->max -= lb * Aik;
            }
        }
    }
}

static PresolveStatus process_colston_eq(RowView *row, ColView *col, Objective *obj,
                                         double Aik, bool impl_free_from_above,
                                         bool impl_free_from_below, Activity *act,
                                         Lock *locks, iVec *rows_to_delete,
                                         int *A_nnz, const Bound *bounds,
                                         iVec *ston_rows, iVec *dton_rows,
                                         PostsolveInfo *postsolve_info)
{

    assert(*row->lhs == *row->rhs);
    assert(!HAS_TAG(*row->tag, (R_TAG_INACTIVE | R_TAG_LHS_INF | R_TAG_RHS_INF)));

    // No reduction is possible for a column singleton that is neither
    // implied free from below nor above, so we skip these.
    if (!impl_free_from_above && !impl_free_from_below)
    {
        return UNCHANGED;
    }

    // must store for dual postsolve
    double ck = obj->c[col->k];

    // substitute variable in objective and mark it as substituted
    // (note that we don't push it to list of substituted cols)
    sub_var_in_obj(obj, row->vals, row->cols, *row->len, col->k, Aik, *row->rhs);
    UPDATE_TAG(*col->tag, C_TAG_SUBSTITUTED);
    *col->len = SIZE_INACTIVE_COL;
    col->range->end = col->range->start;

    // -------------------------------------------------------------------------
    // if var is implied free, the constraint is redundant and we remove it.
    // -------------------------------------------------------------------------
    if (impl_free_from_above && impl_free_from_below)
    {
        save_retrieval_sub_col(postsolve_info, col->k, row->cols, row->vals,
                               *row->len, *row->rhs, row->i, ck);
        set_row_to_inactive(row->i, row->tag, rows_to_delete, postsolve_info,
                            ck / Aik);
        return REDUCED;
    }
    // -------------------------------------------------------------------------
    // If the variable is only implied free from above we must modify the
    // constraint when removing it.
    // -------------------------------------------------------------------------
    else if (impl_free_from_above)
    {
        handle_impl_free_from_above_eq(row, col->k, Aik, *col->lb, *col->ub,
                                       col->tag, act, locks, A_nnz, bounds,
                                       postsolve_info, ck);
    }
    // --------------------------------------------------------------------------
    // If the variable is only implied free from below we must modify the
    // constraint when removing it.
    // --------------------------------------------------------------------------
    else
    {
        assert(impl_free_from_below);
        handle_impl_free_from_below_eq(row, col->k, Aik, *col->lb, *col->ub,
                                       col->tag, act, locks, A_nnz, bounds,
                                       postsolve_info, ck);
    }

    if (*row->len == 1)
    {
        assert(!iVec_contains(ston_rows, row->i));
        iVec_append(ston_rows, row->i);
    }
    else if (*row->len == 2 && HAS_TAG(*row->tag, R_TAG_EQ))
    {
        assert(!iVec_contains(dton_rows, row->i));
        iVec_append(dton_rows, row->i);
    }

    assert(*row->len != 0);

    return REDUCED;
}

/* Fixes a variable that is a column singleton.
    Parameters:
    * val: value that variable should be fixed to
    * Aik: coefficient in the matrix A of the fixed variable
    * (lhs, rhs) belong to the row that the fixed variable appears in
    * (vals, cols, len, row_range) correspond to the row of the col_ston
    * (col_range, col_tag, col_size) are attributes of the variable that
                                      is fixed
    * obj: objective function
    * k: index of the fixed variable
    * A_nnz: number of nonzeros in the matrix A (OBS: this is modified)
    * row_tag: tag of the row that the fixed variable appears in
    * act: activity of the row that the fixed variable appears in
    * bounds: variable bounds

    @Note: we don't have to check for infeasibility here, because
            we only call this function with val equal to one of the
            bounds.
 */
static void fix_col_ston(double val, double Aik, RowView *row, ColView *col,
                         Activity *act, Objective *obj, int *A_nnz,
                         const Bound *bounds, PostsolveInfo *postsolve_info)
{
    // remove coefficient from A (which removes the variable from A)
    remove_coeff(row, col->k);
    (*A_nnz)--;

    // set row k from AT to inactive
    col->range->end = col->range->start;
    *col->len = SIZE_INACTIVE_COL;
    UPDATE_TAG(*col->tag, C_TAG_FIXED);

    // update lhs
    if (!HAS_TAG(*row->tag, R_TAG_LHS_INF))
    {
        *row->lhs -= Aik * val;
    }

    // update rhs
    if (!HAS_TAG(*row->tag, R_TAG_RHS_INF))
    {
        *row->rhs -= Aik * val;
    }

    // update min activity
    if (act->n_inf_min == 0)
    {
        act->min -= Aik * val;
    }
    else if ((Aik > 0 && HAS_TAG(*col->tag, C_TAG_LB_INF)) ||
             (Aik < 0 && HAS_TAG(*col->tag, C_TAG_UB_INF)))
    {
        act->n_inf_min--;

        if (act->n_inf_min == 0)
        {
            act->min =
                compute_min_act_no_tags(row->vals, row->cols, *row->len, bounds);
            assert(!IS_ABS_INF(act->min));
        }
    }

    // update max activity
    if (act->n_inf_max == 0)
    {
        act->max -= Aik * val;
    }
    else if ((Aik > 0 && HAS_TAG(*col->tag, C_TAG_UB_INF)) ||
             (Aik < 0 && HAS_TAG(*col->tag, C_TAG_LB_INF)))
    {
        act->n_inf_max--;

        if (act->n_inf_max == 0)
        {
            act->max =
                compute_max_act_no_tags(row->vals, row->cols, *row->len, bounds);
            assert(!IS_ABS_INF(act->max));
        }
    }

    // fix variable in objective
    fix_var_in_obj(obj, col->k, val);
    save_retrieval_fixed_col(postsolve_info, col->k, val, obj->c[col->k], &Aik,
                             &row->i, 1);
}

static inline PresolveStatus
process_colston_ineq(RowView *row, ColView *col, Objective *obj, double Aik,
                     bool impl_free_from_above, bool impl_free_from_below,
                     Activity *act, Lock *locks, iVec *rows_to_delete, int *A_nnz,
                     const Bound *bounds, iVec *ston_rows, iVec *dton_rows,
                     PostsolveInfo *postsolve_info)
{
    bool is_lhs_inf = HAS_TAG(*row->tag, R_TAG_LHS_INF);
    bool is_rhs_inf = HAS_TAG(*row->tag, R_TAG_RHS_INF);
    bool is_ub_inf = HAS_TAG(*col->tag, C_TAG_UB_INF);
    bool is_lb_inf = HAS_TAG(*col->tag, C_TAG_LB_INF);
    double ck = obj->c[col->k];

    // ------------------------------------------------------------------------
    //                    dual fix to lower bound
    // ------------------------------------------------------------------------
    if ((ck > 0 && Aik > 0 && is_lhs_inf) || (ck > 0 && Aik < 0 && is_rhs_inf))
    {
        if (is_lb_inf)
        {
            return UNBNDORINFEAS;
        }
        assert(!IS_ABS_INF(*col->lb));

        fix_col_ston(*col->lb, Aik, row, col, act, obj, A_nnz, bounds,
                     postsolve_info);

        // should not check for dton row here since the constraint is an ineq
        if (*row->len == 1)
        {
            assert(!iVec_contains(ston_rows, row->i));
            iVec_append(ston_rows, row->i);
        }

        assert(*row->len != 0);
        return REDUCED;
    }

    // ------------------------------------------------------------------------
    //                    dual fix to upper bound
    // ------------------------------------------------------------------------
    if ((ck < 0 && Aik < 0 && is_lhs_inf) || (ck < 0 && Aik > 0 && is_rhs_inf))
    {
        if (is_ub_inf)
        {
            return UNBNDORINFEAS;
        }
        assert(!IS_ABS_INF(*col->ub));

        fix_col_ston(*col->ub, Aik, row, col, act, obj, A_nnz, bounds,
                     postsolve_info);

        // should not check for dton row here since the constraint is an ineq
        if (*row->len == 1)
        {
            assert(!iVec_contains(ston_rows, row->i));
            iVec_append(ston_rows, row->i);
        }

        assert(*row->len != 0);
        return REDUCED;
    }

    // ------------------------------------------------------------------------
    //        Reduce column singletons that are (implied) free
    // ------------------------------------------------------------------------
    if (impl_free_from_above && impl_free_from_below)
    {
        double new_side;

        // lhs active at an optimal solution
        if ((ck > 0 && Aik > 0) || (ck < 0 && Aik < 0))
        {
            if (is_lhs_inf)
            {
                return UNBNDORINFEAS;
            }

            new_side = *row->lhs;
        }
        // rhs active at an optimal solution
        else if ((ck > 0 && Aik < 0) || (ck < 0 && Aik > 0))
        {
            if (is_rhs_inf)
            {
                return UNBNDORINFEAS;
            }

            new_side = *row->rhs;
        }
        else
        {
            assert(ck == 0);
            new_side = (is_lhs_inf) ? *row->rhs : *row->lhs;
        }

        // substitute variable in objective and mark it as substituted
        // (note that we don't push it to list of substituted cols)
        sub_var_in_obj(obj, row->vals, row->cols, *row->len, col->k, Aik, new_side);
        col->range->end = col->range->start;
        UPDATE_TAG(*col->tag, C_TAG_SUBSTITUTED);
        *col->len = SIZE_INACTIVE_COL;

        save_retrieval_sub_col(postsolve_info, col->k, row->cols, row->vals,
                               *row->len, new_side, row->i, ck);
        set_row_to_inactive(row->i, row->tag, rows_to_delete, postsolve_info,
                            ck / Aik);

        return REDUCED;
    }

    // ----------------------------------------------------------------
    // try to conclude that the upper part of the constraint is active
    // at an optimal solution
    // ----------------------------------------------------------------
    if ((ck < 0 && Aik > 0 && impl_free_from_above) ||
        (ck > 0 && Aik < 0 && impl_free_from_below))
    {
        // update locks in case the old lhs was infinite
        if (is_lhs_inf)
        {
            update_locks_ineq_to_eq(locks, row->vals, row->cols, *row->len, true);
        }

        assert(!is_rhs_inf);
        // update lhs = rhs
        *row->lhs = *row->rhs;
        RESET_TAG(*row->tag, R_TAG_EQ);

        if (*row->len == 2)
        {
            assert(!iVec_contains(dton_rows, row->i));
            iVec_append(dton_rows, row->i);
        }
        return REDUCED;
    }

    // ----------------------------------------------------------------
    // try to conclude that the lower part of the constraint is active
    // at an optimal solution
    // ----------------------------------------------------------------
    if ((ck > 0 && Aik > 0 && impl_free_from_below) ||
        (ck < 0 && Aik < 0 && impl_free_from_above))
    {
        // update locks in case the old rhs was infinite
        if (is_rhs_inf)
        {
            update_locks_ineq_to_eq(locks, row->vals, row->cols, *row->len, false);
        }

        assert(!is_lhs_inf);
        // update rhs = lhs
        *row->rhs = *row->lhs;
        RESET_TAG(*row->tag, R_TAG_EQ);

        if (*row->len == 2)
        {
            assert(!iVec_contains(dton_rows, row->i));
            iVec_append(dton_rows, row->i);
        }

        return REDUCED;
    }

    return UNCHANGED;
}

/* Refresh the list of column singletons. A variable that used to be a
   column singleton may no longer be a column singleton for one of the
   following reasons:
  1. It has been processed as a column singleton and is therefore inactive.
  2. It has been fixed by some other reduction.
*/
static inline void refresh_ston_cols(const int *col_sizes, iVec *ston_cols)
{
    int n_ston_cols = 0;
    int shift = 0;
    int len = (int) ston_cols->len;

    for (int i = 0; i < len; i++)
    {
        if (col_sizes[ston_cols->data[i]] != 1)
        {
            shift++;
        }
        else
        {
            ston_cols->data[i - shift] = ston_cols->data[i];
            n_ston_cols++;
        }
    }

    ston_cols->len = (size_t) n_ston_cols;
}

PresolveStatus remove_ston_cols__(Problem *prob)
{
    Constraints *constraints = prob->constraints;
    PostsolveInfo *postsolve_info = constraints->state->postsolve_info;
    const Matrix *A = constraints->A;
    int *A_nnz = &(constraints->A->nnz);
    const Matrix *AT = constraints->AT;
    RowTag *row_tags = constraints->row_tags;
    ColTag *col_tags = constraints->col_tags;
    double *lhs = constraints->lhs;
    double *rhs = constraints->rhs;
    Bound *bounds = constraints->bounds;
    int *row_sizes = constraints->state->row_sizes;
    int *col_sizes = constraints->state->col_sizes;
    iVec *ston_cols = constraints->state->ston_cols;
    iVec *ston_rows = constraints->state->ston_rows;
    iVec *dton_rows = constraints->state->dton_rows;
    Activity *acts = constraints->state->activities;
    Lock *locks = constraints->state->col_locks;
    iVec *rows_to_delete = constraints->state->rows_to_delete;
    double Aik;
    int i, k, kk;
    int *row_size;
    int *cols;
    double *vals;
    bool impl_free_from_above, impl_free_from_below;

    PresolveStatus status = UNCHANGED;

    // get up to date information about which columns are singletons
    refresh_ston_cols(col_sizes, ston_cols);

    // ------------------------------------------------------------------------
    // Loop over col_stons. We index the column with k and the row with i.
    // ------------------------------------------------------------------------
    for (kk = 0; kk < ston_cols->len; ++kk)
    {
        k = ston_cols->data[kk];
        i = AT->i[AT->p[k].start];
        Aik = AT->x[AT->p[k].start];

        assert(AT->p[k].end - AT->p[k].start == 1);
        assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));

#ifndef NDEBUG
        if (IS_HUGE(Aik) || IS_HUGE(1 / Aik))
        {
            printf("debug warning: large coefficient when eliminating col ston\n");
        }
#endif

        // If two column singletons appear in the same constraint, the
        // constraint will be inactive when the second column singleton is
        // processed. Also skip column singletons appearing in singleton rows.
        // It is necessary to both check for inactivity and the row size here.
        if (HAS_TAG(row_tags[i], R_TAG_INACTIVE) || row_sizes[i] <= 1)
        {
            continue;
        }

        vals = A->x + A->p[i].start;
        cols = A->i + A->p[i].start;
        row_size = row_sizes + i;

        ConstRowView const_row_view = new_const_rowview(
            vals, cols, row_size, A->p + i, lhs + i, rhs + i, row_tags + i, i);

        impl_free_from_above = implied_free_from_above_by_row(
            Aik, &const_row_view, bounds[k].lb, bounds[k].ub, col_tags[k], acts + i,
            col_tags, bounds);

        impl_free_from_below = implied_free_from_below_by_row(
            Aik, &const_row_view, bounds[k].lb, bounds[k].ub, col_tags[k], acts + i,
            col_tags, bounds);

#ifndef NDEBUG
        acts[i].min = (acts[i].n_inf_min != 0) ? INVALID_ACT_DEBUG : acts[i].min;
        acts[i].max = (acts[i].n_inf_max != 0) ? INVALID_ACT_DEBUG : acts[i].max;
#endif

        RowView row_view = new_rowview(vals, cols, row_size, A->p + i, lhs + i,
                                       rhs + i, row_tags + i, i);

        ColView col_view =
            new_colview(NULL, NULL, col_sizes + k, AT->p + k, &bounds[k].lb,
                        &bounds[k].ub, col_tags + k, k);

        if (HAS_TAG(row_tags[i], R_TAG_EQ))
        {

            status |= process_colston_eq(
                &row_view, &col_view, prob->obj, Aik, impl_free_from_above,
                impl_free_from_below, acts + i, locks, rows_to_delete, A_nnz, bounds,
                ston_rows, dton_rows, postsolve_info);
        }
        else
        {
            status |= process_colston_ineq(
                &row_view, &col_view, prob->obj, Aik, impl_free_from_above,
                impl_free_from_below, acts + i, locks, rows_to_delete, A_nnz, bounds,
                ston_rows, dton_rows, postsolve_info);
        }
    }

    // process the rows that should be deleted
    constraints->AT->nnz = A->nnz;
    delete_inactive_rows(constraints);

    if (HAS_STATUS(status, UNBNDORINFEAS))
    {
        return UNBNDORINFEAS;
    }
    else if (HAS_STATUS(status, REDUCED))
    {
        return REDUCED;
    }
    else
    {
        assert(status == UNCHANGED);
        return UNCHANGED;
    }
}

PresolveStatus remove_ston_cols(Problem *prob)
{
    assert(prob->constraints->state->ston_rows->len == 0);
    assert(prob->constraints->state->empty_rows->len == 0);
    assert(prob->constraints->state->empty_cols->len == 0);
    DEBUG(verify_problem_up_to_date(prob->constraints));
    DEBUG(verify_no_duplicates_sort(prob->constraints->state->ston_cols));

    PresolveStatus status = UNCHANGED;

    do
    {
        status = remove_ston_cols__(prob);
    } while (status == REDUCED);

    DEBUG(verify_problem_up_to_date(prob->constraints));

    // status will always be UNBNDORINFEAS or UNCHANGED (never REDUCED).
    assert(status != REDUCED);
    return status;
}
