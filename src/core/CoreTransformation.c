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

#include "Activity.h"
#include "Bounds.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Locks.h"
#include "Matrix.h"
#include "Numerics.h"
#include "Problem.h"
#include "RowColViews.h"
#include "State.h"

PresolveStatus fix_col(Constraints *constraints, int col, double val, double ck)
{
    State *data = constraints->state;
    const Matrix *A = constraints->A;
    Bound *bounds = constraints->bounds;
    double old_lb = bounds[col].lb;
    double old_ub = bounds[col].ub;
    ColTag *col_tag = constraints->col_tags + col;
    bool is_ub_inf = HAS_TAG(*col_tag, C_TAG_UB_INF);
    bool is_lb_inf = HAS_TAG(*col_tag, C_TAG_LB_INF);
    assert(constraints->state->col_sizes[col] > 0);
    assert(!HAS_TAG(*col_tag, C_TAG_INACTIVE));

    if ((!is_ub_inf && val >= old_ub + FEAS_TOL) ||
        (!is_lb_inf && val <= old_lb - FEAS_TOL))
    {
        return INFEASIBLE;
    }

    const Matrix *AT = constraints->AT;
    const double *vals = AT->x + AT->p[col].start;
    const int *rows = AT->i + AT->p[col].start;
    int len = AT->p[col].end - AT->p[col].start;

    set_col_to_fixed(col, col_tag, data->fixed_cols_to_delete);
    save_retrieval_fixed_col(data->postsolve_info, col, val, ck, vals, rows, len);

    bool ub_update = (is_ub_inf || val != old_ub);
    bool lb_update = (is_lb_inf || val != old_lb);
    bounds[col].ub = val;
    bounds[col].lb = val;

    update_activities_fixed_col(data->activities, A, bounds, vals, rows, len, old_ub,
                                old_lb, val, ub_update, lb_update, is_ub_inf,
                                is_lb_inf, data->updated_activities);

    return REDUCED;
}

void fix_col_to_negative_inf(Constraints *constraints, int col)
{
    assert(!HAS_TAG(constraints->col_tags[col], C_TAG_INACTIVE));
    assert(HAS_TAG(constraints->col_tags[col], C_TAG_LB_INF));
    assert(constraints->state->col_locks[col].down == 0);

    Matrix *AT = constraints->AT;
    RowTag *row_tags = constraints->row_tags;
    PostsolveInfo *postsolve_info = constraints->state->postsolve_info;

    // store postsolve information
    save_retrieval_fixed_col_inf(postsolve_info, col, -1, constraints,
                                 constraints->bounds[col].ub);

    for (int i = AT->p[col].start; i < AT->p[col].end; ++i)
    {
        int row = AT->i[i];
        // if the same row has two variables that we fix to neg inf, it
        // will be marked as inactive when the second variable is fixed
        if (!HAS_TAG(row_tags[row], R_TAG_INACTIVE))
        {
            set_row_to_inactive(row, row_tags + row,
                                constraints->state->rows_to_delete, postsolve_info,
                                0.0);
        }
    }

    // It is not needed to append the column to fixed_cols_to_delete
    // since all rows that the column appears in will be deleted
    UPDATE_TAG(constraints->col_tags[col], C_TAG_FIXED);
    constraints->state->col_sizes[col] = SIZE_INACTIVE_COL;
    AT->p[col].end = AT->p[col].start;
    constraints->bounds[col].lb = -INF;
    constraints->bounds[col].ub = -INF;
}

void fix_col_to_positive_inf(Constraints *constraints, int col)
{
    assert(!HAS_TAG(constraints->col_tags[col], C_TAG_INACTIVE));
    assert(HAS_TAG(constraints->col_tags[col], C_TAG_UB_INF));
    assert(constraints->state->col_locks[col].up == 0);

    Matrix *AT = constraints->AT;
    RowTag *row_tags = constraints->row_tags;
    PostsolveInfo *postsolve_info = constraints->state->postsolve_info;

    // store postsolve information
    save_retrieval_fixed_col_inf(postsolve_info, col, 1, constraints,
                                 constraints->bounds[col].lb);

    for (int i = AT->p[col].start; i < AT->p[col].end; ++i)
    {
        int row = AT->i[i];

        // if the same row has two variables that we fix to pos inf, it
        // will be marked as inactive when the second variable is fixed
        if (!HAS_TAG(row_tags[row], R_TAG_INACTIVE))
        {
            set_row_to_inactive(row, row_tags + row,
                                constraints->state->rows_to_delete, postsolve_info,
                                0.0);
        }
    }

    // It is not needed to append the column to fixed_cols_to_delete
    // since all rows that the column appears in will be deleted
    UPDATE_TAG(constraints->col_tags[col], C_TAG_FIXED);
    constraints->state->col_sizes[col] = SIZE_INACTIVE_COL;
    AT->p[col].end = AT->p[col].start;
    constraints->bounds[col].lb = INF;
    constraints->bounds[col].ub = INF;
}

PresolveStatus update_lb(Constraints *constraints, int col, double new_lb,
                         int row HUGE_BOUND_PARAM)
{
    assert(!IS_HUGE(new_lb) || huge_bound_ok);
    ColTag *col_tag = constraints->col_tags + col;
    Bound *bounds = constraints->bounds;
    double ub = bounds[col].ub;
    double lb = bounds[col].lb;
    assert(!HAS_TAG(constraints->col_tags[col], C_TAG_INACTIVE));

    if (!HAS_TAG(*col_tag, C_TAG_UB_INF) && new_lb >= ub + FEAS_TOL)
    {
        return INFEASIBLE;
    }

    bool is_lb_inf = HAS_TAG(*col_tag, C_TAG_LB_INF);

    // update the lower bound if it is tighter
    if (is_lb_inf || new_lb >= lb + FEAS_TOL)
    {
        const Matrix *AT = constraints->AT;
        const double *vals = AT->x + AT->p[col].start;
        const int *rows = AT->i + AT->p[col].start;
        int len = AT->p[col].end - AT->p[col].start;
        State *data = constraints->state;
        Activity *activities = data->activities;

        bounds[col].lb = new_lb;
        update_activities_bound_change(activities, constraints->A, bounds, vals,
                                       rows, len, lb, new_lb, !is_lb_inf, true,
                                       data->updated_activities HUGE_BOUND_ARG);
        REMOVE_TAG(*col_tag, C_TAG_LB_INF);

        const Matrix *A = constraints->A;
        const double *row_vals = A->x + A->p[row].start;
        const int *row_cols = A->i + A->p[row].start;
        int row_len = A->p[row].end - A->p[row].start;

        save_retrieval_bound_change_no_row(data->postsolve_info, col, new_lb, ub,
                                           false);
        save_retrieval_bound_change_the_row(data->postsolve_info, row, row_cols,
                                            row_vals, row_len, 1);
        return REDUCED;
    }

    return UNCHANGED;
}

PresolveStatus update_ub(Constraints *constraints, int col, double new_ub,
                         int row HUGE_BOUND_PARAM)
{
    assert(!IS_HUGE(new_ub) || huge_bound_ok);
    ColTag *col_tag = constraints->col_tags + col;
    Bound *bounds = constraints->bounds;
    double ub = bounds[col].ub;
    double lb = bounds[col].lb;

    assert(!HAS_TAG(constraints->col_tags[col], C_TAG_INACTIVE));

    if (!HAS_TAG(*col_tag, C_TAG_LB_INF) && new_ub <= lb - FEAS_TOL)
    {
        return INFEASIBLE;
    }

    bool is_ub_inf = HAS_TAG(*col_tag, C_TAG_UB_INF);

    // update the upper bound if it is tighter
    if (is_ub_inf || new_ub <= ub - FEAS_TOL)
    {
        const Matrix *AT = constraints->AT;
        const double *vals = AT->x + AT->p[col].start;
        const int *rows = AT->i + AT->p[col].start;
        int len = AT->p[col].end - AT->p[col].start;
        State *data = constraints->state;
        Activity *activities = data->activities;

        bounds[col].ub = new_ub;
        update_activities_bound_change(activities, constraints->A, bounds, vals,
                                       rows, len, ub, new_ub, !is_ub_inf, false,
                                       data->updated_activities HUGE_BOUND_ARG);
        REMOVE_TAG(*col_tag, C_TAG_UB_INF);

        const Matrix *A = constraints->A;
        const double *row_vals = A->x + A->p[row].start;
        const int *row_cols = A->i + A->p[row].start;
        int row_len = A->p[row].end - A->p[row].start;

        save_retrieval_bound_change_no_row(data->postsolve_info, col, new_ub, lb,
                                           true);

        save_retrieval_bound_change_the_row(data->postsolve_info, row, row_cols,
                                            row_vals, row_len, 1);
    }

    return UNCHANGED;
}

bool implied_free_from_above_by_row(double Aik, const ConstRowView *row, double lb,
                                    double ub, ColTag col_tag, const Activity *act,
                                    const ColTag *col_tags, const Bound *bounds)
{
    if (HAS_TAG(col_tag, C_TAG_UB_INF))
    {
        return true;
    }

    assert(!IS_POS_INF(ub) && Aik != 0);
    double implied_ub = INF;

    if (Aik > 0 && !HAS_TAG(*row->tag, R_TAG_RHS_INF))
    {
        // standard case
        if (act->n_inf_min == 0)
        {
            assert(!IS_ABS_INF(act->min) && !IS_ABS_INF(*row->rhs) &&
                   !IS_ABS_INF(lb));
            implied_ub = lb + (*row->rhs - act->min) / Aik;
        }
        // Gondzio's case (important for implied free column singletons)
        else if (act->n_inf_min == 1 && HAS_TAG(col_tag, C_TAG_LB_INF))
        {
            double min_temp = compute_min_act_tags(row->vals, row->cols, *row->len,
                                                   bounds, col_tags);
            assert(!IS_ABS_INF(min_temp) && !IS_ABS_INF(*row->rhs));
            implied_ub = (*row->rhs - min_temp) / Aik;
        }
    }
    else if (Aik < 0 && !HAS_TAG(*row->tag, R_TAG_LHS_INF))
    {
        if (act->n_inf_max == 0)
        {
            assert(!IS_ABS_INF(act->max) && !IS_ABS_INF(*row->lhs) &&
                   !IS_ABS_INF(lb));
            implied_ub = lb + (*row->lhs - act->max) / Aik;
        }
        // Gondzio's case (important for implied free column singletons)
        else if (act->n_inf_max == 1 && HAS_TAG(col_tag, C_TAG_LB_INF))
        {
            double max_temp = compute_max_act_tags(row->vals, row->cols, *row->len,
                                                   bounds, col_tags);
            assert(!IS_ABS_INF(max_temp) && !IS_ABS_INF(*row->lhs));
            implied_ub = (*row->lhs - max_temp) / Aik;
        }
    }

    assert(!IS_NEG_INF(implied_ub));
    return (implied_ub - ub <= 0);
}

bool implied_free_from_below_by_row(double Aik, const ConstRowView *row, double lb,
                                    double ub, ColTag col_tag, const Activity *act,
                                    const ColTag *col_tags, const Bound *bounds)
{
    if (HAS_TAG(col_tag, C_TAG_LB_INF))
    {
        return true;
    }

    assert(Aik != 0);
    assert(lb != -INF);
    double implied_lb = -INF;

    if (Aik > 0 && !HAS_TAG(*row->tag, R_TAG_LHS_INF))
    {
        // standard case
        if (act->n_inf_max == 0)
        {
            assert(!IS_ABS_INF(act->max) && !IS_ABS_INF(*row->lhs) &&
                   !IS_ABS_INF(ub));
            implied_lb = ub + (*row->lhs - act->max) / Aik;
        }
        // Gondzio's case (important for implied free column singletons)
        else if (act->n_inf_max == 1 && HAS_TAG(col_tag, C_TAG_UB_INF))
        {
            double max_temp = compute_max_act_tags(row->vals, row->cols, *row->len,
                                                   bounds, col_tags);

            assert(!IS_ABS_INF(max_temp) && !IS_ABS_INF(*row->lhs));
            implied_lb = (*row->lhs - max_temp) / Aik;
        }
    }
    else if (Aik < 0 && !HAS_TAG(*row->tag, R_TAG_RHS_INF))
    {
        if (act->n_inf_min == 0)
        {
            assert(!IS_ABS_INF(act->min) && !IS_ABS_INF(*row->rhs) &&
                   !IS_ABS_INF(ub));
            implied_lb = ub + (*row->rhs - act->min) / Aik;
        }
        // Gondzio's case (important for implied free column singletons)
        else if (act->n_inf_min == 1 && HAS_TAG(col_tag, C_TAG_UB_INF))
        {
            double min_temp = compute_min_act_tags(row->vals, row->cols, *row->len,
                                                   bounds, col_tags);
            assert(!IS_ABS_INF(min_temp) && !IS_ABS_INF(*row->rhs));
            implied_lb = (*row->rhs - min_temp) / Aik;
        }
    }

    printf("implied_lb: %f, lb: %f\n", implied_lb, lb);
    assert(!IS_POS_INF(implied_lb));
    return (lb - implied_lb <= 0);
}

PresolveStatus change_rhs_of_ineq(RowView *row, double new_rhs, Lock *locks,
                                  iVec *dton_rows, PostsolveInfo *postsolve_info,
                                  int other_row_idx, double ratio)
{
    assert(!HAS_TAG(*row->tag, (R_TAG_EQ | R_TAG_INACTIVE)));
    assert(!IS_ABS_INF(new_rhs));
    assert(other_row_idx >= 0);

    if (!HAS_TAG(*row->tag, R_TAG_LHS_INF) && *row->lhs >= new_rhs + FEAS_TOL)
    {
        return INFEASIBLE;
    }

    bool is_rhs_inf = HAS_TAG(*row->tag, R_TAG_RHS_INF);
    assert(is_rhs_inf || *row->rhs > new_rhs);

    // -------------------------------------------------------------------------
    //                          update locks
    // -------------------------------------------------------------------------
    int len = *row->len;
    if (is_rhs_inf)
    {
        for (int i = 0; i < len; ++i)
        {
            if (row->vals[i] > 0)
            {
                locks[row->cols[i]].up++;
            }
            else
            {
                locks[row->cols[i]].down++;
            }
        }
    }

    // update row tag and the actual rhs
    *row->rhs = new_rhs;
    REMOVE_TAG(*row->tag, R_TAG_RHS_INF);

    // check if the row becomes an equality constraint
    if (!HAS_TAG(*row->tag, R_TAG_LHS_INF) &&
        IS_EQUAL_FEAS_TOL(*row->lhs, *row->rhs))
    {
        *row->rhs = *row->lhs;
        UPDATE_TAG(*row->tag, R_TAG_EQ);

        if (len == 2)
        {
            iVec_append(dton_rows, row->i);
        }
    }

    // dual postsolve
    save_retrieval_rhs_or_lhs_change(postsolve_info, row->i, row->vals, row->cols,
                                     len, new_rhs, other_row_idx, ratio, false);

    return REDUCED;
}

PresolveStatus change_lhs_of_ineq(RowView *row, double new_lhs, Lock *locks,
                                  iVec *dton_rows, PostsolveInfo *postsolve_info,
                                  int other_row_idx, double ratio)
{
    assert(!HAS_TAG(*row->tag, (R_TAG_EQ | R_TAG_INACTIVE)));
    assert(!IS_ABS_INF(new_lhs));
    assert(other_row_idx >= 0);

    if (!HAS_TAG(*row->tag, R_TAG_RHS_INF) && *row->rhs <= new_lhs - FEAS_TOL)
    {
        return INFEASIBLE;
    }

    bool is_lhs_inf = HAS_TAG(*row->tag, R_TAG_LHS_INF);
    assert(is_lhs_inf || *row->lhs < new_lhs);

    // -------------------------------------------------------------------------
    //                          update locks
    // -------------------------------------------------------------------------
    int len = *row->len;
    if (is_lhs_inf)
    {
        for (int i = 0; i < len; ++i)
        {
            if (row->vals[i] > 0)
            {
                locks[row->cols[i]].down++;
            }
            else
            {
                locks[row->cols[i]].up++;
            }
        }
    }

    // update row tag and the actual lhs
    *row->lhs = new_lhs;
    REMOVE_TAG(*row->tag, R_TAG_LHS_INF);

    // check if the row becomes an equality constraint
    if (!HAS_TAG(*row->tag, R_TAG_RHS_INF) &&
        IS_EQUAL_FEAS_TOL(*row->rhs, *row->lhs))
    {
        *row->lhs = *row->rhs;
        UPDATE_TAG(*row->tag, R_TAG_EQ);

        if (len == 2)
        {
            iVec_append(dton_rows, row->i);
        }
    }

    // dual postsolve
    save_retrieval_rhs_or_lhs_change(postsolve_info, row->i, row->vals, row->cols,
                                     len, new_lhs, other_row_idx, ratio, true);

    return REDUCED;
}

bool is_ub_implied_free(const ConstColView *col, Activity *acts, const double *lhs,
                        const double *rhs, const RowTag *row_tags, const Matrix *A,
                        const ColTag *col_tags, const Bound *bounds)
{
    assert(!HAS_TAG(*col->tag, (C_TAG_INACTIVE | C_TAG_UB_INF)));

    for (int ii = 0; ii < *col->len; ++ii)
    {
        int i = col->rows[ii];
        assert(!HAS_TAG(row_tags[i], R_TAG_INACTIVE));

        int len = A->p[i].end - A->p[i].start;
        ConstRowView row =
            new_const_rowview(A->x + A->p[i].start, A->i + A->p[i].start, &len,
                              A->p + i, lhs + i, rhs + i, row_tags + i, i);

        if (implied_free_from_above_by_row(col->vals[ii], &row, *col->lb, *col->ub,
                                           *col->tag, acts + i, col_tags, bounds))
        {
            DEBUG(acts[i].min =
                      (acts[i].n_inf_min != 0) ? INVALID_ACT_DEBUG : acts[i].min);
            DEBUG(acts[i].max =
                      (acts[i].n_inf_max != 0) ? INVALID_ACT_DEBUG : acts[i].max);
            return true;
        }

        DEBUG(acts[i].min =
                  (acts[i].n_inf_min != 0) ? INVALID_ACT_DEBUG : acts[i].min);
        DEBUG(acts[i].max =
                  (acts[i].n_inf_max != 0) ? INVALID_ACT_DEBUG : acts[i].max);
    }

    return false;
}

bool is_lb_implied_free(const ConstColView *col, Activity *acts, const double *lhs,
                        const double *rhs, const RowTag *row_tags, const Matrix *A,
                        const ColTag *col_tags, const Bound *bounds)

{
    assert(!HAS_TAG(*col->tag, (C_TAG_INACTIVE | C_TAG_LB_INF)));

    for (int ii = 0; ii < *col->len; ++ii)
    {
        int i = col->rows[ii];
        assert(!HAS_TAG(row_tags[i], R_TAG_INACTIVE));

        int len = A->p[i].end - A->p[i].start;
        ConstRowView row =
            new_const_rowview(A->x + A->p[i].start, A->i + A->p[i].start, &len,
                              A->p + i, lhs + i, rhs + i, row_tags + i, i);

        assert(col->vals[ii] != 0);
        if (implied_free_from_below_by_row(col->vals[ii], &row, *col->lb, *col->ub,
                                           *col->tag, acts + i, col_tags, bounds))
        {
            DEBUG(acts[i].min =
                      (acts[i].n_inf_min != 0) ? INVALID_ACT_DEBUG : acts[i].min);
            DEBUG(acts[i].max =
                      (acts[i].n_inf_max != 0) ? INVALID_ACT_DEBUG : acts[i].max);
            return true;
        }

        DEBUG(acts[i].min =
                  (acts[i].n_inf_min != 0) ? INVALID_ACT_DEBUG : acts[i].min);
        DEBUG(acts[i].max =
                  (acts[i].n_inf_max != 0) ? INVALID_ACT_DEBUG : acts[i].max);
    }

    return false;
}
