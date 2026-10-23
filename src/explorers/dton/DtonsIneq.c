/*
 * Copyright 2025-2026 Daniel Cederberg
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
#include "DtonsEq_internal.h"
#include "Locks.h"
#include "Matrix.h"
#include "Problem.h"
#include "RowColViews.h"
#include "State.h"
#include "Tags.h"
#include "glbopts.h"

/* Returns true (for push_col_down == True) when xj may be decreased without
   making the cost worse and when only one row holds it back to be decreased
   infinitely. */
static bool dton_pushable(const Problem *prob, int j, bool push_col_down)
{
    const Lock *lock = prob->constraints->state->col_locks + j;
    double cost = prob->obj->c[j];
    return push_col_down ? (cost >= 0 && lock->down == 1)
                         : (cost <= 0 && lock->up == 1);
}

/* An inequality can, under certain conditions, be made an equality by pushing
    increasing or decreasing the value of a column in the row. These condtions
    are (1) The column must be (implied) free in the direction of the push. (2)
    The push must not raise the objective function. (3) The row must be the only
    row that holds the column back in that direction.

    This function returns the entry, 0 or 1, and its push direction, or -1 when
    neither column qualifies. */
static int dton_dominated_entry(const Problem *prob, const ConstRowView *row,
                                bool rhs_inf, bool *push_col_down)
{
    const Constraints *constraints = prob->constraints;
    const State *state = constraints->state;
    const Bound *bounds = constraints->bounds;
    const ColTag *col_tags = constraints->col_tags;
    const Activity *act = state->activities + row->i;

    for (int ii = 0; ii < 2; ++ii)
    {
        int j = row->cols[ii];
        double a = row->vals[ii];
        if (col_in_P(state, j))
        {
            continue;
        }

        /* if the constraint is infinite on the right and the coefficient is
           positive, we can make the constraint tighter by decreasing the value
           of the column (ie., push it down) */
        *push_col_down = rhs_inf ? a > 0 : a < 0;
        if (*push_col_down)
        {
            if (dton_pushable(prob, j, true) &&
                implied_free_from_below_by_row(a, row, bounds[j].lb, bounds[j].ub,
                                               col_tags[j], act, col_tags, bounds))
            {
                return ii;
            }
        }
        else
        {
            if (dton_pushable(prob, j, false) &&
                implied_free_from_above_by_row(a, row, bounds[j].lb, bounds[j].ub,
                                               col_tags[j], act, col_tags, bounds))
            {
                return ii;
            }
        }
    }
    return -1;
}

bool dton_tighten_dominated_row(Problem *prob, int i)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT;
    State *state = constraints->state;
    Bound *bounds = constraints->bounds;
    ColTag *col_tags = constraints->col_tags;
    RowTag *tag = constraints->row_tags + i;
    bool lhs_inf = HAS_TAG(*tag, R_TAG_LHS_INF);
    bool rhs_inf = HAS_TAG(*tag, R_TAG_RHS_INF);
    const double *vals = A->x + A->p[i].start;
    const int *cols = A->i + A->p[i].start;
    int len = 2;

    assert(lhs_inf != rhs_inf);
    assert(state->row_sizes[i] == 2 && A->p[i].end - A->p[i].start == 2);
    ConstRowView row =
        new_const_rowview(vals, cols, &len, A->p + i, constraints->lhs + i,
                          constraints->rhs + i, tag, i);
    bool push_col_down;
    int ii = dton_dominated_entry(prob, &row, rhs_inf, &push_col_down);

    /* if we don't find a dominated entry, we cannot tighten the row */
    if (ii < 0)
    {
        return false;
    }

    // -------------------------------------------------------------------------
    //                Make the constraint an equality
    // -------------------------------------------------------------------------
    if (rhs_inf)
    {
        constraints->rhs[i] = constraints->lhs[i];
    }
    else
    {
        constraints->lhs[i] = constraints->rhs[i];
    }
    RESET_TAG(*tag, R_TAG_EQ);
    update_locks_ineq_to_eq(state->col_locks, vals, cols, 2, lhs_inf);

    // -------------------------------------------------------------------------
    // Drops the bound of the dominated column on the side the row implies, and
    // keep the activity bookkeeping consistent.
    // -------------------------------------------------------------------------
    int j = cols[ii];
    if (!HAS_TAG(col_tags[j], push_col_down ? C_TAG_LB_INF : C_TAG_UB_INF))
    {
        ConstColView col = new_const_colview(
            AT->x + AT->p[j].start, AT->i + AT->p[j].start, state->col_sizes + j,
            AT->p + j, &bounds[j].lb, &bounds[j].ub, col_tags + j, j);
        if (push_col_down)
        {
            remove_finite_lb_from_activities(&col, state->activities, bounds[j].lb);
            bounds[j].lb = -INF;
        }
        else
        {
            remove_finite_ub_from_activities(&col, state->activities, bounds[j].ub);
            bounds[j].ub = INF;
        }
        UPDATE_TAG(col_tags[j], push_col_down ? C_TAG_LB_INF : C_TAG_UB_INF);
    }

    /* append the row since it is now an equality */
    iVec_append(state->dton_rows, i);
    return true;
}

/* Phase 0 for one column: the one row that holds it back in a direction its
   cost does not object to is made tight when the rule allows. */
static void dton_tighten_column(Problem *prob, int j)
{
    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT;
    const State *state = constraints->state;
    const RowTag *row_tags = constraints->row_tags;

    if (HAS_TAG(constraints->col_tags[j], C_TAG_INACTIVE) || col_in_P(state, j))
    {
        return;
    }
    bool down = dton_pushable(prob, j, true);
    bool up = dton_pushable(prob, j, false);
    if (!down && !up)
    {
        return;
    }

    for (int p = AT->p[j].start; p < AT->p[j].end; ++p)
    {
        int i = AT->i[p];
        bool lhs_inf = HAS_TAG(row_tags[i], R_TAG_LHS_INF);
        bool rhs_inf = HAS_TAG(row_tags[i], R_TAG_RHS_INF);

        /* we just want one-sided doubleton inequalities */
        if (lhs_inf == rhs_inf || state->row_sizes[i] != 2)
        {
            continue;
        }

        /* A >= row holds prevents a column with a positive coefficient to
           decrease infinitely, since the constraint gets tightened  */
        bool push_col_down = rhs_inf ? AT->x[p] > 0 : AT->x[p] < 0;
        if (push_col_down ? down : up)
        {
            dton_tighten_dominated_row(prob, i);
        }
    }
}

void dton_tighten_columns(Problem *prob)
{
    int n = (int) prob->constraints->n;
    for (int j = 0; j < n; ++j)
    {
        dton_tighten_column(prob, j);
    }
}

void dton_tighten_changed_columns(Problem *prob, const DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT;
    const State *state = constraints->state;
    const RowTag *row_tags = constraints->row_tags;

    /* The round changed the entries and bounds of these columns, so only the
       two-entry rows holding one of them can test differently now. */
    for (int ii = 0; ii < dton_work->targets.n_targets; ++ii)
    {
        int j = dton_work->targets.list[ii].col;
        for (int p = AT->p[j].start; p < AT->p[j].end; ++p)
        {
            int i = AT->i[p];
            bool lhs_inf = HAS_TAG(row_tags[i], R_TAG_LHS_INF);
            bool rhs_inf = HAS_TAG(row_tags[i], R_TAG_RHS_INF);
            if (lhs_inf != rhs_inf && state->row_sizes[i] == 2)
            {
                dton_tighten_dominated_row(prob, i);
            }
        }
    }
}
