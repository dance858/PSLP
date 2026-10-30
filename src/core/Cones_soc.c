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

#include "Cones_soc.h"
#include "Bounds.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Matrix.h"
#include "Numerics.h"
#include "Postsolver_cones.h"
#include "Problem.h"
#include "State.h"
#include "glbopts.h"
#include "iVec.h"

#include <assert.h>
#include <math.h>

/* largest kappa = c_j / (a omega) a rotated block may put into P */
#define MAX_KAPPA 1e8

/* A row that still belongs to its block. Test the tag, not the size: a row
   deleted on this visit keeps its size until delete_inactive_rows. */
static bool row_active(const Constraints *constraints, int row)
{
    return !HAS_TAG(constraints->row_tags[row], R_TAG_INACTIVE);
}

/* Marks the cone row 'row' inactive with a zero multiplier. */
static void delete_cone_row(Constraints *constraints, int row)
{
    RESET_TAG(constraints->row_tags[row], R_TAG_INACTIVE);
    iVec_append(constraints->state->rows_to_delete, row);
    save_retrieval_deleted_row(constraints->state->postsolve_info, row, 0.0);
}

/* Deletes the zero tails of the block (empty rows with a zero side) and
   returns the number of tails left. */
static int delete_zero_tails(Constraints *constraints, const ConeBlock *block,
                             int n_heads)
{
    const int *row_sizes = constraints->state->row_sizes;
    const double *rhs = constraints->rhs;
    int n_tails = 0;

    for (int t = block->start + n_heads; t < block->start + block->dim; ++t)
    {
        if (!row_active(constraints, t))
        {
            continue;
        }
        if (row_sizes[t] == 0 && IS_ZERO_FEAS_TOL(rhs[t]))
        {
            delete_cone_row(constraints, t);
        }
        else
        {
            n_tails++;
        }
    }
    return n_tails;
}

/* every head empty or a single entry */
static bool simple_heads(const Constraints *constraints, const ConeBlock *block,
                         int n_heads)
{
    const int *row_sizes = constraints->state->row_sizes;
    for (int h = block->start; h < block->start + n_heads; ++h)
    {
        if (row_sizes[h] > 1)
        {
            return false;
        }
    }
    return true;
}

/* Turns head row h, a x_j - r_h >= 0 (or -r_h >= 0 when empty), into a bound
   on x_j. */
static PresolveStatus head_to_bound(Constraints *constraints, int h)
{
    const Matrix *A = constraints->A;
    double r_h = constraints->rhs[h];
    if (constraints->state->row_sizes[h] == 0)
    {
        return IS_INFEAS_GT(r_h, 0.0) ? INFEASIBLE : UNCHANGED;
    }

    int j = A->i[A->p[h].start];
    double a = A->x[A->p[h].start];
    return a > 0.0 ? update_lb(constraints, j, r_h / a, h HUGE_BOUND_IS_OK)
                   : update_ub(constraints, j, r_h / a, h HUGE_BOUND_IS_OK);
}

/* A block without tails is its heads alone, s_0 >= 0 and s_1 >= 0 (SOC dim 1
   and ROTATED dim 2 are valid blocks, so no tail is kept). Each head becomes
   a bound and goes. */
static PresolveStatus heads_to_bounds(Constraints *constraints,
                                      const ConeBlock *block, int n_heads)
{
    for (int h = block->start; h < block->start + n_heads; ++h)
    {
        RETURN_IF_INFEASIBLE(head_to_bound(constraints, h));
        delete_cone_row(constraints, h);
    }
    return UNCHANGED;
}

/* checks if every active tail of a rotated block is a single entry */
static bool single_entry_tails(const Constraints *constraints,
                               const ConeBlock *block)
{
    const int *row_sizes = constraints->state->row_sizes;
    for (int t = block->start + 2; t < block->start + block->dim; ++t)
    {
        if (row_active(constraints, t) && row_sizes[t] != 1)
        {
            return false;
        }
    }
    return true;
}

/* get x_j of the head row h if the block converts with w as the constant head,
   else -1 */
static int get_head_column(const Problem *prob, int w, int h)
{
    const Constraints *constraints = prob->constraints;
    const State *state = constraints->state;
    const Matrix *A = constraints->A;
    if (state->row_sizes[w] != 0 || state->row_sizes[h] != 1)
    {
        return -1;
    }

    int j = A->i[A->p[h].start];
    double a = A->x[A->p[h].start];
    double ratio = prob->obj->c[j] / a;
    double omega = -constraints->rhs[w];
    if (state->col_sizes[j] != 1 || col_in_P(state, j) || ratio <= 0.0 ||
        omega < FEAS_TOL || ratio / omega > MAX_KAPPA)
    {
        return -1;
    }

    /* the cone gives a x_j >= r_h, so x_j must be free in the other direction
       and its bound on this side no tighter than r_h / a */
    Bound b = constraints->bounds[j];
    double bound = constraints->rhs[h] / a;
    if (a > 0.0)
    {
        if (IS_POS_INF(b.ub) && (IS_NEG_INF(b.lb) || b.lb <= bound))
        {
            return j;
        }
        return -1;
    }
    if (IS_NEG_INF(b.lb) && (IS_POS_INF(b.ub) || b.ub >= bound))
    {
        return j;
    }
    return -1;
}

static void convert_block(Problem *prob, const ConeBlock *block, int w, int h, int j,
                          int n_tails)
{
    Constraints *constraints = prob->constraints;
    State *state = constraints->state;
    const Matrix *A = constraints->A;
    const double *rhs = constraints->rhs;
    double *c = prob->obj->c;
    double a = A->x[A->p[h].start];
    double omega = -rhs[w];
    double kappa = c[j] / (a * omega);
    int end = block->start + block->dim;

    save_retrieval_rotated_to_p(state->postsolve_info, constraints, c[j], j, w, h,
                                block, n_tails);

    double r2 = 0.0;
    for (int t = block->start + 2; t < end; ++t)
    {
        if (!row_active(constraints, t))
        {
            continue;
        }
        int k = A->i[A->p[t].start];
        double f = A->x[A->p[t].start];
        diagonal_P_add(prob, k, kappa * f * f);
        c[k] -= kappa * f * rhs[t];
        r2 += rhs[t] * rhs[t];
    }
    prob->obj->offset += c[j] * rhs[h] / a + 0.5 * kappa * r2;

    UPDATE_TAG(constraints->col_tags[j], C_TAG_SUBSTITUTED);
    state->col_sizes[j] = SIZE_INACTIVE_COL;
    constraints->AT->p[j].end = constraints->AT->p[j].start;
    mark_P_row_inactive(prob, j);

    /* The record covers the block, so its rows go without their own records.
       A row already listed must not be listed twice. */
    for (int row = block->start; row < end; ++row)
    {
        if (row_active(constraints, row))
        {
            RESET_TAG(constraints->row_tags[row], R_TAG_INACTIVE);
            iVec_append(state->rows_to_delete, row);
        }
    }
}

/* Replaces a rotated block 2 s_0 s_1 >= sum_t (f_t x_k - r_t)^2 with constant
   head omega and the singleton column x_j of the other head a x_j - r_h by
   the term (kappa / 2) sum_t (f_t x_k - r_t)^2 of P, kappa = c_j / (a omega). */
static bool rotated_to_P_block(Problem *prob, const ConeBlock *block, int n_tails)
{
    if (!single_entry_tails(prob->constraints, block))
    {
        return false;
    }

    for (int kw = 0; kw < 2; ++kw)
    {
        int w = block->start + kw, h = block->start + 1 - kw;
        int j = get_head_column(prob, w, h);
        if (j >= 0 && (prob->obj->P != NULL || diagonal_P_alloc(prob)))
        {
            /* the heads are simple, so a block without tails went to bounds */
            assert(n_tails >= 1);
            convert_block(prob, block, w, h, j, n_tails);
            return true;
        }
    }
    return false;
}

PresolveStatus process_soc_block(Problem *prob, const ConeBlock *block,
                                 bool rotated_to_P)
{
    Constraints *constraints = prob->constraints;
    int n_heads = block->type == PSLP_CONE_SOC ? 1 : 2;

    int n_tails = delete_zero_tails(constraints, block, n_heads);
    if (n_tails == 0 && simple_heads(constraints, block, n_heads))
    {
        return heads_to_bounds(constraints, block, n_heads);
    }

    if (rotated_to_P && block->type == PSLP_CONE_ROTATED)
    {
        rotated_to_P_block(prob, block, n_tails);
    }
    return UNCHANGED;
}
