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
#include "Problem.h"
#include "State.h"
#include "glbopts.h"
#include "iVec.h"

#include <assert.h>

/* A row that still belongs to its block. Rows deleted in an earlier round and
   rows marked on this visit both carry the tag (the sizes differ). */
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

/* The head row h of a block without tails reads a x_j - r_h >= 0 (or
   -r_h >= 0 when the row is empty). Turns it into a bound on x_j. */
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

PresolveStatus process_soc_block(Problem *prob, const ConeBlock *block)
{
    Constraints *constraints = prob->constraints;
    int n_heads = block->type == PSLP_CONE_SOC ? 1 : 2;

    int n_tails = delete_zero_tails(constraints, block, n_heads);
    if (n_tails == 0 && simple_heads(constraints, block, n_heads))
    {
        return heads_to_bounds(constraints, block, n_heads);
    }
    return UNCHANGED;
}
