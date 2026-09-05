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

/* Maps a primal-dual point of the original problem to the reduced problem
   (the inverse direction of postsolve), e.g. to warm start a solver on the
   reduced problem. The public entry point is 'map_solution_to_reduced' in
   PSLP_API.h; the core replay of the recorded reductions is
   'postsolver_map_to_reduced', declared in Postsolver.h. */

#include "Constraints.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "PSLP_API.h"
#include "PSLP_sol.h"
#include "Postsolver.h"
#include "Problem.h"
#include "State.h"
#include "Workspace.h"
#include "dVec.h"
#include "iVec.h"
#include "u16Vec.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

void postsolver_map_to_reduced(const PostsolveInfo *info, const int *col_map,
                               const int *row_map, size_t n_cols_orig,
                               size_t n_rows_orig, const double *x, const double *y,
                               double *x_work, double *y_work, double *x_red,
                               double *y_red)
{
    int n_reductions = (int) info->type->len;
    const ReductionType *reductions = info->type->data;
    const int *indices = info->indices->data;
    const double *vals = info->vals->data;
    const int *starts = info->starts->data;
    bool map_primal = (x_red != NULL);
    bool map_dual = (y_red != NULL);
    assert(n_reductions == (int) info->starts->len - 1);
    assert(!map_primal || (x != NULL && x_work != NULL));
    assert(!map_dual || (y != NULL && y_work != NULL));

    if (map_primal)
    {
        memcpy(x_work, x, n_cols_orig * sizeof(double));
    }

    if (map_dual)
    {
        memcpy(y_work, y, n_rows_orig * sizeof(double));
    }

    // Each recorded reduction is undone by 'postsolver_run' using an affine
    // update of (x, y). Here we apply the inverse of these updates in the order
    // the reductions were made, but only for reductions that modify a row or
    // column that survives in the reduced problem. Removed rows and columns are
    // dropped by the maps at the end, and z is not mapped at all since the
    // caller recomputes it from the reduced problem data.
    for (int t = 0; t < n_reductions; ++t)
    {
        ReductionType type = reductions[t];
        int start = starts[t];

        if (type == PARALLEL_COL)
        {
            // (xj, xk) was replaced by x_new = xj + ratio * xk
            if (map_primal)
            {
                int j = indices[start];
                int k = indices[start + 1];
                double ratio = vals[start + 4];
                x_work[j] += ratio * x_work[k];
            }
        }
        else if (type == PARALLEL_ROW)
        {
            // parallel row j (ai = ratio * aj) was removed in favour of row i.
            // Since ai^T yi + aj^T yj = ai^T (yi + yj / ratio), the multiplier
            // of row j is transferred to row i.
            if (map_dual)
            {
                int i = indices[start];
                int j = indices[start + 1];
                double ratio = vals[start];
                y_work[i] += y_work[j] / ratio;
                y_work[j] = 0.0;
            }
        }
        else if (type == EQ_TO_INEQ)
        {
            // postsolve does yi += ck / aik
            if (map_dual)
            {
                y_work[indices[start]] -= vals[start];
            }
        }
        else if (type == ADDED_ROW || type == ADDED_ROWS)
        {
            // row i was added to other rows to eliminate a column. Postsolve
            // adjusts yi, but row i itself is always removed afterwards (see
            // DtonsEq.c and SimpleReductions.c), so there is nothing to map.
            assert(row_map[indices[start]] == -1);
        }
        else if (type == LHS_CHANGE || type == RHS_CHANGE)
        {
            // the transfer of the removed row's multiplier is captured by the
            // PARALLEL_ROW record that follows
        }
        else if (type == FIXED_COL || type == FIXED_COL_INF || type == SUB_COL ||
                 type == DELETED_ROW)
        {
            // only affects a removed row or column
        }
        else if (type == BOUND_CHANGE_NO_ROW || type == BOUND_CHANGE_THE_ROW)
        {
            // bound changes do not touch A, so (x, y) carry over unchanged
        }
        else
        {
            assert(false);
        }
    }

    if (map_primal)
    {
        for (size_t i = 0; i < n_cols_orig; ++i)
        {
            if (col_map[i] != -1)
            {
                x_red[col_map[i]] = x_work[i];
            }
        }
    }

    if (map_dual)
    {
        for (size_t i = 0; i < n_rows_orig; ++i)
        {
            if (row_map[i] != -1)
            {
                y_red[row_map[i]] = y_work[i];
            }
        }
    }
}

/* z = c - A^T y for the reduced problem */
static void compute_reduced_dual_slack(double *z, const PresolvedProblem *prob,
                                       const double *y)
{
    for (size_t j = 0; j < prob->n; ++j)
    {
        z[j] = prob->c[j];
    }

    for (size_t i = 0; i < prob->m; ++i)
    {
        for (int p = prob->Ap[i]; p < prob->Ap[i + 1]; ++p)
        {
            z[prob->Ai[p]] -= prob->Ax[p] * y[i];
        }
    }
}

void map_solution_to_reduced(Presolver *presolver, const double *x, const double *y,
                             double *x_red, double *y_red, double *z_red)
{
    Solution *sol = presolver->sol;
    State *data = presolver->prob->constraints->state;
    const Mapping *maps = data->work->mappings;
    const PresolvedProblem *reduced_prob = presolver->reduced_prob;
    size_t n_cols_orig = sol->dim_x;
    size_t n_rows_orig = sol->dim_y;
    double *work;

    if (x_red && !x)
    {
        fprintf(stderr, "PSLP warning: map_solution_to_reduced needs x to compute "
                        "x_red! x_red is left untouched. \n");
        x_red = NULL;
    }

    if (y_red && !y)
    {
        fprintf(stderr, "PSLP warning: map_solution_to_reduced needs y to compute "
                        "y_red! y_red and z_red are left untouched. \n");
        y_red = NULL;
        z_red = NULL;
    }

    if (z_red && !y_red)
    {
        fprintf(stderr, "PSLP warning: map_solution_to_reduced needs y_red to "
                        "compute z_red! z_red is left untouched. \n");
        z_red = NULL;
    }

    if (!x_red && !y_red)
    {
        return;
    }

    // scratch space of the original dimensions (we do not use presolver->sol
    // since it may hold a solution obtained by postsolve)
    work = (double *) ps_malloc(n_cols_orig + n_rows_orig + 1, sizeof(double));
    if (!work)
    {
        fprintf(stderr, "PSLP warning: map_solution_to_reduced failed to allocate "
                        "memory! Outputs are left untouched. \n");
        return;
    }

    postsolver_map_to_reduced(data->postsolve_info, maps->cols, maps->rows,
                              n_cols_orig, n_rows_orig, x, y, work,
                              work + n_cols_orig, x_red, y_red);
    PS_FREE(work);

    // project onto the (possibly tightened) bounds of the reduced problem. The
    // bounds are gone after 'free_presolver_reduced_problem', in which case
    // x_red is returned unprojected.
    if (x_red && reduced_prob->lbs && reduced_prob->ubs)
    {
        for (size_t j = 0; j < reduced_prob->n; ++j)
        {
            x_red[j] =
                MIN(MAX(x_red[j], reduced_prob->lbs[j]), reduced_prob->ubs[j]);
        }
    }

    if (!z_red)
    {
        return;
    }

    if (!reduced_prob->Ax)
    {
        // The reduced problem data was released by
        // 'free_presolver_reduced_problem', so z_red cannot be computed.
        fprintf(stderr, "PSLP warning: map_solution_to_reduced called after "
                        "free_presolver_reduced_problem! z_red is left "
                        "untouched. \n");
        return;
    }

    compute_reduced_dual_slack(z_red, reduced_prob, y_red);
}
