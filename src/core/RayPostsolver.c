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

#include "Bounds.h"
#include "Constraints.h"
#include "Matrix.h"
#include "Numerics.h"
#include "Postsolver.h"
#include "State.h"
#include "dVec.h"
#include "debug_macros.h"
#include "glbopts.h"
#include "iVec.h"
#include "u16Vec.h"
#include <PSLP_warnings.h>
#include <assert.h>

static inline void copy_reduced_primal_ray_to_orginal(Solution *sol, const double *y,
                                                      const double *z,
                                                      const int *col_map,
                                                      const int *row_map)
{
    size_t dim_x = sol->dim_x;
    for (size_t i = 0; i < dim_x; ++i)
    {
        if (col_map[i] == -1)
        {
            sol->z[i] = COL_NOT_RETRIEVED;
            continue;
        }

        sol->z[i] = z[col_map[i]];
    }

    size_t dim_y = sol->dim_y;
    for (size_t i = 0; i < dim_y; ++i)
    {
        if (row_map[i] == -1)
        {
            sol->y[i] = ROW_NOT_RETRIEVED;
            continue;
        }

        sol->y[i] = y[row_map[i]];
    }
}

static void retrieve_fix_col_primal_ray(Solution *sol, int col,
                                        const double *ak_vals, const int *ak_rows,
                                        int len)
{
    assert(sol->z[col] == COL_NOT_RETRIEVED);

    sol->z[col] = 0.0;
    for (int i = 0; i < len; ++i)
    {
        assert(sol->y[ak_rows[i]] != ROW_NOT_RETRIEVED);
        sol->z[col] -= ak_vals[i] * sol->y[ak_rows[i]];
    }
}

static void retrieve_sub_col_primal_ray(Solution *sol, int k, const int *cols,
                                        const double *vals, int len, int i)
{
    assert(sol->y[i] != ROW_NOT_RETRIEVED);
    assert(sol->z[k] == COL_NOT_RETRIEVED);

    double aik = 0.0;
    for (int ii = 0; ii < len; ++ii)
    {
        if (cols[ii] == k)
        {
            aik = vals[ii];
            break;
        }
    }

    sol->z[k] = -aik * sol->y[i];
}

static void retrieve_fix_col_inf_primal_ray(Solution *sol, const int *indices)
{
    int col = indices[1];
    assert(sol->z[col] == COL_NOT_RETRIEVED);
    sol->z[col] = 0.0;
}

static void retrieve_parallel_col_primal_ray(Solution *sol, const int *indices,
                                             const double *vals)
{
    int j = indices[0];
    int k = indices[1];
    double ratio = vals[4];
    assert(indices[4] == DUMMY_VALUE);
    assert(sol->z[j] != COL_NOT_RETRIEVED && sol->z[k] == COL_NOT_RETRIEVED);

    sol->z[k] = ratio * sol->z[j];
}

static void retrieve_bound_change_primal_ray(Solution *sol, int i, int j,
                                             const int *cols, const double *vals,
                                             int len,
                                             int is_original_other_bound_lower_bound)
{
    int k, ii;
    double aij = 1.0;
    bool implied_bound_is_upper = is_original_other_bound_lower_bound;
    assert(sol->y[i] != ROW_NOT_RETRIEVED);
    assert(sol->z[j] != COL_NOT_RETRIEVED);

    // If the ray does not use the implied bound, we do not have to retrieve
    // anything from the row that implied it.
    if ((implied_bound_is_upper && sol->z[j] <= 0.0) ||
        (!implied_bound_is_upper && sol->z[j] >= 0.0))
    {
        return;
    }

    // find aij
    for (ii = 0; ii < len; ++ii)
    {
        if (cols[ii] == j)
        {
            aij = vals[ii];
            break;
        }
    }
    assert(ii != len && aij != 0.0);

    // update yi for row i that was used in the bound change
    sol->y[i] += sol->z[j] / aij;

    // update zk for all variables k appearing in row i
    for (ii = 0; ii < len; ++ii)
    {
        k = cols[ii];
        if (k == j)
        {
            continue;
        }

        // for now we cheat for dual postsolve of primal propagation fixing variables
        if (sol->z[k] == COL_NOT_RETRIEVED)
        {
            continue;
        }

        sol->z[k] -= (vals[ii] / aij) * sol->z[j];
    }

    sol->z[j] = 0.0;
}

// If both sides of row i were tightened from the same deleted row j, the pair
// (i, j) has two reductions and row j may already carry the multiplier that the
// other one transferred. We must therefore leave y[j] alone when the ray does
// not use this side.
static void retrieve_lhs_change_primal_ray(Solution *sol, int i, int j, double ratio)
{
    assert(sol->y[i] != ROW_NOT_RETRIEVED);

    if (sol->y[j] == ROW_NOT_RETRIEVED)
    {
        sol->y[j] = 0.0;
    }

    if (sol->y[i] >= 0.0)
    {
        return;
    }

    sol->y[j] += ratio * sol->y[i];
    sol->y[i] = 0.0;
}

static void retrieve_rhs_change_primal_ray(Solution *sol, int i, int j, double ratio)
{
    assert(sol->y[i] != ROW_NOT_RETRIEVED);

    if (sol->y[j] == ROW_NOT_RETRIEVED)
    {
        sol->y[j] = 0.0;
    }

    if (sol->y[i] <= 0.0)
    {
        return;
    }

    sol->y[j] += ratio * sol->y[i];
    sol->y[i] = 0.0;
}

static inline void copy_reduced_dual_ray_to_orginal(Solution *sol, const double *x,
                                                    const int *col_map)
{
    size_t dim_x = sol->dim_x;
    for (size_t i = 0; i < dim_x; ++i)
    {
        if (col_map[i] == -1)
        {
            sol->x[i] = COL_NOT_RETRIEVED;
            continue;
        }

        sol->x[i] = x[col_map[i]];
    }
}

static void retrieve_fix_col_dual_ray(Solution *sol, int col)
{
    assert(sol->x[col] == COL_NOT_RETRIEVED);
    sol->x[col] = 0.0;
}

static void retrieve_sub_col_dual_ray(Solution *sol, int k, const int *cols,
                                      const double *vals, int len)
{
    assert(sol->x[k] == COL_NOT_RETRIEVED);

    sol->x[k] = 0.0;
    double aik = 0.0;
    for (int ii = 0; ii < len; ++ii)
    {
        if (cols[ii] == k)
        {
            aik = vals[ii];
            continue;
        }

        assert(sol->x[cols[ii]] != COL_NOT_RETRIEVED);
        sol->x[k] -= vals[ii] * sol->x[cols[ii]];
    }

    sol->x[k] /= aik;
}

static void retrieve_fix_col_inf_dual_ray(Solution *sol, const int *indices,
                                          const double *vals)
{
    int i, j, counter, row_len;
    const int *cols;
    double coeff = 0;
    double val, side;
    const double *coeffs;
    int n_rows = (int) vals[0];
    double extreme_val = 0.0;
    bool fix_to_pos_inf = (indices[0] > 0);
    int col = indices[1];
    assert(sol->x[col] == COL_NOT_RETRIEVED);

    counter = 2;
    for (i = 0; i < n_rows; ++i)
    {
        side = 0.0;
        coeffs = vals + counter + 1;
        row_len = indices[counter];
        cols = indices + counter + 1;
        counter += row_len + 1;

        for (j = 0; j < row_len; ++j)
        {
            if (cols[j] == col)
            {
                coeff = coeffs[j];
                continue;
            }

            //  If two columns are fixed to pos inf in the same row, we
            //  pretend one of them is zero while we compute the other one
            if (sol->x[cols[j]] == COL_NOT_RETRIEVED)
            {
                continue;
            }

            assert(sol->x[cols[j]] != COL_NOT_RETRIEVED);
            side -= coeffs[j] * sol->x[cols[j]];
        }

        val = side / coeff;
        if (fix_to_pos_inf)
        {
            extreme_val = MAX(extreme_val, val);
        }
        else
        {
            extreme_val = MIN(extreme_val, val);
        }
    }

    assert(!IS_ABS_INF(extreme_val));
    sol->x[col] = extreme_val;
}

static void retrieve_parallel_col_dual_ray(Solution *sol, const int *indices,
                                           const double *vals)
{
    int j = indices[0];
    int k = indices[1];
    assert(indices[4] == DUMMY_VALUE);
    double lb_j = vals[0];
    double ub_j = vals[1];
    double lb_k = vals[2];
    double ub_k = vals[3];
    double ratio = vals[4];
    assert(sol->x[j] != COL_NOT_RETRIEVED && sol->x[k] == COL_NOT_RETRIEVED);
    double x_new_sol = sol->x[j];
    double xk_val;

    if ((x_new_sol >= -FEAS_TOL && x_new_sol <= FEAS_TOL) ||
        (x_new_sol < -FEAS_TOL && IS_NEG_INF(lb_j)) ||
        (x_new_sol > FEAS_TOL && IS_POS_INF(ub_j)))
    {
        sol->x[j] = x_new_sol;
        sol->x[k] = 0.0;
        return;
    }

    xk_val = x_new_sol / ratio;
    if ((xk_val >= -FEAS_TOL && xk_val <= FEAS_TOL) ||
        (xk_val < -FEAS_TOL && IS_NEG_INF(lb_k)) ||
        (xk_val > FEAS_TOL && IS_POS_INF(ub_k)))
    {
        sol->x[j] = 0.0;
        sol->x[k] = xk_val;
        return;
    }

    assert(false);
}

void postsolver_run_primal_infeas_ray(const PostsolveInfo *info, Solution *sol,
                                      const double *y, const double *z)
{
    const int *col_map = info->col_map;
    const int *row_map = info->row_map;
    int n_reductions = (int) info->type->len;
    ReductionType *reductions = info->type->data;
    const int *indices = info->indices->data;
    const double *vals = info->vals->data;
    const int *starts = info->starts->data;
    assert(n_reductions == info->starts->len - 1);
    ReductionType type;
    int start, len;

    copy_reduced_primal_ray_to_orginal(sol, y, z, col_map, row_map);

    for (int i = n_reductions - 1; i >= 0; --i)
    {
        type = reductions[i];
        start = starts[i];

        if (type == FIXED_COL)
        {
            len = starts[i + 1] - start - 2;
            // we might have fixed an empty column so len can be 0
            assert(len >= 0);
            retrieve_fix_col_primal_ray(sol, indices[start], vals + start + 2,
                                        indices + start + 2, len);
        }
        else if (type == SUB_COL)
        {
            len = starts[i + 1] - start - 2;
            assert(len > 1);
            retrieve_sub_col_primal_ray(sol, indices[start], indices + start + 1,
                                        vals + start + 1, len,
                                        indices[start + len + 1]);
        }
        else if (type == FIXED_COL_INF)
        {
            retrieve_fix_col_inf_primal_ray(sol, indices + start);
        }
        else if (type == PARALLEL_COL)
        {
            assert(starts[i + 1] - start == 5);
            retrieve_parallel_col_primal_ray(sol, indices + start, vals + start);
        }
        else if (type == DELETED_ROW)
        {
            assert(starts[i + 1] - start == 1);
            retrieve_deleted_row(sol, indices[start], 0.0);
        }
        else if (type == ADDED_ROW)
        {
            assert(starts[i + 1] - start == 2);
            retrieve_added_row(sol, indices + start, vals + start);
        }
        else if (type == ADDED_ROWS)
        {
            len = starts[i + 1] - start - 1;
            assert(len >= 1);
            retrieve_added_rows(sol, indices[start], indices + start + 1,
                                vals + start + 1, len, vals[start]);
        }
        else if (type == BOUND_CHANGE_THE_ROW)
        {
            // get the row that was used to derive the bound changes
            len = starts[i + 1] - start - 1;
            assert(len > 0);
            int num_of_bound_changes = (int) vals[start];
            int row = indices[start];
            const int *row_cols = indices + start + 1;
            const double *row_vals = vals + start + 1;
            int bound_changes_processed = 0;
            int j = i - 1;

            while (bound_changes_processed < num_of_bound_changes)
            {
                type = reductions[j];
                start = starts[j];
                assert(type == BOUND_CHANGE_NO_ROW || type == FIXED_COL);

                if (type == FIXED_COL)
                {
                    retrieve_fix_col_primal_ray(
                        sol, indices[start], vals + start + 2, indices + start + 2,
                        starts[j + 1] - start - 2);
                    assert(reductions[j - 1] == BOUND_CHANGE_NO_ROW);
                }
                else
                {
                    bound_changes_processed += 1;
                    retrieve_bound_change_primal_ray(sol, row, indices[start],
                                                     row_cols, row_vals, len,
                                                     indices[start + 1]);
                }

                j -= 1;
            }

            i = j + 1;
            assert(i >= 0);
            assert(i == 0 || reductions[i - 1] != BOUND_CHANGE_NO_ROW);
        }

        else if (type == LHS_CHANGE)
        {
            len = starts[i + 1] - start - 2;
            assert(len > 1);
            retrieve_lhs_change_primal_ray(sol, indices[start], indices[start + 1],
                                           vals[start + 1]);
        }
        else if (type == RHS_CHANGE)
        {
            len = starts[i + 1] - start - 2;
            assert(len > 1);
            retrieve_rhs_change_primal_ray(sol, indices[start], indices[start + 1],
                                           vals[start + 1]);
        }
        else if (type == EQ_TO_INEQ)
        {
            assert(starts[i + 1] - start == 1);
        }
        else
        {
            assert(type != BOUND_CHANGE_NO_ROW);
            assert(false);
        }
    }

#ifndef NDEBUG
    for (int i = 0; i < sol->dim_x; ++i)
    {
        if (sol->z[i] == COL_NOT_RETRIEVED)
        {
            printf("col %d not fully retrieved \n", i);
        }

        assert(sol->z[i] != COL_NOT_RETRIEVED);
    }

    for (int i = 0; i < sol->dim_y; ++i)
    {
        assert(sol->y[i] != ROW_NOT_RETRIEVED);
    }
#endif
}

void postsolver_run_dual_infeas_ray(const PostsolveInfo *info, Solution *sol,
                                    const double *x)
{
    const int *col_map = info->col_map;
    int n_reductions = (int) info->type->len;
    ReductionType *reductions = info->type->data;
    const int *indices = info->indices->data;
    const double *vals = info->vals->data;
    const int *starts = info->starts->data;
    assert(n_reductions == info->starts->len - 1);
    ReductionType type;
    int start, len;

    copy_reduced_dual_ray_to_orginal(sol, x, col_map);

    for (int i = n_reductions - 1; i >= 0; --i)
    {
        type = reductions[i];
        start = starts[i];

        if (type == FIXED_COL)
        {
            len = starts[i + 1] - start - 2;
            // we might have fixed an empty column so len can be 0
            assert(len >= 0);
            retrieve_fix_col_dual_ray(sol, indices[start]);
        }
        else if (type == SUB_COL)
        {
            len = starts[i + 1] - start - 2;
            assert(len > 1);
            retrieve_sub_col_dual_ray(sol, indices[start], indices + start + 1,
                                      vals + start + 1, len);
        }
        else if (type == FIXED_COL_INF)
        {
            retrieve_fix_col_inf_dual_ray(sol, indices + start, vals + start);
        }
        else if (type == PARALLEL_COL)
        {
            assert(starts[i + 1] - start == 5);
            retrieve_parallel_col_dual_ray(sol, indices + start, vals + start);
        }
        else if (type == DELETED_ROW)
        {
            assert(starts[i + 1] - start == 1);
        }
        else if (type == ADDED_ROW)
        {
            assert(starts[i + 1] - start == 2);
        }
        else if (type == ADDED_ROWS)
        {
            len = starts[i + 1] - start - 1;
            assert(len >= 1);
        }
        else if (type == BOUND_CHANGE_THE_ROW)
        {
            // get the row that was used to derive the bound changes
            len = starts[i + 1] - start - 1;
            assert(len > 0);
            int num_of_bound_changes = (int) vals[start];
            int bound_changes_processed = 0;
            int j = i - 1;

            while (bound_changes_processed < num_of_bound_changes)
            {
                type = reductions[j];
                start = starts[j];
                assert(type == BOUND_CHANGE_NO_ROW || type == FIXED_COL);

                if (type == FIXED_COL)
                {
                    retrieve_fix_col_dual_ray(sol, indices[start]);
                    assert(reductions[j - 1] == BOUND_CHANGE_NO_ROW);
                }
                else
                {
                    bound_changes_processed += 1;
                }

                j -= 1;
            }

            i = j + 1;
            assert(i >= 0);
            assert(i == 0 || reductions[i - 1] != BOUND_CHANGE_NO_ROW);
        }

        else if (type == LHS_CHANGE)
        {
            len = starts[i + 1] - start - 2;
            assert(len > 1);
        }
        else if (type == RHS_CHANGE)
        {
            len = starts[i + 1] - start - 2;
            assert(len > 1);
        }
        else if (type == EQ_TO_INEQ)
        {
            assert(starts[i + 1] - start == 1);
        }
        else
        {
            assert(type != BOUND_CHANGE_NO_ROW);
            assert(false);
        }
    }

#ifndef NDEBUG
    for (int i = 0; i < sol->dim_x; ++i)
    {
        if (sol->x[i] == COL_NOT_RETRIEVED)
        {
            printf("col %d not fully retrieved \n", i);
        }

        assert(sol->x[i] != COL_NOT_RETRIEVED);
    }
#endif
}
