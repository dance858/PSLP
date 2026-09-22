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

static void retrieve_fix_col_inf_primal_ray(Solution *sol,
                                            const FixedColInfRecord *r)
{
    int col = r->col;
    assert(sol->z[col] == COL_NOT_RETRIEVED);
    sol->z[col] = 0.0;
}

static void retrieve_parallel_col_primal_ray(Solution *sol,
                                             const ParallelColRecord *r)
{
    int j = r->j;
    int k = r->k;
    double ratio = r->ratio;
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

static void retrieve_fix_col_inf_dual_ray(Solution *sol, const FixedColInfRecord *r)
{
    int i, j, pos;
    double coeff = 0;
    double val, side;
    int n_rows = r->n_rows;
    double extreme_val = 0.0;
    bool fix_to_pos_inf = r->pos_inf;
    int col = r->col;
    assert(sol->x[col] == COL_NOT_RETRIEVED);

    pos = 0;
    for (i = 0; i < n_rows; ++i)
    {
        // a ray has no right-hand side, so row.side is not used
        FixedColInfRow row = fixed_col_inf_row(r, &pos);
        side = 0.0;

        for (j = 0; j < row.len; ++j)
        {
            if (row.cols[j] == col)
            {
                coeff = row.coeffs[j];
                continue;
            }

            //  If two columns are fixed to pos inf in the same row, we
            //  pretend one of them is zero while we compute the other one
            if (sol->x[row.cols[j]] == COL_NOT_RETRIEVED)
            {
                continue;
            }

            assert(sol->x[row.cols[j]] != COL_NOT_RETRIEVED);
            side -= row.coeffs[j] * sol->x[row.cols[j]];
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

static void retrieve_parallel_col_dual_ray(Solution *sol, const ParallelColRecord *r)
{
    int j = r->j;
    int k = r->k;
    double lb_j = r->lb_j;
    double ub_j = r->ub_j;
    double lb_k = r->lb_k;
    double ub_k = r->ub_k;
    double ratio = r->ratio;
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
        len = starts[i + 1] - start;

        if (type == FIXED_COL)
        {
            FixedColRecord r = decode_fixed_col(indices + start, vals + start, len);
            retrieve_fix_col_primal_ray(sol, r.col, r.vals, r.rows, r.len);
        }
        else if (type == SUB_COL)
        {
            SubColRecord r = decode_sub_col(indices + start, vals + start, len);
            retrieve_sub_col_primal_ray(sol, r.k, r.cols, r.vals, r.len, r.row);
        }
        else if (type == FIXED_COL_INF)
        {
            FixedColInfRecord r =
                decode_fixed_col_inf(indices + start, vals + start, len);
            retrieve_fix_col_inf_primal_ray(sol, &r);
        }
        else if (type == PARALLEL_COL)
        {
            ParallelColRecord r =
                decode_parallel_col(indices + start, vals + start, len);
            retrieve_parallel_col_primal_ray(sol, &r);
        }
        else if (type == DELETED_ROW)
        {
            // the stored multiplier is not used for a ray
            DeletedRowRecord r =
                decode_deleted_row(indices + start, vals + start, len);
            retrieve_deleted_row(sol, r.row, 0.0);
        }
        else if (type == ADDED_ROW)
        {
            AddedRowRecord r = decode_added_row(indices + start, vals + start, len);
            retrieve_added_row(sol, r.i, r.j, r.ratio);
        }
        else if (type == ADDED_ROWS)
        {
            AddedRowsRecord r =
                decode_added_rows(indices + start, vals + start, len);
            retrieve_added_rows(sol, r.i, r.rows, r.vals, r.len, r.aik);
        }
        else if (type == BOUND_CHANGE_THE_ROW)
        {
            // get the row that was used to derive the bound changes; the
            // bound changes themselves are the preceding records
            BoundChangeTheRowRecord row =
                decode_bound_change_the_row(indices + start, vals + start, len);
            int bound_changes_processed = 0;
            int j = i - 1;

            while (bound_changes_processed < row.num_of_bound_changes)
            {
                type = reductions[j];
                start = starts[j];
                len = starts[j + 1] - start;
                assert(type == BOUND_CHANGE_NO_ROW || type == FIXED_COL);

                if (type == FIXED_COL)
                {
                    FixedColRecord r =
                        decode_fixed_col(indices + start, vals + start, len);
                    retrieve_fix_col_primal_ray(sol, r.col, r.vals, r.rows, r.len);
                    assert(reductions[j - 1] == BOUND_CHANGE_NO_ROW);
                }
                else
                {
                    BoundChangeNoRowRecord r = decode_bound_change_no_row(
                        indices + start, vals + start, len);
                    bound_changes_processed += 1;
                    retrieve_bound_change_primal_ray(
                        sol, row.i, r.j, row.cols, row.vals, row.len,
                        r.is_original_other_bound_lower_bound);
                }

                j -= 1;
            }

            i = j + 1;
            assert(i >= 0);
            assert(i == 0 || reductions[i - 1] != BOUND_CHANGE_NO_ROW);
        }
        else if (type == LHS_CHANGE)
        {
            SideChangeRecord r =
                decode_side_change(indices + start, vals + start, len);
            retrieve_lhs_change_primal_ray(sol, r.i, r.j, r.ratio);
        }
        else if (type == RHS_CHANGE)
        {
            SideChangeRecord r =
                decode_side_change(indices + start, vals + start, len);
            retrieve_rhs_change_primal_ray(sol, r.i, r.j, r.ratio);
        }
        else if (type == EQ_TO_INEQ)
        {
            // nothing to do for a ray; decoded to check the record
            (void) decode_eq_to_ineq(indices + start, vals + start, len);
        }
        else if (type == PARALLEL_ROW)
        {
            // only used by postsolver_map_to_reduced; decoded to check the record
            (void) decode_parallel_row(indices + start, vals + start, len);
        }
        else if (type == SIDE_RELAXED)
        {
            // only used by postsolver_map_to_reduced; decoded to check the record
            (void) decode_side_relaxed(indices + start, vals + start, len);
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
        len = starts[i + 1] - start;

        if (type == FIXED_COL)
        {
            FixedColRecord r = decode_fixed_col(indices + start, vals + start, len);
            retrieve_fix_col_dual_ray(sol, r.col);
        }
        else if (type == SUB_COL)
        {
            SubColRecord r = decode_sub_col(indices + start, vals + start, len);
            retrieve_sub_col_dual_ray(sol, r.k, r.cols, r.vals, r.len);
        }
        else if (type == FIXED_COL_INF)
        {
            FixedColInfRecord r =
                decode_fixed_col_inf(indices + start, vals + start, len);
            retrieve_fix_col_inf_dual_ray(sol, &r);
        }
        else if (type == PARALLEL_COL)
        {
            ParallelColRecord r =
                decode_parallel_col(indices + start, vals + start, len);
            retrieve_parallel_col_dual_ray(sol, &r);
        }
        else if (type == DELETED_ROW)
        {
            // nothing to do for a dual ray; decoded to check the record
            (void) decode_deleted_row(indices + start, vals + start, len);
        }
        else if (type == ADDED_ROW)
        {
            (void) decode_added_row(indices + start, vals + start, len);
        }
        else if (type == ADDED_ROWS)
        {
            (void) decode_added_rows(indices + start, vals + start, len);
        }
        else if (type == BOUND_CHANGE_THE_ROW)
        {
            // the bound changes themselves are the preceding records; only
            // the fixed columns among them matter for a dual ray
            BoundChangeTheRowRecord row =
                decode_bound_change_the_row(indices + start, vals + start, len);
            int bound_changes_processed = 0;
            int j = i - 1;

            while (bound_changes_processed < row.num_of_bound_changes)
            {
                type = reductions[j];
                start = starts[j];
                len = starts[j + 1] - start;
                assert(type == BOUND_CHANGE_NO_ROW || type == FIXED_COL);

                if (type == FIXED_COL)
                {
                    FixedColRecord r =
                        decode_fixed_col(indices + start, vals + start, len);
                    retrieve_fix_col_dual_ray(sol, r.col);
                    assert(reductions[j - 1] == BOUND_CHANGE_NO_ROW);
                }
                else
                {
                    (void) decode_bound_change_no_row(indices + start, vals + start,
                                                      len);
                    bound_changes_processed += 1;
                }

                j -= 1;
            }

            i = j + 1;
            assert(i >= 0);
            assert(i == 0 || reductions[i - 1] != BOUND_CHANGE_NO_ROW);
        }
        else if (type == LHS_CHANGE || type == RHS_CHANGE)
        {
            (void) decode_side_change(indices + start, vals + start, len);
        }
        else if (type == EQ_TO_INEQ)
        {
            (void) decode_eq_to_ineq(indices + start, vals + start, len);
        }
        else if (type == PARALLEL_ROW)
        {
            (void) decode_parallel_row(indices + start, vals + start, len);
        }
        else if (type == SIDE_RELAXED)
        {
            (void) decode_side_relaxed(indices + start, vals + start, len);
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
