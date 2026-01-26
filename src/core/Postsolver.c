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

#include "Postsolver.h"
#include "Bounds.h"
#include "Constraints.h"
#include "Matrix.h"
#include "Numerics.h"
#include "State.h"
#include "dVec.h"
#include "debug_macros.h"
#include "glbopts.h"
#include "iVec.h"
#include "u16Vec.h"
#include <assert.h>

#define INIT_FRAC_POSTSOLVE 0.3
#define COL_NOT_RETRIEVED INF
#define ROW_NOT_RETRIEVED INF
#define DUMMY_VALUE -382749

PostsolveInfo *postsolve_info_new(int n_rows, int n_cols)
{
    PostsolveInfo *info = (PostsolveInfo *) ps_malloc(1, sizeof(PostsolveInfo));
    info->starts =
        iVec_new((size_t) (MAX(1, INIT_FRAC_POSTSOLVE * (n_rows + n_cols))));
    info->indices =
        iVec_new((size_t) (MAX(1, INIT_FRAC_POSTSOLVE * (n_rows + n_cols))));
    info->vals =
        dVec_new((size_t) (MAX(1, INIT_FRAC_POSTSOLVE * (n_rows + n_cols))));
    info->type =
        u16Vec_new((size_t) (MAX(1, INIT_FRAC_POSTSOLVE * (n_rows + n_cols))));
    if (!info->starts || !info->indices || !info->vals || !info->type)
    {
        PS_FREE(info->starts);
        PS_FREE(info->indices);
        PS_FREE(info->vals);
        PS_FREE(info->type);
        PS_FREE(info);
        return NULL;
    }

    iVec_append(info->starts, 0);
    return info;
}

void postsolve_info_free(PostsolveInfo *info)
{
    if (!info)
    {
        return;
    }

    iVec_free(info->starts);
    iVec_free(info->indices);
    dVec_free(info->vals);
    u16Vec_free(info->type);
    PS_FREE(info);
}

static inline void copy_reduced_sol_to_orginal(Solution *sol, const double *x,
                                               const double *y, const double *z,
                                               const int *col_map,
                                               const int *row_map)
{
    int dim_x = sol->dim_x;
    for (int i = 0; i < dim_x; ++i)
    {
        if (col_map[i] == -1)
        {
            sol->x[i] = COL_NOT_RETRIEVED;
            sol->z[i] = COL_NOT_RETRIEVED;
            continue;
        }

        sol->x[i] = x[col_map[i]];
        sol->z[i] = z[col_map[i]];
    }

    int dim_y = sol->dim_y;
    for (int i = 0; i < dim_y; ++i)
    {
        if (row_map[i] == -1)
        {
            sol->y[i] = ROW_NOT_RETRIEVED;
            continue;
        }

        sol->y[i] = y[row_map[i]];
    }
}

// fixes xk and recovers zk = ck - ak^T y
static void retrieve_fix_col(Solution *sol, int col, double val, double ck,
                             const double *ak_vals, const int *ak_rows, int len)
{
    assert(sol->x[col] == COL_NOT_RETRIEVED);
    assert(sol->z[col] == COL_NOT_RETRIEVED);

    sol->x[col] = val;
    sol->z[col] = ck;
    for (int i = 0; i < len; ++i)
    {
        assert(sol->y[ak_rows[i]] != ROW_NOT_RETRIEVED);
        sol->z[col] -= ak_vals[i] * sol->y[ak_rows[i]];
    }
}

// retrieves xk and zk
static void retrieve_sub_col(Solution *sol, int k, const int *cols,
                             const double *vals, int len, double rhs, int i,
                             double ck)
{
    assert(sol->y[i] != ROW_NOT_RETRIEVED);
    assert(sol->x[k] == COL_NOT_RETRIEVED);
    assert(sol->z[k] == COL_NOT_RETRIEVED);

    sol->x[k] = rhs;
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

    sol->z[k] = ck - aik * sol->y[i];
    sol->x[k] /= aik;
}

static void retrieve_fix_col_inf(Solution *sol, const int *indices,
                                 const double *vals)
{
    int i, j, counter, row_len;
    const int *cols;
    double coeff = 0;
    double val, side;
    const double *coeffs;
    int n_rows = (int) vals[0];
    double extreme_val = vals[1];
    bool fix_to_pos_inf = (indices[0] > 0);
    int col = indices[1];
    assert(sol->x[col] == COL_NOT_RETRIEVED);
    assert(sol->z[col] == COL_NOT_RETRIEVED);

    counter = 2;
    for (i = 0; i < n_rows; ++i)
    {
        side = vals[counter];
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
    sol->z[col] = 0.0;
}

static void retrieve_parallel_col(Solution *sol, const int *indices,
                                  const double *vals)
{
    int j = indices[0];
    int k = indices[1];
    ColTag cTag_j = (ColTag) indices[2];
    ColTag cTag_k = (ColTag) indices[3];
    assert(indices[4] == DUMMY_VALUE);
    double lb_j = vals[0];
    double ub_j = vals[1];
    double lb_k = vals[2];
    double ub_k = vals[3];
    double ratio = vals[4];
    assert(sol->x[j] != COL_NOT_RETRIEVED && sol->x[k] == COL_NOT_RETRIEVED);
    assert(sol->z[j] != COL_NOT_RETRIEVED && sol->z[k] == COL_NOT_RETRIEVED);
    double x_new_sol = sol->x[j];
    double xj_val, xk_val;

    // ----------------------------------------------------------------
    //          try to set xj equal to one of its bounds
    // ----------------------------------------------------------------
    if (!HAS_TAG(cTag_j, C_TAG_LB_INF))
    {
        xj_val = lb_j;
    }
    else if (!HAS_TAG(cTag_j, C_TAG_UB_INF))
    {
        xj_val = ub_j;
    }
    else
    {
        xj_val = 0;
    }

    // ----------------------------------------------------------------
    //     check if the corresponding value on xk is feasible
    // ----------------------------------------------------------------
    xk_val = (x_new_sol - xj_val) / ratio;

    if (!HAS_TAG(cTag_k, C_TAG_LB_INF) && lb_k >= xk_val + FEAS_TOL)
    {
        xk_val = lb_k;
        xj_val = x_new_sol - ratio * xk_val;
    }
    else if (!HAS_TAG(cTag_k, C_TAG_UB_INF) && xk_val >= ub_k + FEAS_TOL)
    {
        xk_val = ub_k;
        xj_val = x_new_sol - ratio * xk_val;
    }

    assert(!IS_ABS_INF(xj_val) && !IS_ABS_INF(xk_val));
    sol->x[j] = xj_val;
    sol->x[k] = xk_val;

    // dual postsolve
    bool is_xj_at_bound =
        (!HAS_TAG(cTag_j, C_TAG_LB_INF) && IS_EQUAL_FEAS_TOL(xj_val, lb_j)) ||
        (!HAS_TAG(cTag_j, C_TAG_UB_INF) && IS_EQUAL_FEAS_TOL(xj_val, ub_j));
    bool is_xk_at_bound =
        (!HAS_TAG(cTag_k, C_TAG_LB_INF) && IS_EQUAL_FEAS_TOL(xk_val, lb_k)) ||
        (!HAS_TAG(cTag_k, C_TAG_UB_INF) && IS_EQUAL_FEAS_TOL(xk_val, ub_k));

    if (is_xj_at_bound && is_xk_at_bound)
    {
        sol->z[k] = ratio * sol->z[j];
    }
    else
    {
        sol->z[k] = 0.0;
    }

#ifndef NDEBUG
    SET_ZERO_IF_SMALL_DUAL_POSTSOLVE(sol->z[j]);

    printf("sol->x[j]=%f, lb_j=%f\n", sol->x[j], lb_j);

    // if the multiplier is positive the variable should be at its lower bound
    if (sol->z[j] > 0)
    {
        printf("j=%d, z[j]=%f, x[j]=%f, lb_j=%f, ub_j=%f\n", j, sol->z[j], sol->x[j],
               lb_j, ub_j);
        printf("HAS_TAG(cTag_j, C_TAG_LB_INF)=%d\n", HAS_TAG(cTag_j, C_TAG_LB_INF));
        printf("sol->x[j]=%f, lb_j=%f\n", sol->x[j], lb_j);
        printf("\n \n");
        // assert(!HAS_TAG(cTag_j, C_TAG_LB_INF) && IS_EQUAL_FEAS_TOL(sol->x[j],
        // lb_j));
    }
    // if the multiplier is negative the variable should be at its upper bound
    else if (sol->z[j] < 0)
    {
        assert(!HAS_TAG(cTag_j, C_TAG_UB_INF) && IS_EQUAL_FEAS_TOL(sol->x[j], ub_j));
    }

    // similar checks for variable k
    if (sol->z[k] > 0)
    {
        assert(!HAS_TAG(cTag_k, C_TAG_LB_INF) && IS_EQUAL_FEAS_TOL(sol->x[k], lb_k));
    }
    else if (sol->z[k] < 0)
    {
        assert(!HAS_TAG(cTag_k, C_TAG_UB_INF) && IS_EQUAL_FEAS_TOL(sol->x[k], ub_k));
    }
#endif
}

void retrieve_deleted_row(Solution *sol, int row, double val)
{
    assert(sol->y[row] == ROW_NOT_RETRIEVED);
    sol->y[row] = val;
}

// yi = yi - (ajk / aik) yj whenever we make the jth constraint sparser by
// zeroing out xk using equality row i
void retrieve_added_row(Solution *sol, const int *rows, const double *vals)
{
    int i = rows[0];
    int j = rows[1];
    assert(sol->y[i] != ROW_NOT_RETRIEVED);
    assert(sol->y[j] != ROW_NOT_RETRIEVED);
    assert(vals[1] == DUMMY_VALUE);
    // it should be a plus sign because the minus sign is included in
    // save_retrieval_added_row
    sol->y[i] += sol->y[j] * vals[0];
}

void retrieve_added_rows(Solution *sol, int i, const int *rows, const double *vals,
                         int len, double aik)
{
    assert(sol->y[i] != ROW_NOT_RETRIEVED);

    for (int j = 0; j < len; ++j)
    {
        if (rows[j] != i)
        {
            assert(sol->y[rows[j]] != ROW_NOT_RETRIEVED);
            sol->y[i] -= (vals[j] / aik) * sol->y[rows[j]];
        }
    }
}

void retrieve_bound_change(Solution *sol, int i, int j, const int *cols,
                           const double *vals, int len, double implied_bound,
                           double original_other_bound,
                           int is_original_other_bound_lower_bound)
{
    int k, ii;
    double aij = 1.0;
    assert(sol->x[j] != COL_NOT_RETRIEVED);
    assert(sol->y[i] != ROW_NOT_RETRIEVED);
    assert(sol->z[j] != COL_NOT_RETRIEVED);

    // if the variable is equal to one of its original bounds and
    // has the right sign on its multiplier we do nothing
    if (IS_EQUAL_FEAS_TOL(sol->x[j], original_other_bound))
    {
        if ((is_original_other_bound_lower_bound && sol->z[j] >= 0) ||
            (!is_original_other_bound_lower_bound && sol->z[j] <= 0))
        {
            return;
        }
    }

    // if the implied bound is not active we do nothing
    if (!IS_EQUAL_FEAS_TOL(sol->x[j], implied_bound))
    {
        // sol->y[i] does not always have to be zero since we possibly got the
        // bound from a doubleton equality row assert(sol->y[i] == 0.0);
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
    //
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

        assert(sol->z[k] != COL_NOT_RETRIEVED);
        sol->z[k] -= (vals[ii] / aij) * sol->z[j];
    }

    sol->z[j] = 0.0;
}

void retrieve_lhs_change(Solution *sol, int i, const double *vals, const int *cols,
                         int len, double new_side, int j, double ratio)
{
    assert(sol->y[i] != ROW_NOT_RETRIEVED);
    assert(sol->y[j] == ROW_NOT_RETRIEVED || sol->y[j] == 0.0);

    // compute the residual of constraint i
    double residual = -new_side;
    for (int k = 0; k < len; ++k)
    {
        assert(sol->x[cols[k]] != COL_NOT_RETRIEVED);
        residual += vals[k] * sol->x[cols[k]];
    }

    if (!IS_ZERO_FEAS_TOL(residual))
    {
        return;
    }

    double cand_val = ratio * sol->y[i];
    if (ratio * cand_val < 0)
    {
        sol->y[j] = 0.0;
        return;
    }

    sol->y[j] = cand_val;
    sol->y[i] = 0.0;
}

void retrieve_rhs_change(Solution *sol, int i, const double *vals, const int *cols,
                         int len, double new_side, int j, double ratio)
{
    assert(sol->y[i] != ROW_NOT_RETRIEVED);
    assert(sol->y[j] == ROW_NOT_RETRIEVED || sol->y[j] == 0.0);

    // compute the residual of constraint i
    double residual = -new_side;
    for (int k = 0; k < len; ++k)
    {
        assert(sol->x[cols[k]] != COL_NOT_RETRIEVED);
        residual += vals[k] * sol->x[cols[k]];
    }

    if (!IS_ZERO_FEAS_TOL(residual))
    {
        return;
    }

    double cand_val = ratio * sol->y[i];
    if (ratio * cand_val > 0)
    {
        sol->y[j] = 0.0;
        return;
    }

    sol->y[j] = cand_val;
    sol->y[i] = 0.0;
}

void retrieve_eq_to_ineq(Solution *sol, int row, double val)
{
    assert(sol->y[row] != ROW_NOT_RETRIEVED);
    sol->y[row] += val;
}

void postsolver_update(PostsolveInfo *info, int n_cols_reduced, int n_rows_reduced,
                       const int *col_map, const int *row_map)
{
    info->n_cols_reduced = n_cols_reduced;
    info->n_rows_reduced = n_rows_reduced;
    info->col_map = col_map;
    info->row_map = row_map;
}

void polish_z(Solution *sol, const double *lbs, const double *ubs)
{
    int n_cols = sol->dim_x;
    double *x = sol->x;
    double *z = sol->z;
    for (int k = 0; k < n_cols; ++k)
    {
        bool is_xk_equal_to_lb = IS_EQUAL_FEAS_TOL(x[k], lbs[k]);
        bool is_xk_equal_to_ub = IS_EQUAL_FEAS_TOL(x[k], ubs[k]);

        // if wrong sign we set it to zero
        if (is_xk_equal_to_lb && !is_xk_equal_to_ub && z[k] < 0)
        {
            z[k] = 0.0;
        }
        else if (!is_xk_equal_to_lb && is_xk_equal_to_ub && z[k] > 0)
        {
            z[k] = 0.0;
        }
    }
}

void postsolver_run(const PostsolveInfo *info, Solution *sol, const double *x,
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

    copy_reduced_sol_to_orginal(sol, x, y, z, col_map, row_map);

    for (int i = n_reductions - 1; i >= 0; --i)
    {
        type = reductions[i];
        start = starts[i];

        if (type == FIXED_COL)
        {
            len = starts[i + 1] - start - 2;
            // we might have fixed an empty column so len can be 0
            assert(len >= 0);
            retrieve_fix_col(sol, indices[start], vals[start], vals[start + 1],
                             vals + start + 2, indices + start + 2, len);
        }
        else if (type == SUB_COL)
        {
            len = starts[i + 1] - start - 2;
            assert(len > 1);
            retrieve_sub_col(sol, indices[start], indices + start + 1,
                             vals + start + 1, len, vals[start],
                             indices[start + len + 1], vals[start + len + 1]);
        }
        else if (type == FIXED_COL_INF)
        {
            retrieve_fix_col_inf(sol, indices + start, vals + start);
        }
        else if (type == PARALLEL_COL)
        {
            assert(starts[i + 1] - start == 5);
            retrieve_parallel_col(sol, indices + start, vals + start);
        }
        else if (type == DELETED_ROW)
        {
            assert(starts[i + 1] - start == 1);
            retrieve_deleted_row(sol, indices[start], vals[start]);
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
                    retrieve_fix_col(sol, indices[start], vals[start],
                                     vals[start + 1], vals + start + 2,
                                     indices + start + 2, starts[j + 1] - start - 2);
                    assert(reductions[j - 1] == BOUND_CHANGE_NO_ROW);
                }
                else
                {
                    bound_changes_processed += 1;
                    retrieve_bound_change(sol, row, indices[start], row_cols,
                                          row_vals, len, vals[start],
                                          vals[start + 1], indices[start + 1]);
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
            retrieve_lhs_change(sol, indices[start], vals + start + 2,
                                indices + start + 2, len, vals[start],
                                indices[start + 1], vals[start + 1]);
        }
        else if (type == RHS_CHANGE)
        {
            len = starts[i + 1] - start - 2;
            assert(len > 1);
            retrieve_rhs_change(sol, indices[start], vals + start + 2,
                                indices + start + 2, len, vals[start],
                                indices[start + 1], vals[start + 1]);
        }
        else if (type == EQ_TO_INEQ)
        {
            assert(starts[i + 1] - start == 1);
            retrieve_eq_to_ineq(sol, indices[start], vals[start]);
        }
        else
        {
            assert(type != BOUND_CHANGE_NO_ROW);
            assert(false);
        }
    }

    // polish_z(sol, lbs, ubs);

#ifndef NDEBUG
    for (int i = 0; i < sol->dim_x; ++i)
    {
        if (sol->x[i] == COL_NOT_RETRIEVED || sol->z[i] == COL_NOT_RETRIEVED)
        {
            printf("col %d not fully retrieved \n", i);
        }

        assert(sol->x[i] != COL_NOT_RETRIEVED);
        assert(sol->z[i] != COL_NOT_RETRIEVED);
    }

    for (int i = 0; i < sol->dim_y; ++i)
    {
        assert(sol->y[i] != ROW_NOT_RETRIEVED);
    }
#endif
}

// here we should probably store information so zk = ck - ak^T y
// so must save rows and vals
void save_retrieval_fixed_col(PostsolveInfo *info, int col, double val, double ck,
                              const double *vals, const int *rows, int len)
{
    u16Vec_append(info->type, FIXED_COL);
    iVec_append(info->indices, col);
    iVec_append(info->indices, DUMMY_VALUE);
    iVec_append_array(info->indices, rows, (size_t) len);
    dVec_append(info->vals, val);
    dVec_append(info->vals, ck);
    dVec_append_array(info->vals, vals, (size_t) len);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

void save_retrieval_fixed_col_inf(PostsolveInfo *info, int col, int pos_inf,
                                  const Constraints *constraints, double bound)
{
    assert(pos_inf == 1 || pos_inf == -1);
    int i, row;
    const Matrix *AT = constraints->AT;
    const Matrix *A = constraints->A;
    int *row_sizes = constraints->state->row_sizes;
    double *lhs = constraints->lhs;
    double *rhs = constraints->rhs;
    RowTag *row_tags = constraints->row_tags;
    double side;

    int *rows = AT->i + AT->p[col].start;
    int n_rows = AT->p[col].end - constraints->AT->p[col].start;

    u16Vec_append(info->type, FIXED_COL_INF);
    iVec_append(info->indices, pos_inf);
    iVec_append(info->indices, col);
    dVec_append(info->vals, (double) n_rows);
    dVec_append(info->vals, bound);

    for (i = 0; i < n_rows; ++i)
    {
        row = rows[i];
        side = (HAS_TAG(row_tags[row], R_TAG_LHS_INF)) ? rhs[row] : lhs[row];
        dVec_append(info->vals, side);
        dVec_append_array(info->vals, A->x + A->p[row].start,
                          (size_t) row_sizes[row]);
        iVec_append(info->indices, row_sizes[row]);
        iVec_append_array(info->indices, A->i + A->p[row].start,
                          (size_t) row_sizes[row]);

        assert(!(HAS_TAG(row_tags[row], R_TAG_LHS_INF) &&
                 HAS_TAG(row_tags[row], R_TAG_RHS_INF)));
        assert(HAS_TAG(row_tags[row], (R_TAG_LHS_INF | R_TAG_RHS_INF)));
    }

    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

void save_retrieval_sub_col(PostsolveInfo *info, int col, int *cols, double *coeffs,
                            int len, double rhs, int i, double ck)
{
    u16Vec_append(info->type, SUB_COL);
    iVec_append(info->indices, col);
    iVec_append_array(info->indices, cols, (size_t) len);
    iVec_append(info->indices, i);
    dVec_append(info->vals, rhs);
    dVec_append_array(info->vals, coeffs, (size_t) len);
    dVec_append(info->vals, ck);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->indices->len == info->vals->len);
}

void save_retrieval_parallel_col(PostsolveInfo *info, double ub_j, double lb_j,
                                 double lb_k, double ub_k, double ratio, int j,
                                 int k, ColTag cTag_j, ColTag cTag_k)
{
    u16Vec_append(info->type, PARALLEL_COL);
    iVec_append(info->indices, j);
    iVec_append(info->indices, k);
    iVec_append(info->indices, (int) (cTag_j));
    iVec_append(info->indices, (int) (cTag_k));
    iVec_append(info->indices, DUMMY_VALUE);
    dVec_append(info->vals, lb_j);
    dVec_append(info->vals, ub_j);
    dVec_append(info->vals, lb_k);
    dVec_append(info->vals, ub_k);
    dVec_append(info->vals, ratio);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

void save_retrieval_deleted_row(PostsolveInfo *info, int row, double val)
{
    u16Vec_append(info->type, DELETED_ROW);
    iVec_append(info->indices, row);
    dVec_append(info->vals, val);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

void save_retrieval_added_row(PostsolveInfo *info, int i, int j, double ratio)
{
    u16Vec_append(info->type, ADDED_ROW);
    iVec_append(info->indices, i);
    iVec_append(info->indices, j);
    dVec_append(info->vals, ratio);
    dVec_append(info->vals, DUMMY_VALUE);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

void save_retrieval_added_rows(PostsolveInfo *info, int i, const int *rows,
                               const double *vals, int len, double aik)
{
    u16Vec_append(info->type, ADDED_ROWS);
    iVec_append(info->indices, i);
    iVec_append_array(info->indices, rows, (size_t) len);
    dVec_append(info->vals, aik);
    dVec_append_array(info->vals, vals, (size_t) len);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

// original bound can be infinite
void save_retrieval_bound_change_no_row(PostsolveInfo *info, int j,
                                        double implied_bound,
                                        double original_other_bound,
                                        int is_original_other_bound_lower_bound)
{
    u16Vec_append(info->type, BOUND_CHANGE_NO_ROW);
    iVec_append(info->indices, j);
    iVec_append(info->indices, is_original_other_bound_lower_bound);
    dVec_append(info->vals, implied_bound);
    dVec_append(info->vals, original_other_bound);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

void save_retrieval_bound_change_the_row(PostsolveInfo *info, int i, const int *cols,
                                         const double *vals, int len,
                                         int num_of_bound_changes)
{
    u16Vec_append(info->type, BOUND_CHANGE_THE_ROW);
    iVec_append(info->indices, i);
    iVec_append_array(info->indices, cols, (size_t) len);
    dVec_append(info->vals, (double) num_of_bound_changes);
    dVec_append_array(info->vals, vals, (size_t) len);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

void save_retrieval_rhs_or_lhs_change(PostsolveInfo *info, int i, const double *vals,
                                      const int *cols, int len, double new_side,
                                      int j, double ratio, bool is_lhs_change)
{
    if (is_lhs_change)
    {
        u16Vec_append(info->type, LHS_CHANGE);
    }
    else
    {
        u16Vec_append(info->type, RHS_CHANGE);
    }

    iVec_append(info->indices, i);
    iVec_append(info->indices, j);
    iVec_append_array(info->indices, cols, (size_t) len);
    dVec_append(info->vals, new_side);
    dVec_append(info->vals, ratio);
    dVec_append_array(info->vals, vals, (size_t) len);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}

void save_retrieval_eq_to_ineq(PostsolveInfo *info, int row, double val)
{
    u16Vec_append(info->type, EQ_TO_INEQ);
    iVec_append(info->indices, row);
    dVec_append(info->vals, val);
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->vals->len == info->indices->len);
}
