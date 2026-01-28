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

#include "DTonsEq.h"
#include "Activity.h"
#include "Bounds.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "Locks.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "Postsolver.h"
#include "Problem.h"
#include "RowColViews.h"
#include "SimpleReductions.h"
#include "State.h"
#include "Timer.h"
#include "Workspace.h"
#include <string.h>

/* Checks if there is sufficient space in AT for the substitution.
   This includes functionality for shifting rows to create sufficient
   space. Column j stays, column k is substituted */
#ifdef NDEBUG
static inline bool sufficient_space_AT(Matrix *AT, int j, int k, int max_shift)
#else
static inline bool sufficient_space_AT(Matrix *AT, int j, int k, int max_shift,
                                       const RowTag *row_tags)
#endif
{
    int len_col_j, len_col_k, fill_in, jj, kk, free_space;
    const int *rows_col_j, *rows_col_k;

    rows_col_j = AT->i + AT->p[j].start;
    len_col_j = AT->p[j].end - AT->p[j].start;
    rows_col_k = AT->i + AT->p[k].start;
    len_col_k = AT->p[k].end - AT->p[k].start;

    // check that there are no inactive rows in column j and k
    DEBUG(ASSERT_NO_INACTIVE_ROWS(rows_col_j, row_tags, len_col_j););
    DEBUG(ASSERT_NO_INACTIVE_ROWS(rows_col_k, row_tags, len_col_k););

    // -----------------------------------------------------------------
    // compute the fill-in in column j caused by substituting column k
    // -----------------------------------------------------------------
    fill_in = -1;
    jj = 0;
    kk = 0;
    while (jj < len_col_j && kk < len_col_k)
    {
        if (rows_col_j[jj] == rows_col_k[kk])
        {
            jj += 1;
            kk += 1;
        }
        else if (rows_col_k[kk] < rows_col_j[jj])
        {
            kk += 1;
            fill_in += 1;
        }
        else
        {
            jj += 1;
        }
    }
    fill_in += len_col_k - kk;

    //-------------------------------------------------------------------
    // Check if there is sufficient space in AT for the substitution.
    // If not, try to shift rows to create space.
    //-------------------------------------------------------------------
    free_space = AT->p[j + 1].start - AT->p[j].end;
    return (free_space >= fill_in || shift_row(AT, j, fill_in, max_shift));
}

/* This function checks if a doubleton row should be eliminated. It modifies
   the stay and subst variables. If no substitution is possible, stay and subst
   are set to -1. If the substitution is possible, stay and subst are set to
   the *column indices* of the variables that will stay and be substituted.
   The bounds of the variable that will be substituted are NOT modified.
*/
#ifdef NDEBUG
static inline void can_dton_be_eliminated(Matrix *AT, const double *row_vals,
                                          const int *row_cols, int *stay, int *subst,
                                          int col_size0, int col_size1,
                                          int max_shift)
#else
static inline void can_dton_be_eliminated(Matrix *AT, const double *row_vals,
                                          const int *row_cols, int *stay, int *subst,
                                          int col_size0, int col_size1,
                                          int max_shift, const RowTag *row_tags)
#endif
{
    double a_abs = ABS(row_vals[0]);
    double b_abs = ABS(row_vals[1]);
    bool is_integral0 = IS_INTEGRAL(a_abs / b_abs);
    bool is_integral1 = IS_INTEGRAL(b_abs / a_abs);

    // ------------------------------------------------------------------------
    // Decide which variable to substitute. For a row ax + by = c we can do
    // substitutions x = -(b/a)y + c/a or y = -(a/b)x + c/b. If b/a is an
    // integer but not a/b we substitute x. If a/b is an integer but not b/a
    // we substitute y.
    //
    // We might want to add additional criteria here? Do we first want to
    // check if one of the columns has far more nnz than the other one?
    // ------------------------------------------------------------------------
    if (col_size0 == 1 && col_size1 != 1)
    {
        *subst = 0;
    }
    else if (col_size0 != 1 && col_size1 == 1)
    {
        *subst = 1;
    }
    else if (is_integral0 && !is_integral1)
    {
        *subst = 1;
    }
    else if (is_integral1 && !is_integral0)
    {
        *subst = 0;
    }
    else
    {
        if (col_size0 < col_size1)
        {
            *subst = 0;
        }
        else
        {
            *subst = 1;
        }
    }

    *stay = 1 - *subst;

    // ------------------------------------------------------------------------
    // Simple test to reject large or small pivots. We might want to add a more
    // sophisticated check here. We could potentially move this in the code.
    // ------------------------------------------------------------------------
    double abs_pivot = row_vals[*stay] / row_vals[*subst];
    abs_pivot = ABS(abs_pivot);
    if (abs_pivot > MAX_RATIO_PIVOT || abs_pivot < 1 / MAX_RATIO_PIVOT)
    {
        assert("Is this ever triggered now with this small threshold?");
        *stay = -1;
        *subst = -1;
        return;
    }

    // ------------------------------------------------------------------------
    // Check if there is sufficient space in AT for the substitution.
    // ------------------------------------------------------------------------
#ifdef NDEBUG
    if (!sufficient_space_AT(AT, row_cols[*stay], row_cols[*subst], max_shift))
#else
    if (!sufficient_space_AT(AT, row_cols[*stay], row_cols[*subst], max_shift,
                             row_tags))
#endif
    {
        *stay = -1;
        *subst = -1;
        return;
    }

    *stay = row_cols[*stay];
    *subst = row_cols[*subst];
}

// lb and ub are bounds on variable that gets substituted
static inline void modify_bounds(Constraints *constraints, int i, double aij,
                                 double aik, double rhs, double lb_subst,
                                 double ub_subst, int stay, ColTag col_tag_subst)
{
    assert(aij != 0 && aik != 0);
    bool same_sign = (aik * aij > 0.0);

    if (same_sign)
    {
        if (!HAS_TAG(col_tag_subst, C_TAG_LB_INF))
        {
            double new_ub_cand = (rhs - aik * lb_subst) / aij;
            update_ub(constraints, stay, new_ub_cand, i HUGE_BOUND_IS_NOT_OK);

            if (HAS_TAG(constraints->col_tags[stay], C_TAG_FIXED))
            {
                return;
            }
        }

        if (!HAS_TAG(col_tag_subst, C_TAG_UB_INF))
        {
            double new_lb_cand = (rhs - aik * ub_subst) / aij;
            update_lb(constraints, stay, new_lb_cand, i HUGE_BOUND_IS_NOT_OK);
        }
    }
    else
    {
        if (!HAS_TAG(col_tag_subst, C_TAG_LB_INF))
        {
            double new_lb_cand = (rhs - aik * lb_subst) / aij;
            update_lb(constraints, stay, new_lb_cand, i HUGE_BOUND_IS_NOT_OK);

            if (HAS_TAG(constraints->col_tags[stay], C_TAG_FIXED))
            {
                return;
            }
        }

        if (!HAS_TAG(col_tag_subst, C_TAG_UB_INF))
        {
            double new_ub_cand = (rhs - aik * ub_subst) / aij;
            update_ub(constraints, stay, new_ub_cand, i HUGE_BOUND_IS_NOT_OK);
        }
    }
}

#ifndef TESTING
static inline
#endif
    Old_and_new_coeff update_row_A_dton(Matrix *A, int i, int q, int j, int k,
                                        double aij, double aik, int *row_size,
                                        PostsolveInfo *postsolve_info)
{
    int ii, start, end, insertion;
    double old_val = 0.0;
    int subst_idx = -1;
    int diff_row_size = 0;
    start = A->p[q].start;
    end = A->p[q].end;
    insertion = end;

    // -----------------------------------------------------------------
    // find index of variable that is substituted and and the index of
    // insertion of the variable that stays
    // -----------------------------------------------------------------
    for (ii = start; ii < end; ++ii)
    {
        if (A->i[ii] == k)
        {
            subst_idx = ii;
            break;
        }
    }

    assert(subst_idx != -1);

    for (ii = start; ii < end; ++ii)
    {
        if (A->i[ii] >= j)
        {
            insertion = ii;
            break;
        }
    }

    // -----------------------------------------------------------------
    // Compute new coefficient of the variable that stays.
    // TODO: should we have rejected the reduction if the new value is
    // very small? (don't delete todo here)
    // -----------------------------------------------------------------
    if (insertion != end && A->i[insertion] == j)
    {
        old_val = A->x[insertion];
    }

    double new_val = old_val - (aij / aik) * A->x[subst_idx];
    Old_and_new_coeff old_and_new_coeff;
    old_and_new_coeff.old_coeff_stay = old_val;
    old_and_new_coeff.old_coeff_subst = A->x[subst_idx];
    old_and_new_coeff.new_coeff_stay = new_val;

    save_retrieval_added_row(postsolve_info, i, q, -A->x[subst_idx] / aik);

#ifndef NDEBUG
    if (ABS(A->x[subst_idx]) > 1e-6)
    {
        assert(ABS(new_val) > 1e-6 || new_val == 0.0);
    }

    if (ABS(new_val) < 1e-6 && new_val != 0.0)
    {
        // when we eliminate row 32833 on momentum and substitute col x20
        // into row 30698 this is triggered, because x20 has a very small
        // coefficient in row 30698 BEFORE the substitution
        printf("warning debug: small coefficient dtonrow \n");
    }
#endif

    // -------------------------------------------------------------------
    //          remove variable that is substituted
    // -------------------------------------------------------------------
    size_t len = (size_t) (end - subst_idx - 1);
    memmove(A->x + subst_idx, A->x + subst_idx + 1, len * sizeof(double));
    memmove(A->i + subst_idx, A->i + subst_idx + 1, len * sizeof(int));
    end -= 1;
    diff_row_size += 1;
    insertion -= (subst_idx < insertion);

    // -------------------------------------------------------------------
    // if we get unexpected cancellation we should remove the coefficient,
    // otherwise we insert it
    // -------------------------------------------------------------------
    if (ABS(new_val) <= ZERO_TOL)
    {
        size_t len = (size_t) (end - insertion - 1);
        memmove(A->x + insertion, A->x + insertion + 1, len * sizeof(double));
        memmove(A->i + insertion, A->i + insertion + 1, len * sizeof(int));
        end -= 1;
        diff_row_size += 1;
    }
    // -----------------------------------------------------------------
    // Insert the new value. If it exists or should be inserted in
    // the end, we don't need to shift values.
    // -----------------------------------------------------------------
    else
    {
        if (insertion == end)
        {
            A->x[insertion] = new_val;
            A->i[insertion] = j;
            end += 1;
            diff_row_size -= 1;
        }
        else if (A->i[insertion] == j)
        {
            A->x[insertion] = new_val;
        }
        else
        {
            size_t len = (size_t) (end - insertion);
            memmove(A->x + insertion + 1, A->x + insertion, len * sizeof(double));
            memmove(A->i + insertion + 1, A->i + insertion, len * sizeof(int));
            A->x[insertion] = new_val;
            A->i[insertion] = j;
            end += 1;
            diff_row_size -= 1;
        }
    }

    A->p[q].end = end;
    *row_size -= diff_row_size;
    A->nnz -= (size_t) diff_row_size;

    assert(*row_size == end - start);
    DEBUG(ASSERT_INCREASING_I(A->i + start, (size_t) *row_size););
    DEBUG(ASSERT_NO_ZEROS_D(A->x + start, (size_t) *row_size););
    assert(A->p[q].end <= A->p[q + 1].start);
    return old_and_new_coeff;
}

static inline void subst_col_in_obj_dton(Objective *obj, int stay, int subst,
                                         double aik, double aij, double rhs)
{
    obj->c[stay] -= (aij / aik) * obj->c[subst];
    obj->offset += (rhs / aik) * obj->c[subst];
}

// lhs and rhs should point to the row whose lhs and rhs should be modified
// rTag = row_tags[rows_subst[j]], lhs = lhs + rows_subst[j],
// rhs = rhs + rows_subst[j], rhs_dton
static inline void update_lhs_and_rhs(double *lhs, double *rhs, double aik,
                                      RowTag rTag, double rhs_dton,
                                      double coeff_subst)
{
    if (rhs_dton != 0.0)
    {
        double change = coeff_subst / aik * rhs_dton;

        if (!HAS_TAG(rTag, R_TAG_LHS_INF))
        {
            *lhs -= change;
        }

        if (!HAS_TAG(rTag, R_TAG_RHS_INF))
        {
            *rhs -= change;
        }
    }
}

/* Decides which variable that should be substituted. If the substitution
   is rejected for whatever reason, it sets j = -1. */
#ifdef NDEBUG
static void find_substitution(const RowView *row, const ColTag *col_tags, int *j,
                              int *k, double *aik, double *aij, Matrix *AT,
                              int *col_sizes, int max_shift_per_row,
                              int *iwork_n_rows, int *n_new_dton_rows)
#else
static void find_substitution(const RowView *row, const ColTag *col_tags, int *j,
                              int *k, double *aik, double *aij, Matrix *AT,
                              int *col_sizes, int max_shift_per_row,
                              const RowTag *row_tags, int *iwork_n_rows,
                              int *n_new_dton_rows)
#endif
{

    // a previous dtonrow might have changed this dtonrows sparsity pattern,
    // or a row that used to be a dtonrow was pushed to the list of dtonrows
    // but since then one or both of its variables have been fixed because
    // of singleton rows
    if (*row->len < 2)
    {
        *j = -1;
        return;
    }

    assert(HAS_TAG(*row->tag, R_TAG_EQ) && !HAS_TAG(*row->tag, R_TAG_INACTIVE));
    assert(*row->lhs == *row->rhs && !IS_ABS_INF(*row->rhs));

    // One of the variables might have been fixed when eliminating
    // another doubleton row. However, the variable has not been
    // removed from the problem yet so the row size check above
    // is not sufficient, and we also need this check
    if (HAS_TAG((col_tags[row->cols[0]] | col_tags[row->cols[1]]), C_TAG_INACTIVE))
    {
        *j = -1;
        return;
    }

    // -------------------------------------------------------------------
    // check if there is sufficient space to eliminate the doubleton row,
    // if stay == -1 it means that the elimination was rejected
    // -------------------------------------------------------------------
#ifdef NDEBUG
    can_dton_be_eliminated(AT, row->vals, row->cols, j, k, col_sizes[row->cols[0]],
                           col_sizes[row->cols[1]], max_shift_per_row);
#else
    can_dton_be_eliminated(AT, row->vals, row->cols, j, k, col_sizes[row->cols[0]],
                           col_sizes[row->cols[1]], max_shift_per_row, row_tags);
#endif

    if (*j == -1)
    {
        // append the row so we try it again in next round
        iwork_n_rows[*n_new_dton_rows] = row->i;
        *n_new_dton_rows += 1;
        return;
    }

    assert(*j == row->cols[0] || *j == row->cols[1]);
    assert(*k == row->cols[0] || *k == row->cols[1]);
    assert(!HAS_TAG(col_tags[*j], C_TAG_INACTIVE));
    assert(!HAS_TAG(col_tags[*k], C_TAG_INACTIVE));

    // Our current implementation of singleton columns does not take
    // into account that a column singleton in a dtonrow can always be
    // be eliminated. Hence, we might eliminate dtonrows here that
    // contains a singleton column.
    assert(col_sizes[*j] != 1 || (col_sizes[*j] == 1 && col_sizes[*k] == 1));

    // -------------------------------------------------------------------
    // Find coefficient that stays. i is the row index, column j
    // stays, and column k gets substituted
    // -------------------------------------------------------------------
    if (*k == row->cols[0])
    {
        assert(*j == row->cols[1]);
        *aik = row->vals[0];
        *aij = row->vals[1];
    }
    else
    {
        assert(*j == row->cols[0] && *k == row->cols[1]);
        *aik = row->vals[1];
        *aij = row->vals[0];
    }
}

static void execute_substitution(Constraints *constraints, RowView *row, int j,
                                 int k, double aij, double aik, int *n_new_dton_rows,
                                 double ck)
{
    Matrix *A = constraints->A;
    Matrix *AT = constraints->AT;
    double *lhs = constraints->lhs;
    double *rhs = constraints->rhs;
    int *col_sizes = constraints->state->col_sizes;
    int *row_sizes = constraints->state->row_sizes;
    Activity *activities = constraints->state->activities;
    PostsolveInfo *postsolve_info = constraints->state->postsolve_info;
    Bound *bounds = constraints->bounds;
    RowTag *row_tags = constraints->row_tags;
    ColTag *col_tags = constraints->col_tags;

    iVec *empty_rows = constraints->state->empty_rows;
    iVec *ston_rows = constraints->state->ston_rows;
    int *iwork_n_rows = constraints->state->work->iwork_n_rows;

    Altered_Activity altered1, altered2;
    int jj, q, len, old_len;
    int *rows_subst;

    rows_subst = AT->i + AT->p[k].start;
    len = AT->p[k].end - AT->p[k].start;
    RowView rowj_AT = new_rowview(AT->x + AT->p[j].start, AT->i + AT->p[j].start,
                                  col_sizes + j, AT->p + j, NULL, NULL, NULL, j);

    // remove contribution of xk in AT (we do it before the loop to free space)
    remove_coeff(&rowj_AT, row->i);
    AT->p[k].end = AT->p[k].start;

    // we substitute variable xk from all rows q
    for (jj = 0; jj < len; ++jj)
    {
        q = rows_subst[jj];
        if (HAS_TAG(row_tags[q], R_TAG_INACTIVE))
        {
            continue;
        }

        RowView sub_row =
            new_rowview(A->x + A->p[q].start, A->i + A->p[q].start, row_sizes + q,
                        A->p + q, lhs + q, rhs + q, row_tags + q, q);

        assert(sub_row.i != row->i);

        // update constraint matrix A (updates row size)
        old_len = *sub_row.len;
        Old_and_new_coeff coeffs = update_row_A_dton(
            A, row->i, sub_row.i, j, k, aij, aik, sub_row.len, postsolve_info);

        // Check row size and push to list of empty rows, stonrows and
        // dtonrows. The elimination of one dtonrow might lead to another
        // dtonrow, so we must be careful with how we append. We use a
        // temporary buffer to handle this.
        switch (*sub_row.len)
        {
            case 0:
                assert(!iVec_contains(empty_rows, sub_row.i));
                iVec_append(empty_rows, sub_row.i);
                assert(!HAS_TAG(*sub_row.tag, R_TAG_INACTIVE));
                assert(sub_row.range->start == sub_row.range->end);
                break;
            case 1:
                if (old_len != 1)
                {
                    assert(!iVec_contains(ston_rows, sub_row.i));
                    iVec_append(ston_rows, sub_row.i);
                    assert(!HAS_TAG(*sub_row.tag, R_TAG_INACTIVE));
                }
                break;
            case 2:
                if (HAS_TAG(*sub_row.tag, R_TAG_EQ) && old_len != 2)
                {
                    iwork_n_rows[*n_new_dton_rows] = sub_row.i;
                    *n_new_dton_rows += 1;
                    assert(!HAS_TAG(*sub_row.tag, R_TAG_INACTIVE));
                }
                break;
        }

        // update constraint matrix AT (updates column size)
        insert_or_update_coeff(AT, j, sub_row.i, coeffs.new_coeff_stay,
                               col_sizes + j);

        // update lhs and rhs
        update_lhs_and_rhs(sub_row.lhs, sub_row.rhs, aik, *sub_row.tag, rhs[row->i],
                           coeffs.old_coeff_subst);

        // update activity due to the coefficient change of variable that stays
        altered1 = update_activity_coeff_change(activities + q, bounds[j].lb,
                                                bounds[j].ub, coeffs.old_coeff_stay,
                                                coeffs.new_coeff_stay, col_tags[j]);

        // update activity due to the coefficient change of variable that is
        // substituted
        altered2 =
            update_activity_coeff_change(activities + q, bounds[k].lb, bounds[k].ub,
                                         coeffs.old_coeff_subst, 0.0, col_tags[k]);

        if ((altered1 | altered2) & MIN_ALTERED_RECOMPUTE)
        {
            assert(*sub_row.len == A->p[q].end - A->p[q].start);
            activities[sub_row.i].min = compute_min_act_no_tags(
                sub_row.vals, sub_row.cols, *sub_row.len, bounds);

            if (activities[sub_row.i].status == NOT_ADDED)
            {
                activities[sub_row.i].status = ADDED;
                iVec_append(constraints->state->updated_activities, sub_row.i);
            }
        }

        if ((altered1 | altered2) & MAX_ALTERED_RECOMPUTE)
        {
            // if this assertion fails, replace *sub_row.len with the rhs
            assert(*sub_row.len == A->p[q].end - A->p[q].start);
            activities[sub_row.i].max = compute_max_act_no_tags(
                sub_row.vals, sub_row.cols, *sub_row.len, bounds);

            if (activities[sub_row.i].status == NOT_ADDED)
            {
                activities[sub_row.i].status = ADDED;
                iVec_append(constraints->state->updated_activities, sub_row.i);
            }
        }
    }

    // update list of empty and singleton columns
    switch (col_sizes[j])
    {
        case 0:
            assert(!HAS_TAG(col_tags[j], C_TAG_INACTIVE));
            assert(!iVec_contains(constraints->state->empty_cols, j));
            iVec_append(constraints->state->empty_cols, j);
            assert(AT->p[j].start == AT->p[j].end);
            break;
        case 1:
            assert(!HAS_TAG(col_tags[j], C_TAG_INACTIVE));
            assert(!iVec_contains(constraints->state->ston_cols, j));
            iVec_append(constraints->state->ston_cols, j);
            break;
        default:
            break;
    }

    // mark row of AT as inactive (must do it after loop due to
    // assertion in update activities)
    col_tags[k] = C_TAG_INACTIVE;
    col_sizes[k] = SIZE_INACTIVE_COL;

    assert(AT->p[j].end - AT->p[j].start == col_sizes[j]);
    count_locks_one_column(&rowj_AT, constraints->state->col_locks + j, row_tags);

    // postsolve information
    save_retrieval_sub_col(postsolve_info, k, row->cols, row->vals, 2, rhs[row->i],
                           row->i, ck);
    save_retrieval_deleted_row(postsolve_info, row->i, ck / aik);
}

/* We assume the row i is a doubleton equality row with variables
xj and xk, where xk is substituted and xj stays. We eliminate
xk from every row q that contains xk. */
static PresolveStatus remove_dton_eq_rows__(Problem *prob, int max_shift_per_row)
{
    int ii, i, j, k;
    int n_new_dton_rows = 0;
    double aik, aij;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    Matrix *AT = constraints->AT;
    double *rhs = constraints->rhs;
    double *lhs = constraints->lhs;
    RowRange *row_r = A->p;
    RowTag *row_tags = constraints->row_tags;
    ColTag *col_tags = constraints->col_tags;
    const iVec *dton_rows = constraints->state->dton_rows;
    int *col_sizes = constraints->state->col_sizes;
    int *row_sizes = constraints->state->row_sizes;

    const Bound *bounds = constraints->bounds;
    PresolveStatus status = UNCHANGED;
    int *iwork_n_rows = constraints->state->work->iwork_n_rows;

    DEBUG(verify_no_duplicates_sort(dton_rows));

    for (ii = 0; ii < dton_rows->len; ++ii)
    {
        i = dton_rows->data[ii];
        RowView row =
            new_rowview(A->x + row_r[i].start, A->i + row_r[i].start, row_sizes + i,
                        A->p + i, lhs + i, rhs + i, row_tags + i, i);

#ifdef NDEBUG
        //  find which variable that should be substituted (j stays, k goes)
        find_substitution(&row, col_tags, &j, &k, &aik, &aij, AT, col_sizes,
                          max_shift_per_row, iwork_n_rows, &n_new_dton_rows);
#else
        find_substitution(&row, col_tags, &j, &k, &aik, &aij, AT, col_sizes,
                          max_shift_per_row, row_tags, iwork_n_rows,
                          &n_new_dton_rows);
#endif

        // no variable could be substituted
        if (j == -1)
        {
            // todo: here we could run simple primal sparsification?
            continue;
        }

        status = REDUCED;

        // transfer bounds to variable that stays and skip the reduction
        // if the variable that stays was fixed because of the bound change
        modify_bounds(constraints, i, aij, aik, *row.rhs, bounds[k].lb, bounds[k].ub,
                      j, col_tags[k]);

        if (HAS_TAG(col_tags[j], C_TAG_INACTIVE))
        {
            continue;
        }

        //  substitute variable in objective and mark the eliminated row as
        //  inactive
        subst_col_in_obj_dton(prob->obj, j, k, aik, aij, *row.rhs);
        RESET_TAG(*row.tag, R_TAG_INACTIVE);
        *row.len = SIZE_INACTIVE_ROW;
        row.range->end = row.range->start;
        A->nnz -= 2;

        // Substitute the variable in the constraints, using special
        // functionality
        execute_substitution(constraints, &row, j, k, aij, aik, &n_new_dton_rows,
                             prob->obj->c[k]);
    }

    AT->nnz = A->nnz;

    // clear the list of processed dton rows
    iVec_clear_no_resize(constraints->state->dton_rows);

    // append the new dton rows
    if (n_new_dton_rows > 0)
    {
        DEBUG(verify_no_duplicates_sort_ptr(iwork_n_rows, (size_t) n_new_dton_rows));
        iVec_append_array(constraints->state->dton_rows, iwork_n_rows,
                          (size_t) n_new_dton_rows);
    }

    return status;
}

PresolveStatus remove_dton_eq_rows(Problem *prob, int max_shift_per_row)
{
    assert(prob->constraints->state->ston_rows->len == 0);
    assert(prob->constraints->state->empty_rows->len == 0);
    assert(prob->constraints->state->empty_cols->len == 0);
    DEBUG(verify_problem_up_to_date(prob->constraints));

    // important to have double ptr in case of realloc when appending
    iVec **dton_rows = &(prob->constraints->state->dton_rows);
    while ((*dton_rows)->len > 0)
    {
        PresolveStatus temp = remove_dton_eq_rows__(prob, max_shift_per_row);

        if (temp != REDUCED)
        {
            break;
        }
    }

    // some variables might have been fixed because of the bound change
    if (prob->constraints->state->fixed_cols_to_delete->len > 0)
    {
        delete_fixed_cols_from_problem(prob);
        delete_inactive_cols_from_A_and_AT(prob->constraints);
    }

    DEBUG(verify_problem_up_to_date(prob->constraints));

    // always return UNCHANGED because this function can't detect infeas/unbnd.
    return UNCHANGED;
}
