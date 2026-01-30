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

#include "Parallel_rows.h"
#include "Bounds.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "Matrix.h"
#include "Numerics.h"
#include "RowColViews.h"
#include "State.h"
#include "Workspace.h"
#include "iVec.h"
#include "limits.h" // for INT_MAX
#include "utils.h"
#include <PSLP_warnings.h>
#include <math.h> // For round()
#include <stdint.h>

// global variables needed for qsort
static int *global_sparsity_IDs;
static int *global_coeff_hashes;

// djb2 hash function
static inline uint32_t hash_int_array(const int *arr, int size)
{
    /* disable sign-conversion compiler warning*/
    PSLP_DIAG_PUSH();
    PSLP_DIAG_IGNORE_SIGN_CONVERSION();

    uint32_t hash = 5381;
    for (int i = 0; i < size; i++)
    {
        hash = ((hash << 5) + hash) + arr[i];
    }

    /* enable sign-conversion compiler warnings */
    PSLP_DIAG_POP();
    return hash;
}

#define INV_PRECISION 1e6

// djb2 hash function, with normalization of values and rounding of doubles
static inline uint32_t hash_double_array_with_scale(const double *arr, int size)
{
    // find max value of row
    double max = ABS(arr[0]);
    for (int i = 1; i < size; i++)
    {
        if (ABS(arr[i]) > max)
        {
            max = ABS(arr[i]);
        }
    }
    double scale = (arr[0] > 0) ? 1 / max : -1 / max;

    // compute has value
    uint32_t hash = 5381;
    for (int i = 0; i < size; i++)
    {
        uint32_t scaled = (uint32_t) (round((arr[i] * scale) * INV_PRECISION));
        hash = ((hash << 5) + hash) + scaled;
    }

    return hash;
}

#ifndef TESTING
static inline
#endif
    void compute_supp_and_coeff_hash(const Matrix *A, const RowTag *rowtags,
                                     int *sparsity_IDs, int *coeff_hashes,
                                     RowTag INACTIVE_TAG)
{

    for (int i = 0; i < A->m; i++)
    {
        if (HAS_TAG(rowtags[i], INACTIVE_TAG))
        {
            // important to set support_ID for inactive rows to place them
            // last in the sort
            sparsity_IDs[i] = INT_MAX;
            continue;
        }

        // this check makes sense since we also use this function to
        // find parallel columns
        assert(!HAS_TAG(rowtags[i], C_TAG_INACTIVE));

        int len = A->p[i].end - A->p[i].start;
        sparsity_IDs[i] = (int) hash_int_array(A->i + A->p[i].start, len);
        coeff_hashes[i] =
            (int) hash_double_array_with_scale(A->x + A->p[i].start, len);
    }
}

int comparator(const void *a, const void *b)
{
    int idxA = *(const int *) a;
    int idxB = *(const int *) b;

    if (global_sparsity_IDs[idxA] < global_sparsity_IDs[idxB])
    {
        return -1;
    }

    if (global_sparsity_IDs[idxA] > global_sparsity_IDs[idxB])
    {
        return 1;
    }

    if (global_coeff_hashes[idxA] < global_coeff_hashes[idxB])
    {
        return -1;
    }

    if (global_coeff_hashes[idxA] > global_coeff_hashes[idxB])
    {
        return 1;
    }

    return 0;
}

// Sort function
static inline void sort_rows(int *rows, PSLP_uint n_rows, int *sparsity_IDs,
                             int *coeff_hashes)
{
    global_sparsity_IDs = sparsity_IDs;
    global_coeff_hashes = coeff_hashes;
    qsort(rows, n_rows, sizeof(int), comparator);
}

static inline int get_bin_size(int start, PSLP_uint n_rows, const int *rows,
                               const int *sparsity_IDs, const int *coeff_hashes)
{
    int sparsity_ID = sparsity_IDs[rows[start]];
    int coeff_hash = coeff_hashes[rows[start]];
    int i;
    for (i = start + 1; i < n_rows; ++i)
    {
        if (sparsity_IDs[rows[i]] != sparsity_ID ||
            coeff_hashes[rows[i]] != coeff_hash)
        {
            break;
        }
    }
    return i - start;
}

#ifndef NDEBUG
// make sure that the rows in the same bin have the same coefficient
// hash value and the same sparsity ID
void ASSERT_BIN_CORRECT(const int *coeff_hashes, const int *sparsity_IDs,
                        const int *rows, int bin_start, int bin_size,
                        const RowTag *row_tags, RowTag INACTIVE_TAG)
{
    const int *bin_temp = rows + bin_start;
    int coeff_hash_start = coeff_hashes[bin_temp[0]];
    int sparsity_ID_start = sparsity_IDs[bin_temp[0]];

    for (int k = 1; k < bin_size; ++k)
    {
        assert(coeff_hashes[bin_temp[k]] == coeff_hash_start);
        assert(sparsity_IDs[bin_temp[k]] == sparsity_ID_start);
        assert(!HAS_TAG(row_tags[bin_temp[k]], INACTIVE_TAG));
    }
}

void VERIFY_PARALLEL_ROWS(const Matrix *A, const RowTag *rows_tags,
                          const int *parallel_rows, const iVec *group_starts)
{
    int j, k, start, end, n_rows_this_group, len1;
    PSLP_uint i;
    assert(group_starts->len > 0);
    // if (group_starts->len == 0)
    //{
    //     return;
    // }

    PSLP_uint n_groups = group_starts->len - 1;

    for (i = 0; i < n_groups; ++i)
    {
        start = group_starts->data[i];
        end = group_starts->data[i + 1];
        n_rows_this_group = end - start;

        // check that the group has at least two rows
        assert(n_rows_this_group > 1);

        // check that all rows in the group are active
        const int *temp = parallel_rows + start;
        ASSERT_NO_INACTIVE_ROWS(temp, rows_tags, n_rows_this_group);

        // check that all rows have the same support
        len1 = A->p[parallel_rows[start]].end - A->p[parallel_rows[start]].start;
        int *cols1 = A->i + A->p[parallel_rows[start]].start;
        for (j = start + 1; j < end; ++j)
        {
            assert(len1 ==
                   A->p[parallel_rows[j]].end - A->p[parallel_rows[j]].start);
            ARRAYS_EQUAL_INT(cols1, A->i + A->p[parallel_rows[j]].start, len1);
        }

        // check that all rows have same normalized coefficients
        double *vals1 = A->x + A->p[parallel_rows[start]].start;
        for (j = start + 1; j < end; ++j)
        {
            double *vals2 = A->x + A->p[parallel_rows[j]].start;
            double ratio = vals1[0] / vals2[0];
            for (k = 1; k < len1; ++k)
            {
                assert(IS_ZERO_FEAS_TOL(vals1[k] - ratio * vals2[k]));
            }
        }
    }
}
#endif

/*  In one bin there can theoretically be several groups of parallel rows,
    but hopefully it is unlikely that there are more than one group of
    parallel rows. We therefore only try to catch one of the groups.

    Due to our choice of hash function, we must verify that the rows
    have both the same support and the same coefficients (up to a multiple).

    * Returns the number of parallel rows found in the bin.
    * bin: pointer pointing to the FIRST row in the bin
    * bin_size: number of rows in the bin
    * parallel_rows: parallel rows are stored here, starting from
   parallel_rows[0]
*/
static inline int find_parallel_rows_in_bin(const Matrix *A, const int *bin,
                                            int bin_size, int *parallel_rows,
                                            iVec *group_starts)
{
    assert(bin_size > 1);
    int i, j, first_row_in_bin = bin[0];
    int n_new_parallel_rows = 0;

    const int *cols1 = A->i + A->p[bin[0]].start;
    const double *vals1 = A->x + A->p[bin[0]].start;
    int len1 = A->p[bin[0]].end - A->p[bin[0]].start;

    for (i = 1; i < bin_size; ++i)
    {
        const int *cols2 = A->i + A->p[bin[i]].start;
        const double *vals2 = A->x + A->p[bin[i]].start;
        int len2 = A->p[bin[i]].end - A->p[bin[i]].start;

        // check that the rows have same size
        if (len1 != len2)
        {
            continue;
        }

        // check that the rows have same support and coefficients up to multiple
        // (must start at j = 0 because of support check)
        double ratio = vals1[0] / vals2[0];
        for (j = 0; j < len2; ++j)
        {
            // if we run into numerical issues we might want to do
            // diff = vals1[j] / vals2[j] - ratio
            double diff = vals1[j] - ratio * vals2[j];

            if (!IS_ZERO_FEAS_TOL(diff) || cols1[j] != cols2[j])
            {
                break;
            }
        }

        // if j == len2, we have found a parallel row
        if (j == len2)
        {
            parallel_rows[n_new_parallel_rows] = bin[i];
            n_new_parallel_rows++;
        }
    }

    // add first row in bin
    if (n_new_parallel_rows > 0)
    {
        parallel_rows[n_new_parallel_rows] = first_row_in_bin;
        n_new_parallel_rows++;
        iVec_append(group_starts,
                    group_starts->data[group_starts->len - 1] + n_new_parallel_rows);
    }

    return n_new_parallel_rows;
}

// replaced rows with parallel_rows for less memory usage
void find_parallel_rows(const Matrix *A, const RowTag *r_Tags, iVec *group_starts,
                        int *parallel_rows, int *sparsity_IDs, int *coeff_hashes,
                        RowTag INACTIVE_TAG)
{
    int i, bin_size;
    int n_p_rows_total = 0;
    // ------------------------------------------------------------------------
    //             compute sparsity_IDs and coeff_hashes
    // ------------------------------------------------------------------------
    compute_supp_and_coeff_hash(A, r_Tags, sparsity_IDs, coeff_hashes, INACTIVE_TAG);

    // ------------------------------------------------------------------------
    // Make rows in the same bin appear next to each other. The sort makes sure
    // that rows with the same sparsity pattern and the same coefficient hash
    // value end up next to each other.
    // ------------------------------------------------------------------------
    for (i = 0; i < A->m; i++)
    {
        parallel_rows[i] = i;
    }

    sort_rows(parallel_rows, A->m, sparsity_IDs, coeff_hashes);

#ifndef NDEBUG
    if (INACTIVE_TAG == R_TAG_INACTIVE)
    {
        ASSERT_NO_ACTIVE_STON_ROWS(A, r_Tags);
    }
#endif

    // ------------------------------------------------------------------------
    // Start looping through the bins. Inactive rows have been placed in the
    // end of 'rows'. As soon as we meet an inactive row we can stop, since
    // all rows appearing after it are also inactive
    // ------------------------------------------------------------------------
    iVec_clear_no_resize(group_starts);
    iVec_append(group_starts, 0);
    for (i = 0; i < A->m; i++)
    {
        if (HAS_TAG(r_Tags[parallel_rows[i]], INACTIVE_TAG))
        {
            break;
        }

        bin_size = get_bin_size(i, A->m, parallel_rows, sparsity_IDs, coeff_hashes);

        // find the parallel rows in the bin
        if (bin_size > 1)
        {
            DEBUG(ASSERT_BIN_CORRECT(coeff_hashes, sparsity_IDs, parallel_rows, i,
                                     bin_size, r_Tags, INACTIVE_TAG););

            n_p_rows_total += find_parallel_rows_in_bin(
                A, parallel_rows + i, bin_size, parallel_rows + n_p_rows_total,
                group_starts);

            i += bin_size - 1;
        }
    }
    assert(n_p_rows_total == group_starts->data[group_starts->len - 1]);
    DEBUG(VERIFY_PARALLEL_ROWS(A, r_Tags, parallel_rows, group_starts););
}

static inline PresolveStatus process_single_bin(const Constraints *constraints,
                                                const int *bin, int bin_size)
{
    assert(bin_size > 1);
    const RowTag *row_tags = constraints->row_tags;
    const Matrix *A = constraints->A;
    const double *rhs = constraints->rhs;
    const double *lhs = constraints->lhs;
    double remaining_row_new_rhs, remaining_row_new_lhs, remaining_row_coeff;
    double other_row_coeff, ratio, other_row_rhs_scaled, other_row_lhs_scaled;
    double other_row_lhs_scale_factor = 1.0;
    double other_row_rhs_scale_factor = 1.0;
    int other_row_lhs_index = -1;
    int other_row_rhs_index = -1;
    bool is_rhs_inf_other_row, is_lhs_inf_other_row;

    // ------------------------------------------------------------------------
    //  Initialize quantities for the remaining row. These quantities will be
    //  updated as we process the bin.
    // ------------------------------------------------------------------------
    bool is_remaining_row_eq, is_rhs_inf_remaining_row, is_lhs_inf_remaining_row;
    int remaining_row_idx = bin[0];
    is_remaining_row_eq = HAS_TAG(row_tags[remaining_row_idx], R_TAG_EQ);
    is_rhs_inf_remaining_row = HAS_TAG(row_tags[remaining_row_idx], R_TAG_RHS_INF);
    is_lhs_inf_remaining_row = HAS_TAG(row_tags[remaining_row_idx], R_TAG_LHS_INF);
    remaining_row_new_rhs = rhs[remaining_row_idx];
    remaining_row_new_lhs = lhs[remaining_row_idx];
    remaining_row_coeff = A->x[A->p[remaining_row_idx].start];

    for (int i = 1; i < bin_size; ++i)
    {
        int other_row_idx = bin[i];
        bool is_other_row_eq = HAS_TAG(row_tags[other_row_idx], R_TAG_EQ);

        if (!is_remaining_row_eq && is_other_row_eq)
        {
            // swap(remaining_row, other_row)
            int temp = remaining_row_idx;
            remaining_row_idx = other_row_idx;
            other_row_idx = temp;

            is_remaining_row_eq = true;
            is_other_row_eq = false;
            remaining_row_new_rhs = rhs[remaining_row_idx];
            remaining_row_new_lhs = lhs[remaining_row_idx];
            remaining_row_coeff = A->x[A->p[remaining_row_idx].start];
        }

        assert(!(is_other_row_eq && !is_remaining_row_eq));

        other_row_coeff = A->x[A->p[other_row_idx].start];
        ratio = remaining_row_coeff / other_row_coeff;
        other_row_rhs_scaled = rhs[other_row_idx] * ratio;
        other_row_lhs_scaled = lhs[other_row_idx] * ratio;
        is_rhs_inf_other_row = HAS_TAG(row_tags[other_row_idx], R_TAG_RHS_INF);
        is_lhs_inf_other_row = HAS_TAG(row_tags[other_row_idx], R_TAG_LHS_INF);

        // -----------------------------------------------------------------------
        // if both rows are equalities we check that they are consistent
        // -----------------------------------------------------------------------
        if (is_remaining_row_eq && is_other_row_eq)
        {
            if (!IS_EQUAL_FEAS_TOL(remaining_row_new_rhs, other_row_rhs_scaled))
            {
                return INFEASIBLE;
            }
        }
        // -----------------------------------------------------------------------
        // if only the remaining row is an equality we check for infeasibility
        // -----------------------------------------------------------------------
        else if (is_remaining_row_eq)
        {
            bool infeas_rhs =
                !is_rhs_inf_other_row &&
                ((ratio > 0 && remaining_row_new_rhs > other_row_rhs_scaled) ||
                 (ratio < 0 && remaining_row_new_rhs < other_row_rhs_scaled));
            bool infeas_lhs =
                !is_lhs_inf_other_row &&
                ((ratio > 0 && remaining_row_new_rhs < other_row_lhs_scaled) ||
                 (ratio < 0 && remaining_row_new_rhs > other_row_lhs_scaled));

            if (infeas_rhs || infeas_lhs)
            {
                return INFEASIBLE;
            }
        }
        // ----------------------------------------------------------------------
        // if both rows are inequalities we find the tightest lhs and rhs
        // ----------------------------------------------------------------------
        else
        {
            assert(!is_remaining_row_eq && !is_other_row_eq);
            if (ratio > 0)
            {
                // possibly update rhs of remaining row
                if (!is_rhs_inf_other_row &&
                    (is_rhs_inf_remaining_row ||
                     other_row_rhs_scaled < remaining_row_new_rhs))
                {
                    other_row_rhs_index = other_row_idx;
                    other_row_rhs_scale_factor = ratio;
                    remaining_row_new_rhs = other_row_rhs_scaled;
                    is_rhs_inf_remaining_row = false;
                }

                // possibly update lhs of remaining row
                if (!is_lhs_inf_other_row &&
                    (is_lhs_inf_remaining_row ||
                     other_row_lhs_scaled > remaining_row_new_lhs))
                {
                    other_row_lhs_index = other_row_idx;
                    other_row_lhs_scale_factor = ratio;
                    remaining_row_new_lhs = other_row_lhs_scaled;
                    is_lhs_inf_remaining_row = false;
                }
            }
            else
            {
                // possibly update lhs of remaining row
                if (!is_rhs_inf_other_row &&
                    (is_lhs_inf_remaining_row ||
                     other_row_rhs_scaled > remaining_row_new_lhs))
                {
                    other_row_lhs_index = other_row_idx;
                    other_row_lhs_scale_factor = ratio;
                    remaining_row_new_lhs = other_row_rhs_scaled;
                    is_lhs_inf_remaining_row = false;
                }

                // possibly update rhs of remaining row
                if (!is_lhs_inf_other_row &&
                    (is_rhs_inf_remaining_row ||
                     other_row_lhs_scaled < remaining_row_new_rhs))
                {
                    other_row_rhs_index = other_row_idx;
                    other_row_rhs_scale_factor = ratio;
                    remaining_row_new_rhs = other_row_lhs_scaled;
                    is_rhs_inf_remaining_row = false;
                }
            }
        }
    }

    // ---------------------------------------------------------------------------------
    // if the remaining row is an inequality we may need to change the rhs and
    // lhs
    // ---------------------------------------------------------------------------------
    PostsolveInfo *postsolve_info = constraints->state->postsolve_info;
    if (!is_remaining_row_eq)
    {
        assert(!HAS_TAG(row_tags[remaining_row_idx], R_TAG_EQ));
        int start = A->p[remaining_row_idx].start;
        int len = A->p[remaining_row_idx].end - start;
        RowTag *row_tag = constraints->row_tags + remaining_row_idx;
        RowView remaining_row = new_rowview(A->x + start, A->i + start, &len, NULL,
                                            constraints->lhs + remaining_row_idx,
                                            constraints->rhs + remaining_row_idx,
                                            row_tag, remaining_row_idx);

        if (!is_rhs_inf_remaining_row &&
            remaining_row_new_rhs != rhs[remaining_row_idx])
        {
            assert(!IS_ABS_INF(remaining_row_new_rhs));
            assert(remaining_row_new_rhs < rhs[remaining_row_idx]);

            if (change_rhs_of_ineq(&remaining_row, remaining_row_new_rhs,
                                   constraints->state->col_locks,
                                   constraints->state->dton_rows, postsolve_info,
                                   other_row_rhs_index,
                                   other_row_rhs_scale_factor) == INFEASIBLE)
            {
                return INFEASIBLE;
            }
        }

        if (!is_lhs_inf_remaining_row &&
            remaining_row_new_lhs != lhs[remaining_row_idx])
        {
            assert(!IS_ABS_INF(remaining_row_new_lhs));
            assert(remaining_row_new_lhs > lhs[remaining_row_idx]);

            if (change_lhs_of_ineq(&remaining_row, remaining_row_new_lhs,
                                   constraints->state->col_locks,
                                   constraints->state->dton_rows, postsolve_info,
                                   other_row_lhs_index,
                                   other_row_lhs_scale_factor) == INFEASIBLE)
            {
                return INFEASIBLE;
            }
        }
    }

    // --------------------------------------------------------------------------------
    // When we arrive here we know which row should remain in the problem. We
    // mark all other rows as inactive.
    // --------------------------------------------------------------------------------
    for (int i = 0; i < bin_size; ++i)
    {
        if (bin[i] != remaining_row_idx)
        {
            set_row_to_inactive(bin[i], constraints->row_tags + bin[i],
                                constraints->state->rows_to_delete, postsolve_info,
                                0.0);
        }
    }

    return UNCHANGED;
}

// groups is a vector of size 'n_groups', where 'parallel_rows[i]' is the
// index of the first row in group 'i'
static PresolveStatus process_all_bins(const Constraints *constraints,
                                       const int *parallel_rows, const iVec *groups)
{
    if (groups->len == 0)
    {
        return UNCHANGED;
    }

    PSLP_uint n_groups = groups->len - 1;

    for (PSLP_uint i = 0; i < n_groups; ++i)
    {
        int n_rows_this_group = groups->data[i + 1] - groups->data[i];
        if (process_single_bin(constraints, parallel_rows + groups->data[i],
                               n_rows_this_group) == INFEASIBLE)
        {
            return INFEASIBLE;
        }
    }

    return UNCHANGED;
}

PresolveStatus remove_parallel_rows(Constraints *constraints)
{
    assert(constraints->state->ston_rows->len == 0);
    assert(constraints->state->empty_rows->len == 0);
    assert(constraints->state->empty_cols->len == 0);
    DEBUG(verify_problem_up_to_date(constraints));

    int *parallel_rows = constraints->state->work->iwork_n_rows;
    int *sparsity_IDs = constraints->state->work->iwork1_max_nrows_ncols;
    int *coeff_hashes = constraints->state->work->iwork2_max_nrows_ncols;
    iVec *group_starts = constraints->state->work->int_vec;

    find_parallel_rows(constraints->A, constraints->row_tags, group_starts,
                       parallel_rows, sparsity_IDs, coeff_hashes, R_TAG_INACTIVE);

    PresolveStatus status =
        process_all_bins(constraints, parallel_rows, group_starts);

    delete_inactive_rows(constraints);

    DEBUG(verify_problem_up_to_date(constraints));

    return status;
}
