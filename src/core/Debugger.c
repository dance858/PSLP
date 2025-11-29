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

#include "Debugger.h"
#include "Activity.h"
#include "Bounds.h"
#include "Constraints.h"
#include "Locks.h"
#include "Matrix.h"
#include "Numerics.h"
#include "PSLP_stats.h"
#include "State.h"
#include "glbopts.h"
#include <assert.h>
#include <stdio.h>

#ifdef __has_include
#if __has_include(<valgrind/valgrind.h>)
#include <valgrind/valgrind.h>
#else
#define RUNNING_ON_VALGRIND 0
#endif
#else
#include <valgrind/valgrind.h> // fallback
#endif

void run_debugger(const Constraints *constraints, bool finished)
{
    verify_activities(constraints);
    verify_row_and_col_sizes(constraints);
    verify_empty_rows(constraints);
    verify_ston_rows(constraints);
    verify_doubleton_rows(constraints);
    verify_empty_cols(constraints);
    verify_ston_cols(constraints);
    verify_no_duplicates_lists(constraints->state);
    verify_A_and_AT(constraints, false);
    verify_locks(constraints->AT, constraints->state->col_locks,
                 constraints->col_tags, constraints->row_tags);
    verify_row_tags(constraints);

    if (finished)
    {
        verify_empty_when_finished(constraints);
        verify_A_and_AT(constraints, true);
    }
}

int ASSERT_NO_ZEROS_D(const double *x, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        assert(x[i] != 0);
    }

    return 0;
}

void print_int_array(const int *arr, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

void print_double_array(const double *arr, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        printf("%f ", arr[i]);
    }
    printf("\n");
}

void verify_empty_rows(const Constraints *constraints)
{
    const iVec *empty_rows = constraints->state->empty_rows;
    const int *row_sizes = constraints->state->row_sizes;
    int n_rows = constraints->m;
    char *lookup = ps_calloc(n_rows, sizeof(char));

    // verify that all appended rows are empty rows
    for (int i = 0; i < empty_rows->len; ++i)
    {
        assert(row_sizes[empty_rows->data[i]] == 0);
        lookup[empty_rows->data[i]] = 1;
    }

    // verify that all empty rows have been appended using a look-up
    // table
    for (int i = 0; i < n_rows; i++)
    {
        if (row_sizes[i] == 0)
        {
            assert(lookup[i]);
        }
    }

    PS_FREE(lookup);
}

void verify_empty_cols(const Constraints *constraints)
{
    const iVec *empty_cols = constraints->state->empty_cols;
    const int *col_sizes = constraints->state->col_sizes;
    int n_cols = constraints->n;
    char *lookup = ps_calloc(n_cols, sizeof(char));

    // verify that all appended columns are empty columns
    for (int i = 0; i < empty_cols->len; ++i)
    {
        assert(col_sizes[empty_cols->data[i]] == 0);
        lookup[empty_cols->data[i]] = 1;
    }

    // verify that all empty rows have been appended using a look-up
    // table
    for (int i = 0; i < n_cols; i++)
    {
        if (col_sizes[i] == 0)
        {
            assert(lookup[i]);
        }
    }

    PS_FREE(lookup);
}

void verify_ston_cols(const Constraints *constraints)
{
    const iVec *ston_cols = constraints->state->ston_cols;
    const int *col_sizes = constraints->state->col_sizes;
    int n_cols = constraints->n;
    char *lookup = ps_calloc(n_cols, sizeof(char));

    // verify that all appended columns are ston columns or
    // empty/inactive
    for (int i = 0; i < ston_cols->len; ++i)
    {
        assert(col_sizes[ston_cols->data[i]] <= 1);
        lookup[ston_cols->data[i]] = 1;
    }

    // verify that all ston columns have been appended using a look-up
    // table
    for (int i = 0; i < constraints->n; i++)
    {
        if (col_sizes[i] == 1)
        {
            assert(lookup[i]);
        }
    }

    PS_FREE(lookup);
}

void verify_ston_rows(const Constraints *constraints)
{
    const iVec *ston_rows = constraints->state->ston_rows;
    const int *row_sizes = constraints->state->row_sizes;
    int n_rows = constraints->m;
    char *lookup = ps_calloc(n_rows, sizeof(char));

    // verify that all appended rows are empty rows
    for (int i = 0; i < ston_rows->len; ++i)
    {
        // if (!colTag.compare(ColTag::kFixed))??
        assert(row_sizes[ston_rows->data[i]] == 1);
        lookup[ston_rows->data[i]] = 1;
    }

    // verify that all ston rows have been appended using a look-up
    // table
    for (int i = 0; i < n_rows; i++)
    {
        if (row_sizes[i] == 1)
        {
            assert(lookup[i]);
        }
    }

    PS_FREE(lookup);
}

void verify_doubleton_rows(const Constraints *constraints)
{
    const iVec *doubleton_rows = constraints->state->dton_rows;
    const int *row_sizes = constraints->state->row_sizes;
    RowTag *row_tags = constraints->row_tags;
    int n_rows = constraints->m;
    char *lookup = ps_calloc(n_rows, sizeof(char));

    // verify that all appended rows have either size 2 and are equality rows,
    // or they have been reduced to size 1 or less
    for (int i = 0; i < doubleton_rows->len; ++i)
    {
        lookup[doubleton_rows->data[i]] = 1;
        assert(row_sizes[doubleton_rows->data[i]] <= 2);
        bool condition = HAS_TAG(row_tags[doubleton_rows->data[i]],
                                 (R_TAG_EQ | R_TAG_INACTIVE)) ||
                         row_sizes[doubleton_rows->data[i]] < 2;

        if (!condition)
        {
            printf("row_tags[%d] = %d\n", doubleton_rows->data[i],
                   row_tags[doubleton_rows->data[i]]);
            printf("row_sizes[%d] = %d\n", doubleton_rows->data[i],
                   row_sizes[doubleton_rows->data[i]]);
        }

        assert(condition);
    }

    // verify that all dton rows have been appended using a look-up
    // table
    for (int i = 0; i < n_rows; i++)
    {
        if (row_sizes[i] == 2 && HAS_TAG(row_tags[i], R_TAG_EQ))
        {
            assert(lookup[i]);
        }
    }

    PS_FREE(lookup);
}

// Comparison function for qsort
static int compare_ints(const void *a, const void *b)
{
    return (*(int *) a - *(int *) b);
}

void verify_no_duplicates(const iVec *vec)
{
    for (int i = 0; i < vec->len; ++i)
    {
        for (int j = i + 1; j < vec->len; ++j)
        {
            assert(vec->data[i] != vec->data[j]);
        }
    }
}

void verify_no_duplicates_sort_ptr(const int *data, int len)
{
    int *temp = (int *) ps_malloc(len, sizeof(int));
    memcpy(temp, data, len * sizeof(int));
    qsort(temp, len, sizeof(int), compare_ints);

    for (int i = 0; i < len - 1; ++i)
    {
        assert(temp[i] >= 0);
        assert(temp[i] != temp[i + 1]);
    }

    PS_FREE(temp);
}

void verify_no_duplicates_sort(const iVec *vec)
{
    verify_no_duplicates_sort_ptr(vec->data, vec->len);
}

void verify_no_duplicates_ptr(const int *data, int len)
{
    for (int i = 0; i < len; ++i)
    {
        for (int j = i + 1; j < len; ++j)
        {
            assert(data[i] != data[j]);
        }
    }
}

void verify_nonnegative_iVec(const iVec *vec)
{
    for (int i = 0; i < vec->len; ++i)
    {
        assert(vec->data[i] >= 0);
    }
}

void verify_no_duplicates_lists(const State *data)
{
    verify_nonnegative_iVec(data->updated_activities);
    verify_nonnegative_iVec(data->dton_rows);
    verify_nonnegative_iVec(data->ston_rows);
    verify_nonnegative_iVec(data->fixed_cols_to_delete);
    verify_nonnegative_iVec(data->sub_cols_to_delete);
    verify_nonnegative_iVec(data->rows_to_delete);

    verify_no_duplicates_sort(data->updated_activities);
    verify_no_duplicates_sort(data->dton_rows);
    verify_no_duplicates_sort(data->ston_rows);
    verify_no_duplicates_sort(data->fixed_cols_to_delete);
    verify_no_duplicates_sort(data->sub_cols_to_delete);
    verify_no_duplicates_sort(data->rows_to_delete);
}

void verify_row_tags(const Constraints *constraints)
{
    const RowTag *row_tags = constraints->row_tags;
    const double *lhs = constraints->lhs;
    const double *rhs = constraints->rhs;

    for (int i = 0; i < constraints->m; ++i)
    {
        if (IS_NEG_INF(lhs[i]) && !HAS_TAG(row_tags[i], R_TAG_INACTIVE))
        {
            assert(HAS_TAG(row_tags[i], R_TAG_LHS_INF));
        }

        if (IS_POS_INF(rhs[i]) && !HAS_TAG(row_tags[i], R_TAG_INACTIVE))
        {
            assert(HAS_TAG(row_tags[i], R_TAG_RHS_INF));
        }

        if (HAS_TAG(row_tags[i], R_TAG_LHS_INF))
        {
            assert(IS_NEG_INF(lhs[i]));
        }

        if (HAS_TAG(row_tags[i], R_TAG_RHS_INF))
        {
            assert(IS_POS_INF(rhs[i]));
        }

        if (!HAS_TAG(row_tags[i], R_TAG_LHS_INF) &&
            !HAS_TAG(row_tags[i], R_TAG_RHS_INF) && lhs[i] == rhs[i])
        {
            assert(HAS_TAG(row_tags[i], (R_TAG_EQ | R_TAG_INACTIVE)));
        }
    }
}

static void verify_CSR_matrix(const Matrix *A, bool compressed)
{
    int i, j, nnz;
    assert(A->m >= 0 && A->n >= 0 && A->nnz >= 0);
    assert(A->n_alloc >= A->nnz);

    if (A->nnz > 0)
    {
        assert(A->i != NULL && A->x != NULL && A->p != NULL);
    }
    // verify that row ranges are ascending
    for (i = 0; i < A->m; ++i)
    {
        assert(A->p[i + 1].start >= A->p[i].end && A->p[i].end >= A->p[i].start);
    }

    // verify that columns within a row are sorted
    for (i = 0; i < A->m; ++i)
    {
        for (j = A->p[i].start + 1; j < A->p[i].end; ++j)
        {
            assert(A->i[j] > A->i[j - 1]);
        }
    }

    // verify that no explicit zeros are stored are that nnz count is correct
    nnz = 0;
    for (i = 0; i < A->m; ++i)
    {
        for (j = A->p[i].start; j < A->p[i].end; ++j)
        {
            assert(A->x[j] != 0);
        }

        nnz += A->p[i].end - A->p[i].start;
    }

    assert(nnz == A->nnz);

    if (compressed)
    {
        assert(A->p[A->m].start == A->nnz);
        assert(A->p[A->m].end == A->nnz);
        assert(A->p[0].start == 0);
    }
}

bool verify_A_and_AT_consistency(const Matrix *A, const Matrix *AT)
{
    int *work_n_cols = (int *) ps_malloc(A->n, sizeof(int));
    Matrix *real_AT = transpose(A, work_n_cols);

    // check that nnz and dimensions are consistent
    assert(real_AT->m == AT->m);
    assert(real_AT->n == AT->n);
    assert(real_AT->nnz == AT->nnz);

    int i, j, k, start_AT, end_AT, start_real_AT, end_real_AT;

    // check that values are consistent
    for (i = 0; i < AT->m; ++i)
    {
        start_AT = AT->p[i].start;
        end_AT = AT->p[i].end;
        start_real_AT = real_AT->p[i].start;
        end_real_AT = real_AT->p[i].end;

        assert(end_AT - start_AT == end_real_AT - start_real_AT);

        for (j = start_AT, k = start_real_AT; j < end_AT; ++j, ++k)
        {
            assert(!IS_ABS_INF(AT->x[j]) && !IS_ABS_INF(real_AT->x[k]));
            if (AT->x[j] != real_AT->x[k] || AT->i[j] != real_AT->i[k])
            {
                PS_FREE(work_n_cols);
                free_matrix(real_AT);
                return false;
            }
        }
    }

    PS_FREE(work_n_cols);
    free_matrix(real_AT);
    return true;
}

static void verify_no_inactive_cols(const Matrix *A, ColTag *col_tags)
{
    int i, j, start, end;
    for (i = 0; i < A->m; ++i)
    {
        start = A->p[i].start;
        end = A->p[i].end;

        for (j = start; j < end; ++j)
        {
            assert(!HAS_TAG(col_tags[A->i[j]], C_TAG_INACTIVE));
        }
    }
}

void verify_A_and_AT(const Constraints *constraints, bool compressed)
{
    const Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT;

    // verify that both A and AT are valid CSR matrices
    verify_CSR_matrix(A, compressed);
    verify_CSR_matrix(AT, compressed);

    // verify that A and AT are consistent
    assert(verify_A_and_AT_consistency(A, AT));

    // verify that A does not store the coefficients of inactive columns
    verify_no_inactive_cols(A, constraints->col_tags);

    // verify that AT does not store the coefficients of inactive rows
    verify_no_inactive_cols(AT, constraints->row_tags);
}

void verify_locks(const Matrix *AT, const Lock *locks, const ColTag *col_tags,
                  const RowTag *row_tags)
{
    int up, down, i, col, row, n_cols;
    n_cols = AT->m;
    double Aik;

    for (col = 0; col < n_cols; ++col)
    {
        if (HAS_TAG(col_tags[col], C_TAG_INACTIVE))
        {
            continue;
        }

        up = 0;
        down = 0;

        for (i = AT->p[col].start; i < AT->p[col].end; ++i)
        {
            row = AT->i[i];
            Aik = AT->x[i];

            if ((Aik > 0 && !HAS_TAG(row_tags[row], R_TAG_RHS_INF)) ||
                (Aik < 0 && !HAS_TAG(row_tags[row], R_TAG_LHS_INF)))
            {
                up++;
            }

            if ((Aik > 0 && !HAS_TAG(row_tags[row], R_TAG_LHS_INF)) ||
                (Aik < 0 && !HAS_TAG(row_tags[row], R_TAG_RHS_INF)))
            {
                down++;
            }
        }

        assert(up == locks[col].up);
        assert(down == locks[col].down);
    }
}

void verify_empty_when_finished(const Constraints *constraints)
{
    assert(constraints->state->empty_cols->len == 0);
    assert(constraints->state->empty_rows->len == 0);
    assert(constraints->state->rows_to_delete->len == 0);
    assert(constraints->state->fixed_cols_to_delete->len == 0);
    assert(constraints->state->sub_cols_to_delete->len == 0);

    for (int i = 0; i < constraints->A->m; ++i)
    {
        assert(!HAS_TAG(constraints->row_tags[i], R_TAG_INACTIVE));
    }

    for (int i = 0; i < constraints->A->n; ++i)
    {
        assert(!HAS_TAG(constraints->col_tags[i], C_TAG_INACTIVE));
    }
}

void verify_row_and_col_sizes(const Constraints *constraints)
{
    const int *row_sizes = constraints->state->row_sizes;
    const int *col_sizes = constraints->state->col_sizes;
    const Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT;
    const ColTag *col_tags = constraints->col_tags;
    const RowTag *row_tags = constraints->row_tags;
    int i;

    // check that row sizes are correct
    for (i = 0; i < constraints->A->m; ++i)
    {
        if (HAS_TAG(row_tags[i], R_TAG_INACTIVE))
        {
            assert(row_sizes[i] == SIZE_INACTIVE_ROW);
            assert(A->p[i].start == A->p[i].end);
        }
        else
        {
            assert(row_sizes[i] == A->p[i].end - A->p[i].start);
        }

        if (row_sizes[i] == SIZE_INACTIVE_ROW)
        {
            assert(HAS_TAG(row_tags[i], R_TAG_INACTIVE));
        }
    }

    // check that column sizes are correct
    for (i = 0; i < constraints->A->n; ++i)
    {
        if (HAS_TAG(col_tags[i], C_TAG_INACTIVE))
        {
            if (col_sizes[i] != SIZE_INACTIVE_COL)
            {
                printf("col_sizes[%d] = %d\n", i, col_sizes[i]);
            }

            assert(col_sizes[i] == SIZE_INACTIVE_COL);
            if (AT->p[i].start != AT->p[i].end)
            {
                printf("AT->p[%d].start = %d\n", i, AT->p[i].start);
                printf("AT->p[%d].end = %d\n", i, AT->p[i].end);
            }
            assert(AT->p[i].start == AT->p[i].end);
        }
        else
        {
            assert(col_sizes[i] == AT->p[i].end - AT->p[i].start);
        }

        if (col_sizes[i] == SIZE_INACTIVE_COL)
        {
            assert(HAS_TAG(col_tags[i], C_TAG_INACTIVE));
        }
    }
}

void verify_activity(const ColTag *col_tags, const Bound *bounds, Activity activity,
                     RowTag row_tag, const int *cols, const double *vals, int len)
{
    int n_inf_min = 0;
    int n_inf_max = 0;
    double min = 0.0;
    double max = 0.0;
    int i, col;

    if (HAS_TAG(row_tag, R_TAG_INACTIVE))
    {
        return;
    }

    // count the infinite contributions
    for (i = 0; i < len; ++i)
    {
        col = cols[i];
        if (HAS_TAG(col_tags[col], C_TAG_INACTIVE))
        {
            continue;
        }

        if (vals[i] > 0)
        {
            if (HAS_TAG(col_tags[col], C_TAG_UB_INF))
            {
                n_inf_max++;
            }

            if (HAS_TAG(col_tags[col], C_TAG_LB_INF))
            {
                n_inf_min++;
            }
        }
        else
        {
            if (HAS_TAG(col_tags[col], C_TAG_LB_INF))
            {
                n_inf_max++;
            }

            if (HAS_TAG(col_tags[col], C_TAG_UB_INF))
            {
                n_inf_min++;
            }
        }
    }

    assert(n_inf_max == activity.n_inf_max);
    assert(n_inf_min == activity.n_inf_min);

    // if there is a max infinite contribution, the max activity should be
    // INVALID_ACT_DEBUG
    if (n_inf_max > 0)
    {
        if (activity.max != INVALID_ACT_DEBUG)
        {
            printf("n_inf_max = %d\n", n_inf_max);
            printf("activity.max = %f\n", activity.max);
        }
        assert(activity.max == INVALID_ACT_DEBUG);
    }
    else
    {
        // count max activity
        for (i = 0; i < len; ++i)
        {
            col = cols[i];
            if (HAS_TAG(col_tags[col], C_TAG_INACTIVE))
            {
                continue;
            }

            if (vals[i] > 0)
            {
                assert(!HAS_TAG(col_tags[col], C_TAG_UB_INF) &&
                       !IS_POS_INF(bounds[col].ub));
                max += vals[i] * bounds[col].ub;
            }
            else
            {
                assert(!HAS_TAG(col_tags[col], C_TAG_LB_INF) &&
                       !IS_NEG_INF(bounds[col].lb));
                max += vals[i] * bounds[col].lb;
            }
        }
    }

    // if there is a min infinite contribution, the min activity should be
    // INVALID_ACT_DEBUG
    if (n_inf_min > 0)
    {
        assert(activity.min == INVALID_ACT_DEBUG);
    }
    else
    {
        // count min activity
        for (i = 0; i < len; ++i)
        {
            col = cols[i];
            if (HAS_TAG(col_tags[col], C_TAG_INACTIVE))
            {
                continue;
            }

            if (vals[i] > 0)
            {
                assert(!HAS_TAG(col_tags[col], C_TAG_LB_INF) &&
                       !IS_NEG_INF(bounds[col].lb));
                min += vals[i] * bounds[col].lb;
            }
            else
            {
                assert(!HAS_TAG(col_tags[col], C_TAG_UB_INF) &&
                       !IS_POS_INF(bounds[col].ub));
                min += vals[i] * bounds[col].ub;
            }
        }
    }
}

void verify_activities(const Constraints *constraints)
{
    int n_rows = constraints->m;
    const Activity *activities = constraints->state->activities;
    const RowTag *row_tags = constraints->row_tags;
    const RowRange *row_ranges = constraints->A->p;
    const int *cols = constraints->A->i;
    const double *vals = constraints->A->x;

    for (int i = 0; i < n_rows; ++i)
    {
        verify_activity(constraints->col_tags, constraints->bounds, activities[i],
                        row_tags[i], cols + row_ranges[i].start,
                        vals + row_ranges[i].start,
                        row_ranges[i].end - row_ranges[i].start);
    }
}

void ASSERT_NO_ACTIVE_STON_ROWS(const Matrix *A, const RowTag *row_tags)
{
    for (int i = 0; i < A->m; ++i)
    {
        assert(A->p[i].end - A->p[i].start != 1 ||
               HAS_TAG(row_tags[i], R_TAG_INACTIVE));
    }
}

int ASSERT_INCREASING_I(const int *x, size_t len)
{
    for (size_t i = 1; i < len; ++i)
    {
        assert(x[i] > x[i - 1]);
    }

    return 0;
}

void print_matrix(const Matrix *A)
{
    for (int i = 0; i < A->m; i++)
    {
        int k = A->p[i].start;

        for (int j = 0; j < A->n; j++)
        {
            if (k < A->p[i].end && A->i[k] == j)
            {
                printf("%6.2f ", A->x[k]);
                k++;
            }
            else
            {
                printf("%6.2f ", 0.0);
            }
        }
        printf("\n");
    }
}

void verify_problem_up_to_date(const Constraints *constraints)
{
    assert(constraints->state->rows_to_delete->len == 0);
    assert(constraints->state->fixed_cols_to_delete->len == 0);
    assert(constraints->state->sub_cols_to_delete->len == 0);
}

void verify_row_states(const Activity *acts, const iVec *updated_activities)
{
    for (int i = 0; i < updated_activities->len; ++i)
    {
        int row = updated_activities->data[i];
        assert(acts[row].status == ADDED);
        assert(acts[row].n_inf_min == 0 || acts[row].n_inf_max == 0);
    }
}

void run_debugger_stats_consistency_check(const PresolveStats *stats)
{
    int total_removed = stats->nnz_removed_trivial + stats->nnz_removed_fast +
                        stats->nnz_removed_primal_propagation +
                        stats->nnz_removed_parallel_rows +
                        stats->nnz_removed_parallel_cols;

    assert(stats->nnz_original - stats->nnz_reduced == total_removed);

    // valgrind timing is not reliable
#ifndef RUNNING_ON_VALGRIND
    double time_medium = stats->ps_time_primal_propagation +
                         stats->ps_time_parallel_rows + stats->ps_time_parallel_cols;
    assert((stats->ps_time_medium - time_medium) / MAX(1e-2, time_medium) < 0.05);
#endif
}
