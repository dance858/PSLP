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

#include "Matrix.h"
#include "Debugger.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "RowColViews.h"
#include "glbopts.h"
#include "stdlib.h"
#include "string.h"

Matrix *matrix_new(const double *Ax, const int *Ai, const int *Ap, int n_rows,
                   int n_cols, int nnz)
{
    DEBUG(ASSERT_NO_ZEROS_D(Ax, nnz));
    Matrix *A = matrix_alloc(n_rows, n_cols, nnz);
    RETURN_PTR_IF_NULL(A, NULL);
    int offset, i, row_size, row_alloc;

    offset = 0;
    for (i = 0; i < n_rows; ++i)
    {
        A->p[i].start = Ap[i] + offset;
        memcpy(A->x + A->p[i].start, Ax + Ap[i],
               (Ap[i + 1] - Ap[i]) * sizeof(double));

        memcpy(A->i + A->p[i].start, Ai + Ap[i], (Ap[i + 1] - Ap[i]) * sizeof(int));

        A->p[i].end = Ap[i + 1] + offset;
        row_size = A->p[i].end - A->p[i].start;
        row_alloc = calc_memory_row(row_size, EXTRA_ROW_SPACE, EXTRA_MEMORY_RATIO);
        offset += row_alloc - row_size;
    }

    A->p[n_rows].start = Ap[n_rows] + offset;
    A->p[n_rows].end = A->p[n_rows].start;

    return A;
}

// needed for the transpose function
Matrix *matrix_alloc(int n_rows, int n_cols, int nnz)
{
    Matrix *A = (Matrix *) ps_malloc(1, sizeof(Matrix));
    RETURN_PTR_IF_NULL(A, NULL);

    A->m = n_rows;
    A->n = n_cols;
    A->nnz = nnz;
    A->n_alloc = calc_memory(nnz, n_rows, EXTRA_ROW_SPACE, EXTRA_MEMORY_RATIO);

#ifdef TESTING
    A->i = (int *) ps_calloc(A->n_alloc, sizeof(int));
    A->p = (RowRange *) ps_calloc(n_rows + 1, sizeof(RowRange));
    A->x = (double *) ps_calloc(A->n_alloc, sizeof(double));
#else
    A->i = (int *) ps_malloc(A->n_alloc, sizeof(int));
    A->p = (RowRange *) ps_malloc(n_rows + 1, sizeof(RowRange));
    A->x = (double *) ps_malloc(A->n_alloc, sizeof(double));
#endif

    if (!A->i || !A->p || !A->x)
    {
        free_matrix(A);
        return NULL;
    }

    return A;
}

Matrix *matrix_new_no_extra_space(const double *Ax, const int *Ai, const int *Ap,
                                  int n_rows, int n_cols, int nnz)
{
    Matrix *A = (Matrix *) ps_malloc(1, sizeof(Matrix));
    RETURN_PTR_IF_NULL(A, NULL);

    A->m = n_rows;
    A->n = n_cols;
    A->nnz = nnz;
    A->n_alloc = nnz;
    A->i = (int *) ps_malloc(A->n_alloc, sizeof(int));
    A->p = (RowRange *) ps_malloc(n_rows + 1, sizeof(RowRange));
    A->x = (double *) ps_malloc(A->n_alloc, sizeof(double));

    if (!A->i || !A->p || !A->x)
    {
        free_matrix(A);
        return NULL;
    }

    memcpy(A->x, Ax, nnz * sizeof(double));
    memcpy(A->i, Ai, nnz * sizeof(int));

    for (int i = 0; i <= n_rows; ++i)
    {
        A->p[i].start = Ap[i];
        A->p[i].end = Ap[i + 1];
    }

    return A;
}

Matrix *transpose(const Matrix *A, int *work_n_cols)
{
    Matrix *AT = matrix_alloc(A->n, A->m, A->nnz);
    RETURN_PTR_IF_NULL(AT, NULL);
    int i, j, start;
    int *count = work_n_cols;
    memset(count, 0, A->n * sizeof(int));

    // -------------------------------------------------------------------
    //  compute nnz in each column of A
    // -------------------------------------------------------------------
    for (i = 0; i < A->m; ++i)
    {
        for (j = A->p[i].start; j < A->p[i].end; ++j)
        {
            count[A->i[j]]++;
        }
    }
    // ------------------------------------------------------------------
    //  compute row pointers, taking the extra space into account
    // ------------------------------------------------------------------
    AT->p[0].start = 0;
    for (i = 0; i < A->n; ++i)
    {
        start = AT->p[i].start;
        AT->p[i].end = start + count[i];
        AT->p[i + 1].start =
            start + calc_memory_row(count[i], EXTRA_ROW_SPACE, EXTRA_MEMORY_RATIO);
        count[i] = start;
    }

    AT->p[A->n].start = AT->n_alloc;
    AT->p[A->n].end = AT->n_alloc;

    // ------------------------------------------------------------------
    //  fill transposed matrix (this is a bottleneck)
    // ------------------------------------------------------------------
    for (i = 0; i < A->m; ++i)
    {
        for (j = A->p[i].start; j < A->p[i].end; j++)
        {
            AT->x[count[A->i[j]]] = A->x[j];
            AT->i[count[A->i[j]]] = i;
            count[A->i[j]]++;
        }
    }

    return AT;
}

int calc_memory(int nnz, int n_rows, int extra_row_space, double memory_ratio)
{
    return (int) (nnz * memory_ratio) + n_rows * extra_row_space;
}

int calc_memory_row(int size, int extra_row_space, double memory_ratio)
{
    return (int) (size * memory_ratio) + extra_row_space;
}

void free_matrix(Matrix *A)
{
    if (A)
    {
        PS_FREE(A->i);
        PS_FREE(A->p);
        PS_FREE(A->x);
    }

    PS_FREE(A);
}

void remove_extra_space(Matrix *A, const int *row_sizes, const int *col_sizes,
                        bool remove_all, int *col_idxs_map)
{
    int i, j, start, end, len, row_alloc, curr, n_deleted_rows, col_count;
    double extra_row_space = (remove_all) ? 0.0 : EXTRA_ROW_SPACE;
    double extra_mem_ratio = (remove_all) ? 1.0 : EXTRA_MEMORY_RATIO;
    curr = 0;
    n_deleted_rows = 0;

    // --------------------------------------------------------------------------
    // loop through the rows and remove redundant space, including inactive
    // rows.
    // --------------------------------------------------------------------------
    for (i = 0; i < A->m; ++i)
    {
        if (row_sizes[i] == SIZE_INACTIVE_ROW)
        {
            n_deleted_rows++;
            continue;
        }

        start = A->p[i].start;
        end = A->p[i].end;
        len = end - start;
        row_alloc = calc_memory_row(len, extra_row_space, extra_mem_ratio);
        memmove(A->x + curr, A->x + start, len * sizeof(double));
        memmove(A->i + curr, A->i + start, len * sizeof(int));
        A->p[i - n_deleted_rows].start = curr;
        A->p[i - n_deleted_rows].end = curr + len;
        curr += row_alloc;
    }

    A->m -= n_deleted_rows;
    A->p[A->m].start = curr;
    A->p[A->m].end = curr;

    // shrink size
    A->x = (double *) ps_realloc(A->x, MAX(curr, 1), sizeof(double));
    A->i = (int *) ps_realloc(A->i, MAX(curr, 1), sizeof(int));
    A->p = (RowRange *) ps_realloc(A->p, A->m + 1, sizeof(RowRange));

    // -------------------------------------------------------------------------
    //                      compute new column indices
    // -------------------------------------------------------------------------
    col_count = 0;
    for (i = 0; i < A->n; ++i)
    {
        if (col_sizes[i] == SIZE_INACTIVE_COL)
        {
            col_idxs_map[i] = -1;
        }
        else
        {
            col_idxs_map[i] = (col_count++);
        }
    }
    A->n = col_count;

    // -------------------------------------------------------------------------
    //                        update column indices
    // -------------------------------------------------------------------------
    for (i = 0; i < A->m; ++i)
    {
        for (j = A->p[i].start; j < A->p[i].end; ++j)
        {
            A->i[j] = col_idxs_map[A->i[j]];
        }
    }
}

bool shift_row(Matrix *A, int row, int extra_space, int max_shift)
{
    int left, right, missing_space, remaining_shifts, left_shifts;
    int right_shifts, space_left, space_right, n_move_right, n_move_left;
    int next_start, next_end;
    bool shift_left;
    RowRange *row_r = A->p;
    left = row;
    right = row + 1;
    remaining_shifts = max_shift;
    left_shifts = 0;
    right_shifts = 0;
    missing_space = extra_space - (row_r[right].start - row_r[row].end);

    if (missing_space <= 0)
    {
        return true;
    }

    // ------------------------------------------------------------------------
    // compute the new start index for row 'row' and a lower bound on the start
    // index for the the first active row following row 'row'.
    // ------------------------------------------------------------------------
    while (missing_space > 0)
    {
        if (left == 0 && right == A->m)
        {
            return false;
        }

        // space_left is the number of steps we can shift 'left' row to the
        // left without overwriting row 'left - 1'.
        // space_right is the number of steps we can shift 'right' row to
        // the right without overwriting row 'right + 1'.
        assert(left >= 0 && right <= A->m);
        space_left = (left == 0) ? 0 : row_r[left].start - row_r[left - 1].end;
        space_right =
            (right == A->m) ? 0 : row_r[right + 1].start - row_r[right].end;
        assert(space_left >= 0 && space_right >= 0);

        // number of elements that must be moved left resp. right if we shift
        // in a certain direction
        n_move_right = row_r[right].end - row_r[right].start;
        n_move_left = row_r[left].end - row_r[left].start;

        // decide which direction to shift
        if (left == 0)
        {
            if (right != A->m && n_move_right <= remaining_shifts)
            {
                shift_left = false;
            }
            else
            {
                return false;
            }
        }
        else if (right == A->m)
        {
            if (left != 0 && n_move_left <= remaining_shifts)
            {
                shift_left = true;
            }
            else
            {
                return false;
            }
        }
        else if (n_move_left == 0)
        {
            shift_left = true;
        }
        else if (n_move_right == 0)
        {
            shift_left = false;
        }
        else if (n_move_left <= remaining_shifts &&
                 (space_left / (double) (n_move_left) >=
                  space_right / (double) (n_move_right)))
        {
            shift_left = true;
        }
        else if (n_move_right <= remaining_shifts)
        {
            shift_left = false;
        }
        else
        {
            return false;
        }

        assert(!(shift_left && left == 0) && !(!shift_left && right == A->m));

        if (shift_left)
        {
            left_shifts = MIN(missing_space, space_left);
            missing_space -= left_shifts;
            remaining_shifts -= n_move_left;
            left -= 1;
        }
        else
        {
            right_shifts = MIN(missing_space, space_right);
            missing_space -= right_shifts;
            remaining_shifts -= n_move_right;
            right += 1;
        }
    }
    assert(remaining_shifts >= 0);

    // ------------------------------------------------------------------------
    //                 execute total left shift
    // ------------------------------------------------------------------------
    next_start = row_r[left + 1].start - left_shifts;
    for (; left < row; left++)
    {
        int len = row_r[left + 1].end - row_r[left + 1].start;
        if (len > 0)
        {
            memmove(A->x + next_start, A->x + row_r[left + 1].start,
                    len * sizeof(double));
            memmove(A->i + next_start, A->i + row_r[left + 1].start,
                    len * sizeof(int));
        }
        row_r[left + 1].start = next_start;
        row_r[left + 1].end = next_start + len;
        next_start += len;
    }

    // ------------------------------------------------------------------------
    //                 execute total right shift
    // ------------------------------------------------------------------------
    next_end = row_r[right - 1].end + right_shifts;
    for (; right > row + 1; right--)
    {
        int len = row_r[right - 1].end - row_r[right - 1].start;
        if (len > 0)
        {
            memmove(A->x + next_end - len, A->x + row_r[right - 1].start,
                    len * sizeof(double));
            memmove(A->i + next_end - len, A->i + row_r[right - 1].start,
                    len * sizeof(int));
        }
        row_r[right - 1].start = next_end - len;
        row_r[right - 1].end = next_end;
        next_end = row_r[right - 1].start;
    }

    assert(row_r[row + 1].start - row_r[row].end == extra_space);
    return true;
}

void print_row_starts(const RowRange *row_ranges, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        printf("%d ", row_ranges[i].start);
    }
    printf("\n");
}

double insert_or_update_coeff(Matrix *A, int row, int col, double val, int *row_size)
{
    int i, start, end, insertion;
    double old_val = 0.0;
    start = A->p[row].start;
    end = A->p[row].end;
    insertion = end;

    // -----------------------------------------------------------------
    //             find where it should be inserted
    // -----------------------------------------------------------------
    for (i = start; i < end; ++i)
    {
        if (A->i[i] >= col)
        {
            insertion = i;
            break;
        }
    }

    // -----------------------------------------------------------------
    // Insert the new value if it is nonzero. If it exists or should be
    // inserted in the end, we don't need to shift values.
    // -----------------------------------------------------------------
    if (ABS(val) > ZERO_TOL)
    {
        // assert(!IS_ZERO_FEAS_TOL(val));
        if (insertion == end)
        {
            A->x[insertion] = val;
            A->i[insertion] = col;
            A->p[row].end += 1;
            A->nnz += 1;
            *row_size += 1;
        }
        else if (A->i[insertion] == col)
        {
            old_val = A->x[insertion];
            A->x[insertion] = val;
        }
        else
        {
            memmove(A->x + insertion + 1, A->x + insertion,
                    (end - insertion) * sizeof(double));
            memmove(A->i + insertion + 1, A->i + insertion,
                    (end - insertion) * sizeof(int));

            // insert new value
            A->x[insertion] = val;
            A->i[insertion] = col;
            A->p[row].end += 1;
            A->nnz += 1;
            *row_size += 1;
        }
    }
    // if the new value is zero, we just have to shift
    else
    {
        // we only expect that the new value is zero if the coefficient
        // already exists
        assert(A->i[insertion] == col);

        // we only have to shift values if the zero is not in the end
        if (insertion != end - 1)
        {
            memmove(A->x + insertion, A->x + insertion + 1,
                    (end - insertion - 1) * sizeof(double));
            memmove(A->i + insertion, A->i + insertion + 1,
                    (end - insertion - 1) * sizeof(int));
        }

        A->p[row].end -= 1;
        A->nnz -= 1;
        *row_size -= 1;
    }

    assert(A->p[row].end <= A->p[row + 1].start);
    return old_val;
}

void remove_coeff(RowView *row, int col)
{
    int shift = 0;
    int len = *row->len;
    for (int i = 0; i < len; ++i)
    {
        if (row->cols[i] == col)
        {
            shift = 1;
        }

        row->vals[i] = row->vals[i + shift];
        row->cols[i] = row->cols[i + shift];
    }

    assert(shift != 0);
    (*row->range).end -= 1;
    *row->len -= 1;
}

void count_rows(const Matrix *A, int *row_sizes)
{
    printf("inside count rows\n");
    printf("A->m = %d\n", A->m);
    for (int i = 0; i < A->m; ++i)
    {
        row_sizes[i] = A->p[i].end - A->p[i].start;
    }
    printf("exiting count rows\n");
}

#ifdef TESTING
// Function to create a random CSR matrix
Matrix *random_matrix_new(int n_rows, int n_cols, double density)
{
    // allocate memory
    int n_alloc_nnz = (int) (density * n_rows * n_cols);
    double *Ax = (double *) ps_malloc(n_alloc_nnz, sizeof(double));
    int *Ai = (int *) ps_malloc(n_alloc_nnz, sizeof(int));
    int *Ap = (int *) ps_malloc(n_rows + 1, sizeof(int));
    if (!Ax || !Ai || !Ap)
    {
        PS_FREE(Ax);
        PS_FREE(Ai);
        PS_FREE(Ap);
        return NULL;
    }

    // Initialize random number generator
    srand(1);

    int nnz_count = 0; // Counter for nonzero elements
    Ap[0] = 0;

    for (int i = 0; i < n_rows; ++i)
    {
        int row_nnz = 0;

        // Randomly determine the number of nonzeros in this row
        for (int j = 0; j < n_cols; ++j)
        {
            if ((double) rand() / RAND_MAX < density)
            {
                if (nnz_count >= n_alloc_nnz)
                {
                    break;
                }

                Ax[nnz_count] = ((double) (rand() - rand()) / RAND_MAX) * 20.0;
                Ai[nnz_count] = j;
                ++nnz_count;
                ++row_nnz;
            }
        }
        Ap[i + 1] = Ap[i] + row_nnz;
    }

    // create matrix in modified CSR format
    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz_count);
    PS_FREE(Ax);
    PS_FREE(Ai);
    PS_FREE(Ap);

    return A;
}

void replace_row_A(Matrix *A, int row, double ratio, double *new_vals, int *cols_new,
                   int new_len)
{
    int i, len, start, n_new_elements;
    len = A->p[row].end - A->p[row].start;
    n_new_elements = new_len - len;

    // potentially shift row to get extra space
    if (n_new_elements > 0)
    {
        assert(shift_row(A, row, n_new_elements, 2000));
    }

    // replace the row
    start = A->p[row].start;
    for (i = 0; i < new_len; ++i)
    {
        A->x[start + i] = ratio * new_vals[i];
        A->i[start + i] = cols_new[i];
    }
    A->p[row].end = A->p[row].start + new_len;
    assert(A->p[row].end <= A->p[row + 1].start);
}

#endif
