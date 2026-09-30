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

#include "Matrix.h"
#include "Binary_search.h"
#include "Debugger.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "PSLP_warnings.h"
#include "RowColViews.h"
#include "glbopts.h"
#include "stdlib.h"
#include "string.h"
#include <assert.h>
#include <limits.h>

/* forward declaration */
static inline void remove_explicit_zeros(Matrix *A);

Matrix *matrix_new(const double *Ax, const int *Ai, const int *Ap, size_t n_rows,
                   size_t n_cols, size_t nnz)
{
    Matrix *A = matrix_alloc(n_rows, n_cols, nnz);
    RETURN_PTR_IF_NULL(A, NULL);

    memcpy(A->x, Ax, nnz * sizeof(double));
    memcpy(A->i, Ai, nnz * sizeof(int));

    for (size_t i = 0; i < n_rows; ++i)
    {
        A->p[i].start = Ap[i];
        A->p[i].end = Ap[i + 1];
    }
    A->p[n_rows].start = Ap[n_rows];
    A->p[n_rows].end = Ap[n_rows];

    /* the presolver assumes that only nonzero entries are stored */
    remove_explicit_zeros(A);

    return A;
}

// needed for the transpose function
Matrix *matrix_alloc(size_t n_rows, size_t n_cols, size_t nnz)
{
    Matrix *A = (Matrix *) ps_malloc(1, sizeof(Matrix));
    RETURN_PTR_IF_NULL(A, NULL);

    A->m = n_rows;
    A->n = n_cols;
    A->nnz = nnz;
    A->n_alloc = MAX(nnz, 1);

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

bool matrix_valid_csr_input(const int *Ai, const int *Ap, size_t n_rows,
                            size_t n_cols, size_t nnz)
{
    if (Ap[0] != 0)
    {
        return false;
    }

    for (size_t i = 0; i < n_rows; ++i)
    {
        if (Ap[i + 1] < Ap[i])
        {
            return false;
        }

        for (int j = Ap[i]; j < Ap[i + 1]; ++j)
        {
            if (Ai[j] < 0 || (size_t) Ai[j] >= n_cols ||
                (j > Ap[i] && Ai[j] <= Ai[j - 1]))
            {
                return false;
            }
        }
    }

    return Ap[n_rows] >= 0 && (size_t) Ap[n_rows] == nnz;
}

static inline void remove_explicit_zeros(Matrix *A)
{
    int i, j, shift;
    for (i = 0; i < A->m; ++i)
    {
        shift = 0;
        for (j = A->p[i].start; j < A->p[i].end; ++j)
        {
            if (A->x[j] == 0.0)
            {
                shift++;
            }
            else if (shift > 0)
            {
                A->x[j - shift] = A->x[j];
                A->i[j - shift] = A->i[j];
            }
        }
        A->p[i].end -= shift;
        A->nnz -= (size_t) shift;
    }
}

Matrix *transpose(const Matrix *A, int *work_n_cols, size_t tail)
{
    if (A->nnz + tail > (size_t) INT_MAX)
    {
        return NULL;
    }
    Matrix *AT = matrix_alloc(A->n, A->m, A->nnz + tail);
    RETURN_PTR_IF_NULL(AT, NULL);
    transpose_into(A, AT, work_n_cols);
    return AT;
}

void transpose_into(const Matrix *A, Matrix *AT, int *work_n_cols)
{
    assert(AT->m == A->n && AT->n == A->m && AT->n_alloc >= A->nnz);
    AT->nnz = A->nnz;

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
    //  compute row pointers
    // ------------------------------------------------------------------
    AT->p[0].start = 0;
    for (i = 0; i < A->n; ++i)
    {
        start = AT->p[i].start;
        AT->p[i].end = start + count[i];
        AT->p[i + 1].start = AT->p[i].end;
        count[i] = start;
    }
    AT->p[A->n].end = AT->p[A->n].start; // == nnz: the tail arena begins here

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

void remove_extra_space(Matrix *A, const int *row_sizes, const int *col_idxs_map,
                        size_t new_n_cols)
{
    int j, start, end, len, curr;
    curr = 0;
    size_t i, n_deleted_rows;

    // --------------------------------------------------------------------------
    // loop through the rows and remove redundant space, including inactive
    // rows.
    // --------------------------------------------------------------------------
    n_deleted_rows = 0;
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
        memmove(A->x + curr, A->x + start, (size_t) (len) * sizeof(double));
        memmove(A->i + curr, A->i + start, (size_t) (len) * sizeof(int));
        A->p[i - n_deleted_rows].start = curr;
        A->p[i - n_deleted_rows].end = curr + len;
        curr += len;
    }

    A->m -= n_deleted_rows;
    A->p[A->m].start = curr;
    A->p[A->m].end = curr;

    // shrink size
    A->x = (double *) ps_realloc(A->x, (size_t) MAX(curr, 1), sizeof(double));
    A->i = (int *) ps_realloc(A->i, (size_t) MAX(curr, 1), sizeof(int));
    A->p = (RowRange *) ps_realloc(A->p, (size_t) (A->m + 1), sizeof(RowRange));
    A->n_alloc = (size_t) MAX(curr, 1);

    // -------------------------------------------------------------------------
    //                        update column indices
    // -------------------------------------------------------------------------
    A->n = new_n_cols;
    for (i = 0; i < A->m; ++i)
    {
        for (j = A->p[i].start; j < A->p[i].end; ++j)
        {
            A->i[j] = col_idxs_map[A->i[j]];
        }
    }
}

void print_row_starts(const RowRange *row_ranges, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        printf("%d ", row_ranges[i].start);
    }
    printf("\n");
}

void remove_coeff(RowView *row, int col)
{
    int len = *row->len;
    int i = 0;

    // find the coefficient
    while (i < len && row->cols[i] != col)
    {
        ++i;
    }
    assert(i < len);

    // shift the remaining coefficients one step to the left. The loop stops
    // at len - 1 so that nothing is read past the end of the row.
    for (; i + 1 < len; ++i)
    {
        row->vals[i] = row->vals[i + 1];
        row->cols[i] = row->cols[i + 1];
    }

    (*row->range).end -= 1;
    *row->len -= 1;
}

void count_rows(const Matrix *A, int *row_sizes)
{
    for (int i = 0; i < A->m; ++i)
    {
        row_sizes[i] = A->p[i].end - A->p[i].start;
    }
}

#ifdef TESTING
// Function to create a random CSR matrix
Matrix *random_matrix_new(size_t n_rows, size_t n_cols, double density)
{
    // allocate memory

    /* disable conversion compiler warning*/
    PSLP_DIAG_PUSH();
    PSLP_DIAG_IGNORE_CONVERSION();
    /* intentional truncation */
    size_t n_alloc_nnz = (size_t) (density * n_rows * n_cols);

    /* enable conversion compiler warnings */
    PSLP_DIAG_POP();
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

    size_t nnz_count = 0; // Counter for nonzero elements
    Ap[0] = 0;

    for (size_t i = 0; i < n_rows; ++i)
    {
        size_t row_nnz = 0;

        // Randomly determine the number of nonzeros in this row
        for (size_t j = 0; j < n_cols; ++j)
        {
            if ((double) rand() / RAND_MAX < density)
            {
                if (nnz_count >= n_alloc_nnz)
                {
                    break;
                }

                Ax[nnz_count] = ((double) (rand() - rand()) / RAND_MAX) * 20.0;
                Ai[nnz_count] = (int) j;
                ++nnz_count;
                ++row_nnz;
            }
        }
        Ap[i + 1] = Ap[i] + (int) row_nnz;
    }

    // create matrix in modified CSR format
    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz_count);
    PS_FREE(Ax);
    PS_FREE(Ai);
    PS_FREE(Ap);

    return A;
}

/* Test helper: replaces row 'row' by ratio * new_vals on cols_new, shifting the
   rows after it (n_alloc must have room; the arrays are never reallocated, so
   pointers into earlier rows stay valid). */
void replace_row_A(Matrix *A, int row, double ratio, double *new_vals, int *cols_new,
                   int new_len)
{
    int old_start = A->p[row].start;
    int old_end = A->p[row].end;
    int total = A->p[A->m].start;
    int delta = new_len - (old_end - old_start);
    assert((size_t) (total + delta) <= A->n_alloc);

    memmove(A->x + old_end + delta, A->x + old_end,
            (size_t) (total - old_end) * sizeof(double));
    memmove(A->i + old_end + delta, A->i + old_end,
            (size_t) (total - old_end) * sizeof(int));
    for (int i = 0; i < new_len; ++i)
    {
        A->x[old_start + i] = ratio * new_vals[i];
        A->i[old_start + i] = cols_new[i];
    }
    A->nnz = (size_t) (total + delta);
    A->p[row].end = old_start + new_len;
    for (size_t r = (size_t) row + 1; r <= A->m; ++r)
    {
        A->p[r].start += delta;
        A->p[r].end += delta;
    }
}

#endif
