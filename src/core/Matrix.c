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
#include "thread_pool.h"
#include <assert.h>
#include <limits.h>

/* forward declaration */
static size_t remove_explicit_zeros(Matrix *A, int first, int last);

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
    A->nnz -= remove_explicit_zeros(A, 0, (int) n_rows);

    return A;
}

/* A chunk of matrix_new_parallel is a range of rows with about the same
   number of entries as the others. */
#define MATRIX_NEW_NNZ_PER_CHUNK 65536
#define MATRIX_NEW_MAX_CHUNKS (8 * THREAD_POOL_MAX_THREADS)
typedef struct
{
    const double *Ax;
    const int *Ai;
    const int *Ap;
    Matrix *A;
    size_t n_zeros[MATRIX_NEW_MAX_CHUNKS];
    int row_bounds[MATRIX_NEW_MAX_CHUNKS + 1];
} MatrixNewArgs;

static void matrix_new_chunk(void *context, int chunk, int thread)
{
    (void) thread;
    MatrixNewArgs *args = (MatrixNewArgs *) context;
    const int *Ap = args->Ap;
    Matrix *A = args->A;
    int first = args->row_bounds[chunk];
    int last = args->row_bounds[chunk + 1];
    size_t len = (size_t) (Ap[last] - Ap[first]);
    memcpy(A->x + Ap[first], args->Ax + Ap[first], len * sizeof(double));
    memcpy(A->i + Ap[first], args->Ai + Ap[first], len * sizeof(int));
    for (int row = first; row < last; ++row)
    {
        A->p[row].start = Ap[row];
        A->p[row].end = Ap[row + 1];
    }
    args->n_zeros[chunk] = remove_explicit_zeros(A, first, last);
}

Matrix *matrix_new_parallel(const double *Ax, const int *Ai, const int *Ap,
                            size_t n_rows, size_t n_cols, size_t nnz,
                            ThreadPool *pool)
{
    Matrix *A = matrix_alloc(n_rows, n_cols, nnz);
    RETURN_PTR_IF_NULL(A, NULL);
    MatrixNewArgs args = {Ax, Ai, Ap, A, {0}, {0}};
    int n_chunks = thread_pool_chunk_count(thread_pool_n_threads(pool), (int) nnz,
                                           MATRIX_NEW_NNZ_PER_CHUNK);
    if (n_chunks > MATRIX_NEW_MAX_CHUNKS) n_chunks = MATRIX_NEW_MAX_CHUNKS;
    if ((size_t) n_chunks > n_rows) n_chunks = n_rows > 0 ? (int) n_rows : 1;
    thread_pool_chunk_bounds_by_prefix(Ap, (int) n_rows, n_chunks, args.row_bounds);
    thread_pool_run(pool, matrix_new_chunk, &args, n_chunks);
    A->p[n_rows].start = Ap[n_rows];
    A->p[n_rows].end = Ap[n_rows];
    for (int chunk = 0; chunk < n_chunks; ++chunk)
    {
        A->nnz -= args.n_zeros[chunk];
    }
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

Matrix *matrix_alloc_diagonal(size_t n)
{
    Matrix *M = matrix_alloc(n, n, n);
    RETURN_PTR_IF_NULL(M, NULL);

    M->nnz = 0;
    for (size_t k = 0; k <= n; ++k)
    {
        M->p[k].start = M->p[k].end = (int) k;
    }
    return M;
}

void matrix_add_to_diagonal(Matrix *M, int k, double v)
{
    assert(M->p[k].start == k); // the layout of matrix_alloc_diagonal
    if (M->p[k].end == k)
    {
        M->i[k] = k;
        M->x[k] = 0.0;
        M->p[k].end = k + 1;
        M->nnz++;
    }
    M->x[k] += v;
}

bool matrix_is_symmetric(const Matrix *A)
{
    int j, k, l, rel;

    if (A->m != A->n)
    {
        return false;
    }

    for (k = 0; k < A->m; ++k)
    {
        for (j = A->p[k].start; j < A->p[k].end; ++j)
        {
            l = A->i[j];

            if (l == k)
            {
                continue;
            }

            rel = sorted_find(A->i + A->p[l].start, A->p[l].end - A->p[l].start, k);
            if (rel == -1 || A->x[A->p[l].start + rel] != A->x[j])
            {
                return false;
            }
        }
    }

    return true;
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

/* Drops the explicit zeros of the rows [first, last) and returns how many
   there were. Does not update A->nnz, so that ranges of rows can be handled
   on different threads. */
static size_t remove_explicit_zeros(Matrix *A, int first, int last)
{
    int i, j, shift;
    size_t n_zeros = 0;
    for (i = first; i < last; ++i)
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
        n_zeros += (size_t) shift;
    }
    return n_zeros;
}

Matrix *transpose(const Matrix *A, int *work_n_cols, size_t tail)
{
    if (A->nnz + tail > (size_t) INT_MAX)
    {
        return NULL;
    }
    Matrix *AT = matrix_alloc(A->n, A->m, A->nnz + tail);
    RETURN_PTR_IF_NULL(AT, NULL);
    transpose_into(A, AT, work_n_cols, NULL);
    return AT;
}

static bool transpose_on_pool(const Matrix *A, Matrix *AT, int *work_n_cols,
                              ThreadPool *pool);

void transpose_into(const Matrix *A, Matrix *AT, int *work_n_cols, ThreadPool *pool)
{
    assert(AT->m == A->n && AT->n == A->m && AT->n_alloc >= A->nnz);
    AT->nnz = A->nnz;
    if (transpose_on_pool(A, AT, work_n_cols, pool))
    {
        return;
    }

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

/* below this many nonzeros the transpose runs on the calling thread */
#define PARALLEL_TRANSPOSE_THRESHOLD 1000000

typedef struct
{
    const Matrix *A;
    Matrix *AT;
    int *count;
    const int
        *col_bounds;  /* chunk c owns columns [col_bounds[c], col_bounds[c+1]) */
    int *chunk_start; /* entries in the columns before those of chunk c */
} TransposeArgs;

/* counts the entries of the columns of one chunk */
static void transpose_count_chunk(void *context, int chunk, int thread)
{
    (void) thread;
    const TransposeArgs *args = (const TransposeArgs *) context;
    const Matrix *A = args->A;
    int *count = args->count;
    const unsigned lo = (unsigned) args->col_bounds[chunk];
    const unsigned width = (unsigned) args->col_bounds[chunk + 1] - lo;
    const int *cols = A->i;
    int total = 0;

    memset(count + lo, 0, width * sizeof(int));
    for (size_t i = 0; i < A->m; ++i)
    {
        const int end = A->p[i].end;
        for (int j = A->p[i].start; j < end; ++j)
        {
            if ((unsigned) cols[j] - lo < width)
            {
                count[cols[j]]++;
                total++;
            }
        }
    }
    args->chunk_start[chunk] = total;
}

/* row pointers of AT for the columns of one chunk, from their counts and the
   position of the chunk's first entry; the count of a column becomes its
   start */
static void transpose_pointers_chunk(void *context, int chunk, int thread)
{
    (void) thread;
    const TransposeArgs *args = (const TransposeArgs *) context;
    RowRange *p = args->AT->p;
    int *count = args->count;
    int start = args->chunk_start[chunk];
    const int end_col = args->col_bounds[chunk + 1];
    for (int i = args->col_bounds[chunk]; i < end_col; ++i)
    {
        p[i].start = start;
        p[i].end = start + count[i];
        count[i] = start;
        start = p[i].end;
    }
}

/* moves the entries of the columns of one chunk into AT */
static void transpose_fill_chunk(void *context, int chunk, int thread)
{
    (void) thread;
    const TransposeArgs *args = (const TransposeArgs *) context;
    const Matrix *A = args->A;
    Matrix *AT = args->AT;
    int *count = args->count;
    const unsigned lo = (unsigned) args->col_bounds[chunk];
    const unsigned width = (unsigned) args->col_bounds[chunk + 1] - lo;
    const int *cols = A->i;

    for (size_t i = 0; i < A->m; ++i)
    {
        const int end = A->p[i].end;
        for (int j = A->p[i].start; j < end; ++j)
        {
            const int col = cols[j];
            if ((unsigned) col - lo < width)
            {
                const int pos = count[col]++;
                AT->x[pos] = A->x[j];
                AT->i[pos] = (int) i;
            }
        }
    }
}

/* Transposes a large matrix on the pool. False, with nothing done, when the
   matrix is small or there is no pool with workers. */
static bool transpose_on_pool(const Matrix *A, Matrix *AT, int *work_n_cols,
                              ThreadPool *pool)
{
    int n_chunks = pool ? thread_pool_n_threads(pool) : 1;
    if (n_chunks < 2 || A->nnz < PARALLEL_TRANSPOSE_THRESHOLD ||
        A->n < (size_t) n_chunks)
    {
        return false;
    }

    int col_bounds[THREAD_POOL_MAX_THREADS + 1];
    int chunk_start[THREAD_POOL_MAX_THREADS + 1];
    int *count = work_n_cols;
    TransposeArgs args = {A, AT, count, col_bounds, chunk_start};
    const int n_cols = (int) A->n;

    /* count: even column ranges, as the column sizes are not known yet */
    thread_pool_chunk_bounds(n_cols, n_chunks, col_bounds);
    thread_pool_run(pool, transpose_count_chunk, &args, n_chunks);

    /* row pointers of AT, on the same column ranges */
    int start = 0;
    for (int c = 0; c < n_chunks; ++c)
    {
        int total = chunk_start[c];
        chunk_start[c] = start;
        start += total;
    }
    thread_pool_run(pool, transpose_pointers_chunk, &args, n_chunks);
    AT->p[n_cols].start = start; // == nnz: the tail arena begins here
    AT->p[n_cols].end = start;

    /* fill: column ranges with about the same number of entries each */
    col_bounds[0] = 0;
    for (int c = 1; c < n_chunks; ++c)
    {
        const int target = (int) (A->nnz * (size_t) c / (size_t) n_chunks);
        int lo = col_bounds[c - 1], hi = n_cols;
        while (lo < hi)
        {
            int mid = lo + (hi - lo) / 2;
            if (AT->p[mid].start < target)
            {
                lo = mid + 1;
            }
            else
            {
                hi = mid;
            }
        }
        col_bounds[c] = lo;
    }
    col_bounds[n_chunks] = n_cols;
    thread_pool_run(pool, transpose_fill_chunk, &args, n_chunks);
    return true;
}

Matrix *transpose_parallel(const Matrix *A, int *work_n_cols, size_t tail,
                           ThreadPool *pool)
{
    if (A->nnz + tail > (size_t) INT_MAX)
    {
        return NULL;
    }
    Matrix *AT = matrix_alloc(A->n, A->m, A->nnz + tail);
    RETURN_PTR_IF_NULL(AT, NULL);
    transpose_into(A, AT, work_n_cols, pool);
    return AT;
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
