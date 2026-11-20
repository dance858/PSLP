#ifndef TEST_MATRIX_H
#define TEST_MATRIX_H

#include "Constraints.h"
#include "Debugger.h"
#include "Matrix.h"
#include "glbopts.h"
#include "minunit.h"
#include "test_macros.h"
#include "thread_pool.h"
#include <limits.h>
#include <stdio.h>

int counter_matrix = 0;

/* Summary of tests:
Test 0: matrix_new: allocation and free of a matrix.
Test 1: transpose: transpose a matrix.
Test 2: transpose_into: transpose into an over-allocated matrix, twice.
transpose_parallel: same result as transpose for several pool sizes.
Test 2-14:  shiftRow
Then some adhoc tests.
*/

// test allocation and free
static char *test_0_matrix()
{

    double vals[] = {1, -1, 2, 1, 1, 1, 1, 3, 1, 1, 1};
    int cols[] = {0, 1, 2, 1, 2, 4, 4, 0, 1, 2, 3};
    int row_starts[] = {0, 3, 6, 7, 9, 10, 11};

    int work_n_cols[5];

    Matrix *A = matrix_new(vals, cols, row_starts, 6, 5, 11);
    Matrix *AT = transpose(A, work_n_cols, 0);

    mu_assert("error", A);
    mu_assert("error", AT);

    free_matrix(A);
    free_matrix(AT);

    return 0;
}

// test transpose
static char *test_1_matrix()
{

    double vals[] = {1, -1, 2, 1, 1, 1, 1, 3, 1, 1, 1};
    int cols[] = {0, 1, 2, 1, 2, 4, 4, 0, 1, 2, 3};
    int row_starts[] = {0, 3, 6, 7, 9, 10, 11};
    int work_n_cols[5];

    Matrix *A = matrix_new(vals, cols, row_starts, 6, 5, 11);
    Matrix *AT = transpose(A, work_n_cols, 0);

    // compact to drop nothing: checks the transpose layout
    int row_sizes_AT[6] = {2, 3, 3, 1, 2};
    int col_sizes_AT[6] = {3, 3, 1, 2, 1, 1};
    int col_idxs_map_AT[6] = {0};
    int n_new_cols_AT = update_column_map(col_sizes_AT, col_idxs_map_AT, 6);
    remove_extra_space(AT, row_sizes_AT, col_idxs_map_AT, n_new_cols_AT);

    // correct answer
    double AT_vals_correct[] = {1, 3, -1, 1, 1, 2, 1, 1, 1, 1, 1};
    int AT_cols_correct[] = {0, 3, 0, 1, 3, 0, 1, 4, 5, 1, 2};
    int AT_row_starts_correct[] = {0, 2, 5, 8, 9, 11};

    mu_assert("error, vals not equal",
              ARRAYS_EQUAL_DOUBLE(AT_vals_correct, AT->x, 11));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(AT_cols_correct, AT->i, 11));
    mu_assert("row starts", check_row_starts(AT, AT_row_starts_correct));

    free_matrix(A);
    free_matrix(AT);

    return 0;
}

// transpose_into reuses an over-allocated AT and leaves its tail alone
static char *test_2_matrix()
{
    double vals[] = {1, -1, 2, 1, 1, 1, 1, 3, 1, 1, 1};
    int cols[] = {0, 1, 2, 1, 2, 4, 4, 0, 1, 2, 3};
    int row_starts[] = {0, 3, 6, 7, 9, 10, 11};
    int work_n_cols[5];

    Matrix *A = matrix_new(vals, cols, row_starts, 6, 5, 11);
    Matrix *AT = matrix_alloc(5, 6, 11 + 7); // room for a tail of 7
    for (size_t p = 0; p < AT->n_alloc; ++p)
    {
        AT->i[p] = -7;
        AT->x[p] = -7.0;
    }
    transpose_into(A, AT, work_n_cols, NULL);

    double AT_vals_correct[] = {1, 3, -1, 1, 1, 2, 1, 1, 1, 1, 1};
    int AT_cols_correct[] = {0, 3, 0, 1, 3, 0, 1, 4, 5, 1, 2};
    int AT_row_starts_correct[] = {0, 2, 5, 8, 9, 11};
    mu_assert("vals", ARRAYS_EQUAL_DOUBLE(AT_vals_correct, AT->x, 11));
    mu_assert("cols", ARRAYS_EQUAL_INT(AT_cols_correct, AT->i, 11));
    mu_assert("row starts", check_row_starts(AT, AT_row_starts_correct));
    mu_assert("nnz", AT->nnz == 11);
    mu_assert("sentinel", AT->p[5].start == 11 && AT->p[5].end == 11);
    mu_assert("n_alloc kept", AT->n_alloc == 18);
    for (size_t p = 11; p < 18; ++p)
    {
        mu_assert("tail untouched", AT->i[p] == -7 && AT->x[p] == -7.0);
    }

    // row 0 of A loses its last entry (column 2); transpose again in place
    A->p[0].end -= 1;
    A->nnz -= 1;
    transpose_into(A, AT, work_n_cols, NULL);

    double AT_vals_2[] = {1, 3, -1, 1, 1, 1, 1, 1, 1, 1};
    int AT_cols_2[] = {0, 3, 0, 1, 3, 1, 4, 5, 1, 2};
    int AT_row_starts_2[] = {0, 2, 5, 7, 8, 10};
    mu_assert("vals 2", ARRAYS_EQUAL_DOUBLE(AT_vals_2, AT->x, 10));
    mu_assert("cols 2", ARRAYS_EQUAL_INT(AT_cols_2, AT->i, 10));
    mu_assert("row starts 2", check_row_starts(AT, AT_row_starts_2));
    mu_assert("nnz 2", AT->nnz == 10);
    mu_assert("sentinel 2", AT->p[5].start == 10 && AT->p[5].end == 10);
    mu_assert("n_alloc kept 2", AT->n_alloc == 18);

    free_matrix(A);
    free_matrix(AT);
    return 0;
}

// matrix_valid_csr_input on the 2 x 3 matrix [1 0 2; 0 3 0]
static char *test_20_matrix()
{
    int Ai[] = {0, 2, 1};
    int Ap[] = {0, 2, 3};
    mu_assert("valid", matrix_valid_csr_input(Ai, Ap, 2, 3, 3));

    int Ap_empty[] = {0, 0, 0};
    mu_assert("empty", matrix_valid_csr_input(NULL, Ap_empty, 2, 3, 0));

    int Ap_first[] = {1, 2, 3};
    mu_assert("Ap[0] != 0", !matrix_valid_csr_input(Ai, Ap_first, 2, 3, 3));

    int Ap_decreasing[] = {0, 2, 1};
    mu_assert("decreasing Ap", !matrix_valid_csr_input(Ai, Ap_decreasing, 2, 3, 1));

    mu_assert("Ap[m] != nnz", !matrix_valid_csr_input(Ai, Ap, 2, 3, 2));

    int Ai_large[] = {0, 3, 1};
    mu_assert("index >= n", !matrix_valid_csr_input(Ai_large, Ap, 2, 3, 3));

    int Ai_negative[] = {-1, 2, 1};
    mu_assert("negative index", !matrix_valid_csr_input(Ai_negative, Ap, 2, 3, 3));

    int Ai_unsorted[] = {2, 0, 1};
    mu_assert("unsorted row", !matrix_valid_csr_input(Ai_unsorted, Ap, 2, 3, 3));

    int Ai_duplicate[] = {0, 0, 1};
    mu_assert("duplicate index", !matrix_valid_csr_input(Ai_duplicate, Ap, 2, 3, 3));

    return 0;
}

/* matrix_alloc_diagonal and matrix_add_to_diagonal on a 3 by 3 matrix: the
   rows start empty in their reserved positions; adding to a row creates its
   entry once and accumulates after that, leaving the other rows alone. */
static char *test_diagonal_matrix()
{
    Matrix *M = matrix_alloc_diagonal(3);
    mu_assert("alloc", M != NULL && M->m == 3 && M->n == 3 && M->nnz == 0);
    for (int k = 0; k <= 3; ++k)
    {
        mu_assert("empty rows", M->p[k].start == k && M->p[k].end == k);
    }

    matrix_add_to_diagonal(M, 1, 2.0);
    mu_assert("entry created",
              M->nnz == 1 && M->p[1].end == 2 && M->i[1] == 1 && M->x[1] == 2.0);
    matrix_add_to_diagonal(M, 1, 0.5);
    mu_assert("accumulated", M->nnz == 1 && M->p[1].end == 2 && M->x[1] == 2.5);
    mu_assert("others untouched", M->p[0].end == 0 && M->p[2].end == 2);

    matrix_add_to_diagonal(M, 2, -1.0);
    mu_assert("second row", M->nnz == 2 && M->i[2] == 2 && M->x[2] == -1.0);
    mu_assert("sorted and symmetric", matrix_is_symmetric(M));

    free_matrix(M);
    return 0;
}

/* transpose_parallel gives the same matrix as transpose, for a matrix above the
   parallel threshold and pools of one, two and eight threads, and falls back
   to transpose without a pool. */
static char *test_transpose_parallel()
{
    Matrix *A = random_matrix_new(3000, 5000, 0.1);
    mu_assert("random matrix", A != NULL && A->nnz >= 1000000);
    int *work_n_cols = (int *) ps_malloc(A->n, sizeof(int));
    const size_t tail = 7;
    Matrix *ref = transpose(A, work_n_cols, tail);
    mu_assert("reference transpose", ref != NULL);

    const int n_workers[] = {-1, 0, 1, 7}; /* -1: no pool */
    for (int k = 0; k < 4; ++k)
    {
        ThreadPool *pool = n_workers[k] < 0 ? NULL : thread_pool_new(n_workers[k]);
        mu_assert("pool", n_workers[k] < 0 || pool != NULL);
        Matrix *AT = transpose_parallel(A, work_n_cols, tail, pool);
        mu_assert("parallel transpose", AT != NULL);
        mu_assert("dimensions", AT->m == ref->m && AT->n == ref->n &&
                                    AT->nnz == ref->nnz &&
                                    AT->n_alloc == ref->n_alloc);
        mu_assert("row pointers",
                  memcmp(AT->p, ref->p, (ref->m + 1) * sizeof(RowRange)) == 0);
        mu_assert("column indices",
                  memcmp(AT->i, ref->i, ref->nnz * sizeof(int)) == 0);
        mu_assert("values", memcmp(AT->x, ref->x, ref->nnz * sizeof(double)) == 0);
        free_matrix(AT);
        thread_pool_free(pool);
    }

    free_matrix(ref);
    free_matrix(A);
    PS_FREE(work_n_cols);
    return 0;
}

/* matrix_new_parallel gives the same matrix as matrix_new on input with
   explicit zeros, for pools of one, two and eight threads. */
static char *test_matrix_new_parallel()
{
    Matrix *R = random_matrix_new(1000, 2000, 0.2);
    mu_assert("random matrix", R != NULL && R->nnz >= 300000);
    size_t nnz = R->nnz;
    int *Ap = (int *) ps_malloc(R->m + 1, sizeof(int));
    for (size_t i = 0; i <= R->m; ++i)
    {
        Ap[i] = R->p[i].start;
    }
    mu_assert("contiguous rows", Ap[0] == 0 && (size_t) Ap[R->m] == nnz);
    for (size_t j = 0; j < nnz; j += 7)
    {
        R->x[j] = 0.0;
    }
    /* a row of zeros only */
    for (int j = Ap[500]; j < Ap[501]; ++j)
    {
        R->x[j] = 0.0;
    }

    Matrix *ref = matrix_new(R->x, R->i, Ap, R->m, R->n, nnz);
    mu_assert("reference", ref != NULL && ref->nnz < nnz);

    const int n_workers[] = {0, 1, 7};
    for (int k = 0; k < 3; ++k)
    {
        ThreadPool *pool = thread_pool_new(n_workers[k]);
        mu_assert("pool", pool != NULL);
        Matrix *A = matrix_new_parallel(R->x, R->i, Ap, R->m, R->n, nnz, pool);
        mu_assert("parallel copy", A != NULL);
        mu_assert("dimensions", A->m == ref->m && A->n == ref->n &&
                                    A->nnz == ref->nnz &&
                                    A->n_alloc == ref->n_alloc);
        mu_assert("row pointers",
                  memcmp(A->p, ref->p, (ref->m + 1) * sizeof(RowRange)) == 0);
        for (size_t i = 0; i < ref->m; ++i)
        {
            int start = ref->p[i].start;
            size_t len = (size_t) (ref->p[i].end - start);
            mu_assert("column indices",
                      memcmp(A->i + start, ref->i + start, len * sizeof(int)) == 0);
            mu_assert("values", memcmp(A->x + start, ref->x + start,
                                       len * sizeof(double)) == 0);
        }
        free_matrix(A);
        thread_pool_free(pool);
    }

    free_matrix(ref);
    free_matrix(R);
    PS_FREE(Ap);
    return 0;
}

static const char *all_tests_matrix()
{
    mu_run_test(test_0_matrix, counter_matrix);
    mu_run_test(test_1_matrix, counter_matrix);
    mu_run_test(test_2_matrix, counter_matrix);
    mu_run_test(test_20_matrix, counter_matrix);
    mu_run_test(test_diagonal_matrix, counter_matrix);
    mu_run_test(test_transpose_parallel, counter_matrix);
    mu_run_test(test_matrix_new_parallel, counter_matrix);
    return 0;
}

int test_matrix()
{
    const char *result = all_tests_matrix();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("Matrix: TEST FAILED!\n");
    }
    else
    {
        printf("Matrix: ALL TESTS PASSED\n");
    }
    printf("Matrix: Tests run: %d\n", counter_matrix);
    return result == 0;
}

#endif // TEST_MATRIX_H
