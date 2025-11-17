#ifndef TEST_MATRIX_H
#define TEST_MATRIX_H

#include "Debugger.h"
#include "Matrix.h"
#include "minunit.h"
#include "test_macros.h"
#include <stdio.h>

int counter_matrix = 0;

/* Summary of tests:
Test 0: matrix_new: allocation and free of a matrix.
Test 1: transpose: transpose a matrix.
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
    Matrix *AT = transpose(A, work_n_cols);

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
    Matrix *AT = transpose(A, work_n_cols);

    // remove extra space to simplify test
    int row_sizes_AT[6] = {2, 3, 3, 1, 2};
    int col_sizes_AT[6] = {3, 3, 1, 2, 1, 1};
    int col_idxs_map_AT[6];
    remove_extra_space(AT, row_sizes_AT, col_sizes_AT, true, col_idxs_map_AT);

    // correct answer
    double AT_vals_correct[] = {1, 3, -1, 1, 1, 2, 1, 1, 1, 1, 1};
    int AT_cols_correct[] = {0, 3, 0, 1, 3, 0, 1, 4, 5, 1, 2};
    int AT_row_starts_correct[] = {0, 2, 5, 8, 9, 11};

    mu_assert("error, vals not equal",
              ARRAYS_EQUAL_DOUBLE(AT_vals_correct, AT->x, 11));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(AT_cols_correct, AT->i, 11));
    CHECK_ROW_STARTS(AT, AT_row_starts_correct);

    free_matrix(A);
    free_matrix(AT);

    return 0;
}

/*
A = [1 2  3  0      (2 extra)
     4 0  5  6      (2 extra)
     7 0  8  0      (2 extra)
     9 0 10  0]     (2 extra)

    Failure because we want the first row to have 3 extra spaces but
    max_shift = 2. Then we must shift 4, 5, 6 right with one step.
    But this should be rejected.
*/
static char *test_2_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 2};
    int row_starts[] = {0, 3, 6, 8, 10};
    int nnz = 10;
    int n_rows = 4;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    bool success = shift_row(A, 0, 3, 2);
    mu_assert("error, shift_row did not reject shift", !success);
    free_matrix(A);

    return 0;
}

/*
A = [1 2  3  0      (2 extra)
     4 0  5  6      (2 extra)
     7 0  8  0      (2 extra)
     9 0 10  0]     (2 extra)

     Previous test, but with shift allowed
*/
static char *test_3_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 2};
    int row_starts[] = {0, 3, 6, 8, 10};
    int nnz = 10;
    int n_rows = 4;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    bool success = shift_row(A, 0, 3, 3);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3, 0, 0, 4, 4, 5, 6, 0, 7, 8, 0, 0, 9, 10, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    free_matrix(A);

    return 0;
}

/*
A = [1 2  3  0      (2 extra)
     4 0  5  6      (2 extra)
     7 0  8  0      (2 extra)
     9 0 10  0]     (2 extra)

     Success with maxShift = 5.
*/
static char *test_4_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 2};
    int row_starts[] = {0, 3, 6, 8, 10};
    int nnz = 10;
    int n_rows = 4;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    bool success = shift_row(A, 0, 5, 5);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3, 0, 0, 4, 5, 6, 4, 5, 6, 7, 8, 0, 9, 10, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));

    int row_starts_correct[] = {0, 8, 11, 14, 18};
    int row_ends_correct[] = {3, 11, 13, 16, 18};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
A = [1 2  3  0      (2 extra)
     4 0  5  6      (2 extra)
     7 0  8  0      (2 extra)
     9 0 10  0]     (2 extra)

     Failure when maxShift = 4
*/
static char *test_5_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 2};
    int row_starts[] = {0, 3, 6, 8, 10};
    int nnz = 10;
    int n_rows = 4;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    bool success = shift_row(A, 0, 5, 4);
    mu_assert("error, shift_row did not reject shift", !success);
    free_matrix(A);

    return 0;
}

/*
A = [1 2  3  0      (2 extra)
     4 0  5  6      (2 extra)
     7 0  8  0      (2 extra)
     9 0 10  0]     (2 extra)

     Success when we want the last row to have 5 extra spaces.
*/
static char *test_6_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 2};
    int row_starts[] = {0, 3, 6, 8, 10};
    int nnz = 10;
    int n_rows = 4;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    bool success = shift_row(A, 3, 5, 10);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3, 0, 0, 4, 5, 6, 0, 7, 8, 9, 10, 0, 9, 10, 0, 0};
    int cols_correct[] = {0, 1, 2, 0, 0, 0, 2, 3, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, 18));

    int row_starts_correct[] = {0, 5, 9, 11, 18};
    int row_ends_correct[] = {3, 8, 11, 13, 18};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
A = [1  2   3  0      (2 extra)
     4  0   5  6      (2 extra)
     7  0   8  0      (2 extra)
     9  0   0  0      (2 extra)
     10 0   0  0]     (2 extra)

     Row 2 five extra spaces.
*/
static char *test_7_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 0};
    int row_starts[] = {0, 3, 6, 8, 9, 10};
    int nnz = 10;
    int n_rows = 5;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    bool success = shift_row(A, 2, 5, 10);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3, 0, 0, 4, 5, 6, 0,  0,
                             7, 8, 0, 0, 9, 0, 0, 9, 10, 0};
    int cols_correct[] = {0, 1, 2, 0, 0, 0, 2, 3, 0, 0,
                          0, 2, 0, 0, 0, 0, 0, 0, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, 18));

    int row_starts_correct[] = {0, 5, 10, 17, 18, 20};
    int row_ends_correct[] = {3, 8, 12, 18, 19, 20};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
A = [1  2   3  0      (2 extra)
     4  0   5  6      (2 extra)
     7  0   8  0      (2 extra)
     9  0   0  0      (2 extra)
     10 0   0  0]     (2 extra)

     Row 2 seven extra spaces.
*/
static char *test_8_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 0};
    int row_starts[] = {0, 3, 6, 8, 9, 10};
    int nnz = 10;
    int n_rows = 5;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    bool success = shift_row(A, 2, 7, 10);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3, 0, 0, 4, 5, 6,  0, 7,
                             8, 8, 0, 0, 9, 0, 0, 10, 9, 10};
    int cols_correct[] = {0, 1, 2, 0, 0, 0, 2, 3, 0, 0,
                          2, 2, 0, 0, 0, 0, 0, 0, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, 18));

    int row_starts_correct[] = {0, 5, 9, 18, 19, 20};
    int row_ends_correct[] = {3, 8, 11, 19, 20, 20};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
A = [1  2   3  0      (2 extra)
     4  0   5  6      (2 extra)
     7  0   8  0      (2 extra)
     9  0   0  0      (2 extra)
     10 0   0  0]     (2 extra)

     Row 2 seven extra spaces when row 3 is made empty.
*/
static char *test_9_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 0};
    int row_starts[] = {0, 3, 6, 8, 9, 10};
    int nnz = 10;
    int n_rows = 5;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    A->p[3].end = A->p[3].start;

    bool success = shift_row(A, 2, 7, 10);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3, 0, 0, 4, 5, 6,  0, 0,
                             7, 8, 0, 0, 9, 0, 0, 10, 0,

                             10};
    int cols_correct[] = {0, 1, 2, 0, 0, 0, 2, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,

                          0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, 18));

    int row_starts_correct[] = {0, 5, 10, 19, 19, 20};
    int row_ends_correct[] = {3, 8, 12, 19, 20, 20};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
 A = [1    2    3     0     (2 extra)
      4    0    5     6     (2 extra)
      7    0    8     0     (2 extra)
      9    0    0     0     (2 extra)
      10  11    0     0     (2 extra)
      12  13    14    0     (2 extra)
      0   15    0    16]    (2 extra)

      Several empty rows in right direction
*/
static char *test_10_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 0, 1, 0, 1, 2, 1, 3};
    int row_starts[] = {0, 3, 6, 8, 9, 11, 14, 16};
    int nnz = 16;
    int n_rows = 7;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    A->p[3].end = A->p[3].start;
    A->p[4].end = A->p[4].start;
    A->p[5].end = A->p[5].start;

    bool success = shift_row(A, 2, 11, 10);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3,  0,  0, 4, 5,  6,  0,  0, 7, 8,  0,  0, 9,
                             0, 0, 10, 11, 0, 0, 12, 13, 14, 0, 0, 15, 16, 0, 0};
    int cols_correct[] = {0, 1, 2, 0, 0, 0, 2, 3, 0, 0, 0, 2, 0, 0, 0,
                          0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 3, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, 18));

    int row_starts_correct[] = {0, 5, 10, 23, 23, 23, 26, 30};
    int row_ends_correct[] = {3, 8, 12, 23, 23, 23, 28, 30};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
 A = [1    2    3     0     (2 extra)
      4    0    5     6     (2 extra)
      7    0    8     0     (2 extra)
      9    0    10    0     (2 extra)
      11  12    0     0     (2 extra)
      0   13    14    0]    (2 extra)

      Several empty rows in both directions
*/
static char *test_11_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 2, 0, 1, 1, 2};
    int row_starts[] = {0, 3, 6, 8, 10, 12, 14};
    int nnz = 14;
    int n_rows = 6;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    A->p[1].end = A->p[1].start;
    A->p[2].end = A->p[2].start;
    A->p[4].end = A->p[4].start;
    A->p[5].end = A->p[5].start;

    bool success = shift_row(A, 3, 14, 10);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3,  0, 0, 4,  5,  6, 0, 0,  9,  10, 0,
                             0, 9, 10, 0, 0, 11, 12, 0, 0, 13, 14, 0,  0};
    int cols_correct[] = {0, 1, 2, 0, 0, 0, 2, 3, 0, 0, 0, 2, 0,
                          0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 2, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, 18));

    int row_starts_correct[] = {0, 5, 10, 10, 26, 26, 26};
    int row_ends_correct[] = {3, 5, 10, 12, 26, 26, 26};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
A  = [1  2   3   0      empty
      4  0   0   0      (2 extra)
      7  0   8   0      (2 extra)
      9  0   5   6      (2 extra)
      10 0  11   0]     (2 extra)

     Row 2 six extra spaces when first row is made empty.
*/
static char *test_12_matrix()
{
    double vals[] = {1, 2, 3, 4, 7, 8, 9, 5, 6, 10, 11};
    int cols[] = {0, 1, 2, 0, 0, 2, 0, 2, 3, 0, 2};
    int row_starts[] = {0, 3, 4, 6, 9, 11};
    int nnz = 11;
    int n_rows = 5;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    A->p[0].end = A->p[0].start;

    bool success = shift_row(A, 2, 6, 10);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3, 4, 7, 8, 0,  0,  7, 8, 0,
                             0, 9, 5, 6, 0, 0, 10, 11, 0, 0};
    int cols_correct[] = {0, 1, 2, 0, 0, 2, 0, 0, 0, 2, 0,
                          0, 0, 2, 3, 0, 0, 0, 2, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, 18));

    int row_starts_correct[] = {0, 3, 4, 12, 17, 21};
    int row_ends_correct[] = {0, 4, 6, 15, 19, 21};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
A = [1  2   3  0      (2 extra)
     4  0   5  6      (2 extra)
     7  0   8  0      (2 extra)
     9  0   0  0      (2 extra)
     10 0   0  0]     (2 extra)

     Row 2 six extra spaces when last is made empty.
*/
static char *test_13_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int cols[] = {0, 1, 2, 0, 2, 3, 0, 2, 0, 0};
    int row_starts[] = {0, 3, 6, 8, 9, 10};
    int nnz = 10;
    int n_rows = 5;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    A->p[4].end = A->p[4].start;

    bool success = shift_row(A, 2, 6, 10);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1, 2, 3, 0, 0, 4, 5,  6, 0, 0, 7,
                             8, 0, 0, 9, 0, 0, 10, 9, 0, 0};
    int cols_correct[] = {0, 1, 2, 0, 0, 0, 2, 3, 0, 0, 0,
                          2, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));
    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, 18));

    int row_starts_correct[] = {0, 5, 10, 18, 19, 20};
    int row_ends_correct[] = {3, 8, 12, 19, 19, 20};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);

    return 0;
}

/*
A = [1    2    3    0     0      (0 extra)
     1    2    3    0     0      (0 extra)
     6                           (2 extra)
     9        10                 (2 extra) empty
         11   12   12            (2 extra)
     13            14            (2 extra) empty
     15                  16      (2 extra)
                         17      (2 extra) empty
     18 19    20                 (2 extra) empty
     21       22             ]   (2 extra)


       Row 2 21 extra spaces.
*/
static char *test_14_matrix()
{
    double vals[] = {1,  2,  3,  1,  2,  3,  6,  7,  8,  9,  10, 11,
                     12, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
    int cols[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 2, 1,
                  2, 3, 0, 3, 0, 4, 4, 0, 1, 2, 0, 2};
    int row_starts[] = {0, 3, 6, 7, 9, 11, 14, 16, 18, 19, 22, 24};
    int nnz = 24;
    int n_rows = 11;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    A->p[0].end = 5;
    A->p[1].end = 10;
    A->p[4].end = A->p[4].start;
    A->p[6].end = A->p[6].start;
    A->p[8].end = A->p[8].start;
    A->p[9].end = A->p[9].start;

    bool success = shift_row(A, 2, 21, 20);
    mu_assert("error, shift_row did reject shift", success);

    // correct answer
    double vals_correct[] = {1,  2,  3,  0,  0,  1, 2,  3,  0, 0,  6,  0,
                             0,  7,  8,  0,  0,  9, 10, 0,  0, 11, 12, 12,
                             0,  0,  13, 14, 0,  0, 15, 16, 7, 8,

                             11, 12, 12,

                             15, 16,

                             20, 0,  0,  21, 22, 0, 0};

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, 18));

    int row_starts_correct[] = {0, 5, 10, 32, 34, 34, 37, 37, 39, 39, 42, 46};
    int row_ends_correct[] = {5, 10, 11, 34, 34, 37, 37, 39, 39, 39, 44, 46};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);
    printf("%d", EXTRA_ROW_SPACE);

    return 0;
}

//
// A = [9 8 0 0
//      3 0 5 0
//      1 0 0 0
//      0 0 1 1]
static char *test_15_matrix()
{
    double vals[] = {9, 8, 3, 5, 1, 1, 1};
    int cols[] = {0, 1, 0, 2, 0, 2, 3};
    int row_starts[] = {0, 2, 4, 5, 7};
    int nnz = 7;
    int n_rows = 4;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);
    int not_used = 0;

    double old_val;
    old_val = insert_or_update_coeff(A, 1, 0, 2,
                                     &not_used); // first in a row, in-place
    mu_assert("error old_val", old_val == 3);
    old_val = insert_or_update_coeff(A, 3, 1, 4,
                                     &not_used); // first in a row, new space
    mu_assert("error old_val", old_val == 0);
    old_val =
        insert_or_update_coeff(A, 1, 2, 9, &not_used); // end of a row, in-place
    mu_assert("error old_val", old_val == 5);
    old_val = insert_or_update_coeff(A, 0, 2, 1,
                                     &not_used); // end of a row, new space
    mu_assert("error old_val", old_val == 0);
    old_val =
        insert_or_update_coeff(A, 2, 0, 7, &not_used); // end of a row, in-place
    mu_assert("error old_val", old_val == 1);
    old_val = insert_or_update_coeff(A, 2, 2, 3,
                                     &not_used); // end of a row, new space
    mu_assert("error old_val", old_val == 0);

    int row_sizes[] = {3, 2, 2, 3};
    int col_sizes[] = {4, 1, 4, 1};
    int map[4];
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    // correct answer
    double vals_correct[] = {9, 8, 1, 2, 9, 7, 3, 4, 1, 1};
    int cols_correct[] = {0, 1, 2, 0, 2, 0, 2, 1, 2, 3};
    mu_assert("error, vals not equal",
              ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, A->nnz));

    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, A->nnz));

    int row_starts_correct[] = {0, 3, 5, 7, 10};
    int row_ends_correct[] = {3, 5, 7, 10, 10};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);
    return 0;
}

//
// A = [9 8 4 0
//      3 0 5 0
//      1 0 0 0
//      0 0 1 1]
static char *test_16_matrix()
{
    double vals[] = {9, 8, 4, 3, 5, 1, 1, 1};
    int cols[] = {0, 1, 2, 0, 2, 0, 2, 3};
    int row_starts[] = {0, 3, 5, 6, 8};
    int nnz = 8;
    int n_rows = 4;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);

    int not_used = 0;
    double old_val;
    old_val = insert_or_update_coeff(A, 1, 1, 7,
                                     &not_used); // middle of a row, new space
    mu_assert("error old_val", old_val == 0);
    old_val = insert_or_update_coeff(A, 0, 1, 5,
                                     &not_used); // middle of a row, existing
    mu_assert("error old_val", old_val == 8);

    int row_sizes[] = {3, 3, 1, 2};
    int col_sizes[] = {4, 2, 3, 1};
    int map[4];
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    // correct answer
    double vals_correct[] = {9, 5, 4, 3, 7, 5, 1, 1, 1};
    int cols_correct[] = {0, 1, 2, 0, 1, 2, 0, 2, 3};
    mu_assert("error, vals not equal",
              ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, A->nnz));

    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, A->nnz));

    int row_starts_correct[] = {0, 3, 6, 7, 9};
    int row_ends_correct[] = {3, 6, 7, 9, 9};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);
    return 0;
}

/* insert a zero for a variable that exists
 A = [1, 0, 2, 3,
      4, 5, 6, 0,
      3, 0, 0, 1]
*/
static char *test_17_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 3, 1};
    int cols[] = {0, 2, 3, 0, 1, 2, 0, 3};
    int row_starts[] = {0, 3, 6, 8};
    int nnz = 8;
    int n_rows = 3;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);

    int row_size = 3;
    double old_val;
    old_val = insert_or_update_coeff(A, 1, 2, 0,
                                     &row_size); // middle of a row, new space

    mu_assert("error row_size", row_size == 2);

    int row_sizes[] = {3, 2, 2};
    int col_sizes[] = {3, 1, 1, 2};
    int map[4];
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    double vals_correct[] = {1, 2, 3, 4, 5, 3, 1};
    int cols_correct[] = {0, 2, 3, 0, 1, 0, 3};
    mu_assert("error, vals not equal",
              ARRAYS_EQUAL_DOUBLE(vals_correct, A->x, A->nnz));

    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols_correct, A->i, A->nnz));

    int row_starts_correct[] = {0, 3, 5, 7};
    int row_ends_correct[] = {3, 5, 7, 7};

    CHECK_ROW_STARTS(A, row_starts_correct);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);
    return 0;
}

/* insert a zero for a variable that doesn't exist
 A = [1, 0, 2, 3,
      4, 5, 6, 0,
      3, 0, 0, 1]

    OBS: This test should cause an assertion inside insert_or_update_coeff,
    because in the code we only expect to call 'insert_or_update_coeff'
    with val = 0 when the variable already exists.
*/
static char *test_18_matrix()
{
    double vals[] = {1, 2, 3, 4, 5, 6, 3, 1};
    int cols[] = {0, 2, 3, 0, 1, 2, 0, 3};
    int row_starts[] = {0, 3, 6, 8};
    int nnz = 8;
    int n_rows = 3;
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, 4, nnz);

    int row_size = 2;
    double old_val;
    old_val = insert_or_update_coeff(A, 2, 1, 0,
                                     &row_size); // middle of a row, new space

    mu_assert("error row_size", row_size == 2);

    int row_sizes[] = {3, 3, 2};
    int col_sizes[] = {3, 1, 2, 2};
    int map[4];
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    mu_assert("error, vals not equal", ARRAYS_EQUAL_DOUBLE(vals, A->x, A->nnz));

    mu_assert("error, cols not equal", ARRAYS_EQUAL_INT(cols, A->i, A->nnz));

    int row_ends_correct[] = {3, 6, 8, 8};

    CHECK_ROW_STARTS(A, row_starts);
    CHECK_ROW_ENDS(A, row_ends_correct);

    free_matrix(A);
    return 0;
}

static const char *all_tests_matrix()
{
    mu_run_test(test_0_matrix, counter_matrix);
    mu_run_test(test_1_matrix, counter_matrix);
    mu_run_test(test_2_matrix, counter_matrix);
    mu_run_test(test_3_matrix, counter_matrix);
    mu_run_test(test_4_matrix, counter_matrix);
    mu_run_test(test_5_matrix, counter_matrix);
    mu_run_test(test_6_matrix, counter_matrix);
    mu_run_test(test_7_matrix, counter_matrix);
    mu_run_test(test_8_matrix, counter_matrix);
    mu_run_test(test_9_matrix, counter_matrix);
    mu_run_test(test_10_matrix, counter_matrix);
    mu_run_test(test_11_matrix, counter_matrix);
    mu_run_test(test_12_matrix, counter_matrix);
    mu_run_test(test_13_matrix, counter_matrix);
    mu_run_test(test_14_matrix, counter_matrix);
    mu_run_test(test_15_matrix, counter_matrix);
    mu_run_test(test_16_matrix, counter_matrix);
    mu_run_test(test_17_matrix, counter_matrix);
    // mu_run_test(test_18_matrix, counter_matrix); // we don't run this
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
