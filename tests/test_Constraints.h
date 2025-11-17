#ifndef TEST_CONSTRAINTS_H
#define TEST_CONSTRAINTS_H

#include "Activity.h"
#include "Constraints.h"
#include "Locks.h"
#include "PSLP_API.h"
#include "State.h"
#include "glbopts.h"
#include "minunit.h"
#include "test_macros.h"
#include <stdio.h>

// Test 1 - 2: delete_inactive_rows
// Test 3 - 4: delete_inactive_cols_from_A_and_AT

static int counter_constraints = 0;

static char *test_1_constraints()
{

    // build constraints object
    double vals[] = {1, -1, 1, 2, -1, 1, 1, 1, 1};
    int cols[] = {0, 1, 2, 0, 3, 0, 1, 2, 3};
    int row_starts[] = {0, 3, 5, 9};
    int n_rows = 3;
    int n_cols = 4;
    int nnz = 9;
    int work_n_cols[n_cols];
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, n_cols, nnz);
    Matrix *AT = transpose(A, work_n_cols);
    double lhs[3] = {4, 2, -INF};
    double rhs[3] = {4, 2, 1};
    Bound bounds[n_cols];
    INIT_BOUNDS(bounds, -10, 10, n_cols);

    Lock locks[] = {{.up = 3, .down = 2},
                    {.up = 2, .down = 1},
                    {.up = 2, .down = 1},
                    {.up = 2, .down = 1}};

    int row_sizes[] = {3, 2, 4};
    int col_sizes[] = {3, 2, 2, 2};
    Activity activities[3] = {0};
    RowTag *row_tags = new_rowtags(lhs, rhs, n_rows);
    ColTag col_tags[4] = {0};
    Settings *stgs = default_settings();
    State *data = new_state(row_sizes, col_sizes, locks, n_rows, n_cols, activities,
                            NULL, row_tags);
    Constraints *constraints =
        constraints_new(A, AT, lhs, rhs, bounds, data, row_tags, col_tags);

    // remove row 0
    iVec_append(data->rows_to_delete, 0);
    delete_inactive_rows(constraints);

    // test for correctness
    int correct_col_sizes[] = {2, 1, 1, 2};
    int correct_row_sizes[] = {SIZE_INACTIVE_ROW, 2, 4};
    mu_assert("error col_sizes", ARRAYS_EQUAL_INT(col_sizes, correct_col_sizes, 4));
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes, correct_row_sizes, 3));

    int work_map[4] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, work_map);
    remove_extra_space(AT, col_sizes, row_sizes, true, work_map);

    double A_vals_correct[] = {2, -1, 1, 1, 1, 1};
    int A_cols_correct[] = {0, 3, 0, 1, 2, 3};
    int A_row_starts_correct[] = {0, 2, 6};
    mu_assert("error A", ARRAYS_EQUAL_DOUBLE(A->x, A_vals_correct, 6));
    mu_assert("error A", ARRAYS_EQUAL_INT(A->i, A_cols_correct, 6));
    CHECK_ROW_STARTS(A, A_row_starts_correct);

    double AT_vals_correct[] = {2, 1, 1, 1, -1, 1};
    int AT_cols_correct[] = {0, 1, 1, 1, 0, 1};
    int AT_row_starts_correct[] = {0, 2, 3, 4, 6};
    mu_assert("error AT_vals", ARRAYS_EQUAL_DOUBLE(AT->x, AT_vals_correct, 6));
    mu_assert("error AT_cols", ARRAYS_EQUAL_INT(AT->i, AT_cols_correct, 6));
    CHECK_ROW_STARTS(AT, AT_row_starts_correct);

    int correct_up_locks[] = {2, 1, 1, 2};
    int correct_down_locks[] = {1, 0, 0, 1};

    CHECK_LOCKS(locks, correct_up_locks, correct_down_locks, 4);

    // deallocate memory
    free_matrix(A);
    free_matrix(AT);
    PS_FREE(constraints);
    PS_FREE(stgs);
    free_state(data);
    PS_FREE(row_tags);

    return 0;
}

static char *test_2_constraints()
{

    // build constraints object
    double vals[] = {1.3,   -0.6,  0.26,  -0.42, -0.59, -0.9,  0.53,  1.42,
                     1.31,  0.47,  -1.47, 0.63,  -0.93, 1.54,  -1.07, -1.46,
                     -2.37, 0.25,  -1.45, 0.71,  -1.44, -1.57, -0.44, -0.37,
                     -0.5,  0.61,  1.38,  0.17,  0.77,  0.41,  -0.21, -0.73,
                     0.55,  -0.18, 1.14,  0.57,  -0.22, -1.4,  1.02,  0.74};
    int cols[] = {0, 3, 4, 9, 1, 3, 6, 9, 2, 3, 5, 6, 8, 1, 9, 1, 2, 3, 5, 0,
                  2, 4, 5, 7, 4, 6, 9, 1, 4, 6, 7, 8, 1, 3, 5, 7, 1, 2, 4, 8};
    int row_starts[] = {0, 4, 8, 13, 15, 19, 24, 27, 32, 36, 40};
    int n_rows = 10;
    int n_cols = 10;
    int nnz = 40;
    int work_n_cols[n_cols];
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, n_cols, nnz);
    Matrix *AT = transpose(A, work_n_cols);
    double lhs[10] = {0.0};
    double rhs[10] = {0.0};
    Bound bounds[n_cols];
    INIT_BOUNDS(bounds, -10, 10, n_cols);

    Lock locks[n_cols];
    int row_sizes[] = {4, 4, 5, 2, 4, 5, 3, 5, 4, 4};
    int col_sizes[] = {2, 6, 4, 5, 5, 4, 4, 3, 3, 4};
    RowTag *row_tags = new_rowtags(lhs, rhs, n_rows);
    ColTag col_tags[10] = {0};
    Settings *stgs = default_settings();
    Activity activities[10] = {0};
    State *data = new_state(row_sizes, col_sizes, locks, n_rows, n_cols, activities,
                            NULL, row_tags);
    Constraints *constraints =
        constraints_new(A, AT, lhs, rhs, bounds, data, row_tags, col_tags);

    // remove rows 1, 3, 5
    iVec_append_array(data->rows_to_delete, (int[]){1, 3, 5}, 3);
    delete_inactive_rows(constraints);

    // test for correctness
    int correct_col_sizes[] = {1, 4, 3, 4, 4, 3, 3, 2, 3, 2};
    int correct_row_sizes[] = {
        4, SIZE_INACTIVE_ROW, 5, SIZE_INACTIVE_ROW, 4, SIZE_INACTIVE_ROW, 3, 5, 4,
        4};
    mu_assert("error col_sizes", ARRAYS_EQUAL_INT(col_sizes, correct_col_sizes, 10));
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes, correct_row_sizes, 10));

    int work_map[10] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, work_map);
    remove_extra_space(AT, col_sizes, row_sizes, true, work_map);

    double A_vals_correct[] = {1.3,   -0.6,  0.26,  -0.42, 1.31,  0.47, -1.47, 0.63,
                               -0.93, -1.46, -2.37, 0.25,  -1.45, -0.5, 0.61,  1.38,
                               0.17,  0.77,  0.41,  -0.21, -0.73, 0.55, -0.18, 1.14,
                               0.57,  -0.22, -1.4,  1.02,  0.74};
    int A_cols_correct[] = {0, 3, 4, 9, 2, 3, 5, 6, 8, 1, 2, 3, 5, 4, 6,
                            9, 1, 4, 6, 7, 8, 1, 3, 5, 7, 1, 2, 4, 8};
    int A_row_starts_correct[] = {0, 4, 9, 13, 16, 21, 25, 29};

    mu_assert("error A_vals", ARRAYS_EQUAL_DOUBLE(A->x, A_vals_correct, 29));
    mu_assert("error A_cols", ARRAYS_EQUAL_INT(A->i, A_cols_correct, 29));
    CHECK_ROW_STARTS(A, A_row_starts_correct);

    double AT_vals_correct[] = {1.3,   -1.46, 0.17, 0.55,  -0.22, 1.31, -2.37, -1.4,
                                -0.6,  0.47,  0.25, -0.18, 0.26,  -0.5, 0.77,  1.02,
                                -1.47, -1.45, 1.14, 0.63,  0.61,  0.41, -0.21, 0.57,
                                -0.93, -0.73, 0.74, -0.42, 1.38};
    int AT_cols_correct[] = {0, 2, 4, 5, 6, 1, 2, 6, 0, 1, 2, 5, 0, 3, 4,
                             6, 1, 2, 5, 1, 3, 4, 4, 5, 1, 4, 6, 0, 3};
    int AT_row_starts_correct[] = {0, 1, 5, 8, 12, 16, 19, 22, 24, 27, 29};
    mu_assert("error AT_vals", ARRAYS_EQUAL_DOUBLE(AT->x, AT_vals_correct, 27));
    mu_assert("error AT_cols", ARRAYS_EQUAL_INT(AT->i, AT_cols_correct, 27));
    CHECK_ROW_STARTS(AT, AT_row_starts_correct);

    // deallocate memory
    free_matrix(A);
    free_matrix(AT);
    PS_FREE(stgs);
    PS_FREE(constraints);
    free_state(data);
    PS_FREE(row_tags);

    return 0;
}

static char *test_3_constraints()
{

    // build constraints object
    double vals[] = {1, -1, 1, 2, -1, 1, 1, 1, 1};
    int cols[] = {0, 1, 2, 0, 3, 0, 1, 2, 3};
    int row_starts[] = {0, 3, 5, 9};
    int n_rows = 3;
    int n_cols = 4;
    int nnz = 9;
    int work_n_cols[n_cols];
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, n_cols, nnz);
    Matrix *AT = transpose(A, work_n_cols);
    double lhs[3] = {4, 2, -INF};
    double rhs[3] = {4, 2, 1};
    Bound bounds[n_cols];
    INIT_BOUNDS(bounds, -10, 10, n_cols);

    Lock locks[4];

    int row_sizes[] = {3, 2, 4};
    int col_sizes[] = {3, 2, 2, 2};
    RowTag *row_tags = new_rowtags(lhs, rhs, n_rows);
    Activity activities[3] = {0};
    ColTag col_tags[4] = {0};
    Settings *stgs = default_settings();
    State *data = new_state(row_sizes, col_sizes, locks, n_rows, n_cols, activities,
                            NULL, row_tags);
    Constraints *constraints =
        constraints_new(A, AT, lhs, rhs, bounds, data, row_tags, col_tags);

    // fix variable 1
    iVec_append(data->fixed_cols_to_delete, 1);
    delete_inactive_cols_from_A_and_AT(constraints);

    // test for correctness
    int correct_col_sizes[] = {3, SIZE_INACTIVE_COL, 2, 2};
    int correct_row_sizes[] = {2, 2, 3};
    mu_assert("error col_sizes", ARRAYS_EQUAL_INT(col_sizes, correct_col_sizes, 4));
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes, correct_row_sizes, 3));

    int work_map[4] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, work_map);
    remove_extra_space(AT, col_sizes, row_sizes, true, work_map);

    double A_vals_correct[] = {1, 1, 2, -1, 1, 1, 1};
    int A_cols_correct[] = {0, 1, 0, 2, 0, 1, 2};
    int A_row_starts_correct[] = {0, 2, 4, 7};
    mu_assert("error A", ARRAYS_EQUAL_DOUBLE(A->x, A_vals_correct, 7));
    mu_assert("error A", ARRAYS_EQUAL_INT(A->i, A_cols_correct, 7));
    CHECK_ROW_STARTS(A, A_row_starts_correct);

    double AT_vals_correct[] = {1, 2, 1, 1, 1, -1, 1};
    int AT_cols_correct[] = {0, 1, 2, 0, 2, 1, 2};
    int AT_row_starts_correct[] = {0, 3, 5, 7};
    mu_assert("error AT_vals", ARRAYS_EQUAL_DOUBLE(AT->x, AT_vals_correct, 7));
    mu_assert("error AT_cols", ARRAYS_EQUAL_INT(AT->i, AT_cols_correct, 7));
    CHECK_ROW_STARTS(AT, AT_row_starts_correct);

    // deallocate memory
    free_matrix(A);
    free_matrix(AT);
    PS_FREE(stgs);
    PS_FREE(constraints);
    free_state(data);
    PS_FREE(row_tags);

    return 0;
}

static char *test_4_constraints()
{

    // build constraints object
    double vals[] = {1.12, -0.29, -1.41, 0.67,  -0.35, -1.15, -0.72, -0.03,
                     0.82, 0.32,  1.38,  2.23,  -0.44, 0.05,  1.17,  0.21,
                     0.68, 1.38,  -0.56, -1.03, -0.33, -0.1,  0.9,   2.74,
                     0.76, 1.17,  0.04,  -2.16, 0.65,  -0.02, -0.29, -0.07,
                     -0.9, 0.59,  0.45,  1.47,  -1.75, 0.57,  0.27,  1.09};
    int cols[] = {0, 4, 5, 1, 9, 3, 4, 5, 6, 7, 8, 3, 4, 6, 9, 9, 7, 0, 2, 3,
                  4, 6, 7, 2, 3, 5, 6, 8, 9, 1, 2, 3, 5, 6, 8, 9, 2, 6, 7, 9};
    int row_starts[] = {0, 3, 5, 11, 15, 16, 17, 23, 29, 36, 40};
    int n_rows = 10;
    int n_cols = 10;
    int nnz = 40;
    int work_n_cols[n_cols];
    Matrix *A = matrix_new(vals, cols, row_starts, n_rows, n_cols, nnz);
    Matrix *AT = transpose(A, work_n_cols);
    double lhs[10] = {0.0};
    double rhs[10] = {0.0};
    Bound bounds[n_cols];
    INIT_BOUNDS(bounds, -10, 10, n_cols);

    Lock locks[n_cols];

    int row_sizes[] = {3, 2, 6, 4, 1, 1, 6, 6, 7, 4};
    int col_sizes[] = {2, 2, 4, 5, 4, 4, 6, 4, 3, 6};
    RowTag *row_tags = new_rowtags(lhs, rhs, n_rows);
    Activity activities[10] = {0};
    ColTag col_tags[10] = {0};
    Settings *stgs = default_settings();
    State *data = new_state(row_sizes, col_sizes, locks, n_rows, n_cols, activities,
                            NULL, row_tags);
    Constraints *constraints =
        constraints_new(A, AT, lhs, rhs, bounds, data, row_tags, col_tags);

    // remove cols 1, 3, 5
    iVec_append_array(data->fixed_cols_to_delete, (int[]){1, 3, 5}, 3);
    delete_inactive_cols_from_A_and_AT(constraints);

    // test for correctness
    int correct_col_sizes[] = {
        2, SIZE_INACTIVE_COL, 4, SIZE_INACTIVE_COL, 4, SIZE_INACTIVE_COL, 6, 4, 3,
        6};
    int correct_row_sizes[] = {2, 1, 4, 3, 1, 1, 5, 4, 4, 4};
    mu_assert("error col_sizes", ARRAYS_EQUAL_INT(col_sizes, correct_col_sizes, 10));
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes, correct_row_sizes, 10));

    int work_map[10] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, work_map);
    remove_extra_space(AT, col_sizes, row_sizes, true, work_map);

    double A_vals_correct[] = {1.12, -0.29, -0.35, -0.72, 0.82, 0.32,  1.38,  -0.44,
                               0.05, 1.17,  0.21,  0.68,  1.38, -0.56, -0.33, -0.1,
                               0.9,  2.74,  0.04,  -2.16, 0.65, -0.29, 0.59,  0.45,
                               1.47, -1.75, 0.57,  0.27,  1.09};
    int A_cols_correct[] = {0, 2, 6, 2, 3, 4, 5, 2, 3, 6, 6, 4, 0, 1, 2,
                            3, 4, 1, 3, 5, 6, 1, 3, 5, 6, 1, 3, 4, 6};
    int A_row_starts_correct[] = {0, 2, 3, 7, 10, 11, 12, 17, 21, 25, 29};

    mu_assert("error A_vals", ARRAYS_EQUAL_DOUBLE(A->x, A_vals_correct, 29));
    mu_assert("error A_cols", ARRAYS_EQUAL_INT(A->i, A_cols_correct, 29));
    CHECK_ROW_STARTS(A, A_row_starts_correct);

    double AT_vals_correct[] = {
        1.12, 1.38,  -0.56, 2.74,  -0.29, -1.75, -0.29, -0.72, -0.44, -0.33,
        0.82, 0.05,  -0.1,  0.04,  0.59,  0.57,  0.32,  0.68,  0.9,   0.27,
        1.38, -2.16, 0.45,  -0.35, 1.17,  0.21,  0.65,  1.47,  1.09};
    int AT_cols_correct[] = {0, 6, 6, 7, 8, 9, 0, 2, 3, 6, 2, 3, 6, 7, 8,
                             9, 2, 5, 6, 9, 2, 7, 8, 1, 3, 4, 7, 8, 9};
    int AT_row_starts_correct[] = {0, 2, 6, 10, 16, 20, 23, 29};
    mu_assert("error AT_vals", ARRAYS_EQUAL_DOUBLE(AT->x, AT_vals_correct, 29));
    mu_assert("error AT_cols", ARRAYS_EQUAL_INT(AT->i, AT_cols_correct, 29));
    CHECK_ROW_STARTS(AT, AT_row_starts_correct);

    // deallocate memory
    free_matrix(A);
    free_matrix(AT);
    PS_FREE(stgs);
    PS_FREE(constraints);
    free_state(data);
    PS_FREE(row_tags);

    return 0;
}

static const char *all_tests_constraints()
{
    mu_run_test(test_1_constraints, counter_constraints);
    mu_run_test(test_2_constraints, counter_constraints);
    mu_run_test(test_3_constraints, counter_constraints);
    mu_run_test(test_4_constraints, counter_constraints);

    return 0;
}

int test_constraints()
{
    const char *result = all_tests_constraints();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("Constraints: TEST FAILED!\n");
    }
    else
    {
        printf("Constraints: ALL TESTS PASSED\n");
    }
    printf("Constraints: Tests run: %d\n", counter_constraints);
    return result == 0;
}

#endif // TEST_CONSTRAINTS_H
