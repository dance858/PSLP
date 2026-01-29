#ifndef TEST_PARALLEL_ROWS_H
#define TEST_PARALLEL_ROWS_H

#include "Debugger.h"
#include "Matrix.h"
#include "PSLP_API.h"
#include "Parallel_rows.h"

#include "math.h"
#include "minunit.h"
#include "test_macros.h"
#include <stdio.h>

static int counter_parallel_rows = 0;

/* Test that rows with the same sparsity pattern are given the same sparsity ID.
    A_in = [1, 0, 0, 1, 0]
           [0, 1, 0, 1, 0]
           [1, 1, 0, 1, 0]
           [0, 1, 0, 1, 0]
           [1, 0, 0, 1, 0]
 */
static char *test_1_parallel_rows()
{
    double Ax[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 3, 1, 3, 0, 1, 3, 1, 3, 0, 3};
    int Ap[] = {0, 2, 4, 7, 9, 11};
    int nnz = 11;
    int n_rows = 5;
    int n_cols = 5;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);
    RowTag row_tags[] = {R_TAG_NONE, R_TAG_NONE, R_TAG_NONE, R_TAG_NONE, R_TAG_NONE};
    int sparsity_IDs[5];
    int coeff_hashes[5];
    compute_supp_and_coeff_hash(A, row_tags, sparsity_IDs, coeff_hashes,
                                R_TAG_INACTIVE);

    mu_assert("error", (sparsity_IDs[0] == sparsity_IDs[4] &&
                        sparsity_IDs[1] == sparsity_IDs[3]));

    free_matrix(A);
    return 0;
}

/* Test that rows with the same NORMALIZED coefficients in the same order
        are given the same hash value.
    A_in = [   1,   0,   0,  2,    0]
           [   0,   1,   0, -1,    0]
           [-0.5,   0,  -1,  0,    0]
           [ 0.7,   0,   0,  0,  1.4]
           [  -3,   3,   0,  0,    0]
           [   1,   0,   0,  3,    0]
 */
static char *test_2_parallel_rows()
{
    double Ax[] = {1, 2, 1, -1, -0.5, -1, 0.7, 1.4, -3, 3, 1, 3};
    int Ai[] = {0, 3, 1, 3, 0, 2, 0, 4, 0, 1, 0, 3};
    int Ap[] = {0, 2, 4, 6, 8, 10, 12};
    int nnz = 12;
    int n_rows = 6;
    int n_cols = 5;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);
    RowTag row_tags[] = {R_TAG_NONE, R_TAG_NONE, R_TAG_NONE,
                         R_TAG_NONE, R_TAG_NONE, R_TAG_NONE};
    int sparsity_IDs[6];
    int coeff_hashes[6];
    compute_supp_and_coeff_hash(A, row_tags, sparsity_IDs, coeff_hashes,
                                R_TAG_INACTIVE);

    mu_assert("error", (coeff_hashes[0] == coeff_hashes[2] &&
                        coeff_hashes[0] == coeff_hashes[3] &&
                        coeff_hashes[1] == coeff_hashes[4]));
    free_matrix(A);
    return 0;
}

/* Test that two groups of parallel rows are found.
    A_in = [   5/8,           0,      0,            15/8,    0]
           [   0,      sqrt(50),      0,    2 * sqrt(50),    0]
           [   0,             1,      0,               3,    0]
           [0,          -sqrt(5),     0,     -2 * sqrt(5),   0]
           [-3/4,             0,      0,             -9/4,   0]
 */
static char *test_3_parallel_rows()
{
    double Ax[] = {5.0 / 8, 15.0 / 8, sqrt(50),     2 * sqrt(50), 1,
                   3,       -sqrt(5), -2 * sqrt(5), -3.0 / 4,     -9.0 / 4};
    int Ai[] = {0, 3, 1, 3, 1, 3, 1, 3, 0, 3};
    int Ap[] = {0, 2, 4, 6, 8, 10};
    int nnz = 10;
    int n_rows = 5;
    int n_cols = 5;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);
    RowTag row_tags[] = {R_TAG_NONE, R_TAG_NONE, R_TAG_NONE, R_TAG_NONE, R_TAG_NONE};

    iVec *group_start = iVec_new(10);
    int parallel_rows[5];
    int sparsity_IDs[5];
    int coeff_hashes[5];

    find_parallel_rows(A, row_tags, group_start, parallel_rows, sparsity_IDs,
                       coeff_hashes, R_TAG_INACTIVE);

#ifdef _WIN32
    int parallel_rows_correct[] = {0, 4, 1, 3};
#else
    int parallel_rows_correct[] = {4, 0, 3, 1};
#endif

    int group_start_correct[] = {0, 2, 4};
    printf("parallel rows found:\n");
    print_int_array(parallel_rows, 4);
    fflush(stdout);
    mu_assert("error", ARRAYS_EQUAL_INT(parallel_rows_correct, parallel_rows, 4));
    mu_assert("error", group_start->len == 3);
    mu_assert("error", ARRAYS_EQUAL_INT(group_start_correct, group_start->data, 3));

    free_matrix(A);
    iVec_free(group_start);
    return 0;
}

int compare_ints(const void *a, const void *b)
{
    int arg1 = *(const int *) a;
    int arg2 = *(const int *) b;

    if (arg1 < arg2)
    {
        return -1;
    }

    if (arg1 > arg2)
    {
        return 1;
    }
    return 0;
}

/* Test that several groups of parallel rows are found */
static char *test_4_parallel_rows()
{
    int n_rows = 1000;
    int n_cols = 2000;
    double density = 0.1;
    Matrix *A = random_matrix_new(n_rows, n_cols, density);
    int jump = 50;

    // allocate a bunch of extra memory so we can replace rows as we wish
    int new_alloc = 2 * A->n_alloc;
    A->x = (double *) ps_realloc(A->x, new_alloc, sizeof(double));
    A->i = (int *) ps_realloc(A->i, new_alloc, sizeof(int));
    A->n_alloc = new_alloc;
    srand(1337);

    // -------------------------------------------------------------------
    //               create first group of parallel rows
    // -------------------------------------------------------------------
    int start1 = 0;
    int *cols1 = A->i + A->p[start1].start;
    double *vals1 = A->x + A->p[start1].start;
    int len1 = A->p[start1].end - A->p[start1].start;

    for (int i = 1; i < 8; ++i)
    {
        double ratio = (double) (rand() - rand()) / RAND_MAX;
        replace_row_A(A, start1 + i * jump, ratio, vals1, cols1, len1);
    }

    // -------------------------------------------------------------------
    //               create second group of parallel rows
    // -------------------------------------------------------------------
    int start2 = 10;
    int *cols2 = A->i + A->p[start2].start;
    double *vals2 = A->x + A->p[start2].start;
    int len2 = A->p[start2].end - A->p[start2].start;

    for (int i = 1; i < 8; ++i)
    {
        double ratio = (double) (rand() - rand()) / RAND_MAX;
        replace_row_A(A, start2 + i * jump, ratio, vals2, cols2, len2);
    }

    // -------------------------------------------------------------------
    //               create third group of parallel rows
    // -------------------------------------------------------------------
    int start3 = 20;
    int *cols3 = A->i + A->p[start3].start;
    double *vals3 = A->x + A->p[start3].start;
    int len3 = A->p[start3].end - A->p[start3].start;

    for (int i = 1; i < 8; ++i)
    {
        double ratio = (double) (rand() - rand()) / RAND_MAX;
        replace_row_A(A, start3 + i * jump, ratio, vals3, cols3, len3);
    }

    RowTag row_tags[1000];
    for (int i = 0; i < 1000; ++i)
    {
        row_tags[i] = R_TAG_NONE;
    }
    iVec *group_start = iVec_new(10);
    int parallel_rows[1000];
    int sparsity_IDs[1000];
    int coeff_hashes[1000];

    find_parallel_rows(A, row_tags, group_start, parallel_rows, sparsity_IDs,
                       coeff_hashes, R_TAG_INACTIVE);

    mu_assert("error group_starts len", group_start->len == 4);
    int group_start_correct[] = {0, 8, 16, 24};
    mu_assert("error group_starts data",
              ARRAYS_EQUAL_INT(group_start_correct, group_start->data, 4));

    // Below we verify that the answer is correct. We need to sort and do some
    // hacking because the rows appears in a different order on mac compared to
    // linux.)

    int first_group_correct[] = {0, 50, 100, 150, 200, 250, 300, 350};
    int second_group_correct[] = {10, 60, 110, 160, 210, 260, 310, 360};
    int third_group_correct[] = {20, 70, 120, 170, 220, 270, 320, 370};
    int sorted[8];

    // check that first 8 entries of parallell rows corresponds to one of the
    // groups above
    for (int i = 0; i < 8; ++i)
    {
        sorted[i] = parallel_rows[i];
    }
    qsort(sorted, 8, sizeof(int), compare_ints);

    if (sorted[0] == 0)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, first_group_correct, 8));
    }
    else if (sorted[0] == 10)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, second_group_correct, 8));
    }
    else
    {
        assert(sorted[0] == 20);
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, third_group_correct, 8));
    }

    // check that next 8 entries of parallell rows corresponds to one of the
    // groups above
    for (int i = 0; i < 8; ++i)
    {
        sorted[i] = parallel_rows[i + 8];
    }
    qsort(sorted, 8, sizeof(int), compare_ints);

    if (sorted[0] == 0)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, first_group_correct, 8));
    }
    else if (sorted[0] == 10)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, second_group_correct, 8));
    }
    else
    {
        assert(sorted[0] == 20);
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, third_group_correct, 8));
    }

    // check that next 8 entries of parallell rows corresponds to one of the
    // groups above
    for (int i = 0; i < 8; ++i)
    {
        sorted[i] = parallel_rows[i + 16];
    }
    qsort(sorted, 8, sizeof(int), compare_ints);

    if (sorted[0] == 0)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, first_group_correct, 8));
    }
    else if (sorted[0] == 10)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, second_group_correct, 8));
    }
    else
    {
        assert(sorted[0] == 20);
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, third_group_correct, 8));
    }

    free_matrix(A);
    iVec_free(group_start);
    return 0;
}

/* Test that several groups of parallel rows are found, while some are
   manually marked as inactive */
static char *test_5_parallel_rows()
{
    int n_rows = 1000;
    int n_cols = 2000;
    double density = 0.1;
    Matrix *A = random_matrix_new(n_rows, n_cols, density);
    int jump = 50;

    // allocate a bunch of extra memory so we can replace rows as we wish
    int new_alloc = 2 * A->n_alloc;
    A->x = (double *) ps_realloc(A->x, new_alloc, sizeof(double));
    A->i = (int *) ps_realloc(A->i, new_alloc, sizeof(int));
    A->n_alloc = new_alloc;
    srand(1337);

    // -------------------------------------------------------------------
    //               create first group of parallel rows
    // -------------------------------------------------------------------
    int start1 = 0;
    int *cols1 = A->i + A->p[start1].start;
    double *vals1 = A->x + A->p[start1].start;
    int len1 = A->p[start1].end - A->p[start1].start;

    for (int i = 1; i < 8; ++i)
    {
        double ratio = (double) (rand() - rand()) / RAND_MAX;
        replace_row_A(A, start1 + i * jump, ratio, vals1, cols1, len1);
    }

    // -------------------------------------------------------------------
    //               create second group of parallel rows
    // -------------------------------------------------------------------
    int start2 = 10;
    int *cols2 = A->i + A->p[start2].start;
    double *vals2 = A->x + A->p[start2].start;
    int len2 = A->p[start2].end - A->p[start2].start;

    for (int i = 1; i < 8; ++i)
    {
        double ratio = (double) (rand() - rand()) / RAND_MAX;
        replace_row_A(A, start2 + i * jump, ratio, vals2, cols2, len2);
    }

    // -------------------------------------------------------------------
    //               create third group of parallel rows
    // -------------------------------------------------------------------
    int start3 = 20;
    int *cols3 = A->i + A->p[start3].start;
    double *vals3 = A->x + A->p[start3].start;
    int len3 = A->p[start3].end - A->p[start3].start;

    for (int i = 1; i < 8; ++i)
    {
        double ratio = (double) (rand() - rand()) / RAND_MAX;
        replace_row_A(A, start3 + i * jump, ratio, vals3, cols3, len3);
    }

    RowTag row_tags[1000];
    for (int i = 0; i < 1000; ++i)
    {
        row_tags[i] = R_TAG_NONE;
    }
    iVec *group_start = iVec_new(10);
    int parallel_rows[1000] = {0};
    int sparsity_IDs[1000] = {0};
    int coeff_hashes[1000] = {0};

    row_tags[0] = R_TAG_INACTIVE;
    row_tags[200] = R_TAG_INACTIVE;
    row_tags[360] = R_TAG_INACTIVE;
    row_tags[310] = R_TAG_INACTIVE;

    find_parallel_rows(A, row_tags, group_start, parallel_rows, sparsity_IDs,
                       coeff_hashes, R_TAG_INACTIVE);

    mu_assert("error group_starts len", group_start->len == 4);
    // int group_start_correct[] = {0, 8, 14, 20};
    // mu_assert("error group_starts data",
    //           ARRAYS_EQUAL_INT(group_start_correct, group_start->data, 4));

    int first_group_correct[] = {10, 60, 110, 160, 210, 260};
    int second_group_correct[] = {20, 70, 120, 170, 220, 270, 320, 370};
    int third_group_correct[] = {50, 100, 150, 250, 300, 350};

    // figure out length of the groups
    int len_group1 = group_start->data[1] - group_start->data[0];
    int len_group2 = group_start->data[2] - group_start->data[1];
    int len_group3 = group_start->data[3] - group_start->data[2];

    int sorted[8];

    // check that first entries of parallell rows corresponds to one of the
    // groups above
    for (int i = 0; i < len_group1; ++i)
    {
        sorted[i] = parallel_rows[i];
    }
    qsort(sorted, len_group1, sizeof(int), compare_ints);

    if (sorted[0] == 10)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, first_group_correct, 6));
    }
    else if (sorted[0] == 20)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, second_group_correct, 8));
    }
    else
    {
        assert(sorted[0] == 50);
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, third_group_correct, 6));
    }

    // check that next entries of parallell rows corresponds to one of the
    // groups above
    for (int i = 0; i < len_group2; ++i)
    {
        sorted[i] = parallel_rows[i + len_group1];
    }
    qsort(sorted, len_group2, sizeof(int), compare_ints);

    if (sorted[0] == 10)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, first_group_correct, 6));
    }
    else if (sorted[0] == 20)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, second_group_correct, 8));
    }
    else
    {
        assert(sorted[0] == 50);
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, third_group_correct, 6));
    }

    // check that next entries of parallell rows corresponds to one of the
    // groups above
    for (int i = 0; i < len_group3; ++i)
    {
        sorted[i] = parallel_rows[i + len_group1 + len_group2];
    }
    qsort(sorted, len_group3, sizeof(int), compare_ints);

    if (sorted[0] == 10)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, first_group_correct, 6));
    }
    else if (sorted[0] == 20)
    {
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, second_group_correct, 8));
    }
    else
    {
        assert(sorted[0] == 50);
        mu_assert("error first group",
                  ARRAYS_EQUAL_INT(sorted, third_group_correct, 6));
    }

    free_matrix(A);
    iVec_free(group_start);
    return 0;
}

/*  Infeasible problem   (similar to test10 old code)
    min.   [0 0] x
    s.t.   x1 + x2 = 2
           -3x1 + 4x2 <= 0
           2x1 + 2x2 <= 3
*/
static char *test_6_parallel_rows()
{
    double Ax[] = {1, 1, -3, 4, 2, 2};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {2, -INF, -INF};
    double rhs[] = {2, 0, 3};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status == INFEASIBLE);
    problem_clean(prob, true);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    mu_assert("rows", check_row_starts(A, Ap));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Infeasible problem  (similar to test 11 old code)
    min.   [0 0] x
    s.t.   x1 + x2 = 2
           -3x1 + 4x2 <= 0
           -2x1 - 2x2 <= -5
*/
static char *test_7_parallel_rows()
{
    double Ax[] = {1, 1, -3, 4, -2, -2};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {2, -INF, -INF};
    double rhs[] = {2, 0, -5};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status == INFEASIBLE);
    problem_clean(prob, true);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    mu_assert("rows", check_row_starts(A, Ap));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Infeasible problem  (similar to test 12 old code)
    min.   [0 0] x
    s.t.   x1 + x2 = 2
           -3x1 + 4x2 <= 0
           2x1 + 2x2 >= 5
*/
static char *test_8_parallel_rows()
{
    double Ax[] = {1, 1, -3, 4, 2, 2};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {2, -INF, 5};
    double rhs[] = {2, 0, INF};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status == INFEASIBLE);
    problem_clean(prob, true);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    mu_assert("rows", check_row_starts(A, Ap));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Infeasible problem  (similar to test 13 old code)
    min.   [0 0] x
    s.t.   x1 + x2 = 2
           -3x1 + 4x2 <= 0
           -2x1 - 2x2 >= -2
*/
static char *test_9_parallel_rows()
{
    double Ax[] = {1, 1, -3, 4, -2, -2};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {2, -INF, -2};
    double rhs[] = {2, 0, INF};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status == INFEASIBLE);
    problem_clean(prob, true);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    mu_assert("rows", check_row_starts(A, Ap));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Infeasible problem  (similar to test 14 old code)
    min.   [0 0] x
    s.t.   x1 + x2 <= 1
           -3x1 + 4x2 <= 0
           x1 + x2 = 2
*/
static char *test_10_parallel_rows()
{
    double Ax[] = {1, 1, -3, 4, 1, 1};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {-INF, -INF, 2};
    double rhs[] = {1, 0, 2};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status == INFEASIBLE);
    problem_clean(prob, true);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    mu_assert("rows", check_row_starts(A, Ap));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Infeasible problem  (similar to test 15 old code)
    min.   [0 0] x
    s.t.   -x1 - x2 <= -3
           -3x1 + 4x2 <= 0
           x1 + x2 = 2
*/
static char *test_11_parallel_rows()
{
    double Ax[] = {-1, -1, -3, 4, 1, 1};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {-INF, -INF, 2};
    double rhs[] = {-3, 0, 2};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status == INFEASIBLE);
    problem_clean(prob, true);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    mu_assert("rows", check_row_starts(A, Ap));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Infeasible problem  (similar to test 16 old code)
    min.   [0 0] x
    s.t.   x1 + x2 >= 3
           -3x1 + 4x2 <= 0
           x1 + x2 = 2
*/
static char *test_12_parallel_rows()
{
    double Ax[] = {1, 1, -3, 4, 1, 1};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {3, -INF, 2};
    double rhs[] = {INF, 0, 2};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status == INFEASIBLE);
    problem_clean(prob, true);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    mu_assert("rows", check_row_starts(A, Ap));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Infeasible problem  (similar to test 17 old code)
    min.   [0 0] x
    s.t.   -x1 - x2 >= -1
           -3x1 + 4x2 <= 0
           x1 + x2 = 2
*/
static char *test_13_parallel_rows()
{
    double Ax[] = {-1, -1, -3, 4, 1, 1};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {-1, -INF, 2};
    double rhs[] = {INF, 0, 2};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status == INFEASIBLE);
    problem_clean(prob, true);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    mu_assert("rows", check_row_starts(A, Ap));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Remaining row is an equality  (similar to test 18 old code)
    min. c x
    s.t.            random constraint
          -3q1 b1 <= q1 * a @ x <= 2 q1 b1
                         a @ x  =      b1
                    q2 * a @ x  =   q2 b1
                    q3 * a @ x  =   q3 b1
                    random constraint
*/
static char *test_14_parallel_rows()
{

    // -----------------------------------------------
    //            build A matrix with brute force
    // -----------------------------------------------
    double row0_vals[] = {2, 1, -3, 4, -5};
    int row0_cols[] = {0, 1, 2, 3, 4};

    double ax[3] = {1.3, 2.1, -1.7};
    double ai[3] = {0, 1, 3};
    double b1 = 1.4;

    double row1_vals[3];
    double row3_vals[3];
    double row4_vals[3];
    int row1_cols[3];
    int row3_cols[3];
    int row4_cols[3];

    double q1 = 1.2;
    double q2 = -1.3;
    double q3 = 1.5;

    int i;
    for (i = 0; i < 3; ++i)
    {
        row1_vals[i] = q1 * ax[i];
        row3_vals[i] = q2 * ax[i];
        row4_vals[i] = q3 * ax[i];
        row1_cols[i] = ai[i];
        row3_cols[i] = ai[i];
        row4_cols[i] = ai[i];
    }

    double row5_vals[] = {1, 1.1, 1, 1};
    int row5_cols[] = {0, 1, 2, 4};

    int nnz = 5 + 3 + 3 + 3 + 3 + 4;
    int n_rows = 6;
    int n_cols = 5;

    double Ax[5 + 3 + 3 + 3 + 3 + 4];
    int Ai[5 + 3 + 3 + 3 + 3 + 4];
    int Ap[7];

    Ap[0] = 0;
    for (i = 0; i < 5; ++i)
    {
        Ax[i] = row0_vals[i];
        Ai[i] = row0_cols[i];
    }
    Ap[1] = Ap[0] + 5;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[1] + i] = row1_vals[i];
        Ai[Ap[1] + i] = row1_cols[i];
    }
    Ap[2] = Ap[1] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[2] + i] = ax[i];
        Ai[Ap[2] + i] = ai[i];
    }
    Ap[3] = Ap[2] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[3] + i] = row3_vals[i];
        Ai[Ap[3] + i] = row3_cols[i];
    }
    Ap[4] = Ap[3] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[4] + i] = row4_vals[i];
        Ai[Ap[4] + i] = row4_cols[i];
    }
    Ap[5] = Ap[4] + 3;

    for (i = 0; i < 4; ++i)
    {
        Ax[Ap[5] + i] = row5_vals[i];
        Ai[Ap[5] + i] = row5_cols[i];
    }
    Ap[6] = Ap[5] + 4;

    // -----------------------------------------------
    //               begin test
    // -----------------------------------------------
    double lhs[] = {-2.1, -3 * q1 * b1, b1, q2 * b1, q3 * b1, -INF};
    double rhs[] = {3.1, 2 * q1 * b1, b1, q2 * b1, q3 * b1, 5};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {0, 0, 0, 0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status != INFEASIBLE);
    problem_clean(prob, true);

    // check A matrix
    double Ax_correct[] = {row0_vals[0], row0_vals[1], row0_vals[2], row0_vals[3],
                           row0_vals[4], ax[0],        ax[1],        ax[2],
                           row5_vals[0], row5_vals[1], row5_vals[2], row5_vals[3]};

    int Ai_correct[] = {row0_cols[0], row0_cols[1], row0_cols[2], row0_cols[3],
                        row0_cols[4], ai[0],        ai[1],        ai[2],
                        row5_cols[0], row5_cols[1], row5_cols[2], row5_cols[3]};
    int Ap_correct[] = {0, 5, 8, 12};

    printf("A->x:\n");
    print_double_array(A->x, A->nnz);
    printf("Ax-correct:\n");
    print_double_array(Ax_correct, 12);
    fflush(stdout);

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 12));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 12));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check row tags
    RowTag row_tags_correct[] = {R_TAG_NONE, R_TAG_EQ, R_TAG_LHS_INF};
    mu_assert("error row tags",
              ARRAYS_EQUAL_ROWTAG(row_tags_correct, constraints->row_tags, 3));

    // check lhs and rhs
    double lhs_correct[] = {-2.1, b1, -INF};
    double rhs_correct[] = {3.1, b1, 5};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Remaining row is an inequality  (similar to test 19 old code)
    min. c x
    s.t.            random constraint
                    random constraint
                         a @ x  <= 10
                3  <=  2 * a @ x   <=  22
               -18 <= -2 * a @ x
                    random constraint

    After presolve we should have random constraint, random constraint,
    1.5 <= a @ x <= 9 (or something equivalent), random constraint
*/
static char *test_15_parallel_rows()
{

    // -----------------------------------------------
    //            build A matrix with brute force
    // -----------------------------------------------
    double row0_vals[] = {2, 1, -3, 4, -5};
    int row0_cols[] = {0, 1, 2, 3, 4};

    double ax[3] = {1.3, 2.1, -1.7};
    double ai[3] = {0, 1, 3};

    double row1_vals[3];
    double row3_vals[3];
    double row4_vals[3];
    int row1_cols[3];
    int row3_cols[3];
    int row4_cols[3];

    double q1 = 1.2;
    double q2 = 2;
    double q3 = -2;

    int i;
    for (i = 0; i < 3; ++i)
    {
        row1_vals[i] = q1 * ax[i] + i;
        row3_vals[i] = q2 * ax[i];
        row4_vals[i] = q3 * ax[i];
        row1_cols[i] = ai[i];
        row3_cols[i] = ai[i];
        row4_cols[i] = ai[i];
    }

    double row5_vals[] = {1, 1.1, 1, 1};
    int row5_cols[] = {0, 1, 2, 4};

    int nnz = 5 + 3 + 3 + 3 + 3 + 4;
    int n_rows = 6;
    int n_cols = 5;

    double Ax[5 + 3 + 3 + 3 + 3 + 4];
    int Ai[5 + 3 + 3 + 3 + 3 + 4];
    int Ap[7];

    Ap[0] = 0;
    for (i = 0; i < 5; ++i)
    {
        Ax[i] = row0_vals[i];
        Ai[i] = row0_cols[i];
    }
    Ap[1] = Ap[0] + 5;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[1] + i] = row1_vals[i];
        Ai[Ap[1] + i] = row1_cols[i];
    }
    Ap[2] = Ap[1] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[2] + i] = ax[i];
        Ai[Ap[2] + i] = ai[i];
    }
    Ap[3] = Ap[2] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[3] + i] = row3_vals[i];
        Ai[Ap[3] + i] = row3_cols[i];
    }
    Ap[4] = Ap[3] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[4] + i] = row4_vals[i];
        Ai[Ap[4] + i] = row4_cols[i];
    }
    Ap[5] = Ap[4] + 3;

    for (i = 0; i < 4; ++i)
    {
        Ax[Ap[5] + i] = row5_vals[i];
        Ai[Ap[5] + i] = row5_cols[i];
    }
    Ap[6] = Ap[5] + 4;

    // -----------------------------------------------
    //               begin test
    // -----------------------------------------------
    double lhs[] = {-2.1, -1, -INF, 3, -18, -INF};
    double rhs[] = {3.1, 3, 10, 22, INF, 5};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {0, 0, 0, 0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status != INFEASIBLE);
    problem_clean(prob, true);

    // check A matrix
    double Ax_correct[] = {row0_vals[0], row0_vals[1], row0_vals[2], row0_vals[3],
                           row0_vals[4], row1_vals[0], row1_vals[1], row1_vals[2],
                           2 * ax[0],    2 * ax[1],    2 * ax[2],    row5_vals[0],
                           row5_vals[1], row5_vals[2], row5_vals[3]};

    int Ai_correct[] = {row0_cols[0], row0_cols[1], row0_cols[2], row0_cols[3],
                        row0_cols[4], row1_cols[0], row1_cols[1], row1_cols[2],
                        ai[0],        ai[1],        ai[2],        row5_cols[0],
                        row5_cols[1], row5_cols[2], row5_cols[3]};
    int Ap_correct[] = {0, 5, 8, 11, 15};

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 15));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 15));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check row tags
    RowTag row_tags_correct[] = {R_TAG_NONE, R_TAG_NONE, R_TAG_NONE, R_TAG_LHS_INF};
    mu_assert("error row tags",
              ARRAYS_EQUAL_ROWTAG(row_tags_correct, constraints->row_tags, 4));

    // check lhs and rhs
    double lhs_correct[] = {-2.1, -1, 3, -INF};
    double rhs_correct[] = {3.1, 3, 18, 5};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 4));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 4));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Remaining row is an inequality  (similar to test 20 old code)
    min. c x
    s.t.            random constraint
                    random constraint
                3  <=  2 * a @ x   <=  22
               -18 <= -2 * a @ x
                           a @ x   <= 10
                    random constraint

    After presolve we should have random constraint, random constraint,
    1.5 <= a @ x <= 9 (or something equivalent), random constraint
*/
static char *test_16_parallel_rows()
{

    // -----------------------------------------------
    //            build A matrix with brute force
    // -----------------------------------------------
    double row0_vals[] = {2, 1, -3, 4, -5};
    int row0_cols[] = {0, 1, 2, 3, 4};

    double ax[3] = {1.3, 2.1, -1.7};
    double ai[3] = {0, 1, 3};

    double row1_vals[3];
    double row3_vals[3];
    double row2_vals[3];
    int row1_cols[3];
    int row3_cols[3];
    int row2_cols[3];

    int i;
    for (i = 0; i < 3; ++i)
    {
        row1_vals[i] = 1.2 * ax[i] + i;
        row2_vals[i] = 2 * ax[i];
        row3_vals[i] = -2 * ax[i];
        row1_cols[i] = ai[i];
        row2_cols[i] = ai[i];
        row3_cols[i] = ai[i];
    }

    double row5_vals[] = {1, 1.1, 1, 1};
    int row5_cols[] = {0, 1, 2, 4};

    int nnz = 5 + 3 + 3 + 3 + 3 + 4;
    int n_rows = 6;
    int n_cols = 5;

    double Ax[5 + 3 + 3 + 3 + 3 + 4];
    int Ai[5 + 3 + 3 + 3 + 3 + 4];
    int Ap[7];

    Ap[0] = 0;
    for (i = 0; i < 5; ++i)
    {
        Ax[i] = row0_vals[i];
        Ai[i] = row0_cols[i];
    }
    Ap[1] = Ap[0] + 5;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[1] + i] = row1_vals[i];
        Ai[Ap[1] + i] = row1_cols[i];
    }
    Ap[2] = Ap[1] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[2] + i] = row2_vals[i];
        Ai[Ap[2] + i] = row2_cols[i];
    }
    Ap[3] = Ap[2] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[3] + i] = row3_vals[i];
        Ai[Ap[3] + i] = row3_cols[i];
    }
    Ap[4] = Ap[3] + 3;

    for (i = 0; i < 3; ++i)
    {
        Ax[Ap[4] + i] = ax[i];
        Ai[Ap[4] + i] = ai[i];
    }
    Ap[5] = Ap[4] + 3;

    for (i = 0; i < 4; ++i)
    {
        Ax[Ap[5] + i] = row5_vals[i];
        Ai[Ap[5] + i] = row5_cols[i];
    }
    Ap[6] = Ap[5] + 4;

    // -----------------------------------------------
    //               begin test
    // -----------------------------------------------
    double lhs[] = {-2.1, -1, 3, -18, -INF, -INF};
    double rhs[] = {3.1, 3, 22, INF, 10, 5};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {0, 0, 0, 0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_rows(constraints);
    mu_assert("error", status != INFEASIBLE);
    problem_clean(prob, true);

    // check A matrix
    double Ax_correct[] = {row0_vals[0], row0_vals[1], row0_vals[2], row0_vals[3],
                           row0_vals[4], row1_vals[0], row1_vals[1], row1_vals[2],
                           -2 * ax[0],   -2 * ax[1],   -2 * ax[2],   row5_vals[0],
                           row5_vals[1], row5_vals[2], row5_vals[3]};

    int Ai_correct[] = {row0_cols[0], row0_cols[1], row0_cols[2], row0_cols[3],
                        row0_cols[4], row1_cols[0], row1_cols[1], row1_cols[2],
                        ai[0],        ai[1],        ai[2],        row5_cols[0],
                        row5_cols[1], row5_cols[2], row5_cols[3]};
    int Ap_correct[] = {0, 5, 8, 11, 15};

    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 15));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 15));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check row tags
    RowTag row_tags_correct[] = {R_TAG_NONE, R_TAG_NONE, R_TAG_NONE, R_TAG_LHS_INF};
    mu_assert("error row tags",
              ARRAYS_EQUAL_ROWTAG(row_tags_correct, constraints->row_tags, 4));

    // check lhs and rhs
    double lhs_correct[] = {-2.1, -1, -18, -INF};
    double rhs_correct[] = {3.1, 3, -3, 5};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 4));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 4));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_parallel_rows()
{
    mu_run_test(test_1_parallel_rows, counter_parallel_rows);
    mu_run_test(test_2_parallel_rows, counter_parallel_rows);
    mu_run_test(test_3_parallel_rows, counter_parallel_rows);
    mu_run_test(test_4_parallel_rows, counter_parallel_rows);
    mu_run_test(test_5_parallel_rows, counter_parallel_rows);
    mu_run_test(test_6_parallel_rows, counter_parallel_rows);
    mu_run_test(test_7_parallel_rows, counter_parallel_rows);
    mu_run_test(test_8_parallel_rows, counter_parallel_rows);
    mu_run_test(test_9_parallel_rows, counter_parallel_rows);
    mu_run_test(test_10_parallel_rows, counter_parallel_rows);
    mu_run_test(test_11_parallel_rows, counter_parallel_rows);
    mu_run_test(test_12_parallel_rows, counter_parallel_rows);
    mu_run_test(test_13_parallel_rows, counter_parallel_rows);
    mu_run_test(test_14_parallel_rows, counter_parallel_rows);
    mu_run_test(test_15_parallel_rows, counter_parallel_rows);
    mu_run_test(test_16_parallel_rows, counter_parallel_rows);

    return 0;
}

int test_parallel_rows()
{
    const char *result = all_tests_parallel_rows();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("parallel_rows: TEST FAILED!\n");
    }
    else
    {
        printf("parallel_rows: ALL TESTS PASSED\n");
    }
    printf("parallel_rows: Tests run: %d\n", counter_parallel_rows);
    return result == 0;
}

#endif // TEST_PARALLEL_ROWS_H
