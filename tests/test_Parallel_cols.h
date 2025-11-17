#ifndef TEST_PARALLEL_COLS_H
#define TEST_PARALLEL_COLS_H

#include "API.h"
#include "Debugger.h"
#include "Matrix.h"
#include "Parallel_cols.h"

#include "Problem.h"
#include "math.h"
#include "minunit.h"
#include "test_macros.h"
#include <stdio.h>

static int counter_parallel_cols = 0;

/* Example 8 Gurobi paper (with two columns flipped)  */
static char *test_1_parallel_cols()
{
    double Ax[] = {-2, -1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    int nnz = 3;
    int n_rows = 1;
    int n_cols = 3;

    double lhs[] = {-INF};
    double rhs[] = {-10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {4, 3, 5};
    double c[] = {4, 2, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-1, -1};
    int Ai_correct[] = {0, 1};
    int Ap_correct[] = {0, 2};
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, 2));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, 2));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {11, 5};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that the objective function is correct
    double obj_correct[] = {2, 1};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 2));
    mu_assert("error offset", prob->obj->offset == 0);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Based on example 8 gurobi paper
    min. [4 2 1]x
    s.t. [-2 -1 -1]x <= -10
         0 <= x <= [INF 3 5]
*/
static char *test_2_parallel_cols()
{
    double Ax[] = {-2, -1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    int nnz = 3;
    int n_rows = 1;
    int n_cols = 3;

    double lhs[] = {-INF};
    double rhs[] = {-10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {INF, 3, 5};
    double c[] = {4, 2, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-1, -1};
    int Ai_correct[] = {0, 1};
    int Ap_correct[] = {0, 2};
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, 2));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, 2));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {INF, 5};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that the objective function is correct
    double obj_correct[] = {2, 1};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 2));
    mu_assert("error offset", prob->obj->offset == 0);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Based on example 8 gurobi paper
    min. [4 2 1]x
    s.t. [-2 -1 -1]x <= -10
         0 <= x <= [INF INF 5]
*/
static char *test_3_parallel_cols()
{
    double Ax[] = {-2, -1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    int nnz = 3;
    int n_rows = 1;
    int n_cols = 3;

    double lhs[] = {-INF};
    double rhs[] = {-10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {INF, INF, 5};
    double c[] = {4, 2, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-1, -1};
    int Ai_correct[] = {0, 1};
    int Ap_correct[] = {0, 2};
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, 2));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, 2));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {INF, 5};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that the objective function is correct
    double obj_correct[] = {2, 1};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 2));
    mu_assert("error offset", prob->obj->offset == 0);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Based on example 8 gurobi paper
    min. [4 2 1]x
    s.t. [-2 -1 -1]x <= -10
         [0, 0, 0] <= x <= [4 3 5]
*/
static char *test_4_parallel_cols()
{
    double Ax[] = {-2, -1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    int nnz = 3;
    int n_rows = 1;
    int n_cols = 3;

    double lhs[] = {-INF};
    double rhs[] = {-10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {4, 3, 5};
    double c[] = {4, 2, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-1, -1};
    int Ai_correct[] = {0, 1};
    int Ap_correct[] = {0, 2};
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, 2));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, 2));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {11, 5};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that the objective function is correct
    double obj_correct[] = {2, 1};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 2));
    mu_assert("error offset", prob->obj->offset == 0);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Based on example 8 gurobi paper
    min. [4 2 1]x
    s.t. [-2 -1 -1]x <= -10
         [-INF, 0, 0] <= x <= [4 3 5]
*/
static char *test_5_parallel_cols()
{
    double Ax[] = {-2, -1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    int nnz = 3;
    int n_rows = 1;
    int n_cols = 3;

    double lhs[] = {-INF};
    double rhs[] = {-10};
    double lbs[] = {-INF, 0, 0};
    double ubs[] = {4, 3, 5};
    double c[] = {4, 2, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-1, -1};
    int Ai_correct[] = {0, 1};
    int Ap_correct[] = {0, 2};
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, 2));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, 2));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {-INF, 0};
    double ubs_correct[] = {11, 5};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that the objective function is correct
    double obj_correct[] = {2, 1};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 2));
    mu_assert("error offset", prob->obj->offset == 0);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Based on example 8 gurobi paper
    min. [4 2 1]x
    s.t. [-2 -1 -1]x <= -10
         [-INF, -INF, 0] <= x <= [4 3 5]
*/
static char *test_6_parallel_cols()
{
    double Ax[] = {-2, -1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    int nnz = 3;
    int n_rows = 1;
    int n_cols = 3;

    double lhs[] = {-INF};
    double rhs[] = {-10};
    double lbs[] = {-INF, -INF, 0};
    double ubs[] = {4, 3, 5};
    double c[] = {4, 2, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-1};
    int Ai_correct[] = {0};
    int Ap_correct[] = {0, 1};
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, 1));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, 1));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {-INF};
    double ubs_correct[] = {11};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 1));

    // check that the objective function is correct
    double obj_correct[] = {2};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 1));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Based on example 8 gurobi paper
    min. [4 -2 1]x
    s.t. [2 -1 -1]x <= -10
         [-INF, 0, 0] <= x <= [4 3 5]
*/
static char *test_5_negated_parallel_cols()
{
    double Ax[] = {2, -1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    int nnz = 3;
    int n_rows = 1;
    int n_cols = 3;

    double lhs[] = {-INF};
    double rhs[] = {-10};
    double lbs[] = {-INF, 0, 0};
    double ubs[] = {4, 3, 5};
    double c[] = {4, -2, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-1, -1};
    int Ai_correct[] = {0, 1};
    int Ap_correct[] = {0, 2};
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, 2));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, 2));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {-8, 0};
    double ubs_correct[] = {INF, 5};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that the objective function is correct
    double obj_correct[] = {-2, 1};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 2));
    mu_assert("error offset", prob->obj->offset == 0);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Based on example 8 gurobi paper
    min. [4 -2 1]x
    s.t. [2 -1 -1]x <= -10
         [-INF, -INF, 0] <= x <= [4 3 5]
*/
static char *test_6_negated_parallel_cols()
{
    double Ax[] = {2, -1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    int nnz = 3;
    int n_rows = 1;
    int n_cols = 3;

    double lhs[] = {-INF};
    double rhs[] = {-10};
    double lbs[] = {-INF, -INF, 0};
    double ubs[] = {4, 3, 5};
    double c[] = {4, -2, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-1, -1};
    int Ai_correct[] = {0, 1};
    int Ap_correct[] = {0, 2};
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, 2));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, 2));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {-INF, 0};
    double ubs_correct[] = {INF, 5};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that the objective function is correct
    double obj_correct[] = {-2, 1};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 2));
    mu_assert("error offset", prob->obj->offset == 0);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Big test with all bounds equal to infinity
    min. [2  3  -4  -3  -6   6   3   1]x
    s.t. [1  1  -2   1   2   3  -1   5
          0  0   0   0   0  -3   1   1
         -2  2   4  -2  -4  -6   2   2
          0  0   0   0   0  -3   1   3
          3  3  -6   3   6   9  -3   4
          4  0  -8   4   8  12  -4   5
         -2  4   4  -2  -4  -6   2   6
          0  7   0   0   0   0   0   7
          0  5   0   0   0   0   0   8
          0  6   0   0   0   0   0   9]x <= [10 1 20 1 30 40 20 70 50 60]
*/
static char *test_7_parallel_cols()
{
    double Ax[] = {1,  1,  -2, 1,  2, 3, -1, 5,  -3, 1, 1, -2, 2, 4, -2, -4, -6,
                   2,  2,  -3, 1,  3, 3, 3,  -6, 3,  6, 9, -3, 4, 4, -8, 4,  8,
                   12, -4, 5,  -2, 4, 4, -2, -4, -6, 2, 6, 7,  7, 5, 8,  6,  9};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 5, 6, 7, 0, 1, 2, 3, 4, 5,
                6, 7, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 3, 4,
                5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 1, 7, 1, 7, 1, 7};
    int Ap[] = {0, 8, 11, 19, 22, 30, 37, 45, 47, 49, 51};
    int nnz = 51;
    int n_rows = 10;
    int n_cols = 8;

    double lhs[] = {-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {10, 1, 20, 1, 30, 40, 20, 70, 50, 60};
    double lbs[] = {-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double ubs[] = {INF, INF, INF, INF, INF, INF, INF, INF};
    double c[] = {2, 3, -4, -3, -6, 6, 3, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_parallel_cols(prob);
    mu_assert("error status", status == UNBNDORINFEAS);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Big test with finite bounds
    min. [2  3  -4  -3  -6   6   3   1]x
    s.t. [1  1  -2   1   2   3  -1   5
          0  0   0   0   0  -3   1   1
         -2  2   4  -2  -4  -6   2   2
          0  0   0   0   0  -3   1   3
          3  3  -6   3   6   9  -3   4
          4  0  -8   4   8  12  -4   5
         -2  4   4  -2  -4  -6   2   6
          0  7   0   0   0   0   0   7
          0  5   0   0   0   0   0   8
          0  6   0   0   0   0   0   9]x <= [10 1 20 1 30 40 20 70 50 60]
          -i <= x <= i
*/
static char *test_8_parallel_cols()
{
    double Ax[] = {1,  1,  -2, 1,  2, 3, -1, 5,  -3, 1, 1, -2, 2, 4, -2, -4, -6,
                   2,  2,  -3, 1,  3, 3, 3,  -6, 3,  6, 9, -3, 4, 4, -8, 4,  8,
                   12, -4, 5,  -2, 4, 4, -2, -4, -6, 2, 6, 7,  7, 5, 8,  6,  9};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 5, 6, 7, 0, 1, 2, 3, 4, 5,
                6, 7, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 3, 4,
                5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 1, 7, 1, 7, 1, 7};
    int Ap[] = {0, 8, 11, 19, 22, 30, 37, 45, 47, 49, 51};
    int nnz = 51;
    int n_rows = 10;
    int n_cols = 8;

    double lhs[] = {-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {10, 1, 20, 1, 30, 40, 20, 70, 50, 60};
    double lbs[] = {-1, -2, -3, -4, -5, -6, -7, -8};
    double ubs[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double c[] = {2, 3, -4, -3, -6, 6, 3, 1};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_parallel_cols(prob);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {1, -2, 1, 3,  -1, 5,  -3, 1, 1,  2, 4,  -2, -6, 2,
                           2, -3, 1, 3,  3,  -6, 3,  9, -3, 4, -8, 4,  12, -4,
                           5, 4,  4, -2, -6, 2,  6,  7, 7,  5, 8,  6,  9};
    int Ai_correct[] = {0, 1, 2, 3, 4, 5, 3, 4, 5, 0, 1, 2, 3, 4,
                        5, 3, 4, 5, 0, 1, 2, 3, 4, 5, 1, 2, 3, 4,
                        5, 0, 1, 2, 3, 4, 5, 0, 5, 0, 5, 0, 5};
    int Ap_correct[] = {0, 6, 9, 15, 18, 24, 29, 35, 37, 39, 41};

// on mac the order of the columns is different, so we only run this test on
// linux
#ifdef __linux__
    mu_assert("error Ax", ARRAYS_EQUAL(Ax_correct, A->x, A->nnz));
    mu_assert("error Ai", ARRAYS_EQUAL(Ai_correct, A->i, A->nnz));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check that new variable bounds are correct
    double lbs_correct[] = {-2, -3.5, -14, -6, -7, -8};
    double ubs_correct[] = {2, 3.5, 14, 6, 7, 8};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 6));

    // check that the objective function is correct
    double obj_correct[] = {3, -4, -3, 6, 3, 1};
    mu_assert("error obj", ARRAYS_EQUAL(obj_correct, prob->obj->c, 6));
    mu_assert("error offset", prob->obj->offset == 0);
#endif // __linux__

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_parallel_cols()
{
    mu_run_test(test_1_parallel_cols, counter_parallel_cols);
    mu_run_test(test_2_parallel_cols, counter_parallel_cols);
    mu_run_test(test_3_parallel_cols, counter_parallel_cols);
    mu_run_test(test_4_parallel_cols, counter_parallel_cols);
    mu_run_test(test_5_parallel_cols, counter_parallel_cols);
    mu_run_test(test_6_parallel_cols, counter_parallel_cols);
    mu_run_test(test_5_negated_parallel_cols, counter_parallel_cols);
    mu_run_test(test_6_negated_parallel_cols, counter_parallel_cols);
    mu_run_test(test_7_parallel_cols, counter_parallel_cols);
    mu_run_test(test_8_parallel_cols, counter_parallel_cols);

    return 0;
}

int test_parallel_cols()
{
    const char *result = all_tests_parallel_cols();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("parallel_cols: TEST FAILED!\n");
    }
    else
    {
        printf("parallel_cols: ALL TESTS PASSED\n");
    }
    printf("parallel_cols: Tests run: %d\n", counter_parallel_cols);
    return result == 0;
}

#endif // TEST_PARALLEL_COLS_H
