#ifndef TEST_DOMAINPROPAGATION_H
#define TEST_DOMAINPROPAGATION_H

#include "CoreTransformations.h"
#include "Numerics.h"
#include "Primal_propagation.h"
#include "SimpleReductions.h"
#include "Tags.h"
#include "minunit.h"
#include <stdio.h>

static int counter_domain = 0;

// clang-format off
/*
The following tests are/should be implemented for the DomainPropagation class
(a "✓" indicates that the test has been implemented).
---- Tests based on forcing rows ----
Test 1: Lower constraint that is forcing wrt the initial bounds (✓) 
        (integer arithmetic). 
Test 2: Lower constraint that is forcing wrt the initial bounds (✓)
        (floating point number arithmetic).
Test 3: Upper constraint that is forcing wrt the initial bounds (✓) 
        (integer arithmetic). 
Test 4: Upper constraint that is forcing wrt the initial bounds (✓)
        (floating point number arithmetic).
Test 5: Constraint 1 tightens an infinite upper bound so constraint 2 becomes (✓)
         lower forcing (integer arithmetic and NOT all bounds are tightened).
Test 6: Constraint 1 tightens an infinite lower bound so constraint 2 becomes (✓)
         upper forcing (integer arithmetic and NOT all bounds are tightened).
Test 7: Test 5 two rows swapped (✓)
Test 8: Test 6 two rows swapped (✓)
*/
// clang-format on

/*  Lower constraint that is forcing wrt the initial bounds
        (integer arithmetic).
    min. [-1 -1 -1 3 2]x
    s.t. [ 1  0   1   1   3]   >= 4
         [-1 -2   2   0   0] x >= 2
         [ 1  4   0   3  -1]   <= 5
         [ 2  0   0  -2   2]   >= 6
         [ 0  0   0   2  -3]   >= 0
         x1, x2 >= 0
         x3 <= 1
         x4, x5 >= 0
*/
static char *test_1_domain()
{
    double Ax[] = {1, 1, 1, 3, -1, -2, 2, 1, 4, 3, -1, 2, -2, 2, 2, -3};
    int Ai[] = {0, 2, 3, 4, 0, 1, 2, 0, 1, 3, 4, 0, 3, 4, 3, 4};
    int Ap[] = {0, 4, 7, 11, 14, 16};
    int nnz = 16;
    int n_rows = 5;
    int n_cols = 5;

    double lhs[] = {4, 2, -INF, 6, 0};
    double rhs[] = {INF, INF, 5, INF, INF};
    double lbs[5] = {0, 0, -INF, 0, 0};
    double ubs[] = {INF, INF, 1, INF, INF};
    double c[] = {-1, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    // check that new A is correct
    double Ax_correct[] = {1, 3, 3, -1, -2, 2, 2, -3};
    int Ai_correct[] = {0, 1, 0, 1, 0, 1, 0, 1};
    int Ap_correct[] = {0, 2, 4, 6, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {INF, INF};
    mu_assert("error bounds",
              check_bounds(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that lhs and rhs are correct
    double lhs_correct[] = {3, -INF, 6, 0};
    double rhs_correct[] = {INF, 5, INF, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 4));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 4));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Lower constraint that is forcing wrt the initial bounds
        (floating point number arithmetic). The reduction is very
        sensitive to the size of the perturbation of the second constraint.
    min. [-1 -1 -1 3 2]x
    s.t. [ 1      0      1     1   3]   >= 4
         [-2/13  -1/17   7/3   0   0] x >= -2/(13*3) - 1/(17*7) + 7/(3*9) -
   1e-10 [ 1      4      0     3  -1]   <= 5 [ 2      0      0    -2   2]   >= 6
         x1 >= 1/3
         x2 >= 1/7
         x3 <= 1/9
         x4, x5 >= 0
*/
static char *test_2_domain()
{
    double Ax[] = {1, 1, 1, 3, -2.0 / 13, -1.0 / 17, 7.0 / 3, 1, 4, 3, -1, 2, -2, 2};
    int Ai[] = {0, 2, 3, 4, 0, 1, 2, 0, 1, 3, 4, 0, 3, 4};
    int Ap[] = {0, 4, 7, 11, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {4, -2 / (13 * 3) - 1 / (17 * 7) + 7 / (3 * 9) - 1e-10, -INF, 6};
    double rhs[] = {INF, INF, 5, INF};
    double lbs[5] = {1 / 3, 1 / 7, -INF, 0, 0};
    double ubs[] = {INF, INF, 1 / 9, INF, INF};
    double c[] = {-1, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    // check that new A is correct
    double Ax_correct[] = {1, 3, 3, -1, -2, 2};
    int Ai_correct[] = {0, 1, 0, 1, 0, 1};
    int Ap_correct[] = {0, 2, 4, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {INF, INF};
    mu_assert("error bounds",
              check_bounds(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that objective offset is correct
    mu_assert("error offset", prob->obj->offset == -1 / 3 - 1 / 7 - 1 / 9);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4 - 1 / 9 - 1 / 3, -INF, 6 - 2 / 3};
    double rhs_correct[] = {INF, 5 - 1 / 3 - 4 / 7, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Upper constraint that is forcing wrt the initial bounds
        (integer arithmetic).
    min. [1 1 -1 3 2]x
    s.t. [ 1  0   1   1   3]   >= 4
         [-1 -2   2   0   0] x <= 2
         [ 1  4   0   3  -1]   <= 5
         [ 2  0   0  -2   2]   >= 6
         x1, x2 <= 0
         x3 >= 1
         x4, x5 >= 0
*/
static char *test_3_domain()
{
    double Ax[] = {1, 1, 1, 3, -1, -2, 2, 1, 4, 3, -1, 2, -2, 2};
    int Ai[] = {0, 2, 3, 4, 0, 1, 2, 0, 1, 3, 4, 0, 3, 4};
    int Ap[] = {0, 4, 7, 11, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {4, -INF, -INF, 6};
    double rhs[] = {INF, 2, 5, INF};
    double lbs[5] = {-INF, -INF, 1, 0, 0};
    double ubs[] = {0, 0, INF, INF, INF};
    double c[] = {1, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    // check that new A is correct
    double Ax_correct[] = {1, 3, 3, -1, -2, 2};
    int Ai_correct[] = {0, 1, 0, 1, 0, 1};
    int Ap_correct[] = {0, 2, 4, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {INF, INF};
    mu_assert("error bounds",
              check_bounds(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that objective offset is correct
    mu_assert("error offset", prob->obj->offset == -1);

    // check that lhs and rhs are correct
    double lhs_correct[] = {3, -INF, 6};
    double rhs_correct[] = {INF, 5, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Upper constraint that is forcing wrt the initial bounds
        (floating point number arithmetic). The reduction is very
        sensitive to the size of the perturbation of the second constraint.
    min. [1 1 -1 3 2]x
    s.t. [ 1      0      1     1   3]   >= 4
         [-2/13  -1/17   7/3   0   0] x <= -2/(13*3) - 1/(17*7) + 7/(3*9) +
   1e-10 [ 1      4      0     3  -1]   <= 5 [ 2      0      0    -2   2]   >= 6
         x1 <= 1/3
         x2 <= 1/7
         x3 >= 1/9
         x4, x5 >= 0
*/
static char *test_4_domain()
{
    double Ax[] = {1, 1, 1, 3, -2.0 / 13, -1.0 / 17, 7.0 / 3, 1, 4, 3, -1, 2, -2, 2};
    int Ai[] = {0, 2, 3, 4, 0, 1, 2, 0, 1, 3, 4, 0, 3, 4};
    int Ap[] = {0, 4, 7, 11, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {4, -INF, -INF, 6};
    double rhs[] = {INF, -2 / (13 * 3) - 1 / (17 * 7) + 7 / (3 * 9) + 1e-10, 5, INF};
    double lbs[5] = {-INF, -INF, 1 / 9, 0, 0};
    double ubs[] = {1 / 3, 1 / 7, INF, INF, INF};
    double c[] = {1, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    // check that new A is correct
    double Ax_correct[] = {1, 3, 3, -1, -2, 2};
    int Ai_correct[] = {0, 1, 0, 1, 0, 1};
    int Ap_correct[] = {0, 2, 4, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {INF, INF};
    mu_assert("error bounds",
              check_bounds(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that objective offset is correct
    mu_assert("error offset", prob->obj->offset == 1 / 3 + 1 / 7 - 1 / 9);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4 - 1 / 9 - 1 / 3, -INF, 6 - 2 / 3};
    double rhs_correct[] = {INF, 5 - 1 / 3 - 4 / 7, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Constraint 1 tightens an infinite upper bound so constraint 2 becomes
        lower forcing (integer arithmetic and NOT all bounds are tightened).

    min. [1 1 -1 -3 2]x
    s.t. [ 1  0    1   0   0]   <= 4
         [ 2  1    0   0   0] x >= 10
         [ 1  1   -1   3  -1]   <= 7
         [ 2  0    0  -2   2]   >= 6
         x1, x2, x3, x4, x5 >= 0
         x2 <= 2
*/
static char *test_5_domain()
{
    double Ax[] = {1, 1, 2, 1, 1, 1, -1, 3, -1, 2, -2, 2};
    int Ai[] = {0, 2, 0, 1, 0, 1, 2, 3, 4, 0, 3, 4};
    int Ap[] = {0, 2, 4, 9, 12};
    int nnz = 12;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {-INF, 10, -INF, 6};
    double rhs[] = {4, INF, 7, INF};
    double lbs[5] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, 2, INF, INF, INF};
    double c[] = {1, 1, -1, -3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    // check that new A is correct
    double Ax_correct[] = {1, -1, 3, -1, -2, 2};
    int Ai_correct[] = {0, 0, 1, 2, 1, 2};
    int Ap_correct[] = {0, 1, 4, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0};
    double ubs_correct[] = {4, INF, INF};
    mu_assert("error bounds",
              check_bounds(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that objective offset is correct
    mu_assert("error offset", prob->obj->offset == 6);

    // check that lhs and rhs are correct
    double lhs_correct[] = {-INF, -INF, -2};
    double rhs_correct[] = {0, 1, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Constraint 1 tightens an infinite lower bound so constraint 2 becomes
        upper forcing (integer arithmetic and NOT all bounds are tightened).

    min. [1 1 -1 -3 2]x
    s.t. [ 1  0   -1   0   0]   >= 4
         [ 2  1    0   0   0] x <= 10
         [ 1  1   -1   3  -1]   <= 7
         [ 2  0    0  -2   2]   >= 6
         x2, x4, x5 >= 0
         x3 >= 1
         x1 <= 100
*/
static char *test_6_domain()
{
    double Ax[] = {1, -1, 2, 1, 1, 1, -1, 3, -1, 2, -2, 2};
    int Ai[] = {0, 2, 0, 1, 0, 1, 2, 3, 4, 0, 3, 4};
    int Ap[] = {0, 2, 4, 9, 12};
    int nnz = 12;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {4, -INF, -INF, 6};
    double rhs[] = {INF, 10, 7, INF};
    double lbs[5] = {-INF, 0, 1, 0, 0};
    double ubs[] = {100, INF, INF, INF, INF};
    double c[] = {1, 1, -1, -3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    // check that new A is correct
    double Ax_correct[] = {-1, -1, 3, -1, -2, 2};
    int Ai_correct[] = {0, 0, 1, 2, 1, 2};
    int Ap_correct[] = {0, 1, 4, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new variable bounds are correct
    double lbs_correct[] = {1, 0, 0};
    double ubs_correct[] = {96, INF, INF};
    mu_assert("error bounds",
              check_bounds(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that objective offset is correct
    mu_assert("error offset", prob->obj->offset == 5);

    // check that lhs and rhs are correct
    double lhs_correct[] = {-1, -INF, -4};
    double rhs_correct[] = {INF, 2, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Test 5 with two rows swapped.

    min. [1 1 -1 -3 2]x
    s.t. [ 2  1    0   0   0] x >= 10
         [ 1  0    1   0   0]   <= 4
         [ 1  1   -1   3  -1]   <= 7
         [ 2  0    0  -2   2]   >= 6
         x1, x2, x3, x4, x5 >= 0
         x2 <= 2
*/
static char *test_7_domain()
{
    double Ax[] = {2, 1, 1, 1, 1, 1, -1, 3, -1, 2, -2, 2};
    int Ai[] = {0, 1, 0, 2, 0, 1, 2, 3, 4, 0, 3, 4};
    int Ap[] = {0, 2, 4, 9, 12};
    int nnz = 12;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {10, -INF, -INF, 6};
    double rhs[] = {INF, 4, 7, INF};
    double lbs[5] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, 2, INF, INF, INF};
    double c[] = {1, 1, -1, -3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);

    problem_clean(prob);

    // check that new A is correct
    double Ax_correct[] = {1, 1, 3, -1, -2, 2};
    int Ai_correct[] = {0, 0, 1, 2, 1, 2};
    int Ap_correct[] = {0, 1, 4, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0};
    double ubs_correct[] = {2, INF, INF};
    mu_assert("error bounds",
              check_bounds(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that objective offset is correct
    mu_assert("error offset", prob->obj->offset == 4);

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, -INF, -2};
    double rhs_correct[] = {INF, 3, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Test 6 two rows swapped

    min. [1 1 -1 -3 2]x
    s.t. [ 2  1    0   0   0] x <= 10
         [ 1  0   -1   0   0]   >= 4
         [ 1  1   -1   3  -1]   <= 7
         [ 2  0    0  -2   2]   >= 6
         x2, x4, x5 >= 0
         x3 >= 1
         x1 <= 100
*/
static char *test_8_domain()
{
    double Ax[] = {2, 1, 1, -1, 1, 1, -1, 3, -1, 2, -2, 2};
    int Ai[] = {0, 1, 0, 2, 0, 1, 2, 3, 4, 0, 3, 4};
    int Ap[] = {0, 2, 4, 9, 12};
    int nnz = 12;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {-INF, 4, -INF, 6};
    double rhs[] = {10, INF, 7, INF};
    double lbs[5] = {-INF, 0, 1, 0, 0};
    double ubs[] = {100, INF, INF, INF, INF};
    double c[] = {1, 1, -1, -3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    // check that new A is correct
    double Ax_correct[] = {-1, -1, 3, -1, -2, 2};
    int Ai_correct[] = {0, 0, 1, 2, 1, 2};
    int Ap_correct[] = {0, 1, 4, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new variable bounds are correct
    double lbs_correct[] = {1, 0, 0};
    double ubs_correct[] = {96, INF, INF};
    mu_assert("error bounds",
              check_bounds(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that objective offset is correct
    mu_assert("error offset", prob->obj->offset == 5);

    // check that lhs and rhs are correct
    double lhs_correct[] = {-1, -INF, -4};
    double rhs_correct[] = {INF, 2, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Confirm that we can derive integer bounds and then deduce that the variables
   are implied free*/
static char *test_9_domain_integer()
{
    double Ax[] = {2, 3, 1, 8};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 2;

    double lhs[] = {-INF, -INF};
    double rhs[] = {8, 8};
    double lbs[5] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {1, 1};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    mu_assert("error bound", constraints->bounds[0].ub == 4);
    mu_assert("error bound", constraints->bounds[1].ub == 1);
    mu_assert("error bound", !HAS_TAG(constraints->col_tags[0], C_TAG_UB_INF));
    mu_assert("error bound", !HAS_TAG(constraints->col_tags[1], C_TAG_UB_INF));

    remove_redundant_bounds(constraints);

    mu_assert("error bound", HAS_TAG(constraints->col_tags[0], C_TAG_UB_INF));
    mu_assert("error bound", HAS_TAG(constraints->col_tags[1], C_TAG_UB_INF));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Confirm that we can derive decimal bounds and then deduce that the
   variables are implied free*/
static char *test_9_domain_decimal()
{
    double Ax[] = {3, 3, 1, 7};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 2;

    double lhs[] = {-INF, -INF};
    double rhs[] = {8, 8};
    double lbs[5] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {1, 1};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = propagate_primal(prob, true);
    mu_assert("error status", status != INFEASIBLE);
    problem_clean(prob);

    // mu_assert("error bound", IS_EQUAL_FEAS_TOL(constraints->bounds[0].ub, 8
    // / 3.0)); mu_assert("error bound", IS_EQUAL_FEAS_TOL(constraints->bounds[1].ub,
    // 8 / 7.0));
    mu_assert("error bound", !HAS_TAG(constraints->col_tags[0], C_TAG_UB_INF));
    mu_assert("error bound", !HAS_TAG(constraints->col_tags[1], C_TAG_UB_INF));

    remove_redundant_bounds(constraints);

    mu_assert("error bound", HAS_TAG(constraints->col_tags[0], C_TAG_UB_INF));
    mu_assert("error bound", HAS_TAG(constraints->col_tags[1], C_TAG_UB_INF));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Rounding of large terms. The row
       1000 (y1 + y2 + y3 + y4) + 1e-4 x <= 7169425.2401
   with y_i >= (1421.21481, 2107.28224, 1892.50675, 1748.42144), y_i <= 1e4 and
   x in [1, 100] implies x <= 1 (plus 2e-7). The four products are about 1e6
   and their rounding, about 1e-9, is divided by the coefficient 1e-4 of x: the
   computed implied bound is 1 - 9.3e-6, below the lower bound by more than
   the verdict tolerance. Without a margin sized to the terms, propagation
   returns INFEASIBLE on a feasible LP (x = 1, y = y0). */
static char *test_10_domain_rounding_large_terms()
{
    double Ax[] = {1e3, 1e3, 1e3, 1e3, 1e-4};
    int Ai[] = {0, 1, 2, 3, 4};
    int Ap[] = {0, 5};
    double lhs[] = {-INF};
    double rhs[] = {7169425.2401};
    double lbs[] = {1421.21481, 2107.28224, 1892.50675, 1748.42144, 1.0};
    double ubs[] = {1e4, 1e4, 1e4, 1e4, 100.0};
    double c[] = {-1, -1, -1, -1, -1};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 1, 5, 5, lhs, rhs, lbs, ubs, c, stgs);
    PresolveStatus status = propagate_primal(presolver->prob, true);
    mu_assert("feasible", status != INFEASIBLE);
    mu_assert("x <= 1 not cut off",
              presolver->prob->constraints->bounds[4].ub >= 1.0 - 1e-9);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Rounding under cancellation. The row
       1000 y - 1000 z + 1e-4 x <= 11.1001
   with y in [5000.0372, 1e4], z in [0, 5000.0261] and x in [1, 100] implies
   x <= 1 (plus 3e-6). The two terms are about 5e6 while the activity is 11.1,
   so a margin relative to the activity or to the bound misses the rounding of
   the terms (1e-9, divided by 1e-4). Feasible with x = 1. */
static char *test_11_domain_rounding_cancellation()
{
    double Ax[] = {1e3, -1e3, 1e-4};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    double lhs[] = {-INF};
    double rhs[] = {11.1001};
    double lbs[] = {5000.0372, 0.0, 1.0};
    double ubs[] = {1e4, 5000.0261, 100.0};
    double c[] = {-1, 1, -1};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 1, 3, 3, lhs, rhs, lbs, ubs, c, stgs);
    PresolveStatus status = propagate_primal(presolver->prob, true);
    mu_assert("feasible", status != INFEASIBLE);
    mu_assert("x <= 1 not cut off",
              presolver->prob->constraints->bounds[2].ub >= 1.0 - 1e-9);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* A huge upper bound that propagation tightens twice leaves a rounding residue
   in the cached minimum activity of row 0.
       r0:  x1 + x2 - 24 y <= 0
       r1:  y - w <= 0
       r2:  y - v <= 0
       r3:  v + u <= 5
       r4:  x1 + 2 x2 <= 10
   x1, x2, v >= 0, y in [0, 1e9], w in [0, 0.3], u in [5, 6], min -x1 - x2.
   Feasible with optimal value 0 (everything 0 except u = 5). The upper bound
   of y goes 1e9 -> 0.3 -> 0; with incremental updates the cached minimum
   activity of r0 goes -2.4e10 -> -7.1999988555908203 -> 1.14e-06 instead of
   -7.2 -> 0, and x1 <= -1.14e-06 is derived against x1 >= 0. */
static char *test_12_domain_huge_bound_drift()
{
    double Ax[] = {1, 1, -24, 1, -1, 1, -1, 1, 1, 1, 2};
    int Ai[] = {0, 1, 2, 2, 3, 2, 4, 4, 5, 0, 1};
    int Ap[] = {0, 3, 5, 7, 9, 11};
    double lhs[] = {-INF, -INF, -INF, -INF, -INF};
    double rhs[] = {0, 0, 0, 5, 10};
    double lbs[] = {0, 0, 0, 0, 0, 5};
    double ubs[] = {INF, INF, 1e9, 0.3, INF, 6};
    double c[] = {-1, -1, 0, 0, 0, 0};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 5, 6, 11, lhs, rhs, lbs, ubs, c, stgs);
    PresolveStatus status = run_presolver(presolver);
    mu_assert("feasible", (status & INFEASIBLE) == 0);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_domain()
{
    mu_run_test(test_1_domain, counter_domain);
    mu_run_test(test_2_domain, counter_domain);
    mu_run_test(test_3_domain, counter_domain);
    mu_run_test(test_4_domain, counter_domain);
    mu_run_test(test_5_domain, counter_domain);
    mu_run_test(test_6_domain, counter_domain);
    mu_run_test(test_7_domain, counter_domain);
    mu_run_test(test_8_domain, counter_domain);
    mu_run_test(test_9_domain_integer, counter_domain);
    mu_run_test(test_9_domain_decimal, counter_domain);
    mu_run_test(test_10_domain_rounding_large_terms, counter_domain);
    mu_run_test(test_11_domain_rounding_cancellation, counter_domain);
    mu_run_test(test_12_domain_huge_bound_drift, counter_domain);
    return 0;
}

int test_domain()
{
    const char *result = all_tests_domain();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("domain: TEST FAILED!\n");
    }
    else
    {
        printf("domain: ALL TESTS PASSED\n");
    }
    printf("domain: Tests run: %d\n", counter_domain);
    return result == 0;
}

#endif
