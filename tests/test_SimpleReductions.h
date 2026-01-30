#ifndef TEST_SIMPLEREDUCTIONS_H
#define TEST_SIMPLEREDUCTIONS_H

#include "Debugger.h"
#include "PSLP_API.h"

#include "minunit.h"
#include "test_SimpleReductions.h"
#include "test_macros.h"
#include <stdio.h>

// clang-format off
/*
The following tests are/should be implemented for the TrivialRowReductions class
(a "✓" indicates that the test has been implemented).
---- Tests based on empty / singleton rows ----
Test 1: Removal of two feasible empty rows that exist from the start.
Test 2: Removal of an empty infeasible row that exists from the start.
Test 3: Singleton equality row that is infeasible wrt bounds.
Test 4: Two singleton equality rows with the same variable. The problem is
        feasible.                                                                  (✓)
Test 5: Two singleton equality rows with the same variable. The problem is
        infeasible.
Test 6: One singleton equality row being processed before a singleton
        inequality row with the same variable. The problem is feasible.
Test 7: One singleton equality row being processed before a singleton
        inequality row with the same variable. The problem is infeasible.
Test 8: Three singleton equality rows.
Test 9: Three singleton equality rows followed by an empty feasible row.
Test 10: Removal of a singleton two-sided inequality row
---- Tests based on CHANGED (!) row activities ----
Test 11: Two-sided constraint, upper side declared as redundant                    (✓)
Test 12: Two-sided constraint, upper side declared as redundant. 
         This test should fail (it exists just because we should remember this
         special case) forcing first constraint.
 Test 13: Two-sided constraint, lower side declared as redundant                   (✓)
 Test 14: Two sided constraint, both sides declared as redundant                   (✓)
 Test 15: One sided upper constraint declared as redundant                         (✓)
 Test 16: One sided lower constraint declared as redundant                         (✓)
 Test 17: Constraint declared as infeasible due to min activity                    (✓)
 Test 18: Constraint declared as infeasible due to max activity                    (✓)
*/
// clang-format on

int counter_simple = 0;

/*  Two singleton equality rows with the same variable. The problem is
    feasible.
    min. [-1  4  5  2 -8   2] x
    s.t.    [ 7  5  2  0 -2   4]   = 7
            [ 2 -3  0  -1  3   1]   <= 9
            [ 1  0  0  0  0   0]   = 3
            [ 2  0  0  0  0   0]   = 6
            [ 0  1  0  -1  0  0]   >= -10
        x1, x2, x3, x4, x5, x6 >= 0
*/
static char *test_4_simple()
{
    double Ax[] = {7, 5, 2, -2, 4, 2, -3, -1, 3, 1, 1, 2, 1, -1};
    int Ai[] = {0, 1, 2, 4, 5, 0, 1, 3, 4, 5, 0, 0, 1, 3};
    int Ap[] = {0, 5, 10, 11, 12, 14};
    int nnz = 14;
    int n_rows = 5;
    int n_cols = 6;

    double lhs[] = {7, -INF, 3, 6, -10};
    double rhs[] = {7, 9, 3, 6, INF};
    double lbs[6] = {0};
    double ubs[] = {INF, INF, INF, INF, INF, INF};
    double c[] = {-1, 4, 5, 2, -8, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_ston_rows(prob);
    remove_empty_rows(constraints);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {5, 2, -2, 4, -3, -1, 3, 1, 1, -1};
    int Ai_correct[] = {0, 1, 3, 4, 0, 2, 3, 4, 0, 2};
    int Ap_correct[] = {0, 4, 8, 10};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 10));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 10));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that lhs and rhs are correct
    double lhs_correct[] = {-14, -INF, -10};
    double rhs_correct[] = {-14, 3, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Two-sided constraint, upper side declared as redundant
        min. [-1  4  5  2 -8] x
        s.t.            [ 1  0  0  -1   1]   = 2
                -1 <=   [ 1  2  3   0   0]   <= 6
                        [ 1  1  1  -1  -1]   >= 2
         -1 <= x1, x2, x3, x4, x5 <= 1
*/
static char *test_11_simple()
{
    double Ax[] = {1, -1, 1, 1, 2, 3, 1, 1, 1, -1, -1};
    int Ai[] = {0, 3, 4, 0, 1, 2, 0, 1, 2, 3, 4};
    int Ap[] = {0, 3, 6, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -1, 2};
    double rhs[] = {2, 6, INF};
    double lbs[5] = {-1, -1, -1, -1, -1};
    double ubs[] = {1, 1, 1, 1, 1};
    double c[] = {-1, 4, 5, 2, -8};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = check_activities(prob);
    problem_clean(prob, true);

    // check that new A is correct
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, 11));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, 11));
    mu_assert("rows", check_row_starts(A, Ap));

    // check that new rowtags are correct
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_RHS_INF, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 3));

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, -1, 2};
    double rhs_correct[] = {2, INF, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Two-sided constraint, lower side declared as redundant
        min. [-1  4  5  2 -8] x
        s.t.            [ 1  0  0  -1   1]   = 2
                -6 <=   [ 1  2  3   0   0]   <= 2
                        [ 1  1  1  -1  -1]   >= 2
         -1 <= x1, x2, x3, x4, x5 <= 1
*/
static char *test_13_simple()
{
    double Ax[] = {1, -1, 1, 1, 2, 3, 1, 1, 1, -1, -1};
    int Ai[] = {0, 3, 4, 0, 1, 2, 0, 1, 2, 3, 4};
    int Ap[] = {0, 3, 6, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -6, 2};
    double rhs[] = {2, 2, INF};
    double lbs[5] = {-1, -1, -1, -1, -1};
    double ubs[] = {1, 1, 1, 1, 1};
    double c[] = {-1, 4, 5, 2, -8};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = check_activities(prob);
    problem_clean(prob, true);

    // check that new A is correct
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, 11));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, 11));
    mu_assert("rows", check_row_starts(A, Ap));

    // check that new rowtags are correct
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_LHS_INF, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 3));

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, -INF, 2};
    double rhs_correct[] = {2, 2, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Two-sided constraint, both sides declared as redundant
        min. [-1  4  5  2 -8] x
        s.t.            [ 1  0  0  -1   1]   = 2
                -6 <=   [ 1  2  3   0   0]   <= 6
                        [ 1  1  1  -1  -1]   >= 2
         -1 <= x1, x2, x3, x4, x5 <= 1
*/
static char *test_14_simple()
{
    double Ax[] = {1, -1, 1, 1, 2, 3, 1, 1, 1, -1, -1};
    int Ai[] = {0, 3, 4, 0, 1, 2, 0, 1, 2, 3, 4};
    int Ap[] = {0, 3, 6, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -6, 2};
    double rhs[] = {2, 6, INF};
    double lbs[5] = {-1, -1, -1, -1, -1};
    double ubs[] = {1, 1, 1, 1, 1};
    double c[] = {-1, 4, 5, 2, -8};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = check_activities(prob);
    // delete_inactive_rows(constraints, constraints->state->rows_to_delete);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {1, -1, 1, 1, 1, 1, -1, -1};
    int Ai_correct[] = {0, 3, 4, 0, 1, 2, 3, 4};
    int Ap_correct[] = {0, 3, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new rowtags are correct
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 2));

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, 2};
    double rhs_correct[] = {2, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  One sided upper constraint declared as redundant
        min. [-1  4  5  2 -8] x
        s.t.            [ 1  0  0  -1   1]   = 2
                        [ 1  2  3   0   0]   <= 6
                        [ 1  1  1  -1  -1]   >= 2
         -1 <= x1, x2, x3, x4, x5 <= 1
*/
static char *test_15_simple()
{
    double Ax[] = {1, -1, 1, 1, 2, 3, 1, 1, 1, -1, -1};
    int Ai[] = {0, 3, 4, 0, 1, 2, 0, 1, 2, 3, 4};
    int Ap[] = {0, 3, 6, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -INF, 2};
    double rhs[] = {2, 6, INF};
    double lbs[5] = {-1, -1, -1, -1, -1};
    double ubs[] = {1, 1, 1, 1, 1};
    double c[] = {-1, 4, 5, 2, -8};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = check_activities(prob);
    // delete_inactive_rows(constraints, constraints->state->rows_to_delete);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {1, -1, 1, 1, 1, 1, -1, -1};
    int Ai_correct[] = {0, 3, 4, 0, 1, 2, 3, 4};
    int Ap_correct[] = {0, 3, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new rowtags are correct
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 2));

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, 2};
    double rhs_correct[] = {2, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  One sided lower constraint declared as redundant
        min. [-1  4  5  2 -8] x
        s.t.            [ 1  0  0  -1   1]   = 2
                 -6 <=  [ 1  2  3   0   0]
                        [ 1  1  1  -1  -1]   >= 2
         -1 <= x1, x2, x3, x4, x5 <= 1
*/
static char *test_16_simple()
{
    double Ax[] = {1, -1, 1, 1, 2, 3, 1, 1, 1, -1, -1};
    int Ai[] = {0, 3, 4, 0, 1, 2, 0, 1, 2, 3, 4};
    int Ap[] = {0, 3, 6, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -6, 2};
    double rhs[] = {2, INF, INF};
    double lbs[5] = {-1, -1, -1, -1, -1};
    double ubs[] = {1, 1, 1, 1, 1};
    double c[] = {-1, 4, 5, 2, -8};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = check_activities(prob);
    // delete_inactive_rows(constraints, constraints->state->rows_to_delete);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {1, -1, 1, 1, 1, 1, -1, -1};
    int Ai_correct[] = {0, 3, 4, 0, 1, 2, 3, 4};
    int Ap_correct[] = {0, 3, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    mu_assert("rows", check_row_starts(A, Ap_correct));

    // check that new rowtags are correct
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 2));

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, 2};
    double rhs_correct[] = {2, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Constraint declared as infeasible due to min activity
        min. [-1  4  5  2 -8] x
        s.t.            [ 1  0  0  -1   1]    = 3
                        [ 1  2  3   0   0]   <= 5.99999
                        [ 1  1  1  -1  -1]   >= 2
         1 <= x1, x2, x3, x4, x5 <= 4
*/
static char *test_17_simple()
{
    double Ax[] = {1, -1, 1, 1, 2, 3, 1, 1, 1, -1, -1};
    int Ai[] = {0, 3, 4, 0, 1, 2, 0, 1, 2, 3, 4};
    int Ap[] = {0, 3, 6, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -INF, 2};
    double rhs[] = {2, 5.9999, INF};
    double lbs[5] = {1, 1, 1, 1, 1};
    double ubs[] = {4, 4, 4, 4, 4};
    double c[] = {-1, 4, 5, 2, -8};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = check_activities(prob);
    mu_assert("error status", status == INFEASIBLE);
    // delete_inactive_rows(constraints, constraints->state->rows_to_delete);
    problem_clean(prob, true);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Constraint declared as infeasible due to max activity
        min. [-1  4  5  2 -8] x
        s.t. [ 1  0  0  -1   1]    = 3
             [ 1  2  3   0   0]   >= 6.0001
             [ 1  1  1  -1  -1]   >= 2
         -1 <= x1, x2, x3, x4, x5 <= 1
*/
static char *test_18_simple()
{
    double Ax[] = {1, -1, 1, 1, 2, 3, 1, 1, 1, -1, -1};
    int Ai[] = {0, 3, 4, 0, 1, 2, 0, 1, 2, 3, 4};
    int Ap[] = {0, 3, 6, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, 6.00001, 2};
    double rhs[] = {2, INF, INF};
    double lbs[5] = {-1, -1, -1, -1, -1};
    double ubs[] = {1, 1, 1, 1, 1};
    double c[] = {-1, 4, 5, 2, -8};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = check_activities(prob);
    mu_assert("error status", status == INFEASIBLE);
    // delete_inactive_rows(constraints, constraints->state->rows_to_delete);
    problem_clean(prob, true);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_simple()
{
    mu_run_test(test_4_simple, counter_simple);
    mu_run_test(test_11_simple, counter_simple);
    mu_run_test(test_13_simple, counter_simple);
    mu_run_test(test_14_simple, counter_simple);
    mu_run_test(test_15_simple, counter_simple);
    mu_run_test(test_16_simple, counter_simple);
    mu_run_test(test_17_simple, counter_simple);
    mu_run_test(test_18_simple, counter_simple);

    return 0;
}

int test_simple()
{
    const char *result = all_tests_simple();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("SimpleReductions: TEST FAILED!\n");
    }
    else
    {
        printf("SimpleReductions: ALL TESTS PASSED\n");
    }
    printf("SimpleReductions: Tests run: %d\n", counter_simple);
    return result == 0;
}

#endif // TEST_SIMPLEREDUCTIONS_H
