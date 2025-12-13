#ifndef TEST_ston_H
#define TEST_ston_H

#include "Debugger.h"
#include "PSLP_API.h"

#include "StonCols.h"
#include "debug_macros.h"
#include "minunit.h"
#include <stdio.h>

/*
# The following tests are/should be implemented for the StonCols class
# (a "✓" indicates that the test has been implemented).
# ---- Tests for equality column singletons -----
# Test 1: Free column singleton in equality constraint                       (✓)
# Test 2: Implied free column singleton in equality constraint               (✓)
# Test 3: Two free column singletons in two different equality constraint    (✓)
# Test 4: One free column singleton and one non-free column singleton in     (✓)
#         the same equality constraint (very interesting!)
# Test 5: Column singleton unbounded from below in equality constraint       (✓)
#         (negative coefficient)
# Test 6: Column singleton unbounded from above in equality constraint       (✓)
#         (negative coefficient)
# Test 7: Column singleton unbounded from below in equality constraint       (✓)
#         (positive coefficient)
# Test 8: Column singleton unbounded from above in equality constraint       (✓)
#         (positive coefficient)
# Test 9: Chain of column singletons                                         (✓)
# ---- Tests for inequality column singletons ----
# Test 10: Implied free column singleton in two-sided inequality constraint  (✓)
# Test 12: Dual fix to lower bound                                           (✓)
# Test 13: Dual fix to upper bound                                           (✓)
# Test 14: Upper part of a constraint is active at an optimal solution       (✓)
# Test 15: Lower part of a constraint is active at an optimal solution       (✓)
*/

static int counter_ston = 0;

/* Free column singleton in equality constraint
    min. [1 1 -1 3 2]x
    s.t. [2 1  0 1 0]     [2]
         [1 4 -2 6 8] x = [1]
         [1 4  0 3 0]     [4]
         [2 0  0 0 1]    >= 6
         x1, x2, x5 >= 0
         x3, x4 free
*/
static char *test_01_ston()
{
    double Ax[] = {2, 1, 1, 1, 4, -2, 6, 8, 1, 4, 3, 2, 1};
    int Ai[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 1, 3, 0, 4};
    int Ap[] = {0, 3, 8, 11, 13};
    int nnz = 13;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 1, 4, 6};
    double rhs[] = {2, 1, 4, INF};
    double lbs[] = {0, 0, -INF, -INF, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {1, 1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, 1, 1, 1, 4, 3, 2, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 0, 3};
    int Ap_correct[] = {0, 3, 6, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_EQ, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 3));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF,
                                C_TAG_UB_INF | C_TAG_LB_INF, C_TAG_UB_INF};

    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 4));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, -INF, 0};
    double ubs_correct[] = {INF, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {0.5, -1, 0, -2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == 0.5);

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, 4, 6};
    double rhs_correct[] = {2, 4, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* Implied free column singleton in equality constraint

        min. [1 1 -1 3 2]x
        s.t. [2 1  0 1 1]     [2]
             [1 4 -2 6 8] x = [0]
             [1 4  0 3 1]     [4]
             x1, x2, x3, x4, x5 >= 0
*/
static char *test_02_ston()
{
    double Ax[] = {2, 1, 1, 1, 1, 4, -2, 6, 8, 1, 4, 3, 1};
    int Ai[] = {0, 1, 3, 4, 0, 1, 2, 3, 4, 0, 1, 3, 4};
    int Ap[] = {0, 4, 9, 13};
    int nnz = 13;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, 0, 4};
    double rhs[] = {2, 0, 4};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {1, 1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, 1, 1, 1, 1, 4, 3, 1};
    int Ai_correct[] = {0, 1, 2, 3, 0, 1, 2, 3};
    int Ap_correct[] = {0, 4, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_EQ};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 2));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF,
                                C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 4));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0, 0};
    double ubs_correct[] = {INF, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {0.5, -1, 0, -2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == 0);

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, 4};
    double rhs_correct[] = {2, 4};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* Two free column singletons in different equality constraints

        min. [1 1 -1 3 2]x
        s.t.[-1,  0,  2,  3, 4  ]     [2]
            [ 0, -1,  5,  6, 7  ] x = [4]
            [ 0,  0,  8,  9, 10 ]     [6]
            [ 0,  0, 11, 12, 13 ]     [8]
             x3, x4, x5 >= 0
             x1, x2 free
*/
static char *test_03_ston()
{
    double Ax[] = {-1, 2, 3, 4, -1, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    int Ai[] = {0, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4};
    int Ap[] = {0, 4, 8, 11, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 4, 6, 8};
    double rhs[] = {2, 4, 6, 8};
    double lbs[] = {-INF, -INF, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {1, 1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {8, 9, 10, 11, 12, 13};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2};
    int Ap_correct[] = {0, 3, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_EQ};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 2));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 3));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0};
    double ubs_correct[] = {INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that the objective function is correct
    double obj_correct[] = {6, 12, 13};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 3));
    mu_assert("error obj", prob->obj->offset == -6);

    // check that lhs and rhs are correct
    double lhs_correct[] = {6, 8};
    double rhs_correct[] = {6, 8};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* One free column singleton and one non-free column singleton in
   the same equality constraint (very interesting!)

        min. [1 2 -1 3 2]x
        s.t.[-1, -1,   2,  3, 4  ]     [2]
            [ 0,  0,   5,  6, 7  ] x = [4]
            [ 0,  0,   8,  9, 10 ]     [6]
            [ 0,  0,  11, 12, 13 ]     [8]
             x3, x4, x5 >= 0
             -1 <= x2 <= 5
             x1 free
*/
static char *test_04_ston()
{
    double Ax[] = {-1, -1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    int Ai[] = {0, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4};
    int Ap[] = {0, 5, 8, 11, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 4, 6, 8};
    double rhs[] = {2, 4, 6, 8};
    double lbs[] = {-INF, -1, 0, 0, 0};
    double ubs[] = {INF, 5, INF, INF, INF};
    double c[] = {1, 2, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct (it should containt the empty column)
    double Ax_correct[] = {5, 6, 7, 8, 9, 10, 11, 12, 13};
    int Ai_correct[] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    int Ap_correct[] = {0, 3, 6, 9};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 9));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 9));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_EQ, R_TAG_EQ};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 3));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_NONE, C_TAG_UB_INF, C_TAG_UB_INF,
                                C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 4));

    // check that new variable bounds are correct
    double lbs_correct[] = {-1, 0, 0, 0};
    double ubs_correct[] = {5, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that the objective function is correct
    double obj_correct[] = {1, 1, 6, 6};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == -2);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4, 6, 8};
    double rhs_correct[] = {4, 6, 8};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* Column singleton unbounded from below in equality constraint
        (negative coefficient)

    min. [1 1 -1 3 2]x
    s.t. [2 1  0 1 0]     [2]
         [1 4 -2 6 8] x = [1]
         [1 4  0 3 1]     [4]
         [2 0  0 0 1]    >= 6
         x1, x2, x4, x5 >= 0
         x3 <= 4
*/
static char *test_05_ston()
{
    double Ax[] = {2, 1, 1, 1, 4, -2, 6, 8, 1, 4, 3, 1, 2, 1};
    int Ai[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 1, 3, 4, 0, 4};
    int Ap[] = {0, 3, 8, 12, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 1, 4, 6};
    double rhs[] = {2, 1, 4, INF};
    double lbs[] = {0, 0, -INF, 0, 0};
    double ubs[] = {INF, INF, 4, INF, INF};
    double c[] = {1, 1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, 1, 1, 1, 4, 6, 8, 1, 4, 3, 1, 2, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 3};
    int Ap_correct[] = {0, 3, 7, 11, 13};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 13));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 13));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_LHS_INF, R_TAG_EQ, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 4));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF,
                                C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 4));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0, 0};
    double ubs_correct[] = {INF, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {0.5, -1, 0, -2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == 0.5);

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, -INF, 4, 6};
    double rhs_correct[] = {2, 9, 4, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 4));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 4));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* Column singleton unbounded from above in equality constraint
        (negative coefficient)

    min. [1 1 -1 3 2]x
    s.t. [2 1  0 1 0]     [2]
         [1 4 -2 6 8] x = [1]
         [1 4  0 3 1]     [4]
         [2 0  0 0 1]    >= 6
         x1, x2, x4, x5 >= 0
         x3 >= 4
*/
static char *test_06_ston()
{
    double Ax[] = {2, 1, 1, 1, 4, -2, 6, 8, 1, 4, 3, 1, 2, 1};
    int Ai[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 1, 3, 4, 0, 4};
    int Ap[] = {0, 3, 8, 12, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 1, 4, 6};
    double rhs[] = {2, 1, 4, INF};
    double lbs[] = {0, 0, 4, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {1, 1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, 1, 1, 1, 4, 6, 8, 1, 4, 3, 1, 2, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 3};
    int Ap_correct[] = {0, 3, 7, 11, 13};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 13));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 13));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_RHS_INF, R_TAG_EQ, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 4));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF,
                                C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 4));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0, 0};
    double ubs_correct[] = {INF, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {0.5, -1, 0, -2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == 0.5);

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, 9, 4, 6};
    double rhs_correct[] = {2, INF, 4, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 4));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 4));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* Column singleton unbounded from below in equality constraint
   (positive coefficient)

    min. [1 1 -1 3 2]x
    s.t. [2 1  0 1 0]     [2]
         [1 4  2 6 8] x = [1]
         [1 4  0 3 0]     [4]
         [2 0  0 0 1]    >= 6
         x1, x2, x4 >= 0
         x5 >= -10
         x3 <= 4
*/
static char *test_07_ston()
{
    double Ax[] = {2, 1, 1, 1, 4, 2, 6, 8, 1, 4, 3, 1, 2, 1};
    int Ai[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 1, 3, 4, 0, 4};
    int Ap[] = {0, 3, 8, 12, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 1, 4, 6};
    double rhs[] = {2, 1, 4, INF};
    double lbs[] = {0, 0, -INF, 0, -10};
    double ubs[] = {INF, INF, 4, INF, INF};
    double c[] = {1, 1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, 1, 1, 1, 4, 6, 8, 1, 4, 3, 1, 2, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 3};
    int Ap_correct[] = {0, 3, 7, 11, 13};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 13));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 13));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_RHS_INF, R_TAG_EQ, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 4));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF,
                                C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 4));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0, -10};
    double ubs_correct[] = {INF, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {1.5, 3, 6, 6};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == -0.5);

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, -7, 4, 6};
    double rhs_correct[] = {2, INF, 4, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 4));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 4));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* Column singleton unbounded from above in equality constraint
        (positive coefficient)

    min. [1 1 -1 3 2]x
    s.t. [2 1  0 1  0]     [2]
         [1 4  2 6 -8] x = [1]
         [1 4  0 3  0]     [4]
         [2 0  0 0  1]    >= 2
         x1, x2, x4, x5 >= 0
         x3 >= 4
*/
static char *test_08_ston()
{
    double Ax[] = {2, 1, 1, 1, 4, 2, 6, -8, 1, 4, 3, 1, 2, 1};
    int Ai[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 1, 3, 4, 0, 4};
    int Ap[] = {0, 3, 8, 12, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 1, 4, 2};
    double rhs[] = {2, 1, 4, INF};
    double lbs[] = {0, 0, 4, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {1, 1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, 1, 1, 1, 4, 6, -8, 1, 4, 3, 1, 2, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 3};
    int Ap_correct[] = {0, 3, 7, 11, 13};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 13));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 13));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_LHS_INF, R_TAG_EQ, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 4));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF,
                                C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 4));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0, 0};
    double ubs_correct[] = {INF, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {1.5, 3, 6, -2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == -0.5);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/*  Simple chain of column singletons:

    min. [1 1 1 1 1 1 1 1 1 1 1]x

    s.t  [1 2 2 2 2 2 2 2 2 2 2]     [50]
         [0 1 3 3 3 3 3 3 3 3 3]     [40]
         [0 0 1 4 4 4 4 4 4 4 4] x = [30]
         [0 0 0 1 1 1 1 1 1 1 1]     [20]
         [0 0 0 2 1 1 1 1 1 1 1]     [10]

         x1, x2, x3 free
         x4 up to x10 >= 0,  <= [1, 2, 3, 4, 5, 6, 7]
         x11 <= 8
*/
static char *test_09_ston()
{
    double Ax[] = {1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 3, 3, 3, 3,
                   3, 3, 3, 3, 3, 1, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1,
                   1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 2, 3, 4,  5,  6, 7, 8, 9, 10, 1, 2, 3,  4, 5,
                6, 7, 8, 9, 10, 2,  3, 4, 5, 6, 7,  8, 9, 10, 3, 4,
                5, 6, 7, 8, 9,  10, 3, 4, 5, 6, 7,  8, 9, 10};
    int Ap[] = {0, 11, 21, 30, 38, 46};
    int nnz = 46;
    int n_rows = 5;
    int n_cols = 11;

    double lhs[] = {50, 40, 30, 20, 10};
    double rhs[] = {50, 40, 30, 20, 10};
    double lbs[] = {-INF, -INF, -INF, 0, 0, 0, 0, 0, 0, 0, -INF};
    double ubs[] = {INF, INF, INF, 1, 2, 3, 4, 5, 6, 7, 8};
    double c[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1};
    int Ai_correct[] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7};
    int Ap_correct[] = {0, 8, 16};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 16));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 16));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_EQ, R_TAG_EQ};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 2));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_NONE, C_TAG_NONE, C_TAG_NONE, C_TAG_NONE,
                                C_TAG_NONE, C_TAG_NONE, C_TAG_NONE, C_TAG_LB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 8));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0, 0, 0, 0, 0, -INF};
    double ubs_correct[] = {1, 2, 3, 4, 5, 6, 7, 8};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 8));

    // check that the objective function is correct
    double obj_correct[] = {-6, -6, -6, -6, -6, -6, -6, -6};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 8));
    mu_assert("error obj", prob->obj->offset == 70);

    // check that lhs and rhs are correct
    double lhs_correct[] = {20, 10};
    double rhs_correct[] = {20, 10};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/*  Implied free column singleton in two-sided inequality constraint
        min. x1 + 2x2 - 4x3 - 3x4
        s.t.       [3    0   5   2]   <= [5]
            10 >=  [-1   1   0  -3] x >= [8]
                   [-2   0   -2  1]   >= [6]
                x2 >= 8
                x1, x3, x4 >= 0
*/
static char *test_10_ston()
{
    double Ax[] = {3, 5, 2, -1, 1, -3, -2, -2, 1};
    int Ai[] = {0, 2, 3, 0, 1, 3, 0, 2, 3};
    int Ap[] = {0, 3, 6, 9};
    int nnz = 9;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {-INF, 8, 6};
    double rhs[] = {5, 10, INF};
    double lbs[] = {0, 8, 0, 0};
    double ubs[] = {INF, INF, INF, INF};
    double c[] = {1, 2, -4, -3};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {3, 5, 2, -2, -2, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2};
    int Ap_correct[] = {0, 3, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_LHS_INF, R_TAG_RHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 2));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 3));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0};
    double ubs_correct[] = {INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that the objective function is correct
    double obj_correct[] = {3, -4, 3};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 3));
    mu_assert("error obj", prob->obj->offset == 16);

    // check that lhs and rhs are correct
    double lhs_correct[] = {-INF, 6};
    double rhs_correct[] = {5, INF};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/*  Dual fix to lower bound
    min. 4x1 + x2 - 2x3 + 7x4
    s.t  3x1       - x3  - 6x4  <= 0
         3x1 - 2x2 + 5x3 - x4   >= 0
         4x1       + 3x3 + 4x4  <= 5
        x2 >= 1
        x1, x3, x4 >= 0
*/
static char *test_12_ston()
{
    double Ax[] = {3, -1, -6, 3, -2, 5, -1, 4, 3, 4};
    int Ai[] = {0, 2, 3, 0, 1, 2, 3, 0, 2, 3};
    int Ap[] = {0, 3, 7, 10};
    int nnz = 10;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {-INF, 0, -INF};
    double rhs[] = {0, INF, 5};
    double lbs[] = {0, 1, 0, 0};
    double ubs[] = {INF, INF, INF, INF};
    double c[] = {4, 1, -2, 7};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {3, -1, -6, 3, 5, -1, 4, 3, 4};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    int Ap_correct[] = {0, 3, 6, 9};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 9));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 9));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_LHS_INF, R_TAG_RHS_INF, R_TAG_LHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 3));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 3));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0};
    double ubs_correct[] = {INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that the objective function is correct
    double obj_correct[] = {4, -2, 7};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 3));
    mu_assert("error obj", prob->obj->offset == 1);

    // check that lhs and rhs are correct
    double lhs_correct[] = {-INF, 2, -INF};
    double rhs_correct[] = {0, INF, 5};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/*  Dual fix to upper bound
    min. 4x1 - x2 - 2x3 + 7x4
    s.t  3x1       - x3  - 6x4 <= 0
        -3x1 - 2x2 + 5x3 - x4  <= 0
         4x1       + 3x3 + 4x4 <= 5
        0 <= x2 <= 5
        x1, x3, x4 >= 0
*/
static char *test_13_ston()
{
    double Ax[] = {3, -1, -6, -3, -2, 5, -1, 4, 3, 4};
    int Ai[] = {0, 2, 3, 0, 1, 2, 3, 0, 2, 3};
    int Ap[] = {0, 3, 7, 10};
    int nnz = 10;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {-INF, -INF, -INF};
    double rhs[] = {0, 0, 5};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {INF, 5, INF, INF};
    double c[] = {4, -1, -2, 7};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {3, -1, -6, -3, 5, -1, 4, 3, 4};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    int Ap_correct[] = {0, 3, 6, 9};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 9));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 9));
    CHECK_ROW_STARTS(A, Ap_correct);

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check new rowtags
    RowTag rowtags_correct[] = {R_TAG_LHS_INF, R_TAG_LHS_INF, R_TAG_LHS_INF};
    mu_assert("error rowtags",
              ARRAYS_EQUAL_ROWTAG(rowtags_correct, constraints->row_tags, 3));

    // check new coltags
    ColTag coltags_correct[] = {C_TAG_UB_INF, C_TAG_UB_INF, C_TAG_UB_INF};
    mu_assert("error coltags",
              ARRAYS_EQUAL_COLTAG(coltags_correct, constraints->col_tags, 3));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0};
    double ubs_correct[] = {INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that the objective function is correct
    double obj_correct[] = {4, -2, 7};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 3));
    mu_assert("error obj", prob->obj->offset == -5);

    // check that lhs and rhs are correct
    double lhs_correct[] = {-INF, -INF, -INF};
    double rhs_correct[] = {0, 10, 5};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

static char *test_14_ston()
{
    double Ax[] = {1, 1, 2, 1, 1};
    int Ai[] = {0, 1, 2, 1, 2};
    int Ap[] = {0, 3, 5};
    int nnz = 5;
    int n_rows = 2;
    int n_cols = 3;

    double lhs[] = {1, -INF};
    double rhs[] = {6, 5};
    double lbs[] = {0, 0, 0};
    double ubs[] = {INF, INF, INF};
    double c[] = {-1, 1, 1};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {1, 2, 1, 1};
    int Ai_correct[] = {0, 1, 0, 1};
    int Ap_correct[] = {0, 2, 4};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 4));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 4));
    CHECK_ROW_STARTS(A, Ap_correct);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

static char *test_15_ston()
{
    double Ax[] = {-1, 1, 2, 1, 1};
    int Ai[] = {0, 1, 2, 1, 2};
    int Ap[] = {0, 3, 5};
    int nnz = 5;
    int n_rows = 2;
    int n_cols = 3;

    double lhs[] = {1, -INF};
    double rhs[] = {6, 5};
    double lbs[] = {0, 0, 0};
    double ubs[] = {INF, INF, INF};
    double c[] = {-1, 1, 1};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_ston_cols(prob);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {1, 2, 1, 1};
    int Ai_correct[] = {0, 1, 0, 1};
    int Ap_correct[] = {0, 2, 4};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 4));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 4));
    CHECK_ROW_STARTS(A, Ap_correct);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

static const char *all_tests_ston()
{
    mu_run_test(test_01_ston, counter_ston); // (✓)
    mu_run_test(test_02_ston, counter_ston); // (✓)
    mu_run_test(test_03_ston, counter_ston); // (✓)
    mu_run_test(test_04_ston, counter_ston); // (✓)
    mu_run_test(test_05_ston, counter_ston); // (✓)
    mu_run_test(test_06_ston, counter_ston); // (✓)
    mu_run_test(test_07_ston, counter_ston); // (✓)
    mu_run_test(test_08_ston, counter_ston); // (✓)
    mu_run_test(test_09_ston, counter_ston); // (✓)
    mu_run_test(test_10_ston, counter_ston); // (✓)
    mu_run_test(test_12_ston, counter_ston); // (✓)
    mu_run_test(test_13_ston, counter_ston); // (✓)
    mu_run_test(test_14_ston, counter_ston); // (✓)
    mu_run_test(test_15_ston, counter_ston); // (✓)
    return 0;
}

// not clear what should happen for test 14 and 15

int test_ston()
{
    const char *result = all_tests_ston();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("ston: TEST FAILED!\n");
    }
    else
    {
        printf("ston: ALL TESTS PASSED\n");
    }
    printf("ston: Tests run: %d\n", counter_ston);
    return result == 0;
}

#endif // TEST_ston_H
