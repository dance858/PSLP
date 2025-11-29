#ifndef TEST_CORETRANSFORMATIONS_H
#define TEST_CORETRANSFORMATIONS_H

#include "CoreTransformations.h"
#include "Problem.h"
#include "SimpleReductions.h"
#include "minunit.h"
#include <stdio.h>

static int counter_core = 0;

// test initialization, append, free
static char *test_1_core()
{
    double Ax[] = {1, -1, 1, 2, -1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 2, 0, 3, 0, 1, 2, 3};
    int Ap[] = {0, 3, 5, 9};
    int nnz = 9;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {4, 2, -INF};
    double rhs[] = {4, 2, 1};
    double lbs[] = {-1, -2, -3, -4};
    double ubs[] = {10, 20, 30, 40};
    double c[] = {0, 0, 0, 0};

    Settings *stgs = default_settings();
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Constraints *constraints = presolver->prob->constraints;

    // fix x2 to 1
    fix_col(constraints, 1, 1.0, 0);

    delete_fixed_cols_from_problem(presolver->prob);

    // check that new LB and UB are correct (after call to fixcol)
    mu_assert("error", constraints->bounds[1].lb == 1.0);
    mu_assert("error", constraints->bounds[1].ub == 1.0);

    // check LHS and RHS
    mu_assert("error lhs", constraints->lhs[0] == 5);
    mu_assert("error lhs", constraints->lhs[1] == 2);
    mu_assert("error lhs", IS_NEG_INF(constraints->lhs[2]));
    mu_assert("error rhs", constraints->rhs[0] == 5);
    mu_assert("error rhs", constraints->rhs[1] == 2);
    mu_assert("error rhs", constraints->rhs[2] == 0);

    // check activities
    Activity *act = constraints->state->activities;

    mu_assert("error act", act[0].min == -4 && act[0].max == 40);
    mu_assert("error act", act[1].min == -42 && act[1].max == 24);
    mu_assert("error act", act[2].min == -8 && act[2].max == 80);

    PS_FREE(stgs);
    free_presolver(presolver);

    return 0;
}

static const char *all_tests_core()
{
    mu_run_test(test_1_core, counter_core);
    return 0;
}

int test_core()
{
    const char *result = all_tests_core();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("core: TEST FAILED!\n");
    }
    else
    {
        printf("core: ALL TESTS PASSED\n");
    }
    printf("core: Tests run: %d\n", counter_core);
    return result == 0;
}

#endif // TEST_core_H
