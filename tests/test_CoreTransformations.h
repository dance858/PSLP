#ifndef TEST_CORETRANSFORMATIONS_H
#define TEST_CORETRANSFORMATIONS_H

#include "CoreTransformations.h"
#include "Problem.h"
#include "SimpleReductions.h"
#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>

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
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Constraints *constraints = presolver->prob->constraints;

    // fix x2 to 1
    fix_col(presolver->prob, 1, 1.0);

    delete_fixed_cols_from_problem(presolver->prob);

    // c is zero so the offset must stay zero
    mu_assert("error offset", presolver->prob->obj->offset == 0.0);

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

// the offset is updated when a column is fixed, not when it is flushed
static char *test_2_core()
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
    double c[] = {0, 3, 5, 7};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    Problem *prob = presolver->prob;
    const ColTag *col_tags = prob->constraints->col_tags;

    // fixing x1 to 2 adds 3 * 2 immediately
    mu_assert("error status", fix_col(prob, 1, 2.0) == REDUCED);
    mu_assert("error offset", prob->obj->offset == 6.0);

    // fixing x2 to 0 adds nothing
    mu_assert("error status", fix_col(prob, 2, 0.0) == REDUCED);
    mu_assert("error offset", prob->obj->offset == 6.0);

    // an infeasible fix leaves the objective and the column untouched
    mu_assert("error status", fix_col(prob, 3, 100.0) == INFEASIBLE);
    mu_assert("error offset", prob->obj->offset == 6.0);
    mu_assert("error tag", !HAS_TAG(col_tags[3], C_TAG_FIXED));
    mu_assert("error bounds", prob->constraints->bounds[3].ub == 40);

    // the flush does not add the contributions again
    delete_fixed_cols_from_problem(prob);
    mu_assert("error offset", prob->obj->offset == 6.0);

    PS_FREE(stgs);
    free_presolver(presolver);

    return 0;
}

// an empty column is fixed by remove_empty_cols, which updates the offset
static char *test_3_core()
{
    double Ax[] = {1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 3;

    double lhs[] = {-INF, -INF};
    double rhs[] = {1, 2};
    double lbs[] = {0, 0, 1};
    double ubs[] = {10, 10, 3};
    double c[] = {1, 1, 7};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    Problem *prob = presolver->prob;
    const ColTag *col_tags = prob->constraints->col_tags;

    mu_assert("error empty cols", prob->constraints->state->empty_cols->len == 1);
    mu_assert("error status", remove_empty_cols(prob) == UNCHANGED);

    // c[2] > 0 so x2 is fixed to its lower bound 1
    mu_assert("error offset", prob->obj->offset == 7.0);
    mu_assert("error tag", HAS_TAG(col_tags[2], C_TAG_FIXED));
    mu_assert("error bounds", prob->constraints->bounds[2].ub == 1.0);
    mu_assert("error bounds", prob->constraints->bounds[2].lb == 1.0);
    mu_assert("error empty cols", prob->constraints->state->empty_cols->len == 0);

    PS_FREE(stgs);
    free_presolver(presolver);

    return 0;
}

// substituting x_subst = (rhs - aij x_stay) / aik in the objective
static char *test_4_core()
{
    double *c = (double *) malloc(2 * sizeof(double));
    c[0] = 1;
    c[1] = 2;
    Objective *obj = objective_new(c);

    // stay = 0, subst = 1, aik = 2, aij = 4, rhs = 6:
    // c[0] -= (4 / 2) * 2 = 4, offset += (6 / 2) * 2 = 6
    sub_var_in_obj_dton(obj, 0, 1, 2.0, 4.0, 6.0);
    mu_assert("error c", obj->c[0] == -3.0);
    mu_assert("error c", obj->c[1] == 2.0);
    mu_assert("error offset", obj->offset == 6.0);

    free_objective(obj);

    return 0;
}

static const char *all_tests_core()
{
    mu_run_test(test_1_core, counter_core);
    mu_run_test(test_2_core, counter_core);
    mu_run_test(test_3_core, counter_core);
    mu_run_test(test_4_core, counter_core);
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
