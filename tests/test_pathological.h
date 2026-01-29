#ifndef TEST_PATHOLOGICAL_H
#define TEST_PATHOLOGICAL_H

#include "Debugger.h"
#include "Numerics.h"
#include "PSLP_API.h"
#include "PSLP_sol.h"
#include "glbopts.h"
#include "test_macros.h"

#include "Problem.h"
#include "Workspace.h"
#include "minunit.h"
#include <stdio.h>

static int counter_pathological = 0;

#define POSTSOLVE_TOL_FEAS 1e-6

// Infeasible bounds: x1 >= 2, x1 <= 1
static const char *test_infeasible_bounds()
{
    double Ax[] = {0.1, 0.2, -0.3, -0.1, 0.4,  -0.4, 0.1,  0.3, 1,  1,
                   1,   0.3, 0.5,  -0.6, -0.3, 0.2,  -0.2, 0.4, 0.1};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 0, 1, 2, 3, 4, 5, 6, 7};
    int Ap[] = {0, 8, 9, 11, 19};
    int nnz = 19;
    int n_rows = 4;
    int n_cols = 8;

    double lhs[] = {0.5, 3, 4, -2};
    double rhs[] = {0.5, 3, 4, 2};
    double lbs[] = {2, -10, -10, -10, -10, -10, -10, -10};
    double ubs[] = {1, 10, 10, 10, 4.0 - 8.0 / 5.0, 10, 10, 10};
    double c[] = {-1, 1, -1, 1, -1, 2, -0.2, 0.7};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("infeasible bounds failed", status == INFEASIBLE);
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

// Infeasible bounds after singleton row elimination: x1 >= 2, x1 <= 1
static const char *test_infeasible_bounds_after_stonrow()
{
    double Ax[] = {1};
    int Ai[] = {0};
    int Ap[] = {0, 1};
    int nnz = 1;
    int n_rows = 1;
    int n_cols = 4;

    double lhs[] = {2};
    double rhs[] = {INF};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {1, 4, 4, 4};
    double c[] = {1, -1, 1, -1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("no linear constraints failed", status == INFEASIBLE);
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

// inconsistent parallel rows x1 + x2 <= 1 and -x1 - x2 <= -3
static const char *test_infeasible_parallel_rows()
{
    double Ax[] = {1, 1, -1, -1};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 2;

    double lhs[] = {-INF, -INF};
    double rhs[] = {1, -3};
    double lbs[] = {0, 0};
    double ubs[] = {5, 5};
    double c[] = {1, 1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_rows = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("no linear constraints failed", status == INFEASIBLE);
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

// minimize -2x1 + x2 subject to x1 + x2 <= 1, x1 <= 0, x2 free
static const char *test_unbounded_sol()
{
    double Ax[] = {1, 1};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    int nnz = 2;
    int n_rows = 1;
    int n_cols = 2;

    double lhs[] = {-INF};
    double rhs[] = {1};
    double lbs[] = {-INF, -INF};
    double ubs[] = {0, INF};
    double c[] = {-2, 1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_rows = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("no linear constraints failed", status == UNBNDORINFEAS);
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

static const char *test_no_linear_constraints_infeasible()
{
    double Ax[] = {};
    int Ai[] = {};
    int Ap[] = {};
    int nnz = 0;
    int n_rows = 0;
    int n_cols = 4;

    double lhs[] = {};
    double rhs[] = {};
    double lbs[] = {1, 2, 3, 4};
    double ubs[] = {-1, -2, -3, -4};
    double c[] = {1, -1, 1, -1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("no linear constraints failed", status == INFEASIBLE);
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*
static const char *test_no_linear_constraints_feasible()
{
    double Ax[] = {};
    int Ai[] = {};
    int Ap[] = {};
    int nnz = 0;
    int n_rows = 0;
    int n_cols = 4;

    double lhs[] = {};
    double rhs[] = {};
    double lbs[] = {-1, -2, -3, -4};
    double ubs[] = {1, 2, 3, 4};
    double c[] = {1, -1, 1, -1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("no linear constraints failed", status == REDUCED);

    // check that we can recover solution
    postsolve(presolver, NULL, NULL, NULL);

    double correct_x[] = {-1, 2, -3, 4};
    double correct_y[] = {};
    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, c, n_rows, n_cols,
                                  POSTSOLVE_TOL_FEAS));
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

static const char *test_zedongs_example()
{
    double Ax[] = {};
    int Ai[] = {};
    int Ap[] = {};
    int nnz = 0;
    int n_rows = 0;
    int n_cols = 1;

    double lhs[] = {};
    double rhs[] = {};
    double lbs[] = {-INF};
    double ubs[] = {2};
    double c[] = {-3};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("no linear constraints failed", status == REDUCED);
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

static const char *test_unbounded_Zedongs_example()
{
    double Ax[] = {1, 2, 3, 2};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 2;

    double lhs[] = {5, -INF};
    double rhs[] = {5, 8};
    double lbs[] = {-INF, -INF};
    double ubs[] = {INF, INF};
    double c[] = {1, 1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("no linear constraints failed", status == UNBNDORINFEAS);
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}
    */

static const char *all_tests_pathological()
{
    // mu_run_test(test_zedongs_example, counter_pathological);
    // mu_run_test(test_unbounded_Zedongs_example, counter_pathological);
    // mu_run_test(test_no_linear_constraints_feasible, counter_pathological);
    mu_run_test(test_unbounded_sol, counter_pathological);
    mu_run_test(test_infeasible_parallel_rows, counter_pathological);
    mu_run_test(test_infeasible_bounds_after_stonrow, counter_pathological);
    mu_run_test(test_infeasible_bounds, counter_pathological);
    mu_run_test(test_no_linear_constraints_infeasible, counter_pathological);

    return 0;
}

int test_pathological()
{
    const char *result = all_tests_pathological();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("pathological: TEST FAILED!\n");
    }
    else
    {
        printf("pathological: ALL TESTS PASSED\n");
    }
    printf("pathological: Tests run: %d\n", counter_pathological);
    return result == 0;
}

#endif // TEST_PATHOLOGICAL_H
