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
    double *Ax = NULL;
    int *Ai = NULL;
    int Ap[] = {0};
    int nnz = 0;
    int n_rows = 0;
    int n_cols = 4;

    double *lhs = NULL;
    double *rhs = NULL;
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

static const char *test_no_linear_constraints_feasible()
{
    double *Ax = NULL;
    int *Ai = NULL;
    int Ap[] = {0};
    int nnz = 0;
    int n_rows = 0;
    int n_cols = 4;

    double *lhs = NULL;
    double *rhs = NULL;
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
    double *correct_y = NULL;
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
    double *Ax = NULL;
    int *Ai = NULL;
    int Ap[] = {0};
    int nnz = 0;
    int n_rows = 0;
    int n_cols = 1;

    double *lhs = NULL;
    double *rhs = NULL;
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

/* A feasible LP that used to be declared infeasible. With cap = 1 (row 0) and
   v = 9e-7:
       -v cap + x_k + y_k = 0  (k = 1, 2, 3),   -(y_1 + y_2 + y_3) >= 0,
   and x, y >= 0, so y = 0 and x_k = v. The doubleton rows give x_k in [0, v],
   an interval below the feasibility tolerance. Fixing x_k anywhere but at v
   leaves y_k > 0, and the three shifts add up to 2.7e-6 on the last row. */
static const char *test_derived_narrow_intervals_are_not_fixed()
{
    double v = 9e-7;
    double Ax[] = {1, -v, 1, 1, -v, 1, 1, -v, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 0, 1, 4, 0, 2, 5, 0, 3, 6, 4, 5, 6, 1, 2, 3, 0, 1, 2, 3};
    int Ap[] = {0, 1, 4, 7, 10, 13, 16, 20};
    int nnz = 20;
    int n_rows = 7;
    int n_cols = 7;

    double lhs[] = {1, 0, 0, 0, 0, -INF, -5};
    double rhs[] = {1, 0, 0, 0, INF, 10, INF};
    double lbs[] = {0, 0, 0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF, INF, INF};
    double c[] = {0, 1, 1, 1, 0, 0, 0};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("derived narrow intervals: feasible LP declared infeasible",
              status == REDUCED);
    mu_assert("derived narrow intervals: reduced problem not empty",
              presolver->reduced_prob->n == 0 && presolver->reduced_prob->m == 0);

    postsolve(presolver, NULL, NULL, NULL);
    const double *x = presolver->sol->x;
    mu_assert("derived narrow intervals: cap", x[0] == 1.0);
    for (int k = 1; k <= 3; ++k)
    {
        mu_assert("derived narrow intervals: x_k = v", ABS(x[k] - v) <= 1e-15);
        mu_assert("derived narrow intervals: y_k = 0", ABS(x[k + 3]) <= 1e-15);
    }

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/* The close-bound fixings that must survive: x0 has the input bounds
   [0, 1e-8] and is fixed by the sweep over the input, and the singleton row
   x1 >= 1 makes the bounds of x1 equal, which fixes it through the queue of
   the bound updates. */
static const char *test_close_input_bounds_and_equal_bounds_are_fixed()
{
    double Ax[] = {1, 1, 1, 1, 1, 1, -1};
    int Ai[] = {0, 2, 3, 1, 1, 2, 3};
    int Ap[] = {0, 3, 4, 7};
    int nnz = 7;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {-INF, 1, -INF};
    double rhs[] = {4, INF, 3};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {1e-8, 1, 5, 5};
    double c[] = {1, 1, -1, -2};

    // only the trivial explorers
    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("close bounds: status", status == REDUCED);
    mu_assert("close bounds: x0 and x1 not removed",
              presolver->reduced_prob->n == 2 && presolver->reduced_prob->m == 2);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/* A feasible LP in which x2 lies in an interval of width 1e-8:
       -1e-8 x0 + 1e-8 x1 + x2 = 1,   -x0 + x1 >= 0.25,   0 <= x <= 1,
   with x2 >= lb2. Fixing x2 inside the interval leaves a first row of the
   form 1e-8 (x1 - x0) = shift, which contradicts the second row, although
   x = (0, 0.25, 1 - 2.5e-9) is feasible. */
static PresolveStatus presolve_narrow_interval_lp(double lb2, bool propagation)
{
    double Ax[] = {-1e-8, 1e-8, 1, -1, 1};
    int Ai[] = {0, 1, 2, 0, 1};
    int Ap[] = {0, 3, 5};
    int nnz = 5;
    int n_rows = 2;
    int n_cols = 3;

    double lhs[] = {1, 0.25};
    double rhs[] = {1, INF};
    double lbs[] = {0, 0, lb2};
    double ubs[] = {1, 1, 1};
    double c[] = {1, 2, 0};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = propagation;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    free_presolver(presolver);
    free_settings(stgs);
    return status;
}

/* The interval is implied by the first row (x2 >= 1 - 1e-8); propagation must
   not fix x2 at its upper bound. */
static const char *test_narrow_interval_is_not_fixed()
{
    PresolveStatus status = presolve_narrow_interval_lp(0, true);
    mu_assert("feasible LP declared infeasible or unbounded",
              status == REDUCED || status == UNCHANGED);
    return 0;
}

/* The interval is given by the bounds; the close-bounds rule must not fix x2
   at the midpoint. */
static const char *test_close_bounds_are_not_fixed()
{
    PresolveStatus status = presolve_narrow_interval_lp(1 - 1e-8, false);
    mu_assert("feasible LP declared infeasible or unbounded",
              status == REDUCED || status == UNCHANGED);
    return 0;
}

static const char *all_tests_pathological()
{
    mu_run_test(test_derived_narrow_intervals_are_not_fixed, counter_pathological);
    mu_run_test(test_close_input_bounds_and_equal_bounds_are_fixed,
                counter_pathological);
    mu_run_test(test_narrow_interval_is_not_fixed, counter_pathological);
    mu_run_test(test_close_bounds_are_not_fixed, counter_pathological);
    mu_run_test(test_zedongs_example, counter_pathological);
    mu_run_test(test_unbounded_Zedongs_example, counter_pathological);
    mu_run_test(test_no_linear_constraints_feasible, counter_pathological);
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
