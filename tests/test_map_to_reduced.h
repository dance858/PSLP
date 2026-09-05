#ifndef TEST_MAP_TO_REDUCED_H
#define TEST_MAP_TO_REDUCED_H

#include "Memory_wrapper.h"
#include "Numerics.h"
#include "PSLP_API.h"
#include "PSLP_sol.h"
#include "minunit.h"
#include "test_macros.h"
#include <assert.h>
#include <stdio.h>

#define MAP_TOL 1e-6
#define MAP_MAX_DIM 16
static int counter_map_to_reduced = 0;

/* Maps the primal-dual point (x, y) of the original problem to the reduced
   problem and compares the result with (x_red, y_red, z_red). */
static bool check_map_to_reduced(Presolver *presolver, const double *x,
                                 const double *y, const double *x_red_correct,
                                 const double *y_red_correct,
                                 const double *z_red_correct)
{
    size_t n = presolver->reduced_prob->n;
    size_t m = presolver->reduced_prob->m;
    double x_red[MAP_MAX_DIM], y_red[MAP_MAX_DIM], z_red[MAP_MAX_DIM];
    assert(n <= MAP_MAX_DIM && m <= MAP_MAX_DIM);
    map_solution_to_reduced(presolver, x, y, x_red, y_red, z_red);
    return is_solution_correct(x_red, x_red_correct, y_red, y_red_correct, z_red,
                               z_red_correct, (int) m, (int) n, MAP_TOL);
}

/* Maps (x, y) to the reduced problem, postsolves the result, and checks that
   the original point is recovered. */
static bool check_round_trip(Presolver *presolver, const double *x, const double *y,
                             const double *z, size_t n_rows, size_t n_cols)
{
    size_t n = presolver->reduced_prob->n;
    size_t m = presolver->reduced_prob->m;
    double x_red[MAP_MAX_DIM], y_red[MAP_MAX_DIM], z_red[MAP_MAX_DIM];
    assert(n <= MAP_MAX_DIM && m <= MAP_MAX_DIM);
    map_solution_to_reduced(presolver, x, y, x_red, y_red, z_red);
    postsolve(presolver, x_red, y_red, z_red);
    return is_solution_correct(presolver->sol->x, x, presolver->sol->y, y,
                               presolver->sol->z, z, (int) n_rows, (int) n_cols,
                               MAP_TOL);
}

/*  Singleton equality row, singleton columns and primal propagation
    (same problem as test_0_postsolve). */
static char *test_map_ston_row()
{
    double Ax[] = {0.1, 0.2, -0.3, -0.1, 0.4,  -0.4, 0.1,  0.3, 0,
                   1,   0.3, 0.5,  -0.6, -0.3, 0.2,  -0.2, 0.4, 0.1};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 0, 1, 2, 3, 4, 5, 6, 7};
    int Ap[] = {0, 8, 10, 18};
    int nnz = 18;
    int n_rows = 3;
    int n_cols = 8;

    double lhs[] = {0.5, 3, -2};
    double rhs[] = {0.5, 3, 2};
    double lbs[] = {-10, -10, -10, -10, -10, -10, -10, -10};
    double ubs[] = {10, 10, 10, 10, 10, 10, 10, 10};
    double c[] = {-1, 1, -1, 1, -1, 2, -0.2, 0.7};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);

    // optimal solution to the reduced problem
    double x_red[] = {10., -10., -10., 2.71428571, -10., -6.85714286, -10.};
    double y_red[] = {-2.57142857, 0.14285714};
    double z_red[] = {-0.78571429, 1.44285714, 0.78571429, 0., 1., 0., 1.45714286};

    // optimal solution to the original problem
    double x[] = {10., -10., 3., -10., 2.71428571, -10., -6.85714286, -10.};
    double y[] = {-2.57142857, -1.68571429, 0.14285714};
    double z[] = {-0.78571429, 1.44285714, 0., 0.78571429, 0., 1., 0., 1.45714286};

    mu_assert("map ston row error",
              check_map_to_reduced(presolver, x, y, x_red, y_red, z_red));
    mu_assert("map ston row round trip error",
              check_round_trip(presolver, x, y, z, n_rows, n_cols));

    // mapping must not touch the solution obtained by postsolve above
    double x_only[7], y_only[2];
    map_solution_to_reduced(presolver, x, y, x_only, y_only, NULL);
    mu_assert("map ston row must not modify presolver->sol",
              is_solution_correct(presolver->sol->x, x, presolver->sol->y, y,
                                  presolver->sol->z, z, n_rows, n_cols, MAP_TOL));

    // the primal and dual parts can be mapped independently
    map_solution_to_reduced(presolver, x, NULL, x_only, NULL, NULL);
    map_solution_to_reduced(presolver, NULL, y, NULL, y_only, NULL);
    for (int i = 0; i < 7; ++i)
    {
        mu_assert("map ston row primal only error",
                  ABS(x_only[i] - x_red[i]) <= MAP_TOL);
    }
    for (int i = 0; i < 2; ++i)
    {
        mu_assert("map ston row dual only error",
                  ABS(y_only[i] - y_red[i]) <= MAP_TOL);
    }

    // a missing input leaves the corresponding output untouched
    double untouched[7] = {DUMMY_VALUE, DUMMY_VALUE, DUMMY_VALUE, DUMMY_VALUE,
                           DUMMY_VALUE, DUMMY_VALUE, DUMMY_VALUE};
    map_solution_to_reduced(presolver, NULL, y, untouched, y_only, NULL);
    map_solution_to_reduced(presolver, x, NULL, x_only, untouched, NULL);
    for (int i = 0; i < 7; ++i)
    {
        mu_assert("map ston row missing input error", untouched[i] == DUMMY_VALUE);
    }

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Singleton equality row and doubleton equality row
    (same problem as test_1_postsolve). */
static char *test_map_dton_eq_row()
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
    double lbs[] = {-10, -10, -10, -10, -10, -10, -10, -10};
    double ubs[] = {10, 10, 10, 10, 10, 10, 10, 10};
    double c[] = {-1, 1, -1, 1, -1, 2, -0.2, 0.7};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);

    double x_red[] = {10, -10, 0.5333333, -10, 0.6666666, -10};
    double y_red[] = {-4.66666667, 0.66666667};
    double z_red[] = {-0.73333333, 1.6, 0., 0.26666667, 0., 2.03333333};

    double x[] = {10, -10, 3, 0.5333333, 3.4666666, -10, 0.6666666, -10};
    double y[] = {-4.66666667, -2., 0.73333333, 0.66666667};
    double z[] = {-0.73333333, 1.6, 0., 0., 0., 0.26666667, 0., 2.03333333};

    mu_assert("map dton eq row error",
              check_map_to_reduced(presolver, x, y, x_red, y_red, z_red));
    mu_assert("map dton eq row round trip error",
              check_round_trip(presolver, x, y, z, n_rows, n_cols));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Two free column singletons in different equality constraints
   (same problem as test_2_postsolve). */
static char *test_map_free_col_singletons()
{
    double Ax[] = {-1, 2, 3, 4, -1, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    int Ai[] = {0, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4};
    int Ap[] = {0, 4, 8, 11, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 5, 6, 8};
    double rhs[] = {2, 5, 6, 8};
    double lbs[] = {-INF, -INF, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {1, 1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);

    double x_red[] = {0.33333333, 0.0, 0.33333333};
    double y_red[] = {10.83333333, -7.33333333};
    double z_red[] = {0.0, 2.5, 0.0};

    double x[] = {0.0, -1.0, 0.33333333, 0.0, 0.33333333};
    double y[] = {-1.0, -1.0, 10.83333333, -7.33333333};
    double z[] = {0.0, 0.0, 0.0, 2.5, 0.0};

    mu_assert("map free col singletons error",
              check_map_to_reduced(presolver, x, y, x_red, y_red, z_red));
    mu_assert("map free col singletons round trip error",
              check_round_trip(presolver, x, y, z, n_rows, n_cols));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Implied free column singleton in an inequality row, which turns an equality
   row into an inequality (same problem as
   test_implied_free_col_ston_in_inequality_postsolve). */
static char *test_map_implied_free_col_singleton()
{
    double Ax[] = {3.0, 5, 2, -1, 1, -3, -2, -2, 1};
    int Ai[] = {0, 2, 3, 0, 1, 3, 0, 2, 3};
    int Ap[] = {0, 3, 6, 9};
    int nnz = 9;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {-INF, 8, -INF};
    double rhs[] = {5.0, 10, 6};
    double lbs[] = {0.0, 8, 0.0, 0.0};
    double ubs[] = {INF, INF, INF, INF};
    double c[] = {1.0, 2, -4, -3};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_cols = false;
    stgs->dual_fix = false;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);

    double x_red[] = {0., 1., 0.};
    double y_red[] = {-0.8, 0.};
    double z_red[] = {5.4, 0., 4.6};

    double x[] = {0., 8., 1., 0.};
    double y[] = {-0.8, 2., 0.};
    double z[] = {5.4, 0., 0., 4.6};

    mu_assert("map implied free col singleton error",
              check_map_to_reduced(presolver, x, y, x_red, y_red, z_red));
    mu_assert("map implied free col singleton round trip error",
              check_round_trip(presolver, x, y, z, n_rows, n_cols));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Doubleton row used to substitute a free variable
    (same problem as test_8_postsolve). */
static char *test_map_dton_free_var()
{
    double Ax[] = {2, 1, -1, 3, 2, 1, -2, 1, 4, 3, 1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 3, 4};
    int Ap[] = {0, 5, 7, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -1, 4};
    double rhs[] = {2, -1, 4};
    double lbs[] = {-INF, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->ston_cols = false;
    stgs->primal_propagation = false;
    stgs->dual_fix = false;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);

    double x_red[] = {0.83333333, 0.16666667, 0., 0.};
    double y_red[] = {1., -0.33333333};
    double z_red[] = {0., 0., 1., 0.33333333};

    double x[] = {0.66666667, 0.83333333, 0.16666667, 0., 0.};
    double y[] = {1., 0.33333333, -0.33333333};
    double z[] = {0., 0., 0., 1., 0.33333333};

    mu_assert("map dton free var error",
              check_map_to_reduced(presolver, x, y, x_red, y_red, z_red));
    mu_assert("map dton free var round trip error",
              check_round_trip(presolver, x, y, z, n_rows, n_cols));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Parallel rows. The multiplier of a removed row is transferred to the row
    that remains, whether the removed row tightened a side of the remaining row
    or was simply redundant.
    min. [1 1.1 1.2]
    s.t.     [ 1 1 1] >= 1
             [ 2 2 2] >= 4
             [ -3 -3 -3]  <= -9
                0 <= x1, x2, x3
*/
static char *test_map_parallel_rows()
{
    double Ax[] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0, -3.0, -3.0, -3.0};
    int Ai[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    int Ap[] = {0, 3, 6, 9};
    int nnz = 9;
    int n_rows = 3;
    int n_cols = 3;

    double lhs[] = {1, 4, -INF};
    double rhs[] = {INF, INF, -9};
    double lbs[] = {0, 0, 0};
    double ubs[] = {INF, INF, INF};
    double c[] = {1.0, 1.1, 1.2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_cols = false;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);
    mu_assert("map parallel rows: one row should remain",
              presolver->reduced_prob->m == 1);

    // In the original problem the multiplier can sit on any of the three
    // parallel rows. The reduced problem has a single row a_r x with
    // a_r = ratio_i * a_i, so the mapped multiplier is sum_i y_i / ratio_i.
    double a_r = presolver->reduced_prob->Ax[0];
    double y_on_row_0[] = {1.0, 0.0, 0.0};
    double y_on_row_1[] = {0.0, 0.5, 0.0};
    double y_on_row_2[] = {0.0, 0.0, -1.0 / 3.0};
    double x[] = {3., 0., 0.};
    double x_red[] = {3., 0., 0.};
    double y_red[] = {1.0 / a_r};
    double z_red[] = {0., 0.1, 0.2};

    mu_assert("map parallel rows (row 0) error",
              check_map_to_reduced(presolver, x, y_on_row_0, x_red, y_red, z_red));
    mu_assert("map parallel rows (row 1) error",
              check_map_to_reduced(presolver, x, y_on_row_1, x_red, y_red, z_red));
    mu_assert("map parallel rows (row 2) error",
              check_map_to_reduced(presolver, x, y_on_row_2, x_red, y_red, z_red));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Parallel rows where the remaining row is an equality and the removed
    inequality row is redundant, so no side of the remaining row changes and
    the multiplier transfer relies solely on the PARALLEL_ROW record.
    min. [1 1.1 1.2]
    s.t.     [ 1 1 1]  = 3
             [ 2 2 2] >= 4
                0 <= x1, x2, x3
*/
static char *test_map_parallel_rows_eq_remains()
{
    double Ax[] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    int Ai[] = {0, 1, 2, 0, 1, 2};
    int Ap[] = {0, 3, 6};
    int nnz = 6;
    int n_rows = 2;
    int n_cols = 3;

    double lhs[] = {3, 4};
    double rhs[] = {3, INF};
    double lbs[] = {0, 0, 0};
    double ubs[] = {INF, INF, INF};
    double c[] = {1.0, 1.1, 1.2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->parallel_rows = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);
    mu_assert("map parallel rows eq remains: one row should remain",
              presolver->reduced_prob->m == 1);
    mu_assert("map parallel rows eq remains: the equality row should remain",
              presolver->reduced_prob->Ax[0] == 1.0);

    // the original dual can put the multiplier on either row (or both);
    // the mapped multiplier is y_0 + y_1 / ratio with ratio = 1 / 2
    double x[] = {3., 0., 0.};
    double y_on_eq[] = {1.0, 0.0};
    double y_on_ineq[] = {0.0, 0.5};
    double y_on_both[] = {0.4, 0.3};
    double x_red[] = {3., 0., 0.};
    double y_red[] = {1.0};
    double z_red[] = {0., 0.1, 0.2};

    mu_assert("map parallel rows eq remains (eq) error",
              check_map_to_reduced(presolver, x, y_on_eq, x_red, y_red, z_red));
    mu_assert("map parallel rows eq remains (ineq) error",
              check_map_to_reduced(presolver, x, y_on_ineq, x_red, y_red, z_red));
    mu_assert("map parallel rows eq remains (both) error",
              check_map_to_reduced(presolver, x, y_on_both, x_red, y_red, z_red));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Parallel columns: (x1, x2) are merged into x1 + 2 x2
    (same problem as test_parallel_col_dual_identity). */
static char *test_map_parallel_cols()
{
    double Ax[] = {1, 2, 1, 3, 6, 1};
    int Ai[] = {0, 1, 2, 0, 1, 2};
    int Ap[] = {0, 3, 6};
    int nnz = 6;
    int n_rows = 2;
    int n_cols = 3;

    double lhs[] = {-20, -40};
    double rhs[] = {20, 40};
    double lbs[] = {0, 0, -10};
    double ubs[] = {3, 5, 10};
    double c[] = {1, 2, 0.1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->parallel_cols = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);

    double x_red[] = {7, 10};
    double y_red[] = {0.3, 0.1};
    double z_red[] = {0.4, -0.3};

    double x[] = {1, 3, 10};
    double y[] = {0.3, 0.1};

    mu_assert("map parallel cols error",
              check_map_to_reduced(presolver, x, y, x_red, y_red, z_red));

    // an infeasible point is projected onto the merged bounds [0, 13]
    double x_infeas[] = {10, 3, 10};
    double x_red_proj[] = {13, 10};
    mu_assert(
        "map parallel cols projection error",
        check_map_to_reduced(presolver, x_infeas, y, x_red_proj, y_red, z_red));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Variable fixed to +inf by dual fix; its rows are removed
   (same problem as test_fix_col_inf). */
static char *test_map_fix_col_inf()
{
    double Ax[] = {1.0, 1, 1, 2, 2, 1, 1, 1, 2, 1};
    int Ai[] = {0, 1, 2, 0, 1, 2, 0, 2, 0, 2};
    int Ap[] = {0, 3, 6, 8, 10};
    int nnz = 10;
    int n_rows = 4;
    int n_cols = 3;
    double lhs[] = {1, 3, 1, 2};
    double rhs[] = {INF, INF, INF, INF};
    double lbs[] = {-1.0, -INF, 0.5};
    double ubs[] = {INF, INF, INF};
    double c[] = {1.0, 0, 2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);

    double x_red[] = {0.75, 0.5};
    double y_red[] = {0.0, 0.5};
    double z_red[] = {0.0, 1.5};

    double x[] = {0.75, 0.5, 0.5};
    double y[] = {0.0, 0.0, 0.0, 0.5};
    double z[] = {0.0, 0.0, 1.5};

    mu_assert("map fix col inf error",
              check_map_to_reduced(presolver, x, y, x_red, y_red, z_red));
    mu_assert("map fix col inf round trip error",
              check_round_trip(presolver, x, y, z, n_rows, n_cols));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* After free_presolver_reduced_problem, x and y can still be mapped but z_red
   cannot be computed. */
static char *test_map_after_free_reduced_problem()
{
    double Ax[] = {2, 1, -1, 3, 2, 1, -2, 1, 4, 3, 1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 3, 4};
    int Ap[] = {0, 5, 7, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -1, 4};
    double rhs[] = {2, -1, 4};
    double lbs[] = {-INF, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->ston_cols = false;
    stgs->primal_propagation = false;
    stgs->dual_fix = false;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    run_presolver(presolver);
    free_presolver_reduced_problem(presolver);

    double x_red_correct[] = {0.83333333, 0.16666667, 0., 0.};
    double y_red_correct[] = {1., -0.33333333};
    double x[] = {0.66666667, 0.83333333, 0.16666667, 0., 0.};
    double y[] = {1., 0.33333333, -0.33333333};

    double x_red[4], y_red[2];
    double z_red[4] = {DUMMY_VALUE, DUMMY_VALUE, DUMMY_VALUE, DUMMY_VALUE};
    map_solution_to_reduced(presolver, x, y, x_red, y_red, z_red);

    for (int i = 0; i < 4; ++i)
    {
        mu_assert("map after free: x error",
                  ABS(x_red[i] - x_red_correct[i]) <= MAP_TOL);
        mu_assert("map after free: z should be untouched", z_red[i] == DUMMY_VALUE);
    }
    for (int i = 0; i < 2; ++i)
    {
        mu_assert("map after free: y error",
                  ABS(y_red[i] - y_red_correct[i]) <= MAP_TOL);
    }

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_map_to_reduced()
{
    mu_run_test(test_map_ston_row, counter_map_to_reduced);
    mu_run_test(test_map_dton_eq_row, counter_map_to_reduced);
    mu_run_test(test_map_free_col_singletons, counter_map_to_reduced);
    mu_run_test(test_map_implied_free_col_singleton, counter_map_to_reduced);
    mu_run_test(test_map_dton_free_var, counter_map_to_reduced);
    mu_run_test(test_map_parallel_rows, counter_map_to_reduced);
    mu_run_test(test_map_parallel_rows_eq_remains, counter_map_to_reduced);
    mu_run_test(test_map_parallel_cols, counter_map_to_reduced);
    mu_run_test(test_map_fix_col_inf, counter_map_to_reduced);
    mu_run_test(test_map_after_free_reduced_problem, counter_map_to_reduced);
    return 0;
}

static bool test_map_to_reduced()
{
    const char *result = all_tests_map_to_reduced();
    if (result != 0)
    {
        printf("map_to_reduced: %s\n", result);
    }
    else
    {
        printf("map_to_reduced: ALL TESTS PASSED\n");
    }
    printf("map_to_reduced: Tests run: %d\n", counter_map_to_reduced);

    return result == 0;
}

#endif // TEST_MAP_TO_REDUCED_H
