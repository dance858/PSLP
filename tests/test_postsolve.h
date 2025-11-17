#ifndef TEST_POSTSOLVE_H
#define TEST_POSTSOLVE_H

#include "CoreTransformations.h"
#include "Numerics.h"
#include "PSLP_API.h"

#include "SimpleReductions.h"
#include "minunit.h"
#include <stdio.h>

#define POSTSOLVE_TOL_FEAS 1e-6
static int counter_postsolve = 0;

/*  Singleton equality row.
    min. [-1       1    -1       1    -1      2    -0.2   0.7]
    s.t.     [ 0.1   0.2  -0.3    -0.1   0.4   -0.4     0.1   0.3]  = 0.5
             [   0     0     1       0     0      0       0     0]  = 3
     -2 <=   [  0.3   0.5  -0.6    -0.3   0.2   -0.2     0.4   0.1]  <= 2
                    -10  <=  x <= 10
*/
static char *test_0_postsolve()
{
    double Ax[] = {0.1, 0.2, -0.3, -0.1, 0.4, -0.4, 0.1, 0.3, 1,
                   0.3, 0.5, -0.6, -0.3, 0.2, -0.2, 0.4, 0.1};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 2, 0, 1, 2, 3, 4, 5, 6, 7};
    int Ap[] = {0, 8, 9, 17};
    int nnz = 17;
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
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {10., -10., -10., 2.71428571, -10., -6.85714286, -10.};
    double y[] = {-2.57142857, 0.14285714};
    double z[] = {-0.78571429, 1.44285714, 0.78571429, 0., 1., 0., 1.45714286};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {10., -10., 3., -10., 2.71428571, -10., -6.85714286, -10.};
    double correct_y[] = {-2.57142857, -1.68571429, 0.14285714};
    double correct_z[] = {-0.78571429, 1.44285714, 0., 0.78571429,
                          0.,          1.,         0., 1.45714286};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Singleton equality row and doubleton row.
    min. [-1       1    -1       1    -1      2    -0.2   0.7]
    s.t.     [ 0.1   0.2  -0.3    -0.1   0.4   -0.4     0.1   0.3]  = 0.5
             [   0     0     1       0     0      0       0     0]  = 3
             [   0     0     0       1     1      0       0     0]  = 4
     -2 <=   [  0.3   0.5  -0.6    -0.3   0.2   -0.2     0.4   0.1]  <= 2
                    -10  <=  x <= 10
*/
static char *test_1_postsolve()
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
    // stgs->ston_cols = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {10, -10, 0.5333333, -10, 0.6666666, -10};
    double y[] = {-4.66666667, 0.66666667};
    double z[] = {-0.73333333, 1.6, 0., 0.26666667, 0., 2.03333333};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {10, -10, 3, 0.5333333, 3.4666666, -10, 0.6666666, -10};
    double correct_y[] = {-4.66666667, -2., 0.73333333, 0.66666667};
    double correct_z[] = {-0.73333333, 1.6, 0., 0., 0., 0.26666667, 0., 2.03333333};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Another singleton equality row test. */
static char *test_singleton_eq()
{
    double Ax[] = {1.0, -1, 2, 3, 1, 2, -3, 5, 7, 9, 5, 2, -1, 3, 4, 6};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 2, 0, 1, 2, 3, 4};
    int Ap[] = {0, 5, 10, 11, 16};
    int nnz = 16;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2.0, 3, 1, 4};
    double rhs[] = {2.0, 3, 1, 4};
    double lbs[] = {-5.0, -5, -5, -5, -5};
    double ubs[] = {5.0, 5, 5, 5, 5};
    double c[] = {1.0, 2, -1, 3, -2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {5., -2.525, -1.8875, -0.2625};
    double y[] = {3.0625, -2.8125, 3.375};
    double z[] = {-3.1875, 0., 0., 0., 0.};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {5., -2.525, 0.2, -1.8875, -0.2625};
    double correct_y[] = {3.0625, -2.8125, -0.6375, 3.375};
    double correct_z[] = {-3.1875, 0., 0., 0., 0.};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Singleton equality row and doubleton row with implied bound active.
    min. [-1       1    -1       1    -1      2    -0.2   0.7]
    s.t.     [ 0.1   0.2  -0.3    -0.1   0.4   -0.4     0.1   0.3]  = 0.5
             [   0     0     1       0     0      0       0     0]  = 3
             [   0     0     0       1     1      0       0     0]  = 4
     -2 <=   [  0.3   0.5  -0.6    -0.3   0.2   -0.2     0.4   0.1]  <= 2
                    -10  <=  x <= 10
*/
static char *test_1_postsolve_implied_bound_active()
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
    double ubs[] = {10, 10, 10, 10, 4.0 - 8.0 / 5.0, 10, 10, 10};
    double c[] = {-1, 1, -1, 1, -1, 2, -0.2, 0.7};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {10.0, -10.0, 1.6, -10.0, 6.0, -10.0};
    double y[] = {-2.0, 0.0};
    double z[] = {-0.8, 1.4, 1.0, 1.2, 0.0, 1.3};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {10., -10., 3., 1.6, 2.4, -10., 6., -10.};
    double correct_y[] = {-2., -1.6, 0.8, 0.0};
    double correct_z[] = {-0.8, 1.4, 0., 0., -1., 1.2, 0., 1.3};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Singleton inequality row */
static char *test_singleton_ineq_row()
{
    double Ax[] = {1, 1, 1, 1};
    int Ai[] = {0, 0, 1, 2};
    int Ap[] = {0, 1, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 3;

    double lhs[] = {4.0, 5};
    double rhs[] = {INF, INF};
    double lbs[] = {-1.0, 0, 0};
    double ubs[] = {INF, INF, INF};
    double c[] = {3.0, 1, 2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {4.0, 1.0, 0.0};
    double y[] = {1.0};
    double z[] = {2.0, 0.0, 1.0};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {4.0, 1.0, 0.0};
    double correct_y[] = {2.0, 1.0};
    double correct_z[] = {0.0, 0.0, 1.0};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Two free column singletons in different equality constraints

        min. [1 1 -1 3 2]x
        s.t.[-1,  0,  2,  3, 4  ]     [2]
            [ 0, -1,  5,  6, 7  ] x = [5]
            [ 0,  0,  8,  9, 10 ]     [6]
            [ 0,  0, 11, 12, 13 ]     [8]
             x3, x4, x5 >= 0
             x1, x2 free
*/
static char *test_2_postsolve()
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
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {0.33333333, 0.0, 0.33333333};
    double y[] = {10.83333333, -7.33333333};
    double z[] = {0.0, 2.5, 0.0};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {0.0, -1.0, 0.33333333, 0.0, 0.33333333};
    double correct_y[] = {-1.0, -1.0, 10.83333333, -7.33333333};
    double correct_z[] = {0.0, 0.0, 0.0, 2.5, 0.0};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Example 10 in test_ston, but modified to be feasible */
static char *test_implied_free_col_ston_in_inequality_postsolve()
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
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {0., 1., 0.};
    double y[] = {-0.8, 0.};
    double z[] = {5.4, 0., 4.6};
    double obj = 0.0;

    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {0., 8., 1., 0.};
    double correct_y[] = {-0.8, 2., 0.};
    double correct_z[] = {5.4, 0., 0., 4.6};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Example 12 in test_ston, simple dual fix to lower bound */
static char *test_col_ston_dual_fix()
{
    double Ax[] = {3.0, -1, -6, 3, -2, 5, -1, 4, 3, 4};
    int Ai[] = {0, 2, 3, 0, 1, 2, 3, 0, 2, 3};
    int Ap[] = {0, 3, 7, 10};
    int nnz = 10;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {-INF, 0, -INF};
    double rhs[] = {0.0, INF, 5};
    double lbs[] = {0.0, 1, 0, 0};
    double ubs[] = {INF, INF, INF, INF};
    double c[] = {4.0, 1, -2, 7};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_cols = false;
    stgs->dual_fix = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {0., 1.66666667, 0.};
    double y[] = {0., 0., -0.66666667};
    double z[] = {6.66666667, 0., 9.66666667};
    double obj = 0.0;

    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {0., 1., 1.66666667, 0.};
    double correct_y[] = {0., 0., -0.66666667};
    double correct_z[] = {6.66666667, 1., 0., 9.66666667};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Chain of column singletons, so order in postsolve matters */
static char *test_3_postsolve()
{
    double Ax[] = {1.0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 3, 3, 3, 3,
                   3,   3, 3, 3, 3, 1, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1,
                   1,   1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 2, 3, 4,  5,  6, 7, 8, 9, 10, 1, 2, 3,  4, 5,
                6, 7, 8, 9, 10, 2,  3, 4, 5, 6, 7,  8, 9, 10, 3, 4,
                5, 6, 7, 8, 9,  10, 3, 4, 5, 6, 7,  8, 9, 10};
    int Ap[] = {0, 11, 21, 30, 38, 46};
    int nnz = 46;
    int n_rows = 5;
    int n_cols = 11;

    double lhs[] = {5.0, 4, 3, 1, 2};
    double rhs[] = {5.0, 4, 3, 1, 2};
    double lbs[] = {-INF, -INF, -INF, 0, 0, 0, 0, 0, 0, 0, -INF};
    double ubs[] = {INF, INF, INF, 1, 2, 3, 4, 3, 6, 7, 8};
    double c[] = {1.0, 1, 1.03, -5, 2.3, 3.1, -2, 1.05, -3, 4, 1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_cols = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {1., 0., 0., 4., 0., 6., 0., -10.};
    double y[] = {-0.12, -6.};
    double z[] = {0., 1.3, 2.1, -3., 0.05, -4., 3., 0.};
    double obj = 0.0;

    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {-3., 4., -1., 1., 0., 0., 4., 0., 6., 0., -10.};
    double correct_y[] = {1., -1., 2.03, -0.12, -6.};
    double correct_z[] = {0., 0., 0., 0., 1.3, 2.1, -3., 0.05, -4., 3., 0.};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* One variable fixed to infinity

    min. [1      -2        1        0        3]x

    s.t  [2       3       -4        0        3]      = 7
         [1      -1        1        1        2]      <= 7
         [2       1        3        -1       4]      >= 3
             -10 <=  x1, x2, x3, x5 <= 10
            -INF <=  x4 <= 10
*/
static char *test_4_postsolve()
{
    double Ax[] = {2, 3, -4, 3, 1, -1, 1, 1, 2, 2, 1, 3, -1, 4};
    int Ai[] = {0, 1, 2, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
    int Ap[] = {0, 4, 9, 14};
    int nnz = 14;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {7, -INF, 3};
    double rhs[] = {7, 7, INF};
    double lbs[] = {-10, -10, -10, -INF, -10};
    double ubs[] = {10, 10, 10, 10, 10};
    double c[] = {1, -2, 1, 0, 3};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {-10, 10, -6.75, -10};
    double y[] = {-0.25};
    double z[] = {1.5, -1.25, 0., 3.75};
    double obj = 0.0;

    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {-10.0, 10.0, -6.75, -73.25, -10.0};
    double correct_y[] = {-0.25, 0., 0.};
    double correct_z[] = {1.5, -1.25, 0., 0., 3.75};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Parallel columns and several other reductions

     min. [2  3  -4  -3  -6   6   3   1]x
    s.t. [1  1  -2   1   2   3  -1   5
          0  0   0   0   0   0   0   1
         -2  2   4  -2  -4  -6   2   2
          0  0   0   0   0   0   0   3
          3  3  -6   3   6   9  -3   4
          4  0  -8   4   8  12  -4   5
         -2  4   4  -2  -4  -6   2   6
          0  7   0   0   0   0   0   7
          0  5   0   0   0   0   0   8
          0  6   0   0   0   0   0   9]x <= [10 1 20 1 30 40 20 70 50 60]
          -i <= x <= i
*/
static char *test_6_postsolve()
{
    double Ax[] = {1, 1,  -2, 1, 2,  3,  -1, 5, 1,  -2, 2, 4,  -2, -4, -6, 2,
                   2, 3,  3,  3, -6, 3,  6,  9, -3, 4,  4, -8, 4,  8,  12, -4,
                   5, -2, 4,  4, -2, -4, -6, 2, 6,  7,  7, 5,  8,  6,  9};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 7, 0, 1, 2, 3, 4, 5, 6,
                7, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 3, 4, 5, 6,
                7, 0, 1, 2, 3, 4, 5, 6, 7, 1, 7, 1, 7, 1, 7};
    int Ap[] = {0, 8, 9, 17, 18, 26, 33, 41, 43, 45, 47};
    int nnz = 47;
    int n_rows = 10;
    int n_cols = 8;

    double lhs[] = {-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {10, 1, 20, 1, 30, 40, 20, 70, 50, 60};
    double lbs[] = {-1, -2, -3, -4, -5, -6, -7, -8};
    double ubs[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double c[] = {2, 3, -4, -3, -6, 6, 3, 1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    // with dual reductions and domain propagation enabled we reduce the problem
    // much more. Would be interesting to add this test again with those
    // reductions enabled
    stgs->dual_fix = false;
    stgs->primal_propagation = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    // construct optimal primal solution to reduced problem (computed offline).
    // on linux vs mac the parallel column kept is different, leading to
    // different postsolved solutions
#ifdef __linux__
    double x[] = {-2., 12.5, 21., -8.};
    double y[] = {0., 0., 0., 0., 0., 0., 0., 0.};
    double z[] = {3., -4., -3., 1.};
    double obj = 0.0;
#else
    double x[] = {-2., -8.33333333, -21., -8.};
    double y[] = {0., 0., 0., 0., 0., 0., 0., 0.};
    double z[] = {3., 6., 3., 1.};
    double obj = 0.0;
#endif

    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {-1.0, -2.0, 3.0, 4.0, 5.0, -6.0, -7.0, -8.0};
    double correct_y[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double correct_z[] = {2., 3., -4., -3., -6., 6., 3., 1.};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* same as test 6 but with domain propagation active */
static char *test_7_postsolve()
{
    double Ax[] = {1, 1,  -2, 1, 2,  3,  -1, 5, 1,  -2, 2, 4,  -2, -4, -6, 2,
                   2, 3,  3,  3, -6, 3,  6,  9, -3, 4,  4, -8, 4,  8,  12, -4,
                   5, -2, 4,  4, -2, -4, -6, 2, 6,  7,  7, 5,  8,  6,  9};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 7, 0, 1, 2, 3, 4, 5, 6,
                7, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 3, 4, 5, 6,
                7, 0, 1, 2, 3, 4, 5, 6, 7, 1, 7, 1, 7, 1, 7};
    int Ap[] = {0, 8, 9, 17, 18, 26, 33, 41, 43, 45, 47};
    int nnz = 47;
    int n_rows = 10;
    int n_cols = 8;

    double lhs[] = {-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {10, 1, 20, 1, 30, 40, 20, 70, 50, 60};
    double lbs[] = {-1, -2, -3, -4, -5, -6, -7, -8};
    double ubs[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double c[] = {2, 3, -4, -3, -6, 6, 3, 1};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    // with dual reduction enabled we reduce the problem much more.
    // Would be interesting to add this test again with that reductions enabled
    stgs->dual_fix = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    // construct optimal primal solution to reduced problem (computed offline)
#ifdef __linux__
    double x[] = {-2., 12.5, 21., -8.};
    double y[] = {0., 0., 0., 0., 0.};
    double z[] = {3., -4., -3., 1.};
    double obj = 0.0;
#else
    double x[] = {-25., -2., 21., -8.};
    double y[] = {0., 0., 0., 0., 0.};
    double z[] = {2., 3., -3., 1.};
    double obj = 0.0;

#endif

    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {-1.0, -2.0, 3.0, 4.0, 5.0, -6.0, -7.0, -8.0};
    double correct_y[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double correct_z[] = {2., 3., -4., -3., -6., 6., 3., 1.};

    // on mac another solution seems to be postsolved, because a different parallel
    // column is kept. So for non-linux we just check the objective value (not
    // feasibility)
    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Doubleton row.
    min. [2 -1 -1 3 2]
    s.t.     [ 2     1       -1        3     2]  = 0.5
             [ 1    -2        0        0     0]  = 3
             [ 1,    4,       0,       3,    1]  = 4
                    0  <=  x2, x3, x4, x5
*/
static char *test_8_postsolve()
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
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {0.83333333, 0.16666667, 0., 0.};
    double y[] = {1., -0.33333333};
    double z[] = {0., 0., 1., 0.33333333};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {0.66666667, 0.83333333, 0.16666667, 0., 0.};
    double correct_y[] = {1., 0.33333333, -0.33333333};
    double correct_z[] = {0., 0., 0., 1., 0.33333333};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Parallel rows.
    min. [1 1.1 1.2]
    s.t.     [ 1 1 1] >= 1
             [ 2 2 2] >= 4
             [ -3 -3 -3]  <= -9
                0 <= x1, x2, x3
*/
static char *test_9_postsolve()
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
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {3., 0., 0.};
    double y[] = {0.5};
    double z[] = {0., 0.1, 0.2};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {3., 0., 0.};
    double correct_y[] = {0., 0., -0.33333333};
    double correct_z[] = {0., 0.1, 0.2};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* implied bound active at optimal solution together with original bound*/
static char *test_pathological_ston_one()
{
    double Ax[] = {1.0, 1, 2, 3, 4, -1, -5};
    int Ai[] = {0, 1, 2, 3, 1, 2, 3};
    int Ap[] = {0, 1, 4, 7};
    int nnz = 7;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {-INF, 3, 4};
    double rhs[] = {3.0, 3, 4};
    double lbs[] = {3.0, -5, -5, -5};
    double ubs[] = {INF, 5, 5, 5};
    double c[] = {1.0, 2, 3, -4};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {3.64705882, -5., 3.11764706};
    double y[] = {-0.35294118, 0.58823529};
    double z[] = {0., 4.29411765, 0.};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {3., 3.64705882, -5., 3.11764706};
    double correct_y[] = {0., -0.35294118, 0.58823529};
    double correct_z[] = {1., 0., 4.29411765, 0.};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* implied bound active at optimal solution together with original bound*/
static char *test_pathological_ston_two()
{
    double Ax[] = {1.0, 1, 2, 3, 4, -1, -5};
    int Ai[] = {0, 1, 2, 3, 1, 2, 3};
    int Ap[] = {0, 1, 4, 7};
    int nnz = 7;
    int n_rows = 3;
    int n_cols = 4;

    double lhs[] = {-INF, 3, 4};
    double rhs[] = {3.0, 3, 4};
    double lbs[] = {3.0, -5, -5, -5};
    double ubs[] = {INF, 5, 5, 5};
    double c[] = {-2.0, 2, 3, -4};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    run_presolver(presolver);

    Mapping *maps = prob->constraints->state->work->mappings;
    int *rows_map = maps->rows;
    int *cols_map = maps->cols;

    // construct optimal primal solution to reduced problem (computed offline)
    double x[] = {3.64705882, -5., 3.11764706};
    double y[] = {-0.35294118, 0.58823529};
    double z[] = {0., 4.29411765, 0.};
    double obj = 0.0;
    postsolve(presolver, x, y, z, obj);

    // check that the primal solution to the original problem is correct
    double correct_x[] = {3., 3.64705882, -5., 3.11764706};
    double correct_y[] = {-2., -0.35294118, 0.58823529};
    double correct_z[] = {0., 0., 4.29411765, 0.};

    mu_assert("postsolve error",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, POSTSOLVE_TOL_FEAS));
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_postsolve()
{
    mu_run_test(test_0_postsolve, counter_postsolve);
    mu_run_test(test_1_postsolve, counter_postsolve);
    mu_run_test(test_1_postsolve_implied_bound_active, counter_postsolve);
    mu_run_test(test_singleton_ineq_row, counter_postsolve);
    mu_run_test(test_2_postsolve, counter_postsolve);
    mu_run_test(test_implied_free_col_ston_in_inequality_postsolve,
                counter_postsolve);
    mu_run_test(test_col_ston_dual_fix, counter_postsolve);
    mu_run_test(test_3_postsolve, counter_postsolve);
    mu_run_test(test_4_postsolve, counter_postsolve);
    mu_run_test(test_6_postsolve, counter_postsolve);
    mu_run_test(test_7_postsolve, counter_postsolve);
    mu_run_test(test_8_postsolve, counter_postsolve);
    mu_run_test(test_9_postsolve, counter_postsolve);
    mu_run_test(test_singleton_eq, counter_postsolve);
    mu_run_test(test_pathological_ston_one, counter_postsolve);
    mu_run_test(test_pathological_ston_two, counter_postsolve);
    //        all tests above pass

    return 0;
}

int test_postsolve()
{
    const char *result = all_tests_postsolve();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("postsolve: TEST FAILED!\n");
    }
    else
    {
        printf("postsolve: ALL TESTS PASSED\n");
    }
    printf("postsolve: Tests run: %d\n", counter_postsolve);
    return result == 0;
}

#endif