#ifndef TEST_PRESOLVER_H
#define TEST_PRESOLVER_H

#include "Debugger.h"
#include "PSLP_API.h"

#include "Problem.h"
#include "Workspace.h"
#include "minunit.h"
#include <math.h>
#include <stdio.h>

static int counter_presolver = 0;

/* Test presolver init and free for memory leaks (if valgrind doesn't find any
   leaks the programme should be leak-free when executing regularly without any
   memory allocation failures, because all memory is allocated when initializing
   the presolver) */
static char *test_00_presolver()
{
    double Ax[] = {2, 1, -1, 3, 2, 1, -2, 1, 4, 3, 1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 3, 4};
    int Ap[] = {0, 5, 7, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;
    double lhs[] = {2, -1, 4};
    double rhs[] = {2, -1, 4};
    double lbs[] = {0.4, 0, 0, 0, 0};
    double ubs[] = {0.6, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    mu_assert("Presolver initialization failed", presolver != NULL);
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

// afiro - just test that it runs (we implement this so we can have a similar
// C++ test later)
static char *test_01_presolver()
{
    double Ax[] = {
        -1.,   1.,    1.,    -1.06, 1.,    1.,    -1.,   1.4,   -1.,   -1.,   -1.,
        -1.,   1.,    1.,    -1.06, -1.06, -0.96, -0.86, 1.,    1.,    -1.,   1.,
        -1.,   1.,    -1.,   1.,    -1.,   -1.,   1.,    1.,    1.,    -0.43, 1.,
        1.,    -1.,   1.4,   -0.43, -0.43, -0.39, -0.37, 1.,    1.,    1.,    1.,
        1.,    -1.,   1.,    1.,    1.,    -1.,   1.,    -1.,   1.,    -1.,   1.,
        -1.,   2.364, 2.386, 2.408, 2.429, -1.,   2.191, 2.219, 2.249, 2.279, -1.,
        0.109, -1.,   0.109, 0.108, 0.108, 0.107, 0.301, -1.,   0.301, 0.313, 0.313,
        0.326, -1.,   1.,    1.,    1.,    1.};
    int Ai[] = {0,  1,  2,  0,  3,  0,  1,  12, 4,  5,  6,  7,  12, 13, 4,  5,  6,
                7,  14, 4,  8,  5,  9,  6,  10, 7,  11, 15, 16, 17, 18, 15, 19, 15,
                16, 28, 20, 21, 22, 23, 30, 20, 21, 22, 23, 28, 29, 31, 20, 24, 21,
                25, 22, 26, 23, 27, 8,  9,  10, 11, 18, 24, 25, 26, 27, 2,  15, 13,
                20, 21, 22, 23, 0,  17, 4,  5,  6,  7,  29, 3,  19, 14, 30};
    int Ap[] = {0,  3,  5,  6,  8,  14, 19, 21, 23, 25, 27, 31, 33, 34,
                36, 41, 48, 50, 52, 54, 56, 65, 67, 72, 74, 79, 81, 83};
    int nnz = 83;
    int n_rows = 27;
    int n_cols = 32;
    double lhs[] = {0.,   0.,   -INF, -INF, 0.,   0.,   -INF, -INF, -INF,
                    -INF, 0.,   0.,   -INF, -INF, 0.,   44.,  -INF, -INF,
                    -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {0., 0.,  80.,  0., 0., 0., 80., 0., 0., 0., 0., 0.,   500., 0.,
                    0., 44., 500., 0., 0., 0., 0.,  0., 0., 0., 0., 310., 300.};
    double lbs[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double ubs[] = {INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
                    INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
                    INF, INF, INF, INF, INF, INF, INF, INF, INF, INF};
    double c[] = {0., -0.4,  0., 0., 0., 0.,   0.,    0., 0., 0., 0.,
                  0., -0.32, 0., 0., 0., -0.6, 0.,    0., 0., 0., 0.,
                  0., 0.,    0., 0., 0., 0.,   -0.48, 0., 0., 10.};

    Settings *stgs = default_settings();
    stgs->parallel_rows = false;
    stgs->parallel_cols = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("Presolver initialization failed", presolver != NULL);

    run_presolver(presolver);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* A column that is both empty (no entries in A) and fixed (lb == ub) must
   contribute c * x to the objective offset exactly once. Previously
   remove_variables_with_close_bounds fixed it and remove_empty_cols fixed it
   again, doubling its contribution (Netlib 80bau3b has 107 such columns). */
static char *test_02_presolver()
{
    // min -x0 - x1 + 7 x2  s.t.  x0 + 2 x1 <= 4,  3 x0 + x1 <= 6,
    //                             0 <= x0, x1 <= 10,  x2 = 3 (empty column)
    double Ax[] = {1, 2, 3, 1};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 3;
    double lhs[] = {-INF, -INF};
    double rhs[] = {4, 6};
    double lbs[] = {0, 0, 3};
    double ubs[] = {10, 10, 3};
    double c[] = {-1, -1, 7};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("Presolver initialization failed", presolver != NULL);

    run_presolver(presolver);

    mu_assert("empty fixed column not removed", presolver->reduced_prob->n == 2);
    mu_assert("objective offset of empty fixed column counted more than once",
              fabs(presolver->reduced_prob->obj_offset - 21.0) < 1e-12);

    // the reduced problem's optimum is x = (1.6, 1.2); postsolve must put the
    // fixed column back at its bound
    double x[] = {1.6, 1.2};
    double y[] = {-0.4, -0.2};
    double z[] = {0, 0};
    postsolve(presolver, x, y, z);
    mu_assert("postsolved value of empty fixed column", presolver->sol->x[2] == 3.0);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_presolver()
{
    mu_run_test(test_00_presolver, counter_presolver);
    mu_run_test(test_01_presolver, counter_presolver);
    mu_run_test(test_02_presolver, counter_presolver);
    return 0;
}

int test_presolver()
{
    const char *result = all_tests_presolver();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("presolver: TEST FAILED!\n");
    }
    else
    {
        printf("presolver: ALL TESTS PASSED\n");
    }
    printf("presolver: Tests run: %d\n", counter_presolver);
    return result == 0;
}

#endif // TEST_PRESOLVER_H
