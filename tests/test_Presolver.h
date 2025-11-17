#ifndef TEST_PRESOLVER_H
#define TEST_PRESOLVER_H

#include "API.h"
#include "Debugger.h"

#include "Problem.h"
#include "Workspace.h"
#include "minunit.h"
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
    Presolver *presolver = new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs,
                                         lbs, ubs, c, stgs, true);

    mu_assert("Presolver initialization failed", presolver != NULL);
    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_presolver()
{
    mu_run_test(test_00_presolver, counter_presolver);
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
