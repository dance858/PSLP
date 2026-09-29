#ifndef TEST_SPARSEACCUMULATOR_H
#define TEST_SPARSEACCUMULATOR_H

#include "SparseAccumulator.h"
#include "minunit.h"
#include <stdint.h>
#include <stdio.h>

int counter_sparse_accumulator = 0;

/* Arrays for an accumulator over 5 columns. */
typedef struct AccumulatorArrays
{
    int stamp[5];
    double value[5];
    uint8_t flags[5];
    int touched[5];
} AccumulatorArrays;

static void accumulator_setup(SparseAccumulator *acc, AccumulatorArrays *arr)
{
    sparse_accumulator_init(acc, arr->stamp, arr->value, arr->flags, arr->touched,
                            5);
    sparse_accumulator_clear(acc);
}

/* Two adds to column 3 sum their values, OR their flags and touch it once. */
static char *test_sparse_accumulator_add()
{
    SparseAccumulator acc;
    AccumulatorArrays arr;
    accumulator_setup(&acc, &arr);
    sparse_accumulator_add(&acc, 3, 1.5, 1);
    sparse_accumulator_add(&acc, 3, 2.0, 2);
    mu_assert("touched once", acc.n_touched == 1 && acc.touched[0] == 3);
    mu_assert("value", acc.value[3] == 3.5);
    mu_assert("flags", acc.flags[3] == 3);
    return 0;
}

/* After a clear, a column touched in the previous vector starts from zero. */
static char *test_sparse_accumulator_clear()
{
    SparseAccumulator acc;
    AccumulatorArrays arr;
    accumulator_setup(&acc, &arr);
    sparse_accumulator_add(&acc, 1, 4.0, 1);
    sparse_accumulator_clear(&acc);
    mu_assert("empty", acc.n_touched == 0);
    sparse_accumulator_add(&acc, 1, 2.0, 2);
    mu_assert("touched", acc.n_touched == 1 && acc.touched[0] == 1);
    mu_assert("fresh value", acc.value[1] == 2.0);
    mu_assert("fresh flags", acc.flags[1] == 2);
    return 0;
}

/* The touched columns come out ascending. */
static char *test_sparse_accumulator_sort()
{
    SparseAccumulator acc;
    AccumulatorArrays arr;
    accumulator_setup(&acc, &arr);
    sparse_accumulator_add(&acc, 4, 1.0, 1);
    sparse_accumulator_add(&acc, 0, 1.0, 1);
    sparse_accumulator_add(&acc, 2, 1.0, 1);
    sparse_accumulator_add(&acc, 0, 1.0, 1);
    sparse_accumulator_sort(&acc);
    mu_assert("count", acc.n_touched == 3);
    mu_assert("order",
              acc.touched[0] == 0 && acc.touched[1] == 2 && acc.touched[2] == 4);
    mu_assert("value", acc.value[0] == 2.0);
    return 0;
}

/* Init zeroes stamps left over from an earlier use, so a new accumulator
   treats every column as untouched. */
static char *test_sparse_accumulator_init()
{
    SparseAccumulator acc;
    AccumulatorArrays arr;
    for (int c = 0; c < 5; ++c)
    {
        arr.stamp[c] = 1;
        arr.value[c] = 9.0;
        arr.flags[c] = 7;
    }
    accumulator_setup(&acc, &arr);
    sparse_accumulator_add(&acc, 2, 1.0, 1);
    mu_assert("touched", acc.n_touched == 1 && acc.touched[0] == 2);
    mu_assert("fresh", acc.value[2] == 1.0 && acc.flags[2] == 1);
    return 0;
}

static const char *all_tests_sparse_accumulator()
{
    mu_run_test(test_sparse_accumulator_add, counter_sparse_accumulator);
    mu_run_test(test_sparse_accumulator_clear, counter_sparse_accumulator);
    mu_run_test(test_sparse_accumulator_sort, counter_sparse_accumulator);
    mu_run_test(test_sparse_accumulator_init, counter_sparse_accumulator);
    return 0;
}

int test_sparse_accumulator()
{
    const char *result = all_tests_sparse_accumulator();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("sparse accumulator: TEST FAILED!\n");
    }
    else
    {
        printf("sparse accumulator: ALL TESTS PASSED\n");
    }
    printf("sparse accumulator: Tests run: %d\n", counter_sparse_accumulator);
    return result == 0;
}

#endif // TEST_SPARSEACCUMULATOR_H
