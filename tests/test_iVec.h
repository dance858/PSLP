#ifndef TEST_IVEC_H
#define TEST_IVEC_H

#include "iVec.h"
#include "minunit.h"
#include <stdio.h>

static int counter_iVec = 0;

#define CMP_IVEC(vec, arr)                                                          \
    do                                                                              \
    {                                                                               \
        for (PSLP_uint i = 0; i < (vec)->len; ++i)                                  \
        {                                                                           \
            mu_assert("error, vec->data[i] != arr[i]", (vec)->data[i] == (arr)[i]); \
        }                                                                           \
    } while (0)

// test initialization, append, free
static char *test_1_iVec()
{
    iVec *vec = iVec_new(3);
    iVec_append(vec, 1);
    iVec_append(vec, -1);
    iVec_append(vec, 2);
    mu_assert("error, vec->len != 3", vec->len == 3);
    mu_assert("error, vec->capacity != 3", vec->capacity == 3);

    int correct_data1[] = {1, -1, 2};
    CMP_IVEC(vec, correct_data1);

    iVec_append(vec, 5);
    mu_assert("error, vec->len != 4", vec->len == 4);
    mu_assert("error, vec->capacity != 6", vec->capacity == 6);
    int correct_data2[] = {1, -1, 2, 5};
    CMP_IVEC(vec, correct_data2);

    iVec_free(vec);
    return 0;
}

// test initialization several vectors, append, free
static char *test_2_iVec()
{
    iVec *vec1 = iVec_new(3);
    iVec *vec2 = iVec_new(4);
    iVec_append(vec1, 1);
    iVec_append(vec1, -1);
    iVec_append(vec1, 2);
    iVec_append(vec2, 5);
    iVec_append(vec2, 4);
    iVec_append(vec2, 6);
    iVec_append(vec2, 7);

    int correct_data1[] = {1, -1, 2};
    CMP_IVEC(vec1, correct_data1);

    int correct_data2[] = {5, 4, 6, 7};
    CMP_IVEC(vec2, correct_data2);

    iVec_free(vec1);
    iVec_free(vec2);

    return 0;
}

static const char *all_tests_iVec()
{
    mu_run_test(test_1_iVec, counter_iVec);
    mu_run_test(test_2_iVec, counter_iVec);
    return 0;
}

int test_iVec()
{
    const char *result = all_tests_iVec();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("iVec: TEST FAILED!\n");
    }
    else
    {
        printf("iVec: ALL TESTS PASSED\n");
    }
    printf("iVec: Tests run: %d\n", counter_iVec);
    return result == 0;
}

#endif // TEST_IVEC_H
