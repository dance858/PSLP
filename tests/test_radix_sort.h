#ifndef TEST_RADIX_SORT_H
#define TEST_RADIX_SORT_H

#include "radix_sort.h"

#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int counter_radix_sort = 0;

// helper: check that rows are sorted by (sparsity_IDs[rows[i]],
// coeff_hashes[rows[i]]) treated as unsigned uint32_t, matching radix sort behavior
static int is_sorted_by_keys(const int *rows, int n, const int *sparsity_IDs,
                             const int *coeff_hashes)
{
    for (int i = 1; i < n; i++)
    {
        unsigned int sp_prev = (unsigned int) sparsity_IDs[rows[i - 1]];
        unsigned int sp_curr = (unsigned int) sparsity_IDs[rows[i]];
        if (sp_prev > sp_curr)
        {
            return 0;
        }
        if (sp_prev == sp_curr)
        {
            unsigned int ch_prev = (unsigned int) coeff_hashes[rows[i - 1]];
            unsigned int ch_curr = (unsigned int) coeff_hashes[rows[i]];
            if (ch_prev > ch_curr)
            {
                return 0;
            }
        }
    }
    return 1;
}

// helper: check that sorted array is a permutation of 0..n-1
static int is_permutation(const int *rows, int n)
{
    int *seen = (int *) calloc((size_t) n, sizeof(int));
    if (!seen) return 0;
    for (int i = 0; i < n; i++)
    {
        if (rows[i] < 0 || rows[i] >= n || seen[rows[i]])
        {
            free(seen);
            return 0;
        }
        seen[rows[i]] = 1;
    }
    free(seen);
    return 1;
}

/* Test 1: Empty input (n=0) should not crash */
static char *test_1_radix_sort()
{
    int rows[1] = {0};
    int sparsity_IDs[1] = {0};
    int coeff_hashes[1] = {0};
    int aux[1] = {0};

    radix_sort_rows(rows, 0, sparsity_IDs, coeff_hashes, aux);
    // just check it doesn't crash
    return 0;
}

/* Test 2: Single element */
static char *test_2_radix_sort()
{
    int rows[1] = {0};
    int sparsity_IDs[1] = {42};
    int coeff_hashes[1] = {99};
    int aux[1];

    radix_sort_rows(rows, 1, sparsity_IDs, coeff_hashes, aux);
    mu_assert("single element should remain 0", rows[0] == 0);
    return 0;
}

/* Test 3: Already sorted input (small, uses insertion sort path) */
static char *test_3_radix_sort()
{
    int n = 5;
    int rows[] = {0, 1, 2, 3, 4};
    int sparsity_IDs[] = {1, 2, 3, 4, 5};
    int coeff_hashes[] = {10, 20, 30, 40, 50};
    int aux[5];

    radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));
    return 0;
}

/* Test 4: Reverse sorted input (small, insertion sort path) */
static char *test_4_radix_sort()
{
    int n = 5;
    int rows[] = {0, 1, 2, 3, 4};
    int sparsity_IDs[] = {50, 40, 30, 20, 10};
    int coeff_hashes[] = {0, 0, 0, 0, 0};
    int aux[5];

    radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));

    // row 4 has smallest sparsity_ID (10) so should be first
    mu_assert("first should be row 4", rows[0] == 4);
    mu_assert("last should be row 0", rows[4] == 0);
    return 0;
}

/* Test 5: Secondary key (coeff_hashes) breaks ties in primary key */
static char *test_5_radix_sort()
{
    int n = 4;
    int rows[] = {0, 1, 2, 3};
    int sparsity_IDs[] = {100, 100, 100, 100}; // all same
    int coeff_hashes[] = {40, 10, 30, 20};
    int aux[4];

    radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));

    // sorted by coeff_hash: 10(row1), 20(row3), 30(row2), 40(row0)
    mu_assert("first", rows[0] == 1);
    mu_assert("second", rows[1] == 3);
    mu_assert("third", rows[2] == 2);
    mu_assert("fourth", rows[3] == 0);
    return 0;
}

/* Test 6: Negative sparsity_IDs treated as unsigned in radix sort path (n >= 256).
   High bit set sorts after positive values. */
static char *test_6_radix_sort()
{
    int n = 300;
    int rows[300];
    int sparsity_IDs[300];
    int coeff_hashes[300];
    int aux[300];

    for (int i = 0; i < n; i++)
    {
        rows[i] = i;
        // mix positive and negative sparsity_IDs
        sparsity_IDs[i] = (i % 3 == 0) ? -(i + 1) : (i + 1);
        coeff_hashes[i] = 0;
    }

    radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));

    // positive sparsity_IDs (unsigned < 0x80000000) should come before
    // negative ones (unsigned >= 0x80000000)
    int saw_negative = 0;
    for (int i = 0; i < n; i++)
    {
        int sp = sparsity_IDs[rows[i]];
        if (sp < 0) saw_negative = 1;
        if (saw_negative && sp >= 0)
        {
            mu_assert("positive should not appear after negative in unsigned order",
                      0);
        }
    }
    return 0;
}

/* Test 7: Larger input (n=300) forces the radix sort path (threshold is 256) */
static char *test_7_radix_sort()
{
    int n = 300;
    int rows[300];
    int sparsity_IDs[300];
    int coeff_hashes[300];
    int aux[300];

    srand(12345);
    for (int i = 0; i < n; i++)
    {
        rows[i] = i;
        sparsity_IDs[i] = rand() % 50;
        coeff_hashes[i] = rand() % 100;
    }

    radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));
    return 0;
}

/* Test 8: Large input with random keys including negative values */
static char *test_8_radix_sort()
{
    int n = 1000;
    int *rows = (int *) malloc((size_t) n * sizeof(int));
    int *sparsity_IDs = (int *) malloc((size_t) n * sizeof(int));
    int *coeff_hashes = (int *) malloc((size_t) n * sizeof(int));
    int *aux = (int *) malloc((size_t) n * sizeof(int));

    srand(9999);
    for (int i = 0; i < n; i++)
    {
        rows[i] = i;
        sparsity_IDs[i] = rand() - RAND_MAX / 2;
        coeff_hashes[i] = rand() - RAND_MAX / 2;
    }

    radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));

    free(rows);
    free(sparsity_IDs);
    free(coeff_hashes);
    free(aux);
    return 0;
}

/* Test 9: All identical keys */
static char *test_9_radix_sort()
{
    int n = 300;
    int rows[300];
    int sparsity_IDs[300];
    int coeff_hashes[300];
    int aux[300];

    for (int i = 0; i < n; i++)
    {
        rows[i] = i;
        sparsity_IDs[i] = 42;
        coeff_hashes[i] = 77;
    }

    radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));
    return 0;
}

/* Test 10: parallel_radix_sort_rows on a large array (exercises the parallel path)
 */
static char *test_10_radix_sort()
{
    int n = 200000;
    int *rows = (int *) malloc((size_t) n * sizeof(int));
    int *sparsity_IDs = (int *) malloc((size_t) n * sizeof(int));
    int *coeff_hashes = (int *) malloc((size_t) n * sizeof(int));
    int *aux = (int *) malloc((size_t) n * sizeof(int));

    srand(77777);
    for (int i = 0; i < n; i++)
    {
        rows[i] = i;
        sparsity_IDs[i] = rand() - RAND_MAX / 2;
        coeff_hashes[i] = rand() - RAND_MAX / 2;
    }

    parallel_radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));

    free(rows);
    free(sparsity_IDs);
    free(coeff_hashes);
    free(aux);
    return 0;
}

/* Test 11: parallel_radix_sort_rows with small input (falls back to sequential) */
static char *test_11_radix_sort()
{
    int n = 50;
    int rows[50];
    int sparsity_IDs[50];
    int coeff_hashes[50];
    int aux[50];

    srand(555);
    for (int i = 0; i < n; i++)
    {
        rows[i] = i;
        sparsity_IDs[i] = rand() % 10;
        coeff_hashes[i] = rand() % 10;
    }

    parallel_radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));
    return 0;
}

/* Test 12: Keys with values near INT_MAX and INT_MIN (n >= 256 for radix path) */
static char *test_12_radix_sort()
{
    int n = 300;
    int rows[300];
    int sparsity_IDs[300];
    int coeff_hashes[300];
    int aux[300];

    // fill most with benign values
    for (int i = 0; i < n; i++)
    {
        rows[i] = i;
        sparsity_IDs[i] = i;
        coeff_hashes[i] = 0;
    }

    // set some extreme values
    sparsity_IDs[0] = 0;
    sparsity_IDs[1] = 2147483647;  // INT_MAX = 0x7FFFFFFF
    sparsity_IDs[2] = -2147483648; // INT_MIN = 0x80000000
    sparsity_IDs[3] = -1;          // 0xFFFFFFFF

    radix_sort_rows(rows, (size_t) n, sparsity_IDs, coeff_hashes, aux);
    mu_assert("should be sorted",
              is_sorted_by_keys(rows, n, sparsity_IDs, coeff_hashes));
    mu_assert("should be permutation", is_permutation(rows, n));

    // In unsigned order: row0(sp=0) should come before row1(sp=INT_MAX)
    // which should come before row2(sp=INT_MIN as 0x80000000)
    // which should come before row3(sp=-1 as 0xFFFFFFFF)
    int pos_row0 = -1, pos_row1 = -1, pos_row2 = -1, pos_row3 = -1;
    for (int i = 0; i < n; i++)
    {
        if (rows[i] == 0) pos_row0 = i;
        if (rows[i] == 1) pos_row1 = i;
        if (rows[i] == 2) pos_row2 = i;
        if (rows[i] == 3) pos_row3 = i;
    }
    mu_assert("row0 before row1", pos_row0 < pos_row1);
    mu_assert("row1 before row2", pos_row1 < pos_row2);
    mu_assert("row2 before row3", pos_row2 < pos_row3);
    return 0;
}

static const char *all_tests_radix_sort()
{
    mu_run_test(test_1_radix_sort, counter_radix_sort);
    mu_run_test(test_2_radix_sort, counter_radix_sort);
    mu_run_test(test_3_radix_sort, counter_radix_sort);
    mu_run_test(test_4_radix_sort, counter_radix_sort);
    mu_run_test(test_5_radix_sort, counter_radix_sort);
    mu_run_test(test_6_radix_sort, counter_radix_sort);
    mu_run_test(test_7_radix_sort, counter_radix_sort);
    mu_run_test(test_8_radix_sort, counter_radix_sort);
    mu_run_test(test_9_radix_sort, counter_radix_sort);
    mu_run_test(test_10_radix_sort, counter_radix_sort);
    mu_run_test(test_11_radix_sort, counter_radix_sort);
    mu_run_test(test_12_radix_sort, counter_radix_sort);

    return 0;
}

int test_radix_sort()
{
    const char *result = all_tests_radix_sort();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("radix_sort: TEST FAILED!\n");
    }
    else
    {
        printf("radix_sort: ALL TESTS PASSED\n");
    }
    printf("radix_sort: Tests run: %d\n", counter_radix_sort);
    return result == 0;
}

#endif // TEST_RADIX_SORT_H
