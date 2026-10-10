#ifndef TEST_THREAD_POOL_H
#define TEST_THREAD_POOL_H

#include "minunit.h"
#include "thread_pool.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int counter_thread_pool = 0;

typedef struct SumPass
{
    int n_items;
    int n_chunks;
    int64_t *partials;
} SumPass;

static void sum_chunk(void *context, int chunk, int thread)
{
    (void) thread;
    SumPass *pass = (SumPass *) context;
    int lo = (int) ((int64_t) pass->n_items * chunk / pass->n_chunks);
    int hi = (int) ((int64_t) pass->n_items * (chunk + 1) / pass->n_chunks);
    int64_t sum = 0;
    for (int ii = lo; ii < hi; ++ii)
    {
        sum += ii;
    }
    pass->partials[chunk] = sum;
}

/* Every chunk of every pass runs exactly once, whatever the worker count. */
static const char *test_pool_runs_every_chunk(void)
{
    int n_items = 100000;
    int64_t expected = (int64_t) n_items * (n_items - 1) / 2;
    for (int n_workers = 0; n_workers <= 3; n_workers++)
    {
        ThreadPool *pool = thread_pool_new(n_workers);
        mu_assert("pool", pool != NULL);
        mu_assert("at least the calling thread", thread_pool_n_threads(pool) >= 1);
        for (int pass_index = 0; pass_index < 50; pass_index++)
        {
            SumPass pass = {n_items, 37, (int64_t *) calloc(37, sizeof(int64_t))};
            thread_pool_run(pool, sum_chunk, &pass, pass.n_chunks);
            int64_t total = 0;
            for (int chunk = 0; chunk < pass.n_chunks; chunk++)
            {
                total += pass.partials[chunk];
            }
            free(pass.partials);
            mu_assert("chunk sums add up", total == expected);
        }
        thread_pool_run(pool, sum_chunk, NULL, 0);
        thread_pool_free(pool);
    }
    return 0;
}

static const char *test_pool_chunk_helpers(void)
{
    mu_assert("chunk count by items", thread_pool_chunk_count(4, 10, 4096) == 1);
    mu_assert("chunk count by threads",
              thread_pool_chunk_count(4, 1000000, 4096) == 32);

    int bounds[5];
    thread_pool_chunk_bounds(10, 4, bounds);
    int expected_bounds[5] = {0, 2, 5, 7, 10};
    mu_assert("even bounds", memcmp(bounds, expected_bounds, sizeof(bounds)) == 0);

    /* Item sizes 1, 7, 1, 1: the two chunks split after the heavy item. */
    int prefix[5] = {0, 1, 8, 9, 10};
    thread_pool_chunk_bounds_by_prefix(prefix, 4, 2, bounds);
    mu_assert("bounds by prefix",
              bounds[0] == 0 && bounds[1] == 2 && bounds[2] == 4);

    /* Two chunks, three buckets: counts become write positions, bucket major. */
    int counts[6] = {2, 0, 1, 1, 3, 0};
    int offsets[4];
    int total = thread_pool_bucket_offsets(counts, 2, 3, offsets);
    int expected_counts[6] = {0, 3, 6, 2, 3, 7};
    int expected_offsets[4] = {0, 3, 6, 7};
    mu_assert("bucket total", total == 7);
    mu_assert("bucket positions",
              memcmp(counts, expected_counts, sizeof(counts)) == 0);
    mu_assert("bucket offsets",
              memcmp(offsets, expected_offsets, sizeof(offsets)) == 0);
    return 0;
}

static const char *test_pool_atomics(void)
{
    uint64_t word = 6;
    mu_assert("fetch add returns the old value",
              ps_atomic_fetch_add_u64(&word, 4) == 6);
    mu_assert("fetch add adds", word == 10);
    mu_assert("fetch or returns the old value",
              ps_atomic_fetch_or_u64(&word, 1) == 10);
    mu_assert("fetch or sets", word == 11);
    int flag = 0;
    mu_assert("int fetch add", ps_atomic_fetch_add_int(&flag, 2) == 0 && flag == 2);
    ps_atomic_store_int(&flag, 5);
    mu_assert("load", ps_atomic_load_int(&flag) == 5);
    return 0;
}

static const char *all_tests_thread_pool(void)
{
    mu_run_test(test_pool_runs_every_chunk, counter_thread_pool);
    mu_run_test(test_pool_chunk_helpers, counter_thread_pool);
    mu_run_test(test_pool_atomics, counter_thread_pool);
    return 0;
}

int test_thread_pool(void)
{
    const char *result = all_tests_thread_pool();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("thread pool: TEST FAILED!\n");
    }
    else
    {
        printf("thread pool: ALL TESTS PASSED\n");
    }
    printf("thread pool: Tests run: %d\n", counter_thread_pool);
    return result == 0;
}

#endif // TEST_THREAD_POOL_H
