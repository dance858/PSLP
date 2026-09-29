#ifndef TEST_CHAINS_H
#define TEST_CHAINS_H

#include "Chains.h"
#include "minunit.h"
#include <stdio.h>

int counter_chains = 0;

/* Runs compute_chain_depths on n nodes and checks depth, order and the dropped list
   against the expectations (order and dropped compared up to their lengths). */
static char *check_chain_depths(int n, const int *succ, const int *drop_priority,
                                const int *depth_ok, const int *order_ok,
                                int n_order, const int *dropped_ok, int n_dropped_ok)
{
    int depth[8], order[8], dropped[8], stamp[8], path[8];
    int n_dropped = compute_chain_depths(n, succ, drop_priority, depth, order,
                                         dropped, stamp, path);
    mu_assert("dropped count", n_dropped == n_dropped_ok);
    for (int i = 0; i < n; ++i)
    {
        mu_assert("depth", depth[i] == depth_ok[i]);
    }
    for (int i = 0; i < n_order; ++i)
    {
        mu_assert("order", order[i] == order_ok[i]);
    }
    for (int i = 0; i < n_dropped_ok; ++i)
    {
        mu_assert("dropped", dropped[i] == dropped_ok[i]);
    }
    return 0;
}

/* 0 -> 1 -> 2 -> nothing */
static char *test_chains_path()
{
    int succ[] = {1, 2, -1}, drop_priority[] = {0, 1, 2};
    int depth[] = {2, 1, 0}, order[] = {0, 1, 2};
    return check_chain_depths(3, succ, drop_priority, depth, order, 3, NULL, 0);
}

/* two heads onto one node: 0 -> 2, 1 -> 2, 2 -> 3 -> nothing */
static char *test_chains_tree()
{
    int succ[] = {2, 2, 3, -1}, drop_priority[] = {0, 1, 2, 3};
    int depth[] = {2, 2, 1, 0}, order[] = {0, 1, 2, 3};
    return check_chain_depths(4, succ, drop_priority, depth, order, 4, NULL, 0);
}

/* a node reached after its chain was resolved by an earlier walk: 3 -> 1 */
static char *test_chains_memoised()
{
    int succ[] = {1, 2, -1, 1}, drop_priority[] = {0, 1, 2, 3};
    int depth[] = {2, 1, 0, 2}, order[] = {0, 3, 1, 2};
    return check_chain_depths(4, succ, drop_priority, depth, order, 4, NULL, 0);
}

/* 2-cycle 0 <-> 1: node 1 has the higher drop priority and is dropped, 0 ends
   at it */
static char *test_chains_two_cycle()
{
    int succ[] = {1, 0}, drop_priority[] = {0, 1};
    int depth[] = {0, CHAINS_DROPPED}, order[] = {0}, dropped[] = {1};
    return check_chain_depths(2, succ, drop_priority, depth, order, 1, dropped, 1);
}

/* the highest drop priority sits at the start of the cycle: node 0 is dropped even
   though the walk started there, and node 1 then resolves onto it */
static char *test_chains_two_cycle_drop_start()
{
    int succ[] = {1, 0}, drop_priority[] = {5, 1};
    int depth[] = {CHAINS_DROPPED, 0}, order[] = {1}, dropped[] = {0};
    return check_chain_depths(2, succ, drop_priority, depth, order, 1, dropped, 1);
}

/* a tail into a cycle: 0 -> 1 -> 2 -> 1. Node 2 is dropped, 1 ends at it,
   0 keeps its depth above 1 */
static char *test_chains_tail_into_cycle()
{
    int succ[] = {1, 2, 1}, drop_priority[] = {0, 1, 2};
    int depth[] = {1, 0, CHAINS_DROPPED}, order[] = {0, 1}, dropped[] = {2};
    return check_chain_depths(3, succ, drop_priority, depth, order, 2, dropped, 1);
}

/* every node points to nothing: identity order, depth 0 */
static char *test_chains_roots()
{
    int succ[] = {-1, -1, -1}, drop_priority[] = {0, 1, 2};
    int depth[] = {0, 0, 0}, order[] = {0, 1, 2};
    return check_chain_depths(3, succ, drop_priority, depth, order, 3, NULL, 0);
}

static const char *all_tests_chains()
{
    mu_run_test(test_chains_path, counter_chains);
    mu_run_test(test_chains_tree, counter_chains);
    mu_run_test(test_chains_memoised, counter_chains);
    mu_run_test(test_chains_two_cycle, counter_chains);
    mu_run_test(test_chains_two_cycle_drop_start, counter_chains);
    mu_run_test(test_chains_tail_into_cycle, counter_chains);
    mu_run_test(test_chains_roots, counter_chains);
    return 0;
}

int test_chains()
{
    const char *result = all_tests_chains();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("chains: TEST FAILED!\n");
    }
    else
    {
        printf("chains: ALL TESTS PASSED\n");
    }
    printf("chains: Tests run: %d\n", counter_chains);
    return result == 0;
}

#endif // TEST_CHAINS_H
