#ifndef TEST_CONES_H
#define TEST_CONES_H

#include "Cones.h"
#include "Numerics.h"
#include "PSLP_API.h"
#include "glbopts.h"
#include "kkt.h"
#include "minunit.h"
#include <stdio.h>

static int counter_cones = 0;

/* FEAS_TOL scaled by the block, max(1, ||s||_inf). */
static double cone_tol(const double *s, int dim)
{
    double scale = 1.0;
    for (int ii = 0; ii < dim; ++ii)
    {
        scale = MAX(scale, ABS(s[ii]));
    }
    return FEAS_TOL * scale;
}

static bool in_cone(ConeType type, const double *s, int dim)
{
    return cone_member(type, s, dim, cone_tol(s, dim));
}

static bool in_dual_cone(ConeType type, const double *s, int dim)
{
    return cone_dual_member(type, s, dim, cone_tol(s, dim));
}

static char *test_cone_table()
{
    const ConeTypeInfo *soc = cone_type_info(PSLP_CONE_SOC);
    const ConeTypeInfo *rotated = cone_type_info(PSLP_CONE_ROTATED);
    mu_assert("SOC row", soc != NULL && soc->dim_min == 1);
    mu_assert("ROTATED row", rotated != NULL && rotated->dim_min == 2);
    mu_assert("no EXP row", cone_type_info(PSLP_CONE_EXP) == NULL);
    return 0;
}

static char *test_soc_member()
{
    double inside[] = {6, 3, 4};
    double on[] = {5, 3, 4};
    double outside[] = {5, 3, 4.0001};
    double negative_head[] = {-1, 0, 0};
    double scaled_on[] = {5e6, 3e6, 4e6 + 1};
    double scaled_outside[] = {5e6, 3e6, 4e6 + 10};
    double zero = 0;
    double small_negative = -1e-7;
    double negative = -1;

    mu_assert("inside", in_cone(PSLP_CONE_SOC, inside, 3));
    mu_assert("on", in_cone(PSLP_CONE_SOC, on, 3));
    mu_assert("outside", !in_cone(PSLP_CONE_SOC, outside, 3));
    mu_assert("negative head", !in_cone(PSLP_CONE_SOC, negative_head, 3));
    mu_assert("scaled on", in_cone(PSLP_CONE_SOC, scaled_on, 3));
    mu_assert("scaled outside", !in_cone(PSLP_CONE_SOC, scaled_outside, 3));
    mu_assert("dim 1 zero", in_cone(PSLP_CONE_SOC, &zero, 1));
    mu_assert("dim 1 small negative", in_cone(PSLP_CONE_SOC, &small_negative, 1));
    mu_assert("dim 1 negative", !in_cone(PSLP_CONE_SOC, &negative, 1));

    mu_assert("dual inside", in_dual_cone(PSLP_CONE_SOC, inside, 3));
    mu_assert("dual on", in_dual_cone(PSLP_CONE_SOC, on, 3));
    mu_assert("dual outside", !in_dual_cone(PSLP_CONE_SOC, outside, 3));
    mu_assert("dual negative head", !in_dual_cone(PSLP_CONE_SOC, negative_head, 3));
    mu_assert("EXP has no test", !in_cone(PSLP_CONE_EXP, inside, 3));
    return 0;
}

static char *test_rotated_member()
{
    double inside[] = {2, 2, 2, 0};
    double on[] = {1, 2, 2, 0};
    double outside[] = {1, 2, 2, 0.01};
    double negative_first_head[] = {-1, 2, 0, 0};
    double negative_second_head[] = {1, -2, 0, 0};
    double zero_head[] = {0, 5, 0, 0};
    double zero_heads_tail[] = {0, 0, 1};
    double dim2_zero[] = {0, 0};
    double dim2_inside[] = {1, 1};
    double dim2_negative[] = {1, -1};

    mu_assert("inside", in_cone(PSLP_CONE_ROTATED, inside, 4));
    mu_assert("on", in_cone(PSLP_CONE_ROTATED, on, 4));
    mu_assert("outside", !in_cone(PSLP_CONE_ROTATED, outside, 4));
    mu_assert("negative first head",
              !in_cone(PSLP_CONE_ROTATED, negative_first_head, 4));
    mu_assert("negative second head",
              !in_cone(PSLP_CONE_ROTATED, negative_second_head, 4));
    mu_assert("zero head", in_cone(PSLP_CONE_ROTATED, zero_head, 4));
    mu_assert("zero heads, tail", !in_cone(PSLP_CONE_ROTATED, zero_heads_tail, 3));
    mu_assert("dim 2 zero", in_cone(PSLP_CONE_ROTATED, dim2_zero, 2));
    mu_assert("dim 2 inside", in_cone(PSLP_CONE_ROTATED, dim2_inside, 2));
    mu_assert("dim 2 negative", !in_cone(PSLP_CONE_ROTATED, dim2_negative, 2));

    mu_assert("dual inside", in_dual_cone(PSLP_CONE_ROTATED, inside, 4));
    mu_assert("dual on", in_dual_cone(PSLP_CONE_ROTATED, on, 4));
    mu_assert("dual outside", !in_dual_cone(PSLP_CONE_ROTATED, outside, 4));
    mu_assert("dual negative head",
              !in_dual_cone(PSLP_CONE_ROTATED, negative_first_head, 4));
    return 0;
}

/* Block lists on a problem with 5 rows. */
static char *test_cones_valid()
{
    size_t m = 5;
    double rhs[] = {1, 2, 3, 4, 5};

    ConeBlock one_soc[] = {{PSLP_CONE_SOC, 2, 3}};
    mu_assert("one SOC block", cones_valid(one_soc, 1, m, rhs));

    ConeBlock two_blocks[] = {{PSLP_CONE_ROTATED, 1, 2}, {PSLP_CONE_SOC, 3, 2}};
    mu_assert("two blocks", cones_valid(two_blocks, 2, m, rhs));

    ConeBlock all_rows[] = {{PSLP_CONE_SOC, 0, 5}};
    mu_assert("all rows", cones_valid(all_rows, 1, m, rhs));

    mu_assert("no blocks", cones_valid(NULL, 0, m, rhs));
    mu_assert("NULL with a count", !cones_valid(NULL, 1, m, rhs));

    ConeBlock overlap[] = {{PSLP_CONE_SOC, 2, 2}, {PSLP_CONE_SOC, 3, 2}};
    mu_assert("overlap", !cones_valid(overlap, 2, m, rhs));

    ConeBlock gap[] = {{PSLP_CONE_SOC, 1, 1}, {PSLP_CONE_SOC, 3, 2}};
    mu_assert("gap", !cones_valid(gap, 2, m, rhs));

    ConeBlock not_last[] = {{PSLP_CONE_SOC, 0, 2}};
    mu_assert("not the last rows", !cones_valid(not_last, 1, m, rhs));

    ConeBlock past_end[] = {{PSLP_CONE_SOC, 3, 3}};
    mu_assert("past the last row", !cones_valid(past_end, 1, m, rhs));

    ConeBlock negative_start[] = {{PSLP_CONE_SOC, -1, 6}};
    mu_assert("negative start", !cones_valid(negative_start, 1, m, rhs));

    ConeBlock rotated_dim_1[] = {{PSLP_CONE_ROTATED, 4, 1}};
    mu_assert("ROTATED dim 1", !cones_valid(rotated_dim_1, 1, m, rhs));

    ConeBlock exp_block[] = {{PSLP_CONE_EXP, 2, 3}};
    mu_assert("EXP", !cones_valid(exp_block, 1, m, rhs));

    double infinite_r[] = {1, 2, 3, 4, INF};
    mu_assert("infinite r", !cones_valid(one_soc, 1, m, infinite_r));
    ConeBlock before_infinite[] = {{PSLP_CONE_SOC, 4, 1}};
    double infinite_linear[] = {INF, 2, 3, 4, 5};
    mu_assert("infinite linear rhs",
              cones_valid(before_infinite, 1, m, infinite_linear));
    return 0;
}

/*  min. [1 1] x
    s.t.  [1 1] x <= 4
          [1 -1] x - 1 in SOC (dimension 1)
          0 <= x <= 3
    A problem with a cone block is refused; the same problem without the block
    is presolved. */
static char *test_cones_refused()
{
    double Ax[] = {1, 1, 1, -1};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    double lhs[] = {-INF, 1};
    double rhs[] = {4, 1};
    double lbs[] = {0, 0};
    double ubs[] = {3, 3};
    double c[] = {1, 1};
    ConeBlock cones[] = {{PSLP_CONE_SOC, 1, 1}};
    Settings *stgs = default_settings();
    stgs->verbose = false;

    ProblemData problem =
        problem_from_csr(Ax, Ai, Ap, 2, 2, 4, lhs, rhs, lbs, ubs, c);
    problem.cones = cones;
    problem.n_cones = 1;
    mu_assert("cone block refused",
              new_presolver_from_problem(&problem, stgs) == NULL);

    ConeBlock bad[] = {{PSLP_CONE_SOC, 0, 1}};
    problem.cones = bad;
    mu_assert("bad block refused",
              new_presolver_from_problem(&problem, stgs) == NULL);

    problem.cones = NULL;
    problem.n_cones = 0;
    Presolver *presolver = new_presolver_from_problem(&problem, stgs);
    mu_assert("no cones accepted", presolver != NULL);
    PresolveStatus status = run_presolver(presolver);
    mu_assert("presolved", status == UNCHANGED || status == REDUCED);
    mu_assert("cone stats zero", presolver->stats->n_cones_original == 0 &&
                                     presolver->stats->n_cone_rows_original == 0);
    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

static const char *all_tests_cones()
{
    mu_run_test(test_cone_table, counter_cones);
    mu_run_test(test_soc_member, counter_cones);
    mu_run_test(test_rotated_member, counter_cones);
    mu_run_test(test_cones_valid, counter_cones);
    mu_run_test(test_cones_refused, counter_cones);
    return 0;
}

int test_cones()
{
    const char *result = all_tests_cones();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("cones: TEST FAILED!\n");
    }
    else
    {
        printf("cones: ALL TESTS PASSED\n");
    }
    printf("cones: Tests run: %d\n", counter_cones);
    return result == 0;
}

#endif
