#ifndef TEST_CONES_H
#define TEST_CONES_H

#include "Cones.h"
#include "PSLP_API.h"
#include "glbopts.h"
#include "kkt.h"
#include "minunit.h"
#include "test_kkt.h"
#include <stdio.h>

static int counter_cones = 0;

static char *test_cone_table()
{
    const ConeTypeInfo *soc = cone_type_info(PSLP_CONE_SOC);
    const ConeTypeInfo *rotated = cone_type_info(PSLP_CONE_ROTATED);
    mu_assert("SOC row", soc != NULL && soc->dim_min == 1);
    mu_assert("ROTATED row", rotated != NULL && rotated->dim_min == 2);
    mu_assert("no EXP row", cone_type_info(PSLP_CONE_EXP) == NULL);
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
    The block comes back as the last row of the reduced problem with
    lhs = rhs = r; a malformed block list is refused. */
static char *test_cones_accepted()
{
    double Ax[] = {1, 1, 1, -1};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    double lhs[] = {-INF, -INF};
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
    Presolver *presolver = new_presolver_from_problem(&problem, stgs);
    mu_assert("cone block accepted", presolver != NULL);
    PresolveStatus status = run_presolver(presolver);
    mu_assert("presolved", status == UNCHANGED || status == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("cone block returned",
              reduced->n_cones == 1 && reduced->cones[0].type == PSLP_CONE_SOC &&
                  reduced->cones[0].dim == 1 &&
                  (size_t) reduced->cones[0].start == reduced->m - 1);
    mu_assert("cone row sides", reduced->lhs[reduced->m - 1] == 1 &&
                                    reduced->rhs[reduced->m - 1] == 1);
    mu_assert("cone stats", presolver->stats->n_cones_original == 1 &&
                                presolver->stats->n_cone_rows_original == 1 &&
                                presolver->stats->n_cones_reduced == 1 &&
                                presolver->stats->n_cone_rows_reduced == 1);
    free_presolver(presolver);

    ConeBlock bad[] = {{PSLP_CONE_SOC, 0, 1}};
    problem.cones = bad;
    mu_assert("bad block refused",
              new_presolver_from_problem(&problem, stgs) == NULL);
    free_settings(stgs);
    return 0;
}

/*  The worked example of test_kkt.h through presolve: the doubleton row
    eliminates x0 or x3 into the cone block and x4 is fixed, leaving the SOC
    alone on three columns. The KKT point maps forward, checks on the reduced
    problem, and postsolves back to the original. */
static char *test_cones_soc_flow()
{
    KKT_TEST_SOC
    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("reduced dims", reduced->m == 3 && reduced->n == 3);
    mu_assert("reduced block",
              reduced->n_cones == 1 && reduced->cones[0].type == PSLP_CONE_SOC &&
                  reduced->cones[0].start == 0 && reduced->cones[0].dim == 3);

    double x_red[3], y_red[3], z_red[3];
    map_original_sol_to_reduced(presolver, x, y, x_red, y_red, z_red);
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));

    postsolve(presolver, x_red, y_red, z_red);
    mu_assert("original KKT point",
              is_kkt_point(&prob, presolver->sol->x, presolver->sol->y,
                           presolver->sol->z, 1e-9));
    mu_assert("linear row multiplier", ABS(presolver->sol->y[0] - 1) <= 1e-9);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*  min 2 x0 + x1  s.t.  (x0, x1, x2) in ROTATED,  x2 = 2 (lb = ub).
    Fixing x2 shifts r of the third cone row to -2 and empties the row, which
    stays as a constant component. Optimum x = (1, 2, 2) on the boundary with
    y = (2, 1, -2) in the dual cone and objective 4 (exact in floating point,
    so the KKT checker sees a zero defect). */
static char *test_cones_rotated_flow()
{
    double Ax[] = {1, 1, 1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 1, 2, 3};
    double lhs[] = {-INF, -INF, -INF};
    double rhs[] = {0, 0, 0};
    double lbs[] = {-INF, -INF, 2};
    double ubs[] = {INF, INF, 2};
    double c[] = {2, 1, 0};
    ConeBlock cones[] = {{PSLP_CONE_ROTATED, 0, 3}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 3, 3, 3, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("reduced dims", reduced->m == 3 && reduced->n == 2);
    mu_assert("reduced block",
              reduced->n_cones == 1 && reduced->cones[0].type == PSLP_CONE_ROTATED &&
                  reduced->cones[0].start == 0 && reduced->cones[0].dim == 3);
    mu_assert("empty cone row kept", reduced->Ap[3] == reduced->Ap[2] &&
                                         reduced->lhs[2] == -2 &&
                                         reduced->rhs[2] == -2);

    double x_red[] = {1, 2};
    double y_red[] = {2, 1, -2};
    double z_red[] = {0, 0};
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));

    postsolve(presolver, x_red, y_red, z_red);
    mu_assert("original KKT point",
              is_kkt_point(&prob, presolver->sol->x, presolver->sol->y,
                           presolver->sol->z, 1e-9));

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/* The two ray postsolves refuse a problem that had cones and leave their
   output untouched. */
static char *test_cones_refusals()
{
    KKT_TEST_SOC
    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);

    double y_ray[] = {0, 0, 0};
    double y_orig[] = {7, 7, 7, 7};
    postsolve_primal_infeas_ray(presolver, y_ray, y_orig);
    mu_assert("primal ray untouched",
              y_orig[0] == 7 && y_orig[1] == 7 && y_orig[2] == 7 && y_orig[3] == 7);

    double x_ray[] = {0, 0, 0};
    double x_orig[] = {7, 7, 7, 7, 7};
    postsolve_dual_infeas_ray(presolver, x_ray, x_orig);
    mu_assert("dual ray untouched", x_orig[0] == 7 && x_orig[1] == 7 &&
                                        x_orig[2] == 7 && x_orig[3] == 7 &&
                                        x_orig[4] == 7);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

static const char *all_tests_cones()
{
    mu_run_test(test_cone_table, counter_cones);
    mu_run_test(test_cones_valid, counter_cones);
    mu_run_test(test_cones_accepted, counter_cones);
    mu_run_test(test_cones_soc_flow, counter_cones);
    mu_run_test(test_cones_rotated_flow, counter_cones);
    mu_run_test(test_cones_refusals, counter_cones);
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
