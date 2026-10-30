#ifndef TEST_CONES_H
#define TEST_CONES_H

#include "Cones.h"
#include "Constraints.h"
#include "PSLP_API.h"
#include "Problem.h"
#include "State.h"
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
    const ConeTypeInfo *exp_info = cone_type_info(PSLP_CONE_EXP);
    mu_assert("EXP row",
              exp_info != NULL && exp_info->dim_min == 3 && exp_info->dim_max == 3);
    mu_assert("no row for 0", cone_type_info((ConeType) 0) == NULL);
    mu_assert("no row for 99", cone_type_info((ConeType) 99) == NULL);
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
    mu_assert("EXP", cones_valid(exp_block, 1, m, rhs));

    ConeBlock exp_dim_2[] = {{PSLP_CONE_EXP, 3, 2}};
    mu_assert("EXP dim 2", !cones_valid(exp_dim_2, 1, m, rhs));

    ConeBlock exp_dim_4[] = {{PSLP_CONE_EXP, 1, 4}};
    mu_assert("EXP dim 4", !cones_valid(exp_dim_4, 1, m, rhs));

    ConeBlock exp_and_soc[] = {{PSLP_CONE_EXP, 0, 3}, {PSLP_CONE_SOC, 3, 2}};
    mu_assert("EXP and SOC", cones_valid(exp_and_soc, 2, m, rhs));

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

/*  min x0  s.t.  x1 + x2 >= 25,  (1, x0, 3 x1, 4 x2) in ROTATED (the first
    cone row is empty with r = -1). The block becomes P = diag(9, 16):
    min (1/2)(9 x1^2 + 16 x2^2) s.t. x1 + x2 >= 25, optimum x = (16, 9),
    y = 144. Postsolve gives x0 = 1800 and y = (144, 1800, 1, -48, -36);
    the tails (48, 36) have an exact norm, so the KKT checker sees no
    rounding. */
static char *test_cones_rotated_to_P()
{
    double Ax[] = {1, 1, 1, 3, 4};
    int Ai[] = {1, 2, 0, 1, 2};
    int Ap[] = {0, 2, 2, 3, 4, 5};
    double lhs[] = {25, -INF, -INF, -INF, -INF};
    double rhs[] = {INF, -1, 0, 0, 0};
    double lbs[] = {-INF, -INF, -INF};
    double ubs[] = {INF, INF, INF};
    double c[] = {1, 0, 0};
    ConeBlock cones[] = {{PSLP_CONE_ROTATED, 1, 4}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 5, 3, 5, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("reduced dims", reduced->m == 1 && reduced->n == 2 &&
                                  reduced->n_cones == 0 && reduced->Pnnz == 2);
    mu_assert("reduced P", reduced->Pp[1] == 1 && reduced->Pp[2] == 2 &&
                               reduced->Pi[0] == 0 && reduced->Pi[1] == 1 &&
                               reduced->Px[0] == 9 && reduced->Px[1] == 16);
    mu_assert("reduced objective",
              reduced->c[0] == 0 && reduced->c[1] == 0 && reduced->obj_offset == 0);

    double x_red[] = {16, 9};
    double y_red[] = {144};
    double z_red[] = {0, 0};
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));

    postsolve(presolver, x_red, y_red, z_red);
    const double *x = presolver->sol->x, *y = presolver->sol->y;
    mu_assert("original KKT point",
              is_kkt_point(&prob, x, y, presolver->sol->z, 1e-9));
    mu_assert("postsolved point", x[0] == 1800 && y[0] == 144 && y[1] == 1800 &&
                                      y[2] == 1 && y[3] == -48 && y[4] == -36);

    double x_map[] = {DUMMY_VALUE, DUMMY_VALUE};
    double y_map[] = {DUMMY_VALUE};
    double z_map[] = {DUMMY_VALUE, DUMMY_VALUE};
    double x_orig[] = {1800, 16, 9};
    double y_orig[] = {144, 1800, 1, -48, -36};
    map_original_sol_to_reduced(presolver, x_orig, y_orig, x_map, y_map, z_map);
    mu_assert("map refused", x_map[0] == DUMMY_VALUE && y_map[0] == DUMMY_VALUE &&
                                 z_map[0] == DUMMY_VALUE);
    free_presolver(presolver);

    stgs->rotated_to_P = false;
    presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init without conversion", presolver != NULL);
    run_presolver(presolver);
    mu_assert("block kept", presolver->stats->n_cones_reduced == 1);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*  min x0 - x1  s.t.  x2 + x3 >= 24,
        (1, x0, x2 / 4) in ROTATED              (w first, a = 1, kappa = 1),
        (-x1 - 1, 2, x2 / 2, x3 - 2) in ROTATED (w second, a = -1, kappa = 1/2).
    x2 is a tail of both blocks. Reduced: P = diag(3/16, 1/2), c = (0, -1),
    offset 2, optimum x = (16, 8), y = 3, objective 34. Postsolve gives
    x0 = 8, x1 = -26, y = (3, 8, 1, -4, 1, 25/2, -4, -3); the tails of the
    second block are (8, 6), with an exact norm. */
static char *test_cones_rotated_to_P_signs()
{
    double Ax[] = {1, 1, 1, 0.25, -1, 0.5, 1};
    int Ai[] = {2, 3, 0, 2, 1, 2, 3};
    int Ap[] = {0, 2, 2, 3, 4, 5, 5, 6, 7};
    double lhs[] = {24, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {INF, -1, 0, 0, 1, -2, 0, 2};
    double lbs[] = {-INF, -INF, -INF, -INF};
    double ubs[] = {INF, INF, INF, INF};
    double c[] = {1, -1, 0, 0};
    ConeBlock cones[] = {{PSLP_CONE_ROTATED, 1, 3}, {PSLP_CONE_ROTATED, 4, 4}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 8, 4, 7, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 2;
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("reduced dims", reduced->m == 1 && reduced->n == 2 &&
                                  reduced->n_cones == 0 && reduced->Pnnz == 2);
    mu_assert("reduced P", reduced->Px[0] == 0.1875 && reduced->Px[1] == 0.5);
    mu_assert("reduced objective",
              reduced->c[0] == 0 && reduced->c[1] == -1 && reduced->obj_offset == 2);

    double x_red[] = {16, 8};
    double y_red[] = {3};
    double z_red[] = {0, 0};
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));

    postsolve(presolver, x_red, y_red, z_red);
    const double *x = presolver->sol->x, *y = presolver->sol->y;
    mu_assert("original KKT point",
              is_kkt_point(&prob, x, y, presolver->sol->z, 1e-9));
    mu_assert("postsolved x", x[0] == 8 && x[1] == -26);
    mu_assert("postsolved y", y[1] == 8 && y[2] == 1 && y[3] == -4 && y[4] == 1 &&
                                  y[5] == 12.5 && y[6] == -4 && y[7] == -3);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*  min x0 + x2  s.t.  (x0, 0, x1) in SOC,  (x2, x3, 0) in SOC,  x1 = 2,
    x3 = 1 (lb = ub). The empty rows with a zero side go, one inside its
    block and one last in its block; the fixed columns leave constant rows
    r = -2 and r = -1. Reduced: two blocks of dim 2, x0 >= 2 and x2 >= 1,
    optimum x = (2, 1), y = (1, -1, 1, -1). Postsolve gives
    y = (1, 0, -1, 1, -1, 0). */
static char *test_cones_zero_tails()
{
    double Ax[] = {1, 1, 1, 1};
    int Ai[] = {0, 1, 2, 3};
    int Ap[] = {0, 1, 1, 2, 3, 4, 4};
    double lhs[] = {-INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {0, 0, 0, 0, 0, 0};
    double lbs[] = {-INF, 2, -INF, 1};
    double ubs[] = {INF, 2, INF, 1};
    double c[] = {1, 0, 1, 0};
    ConeBlock cones[] = {{PSLP_CONE_SOC, 0, 3}, {PSLP_CONE_SOC, 3, 3}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 6, 4, 4, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 2;
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("reduced dims",
              reduced->m == 4 && reduced->n == 2 && reduced->n_cones == 2);
    mu_assert("reduced blocks",
              reduced->cones[0].start == 0 && reduced->cones[0].dim == 2 &&
                  reduced->cones[1].start == 2 && reduced->cones[1].dim == 2);
    mu_assert("constant rows kept",
              reduced->Ap[2] == reduced->Ap[1] && reduced->rhs[1] == -2 &&
                  reduced->Ap[4] == reduced->Ap[3] && reduced->rhs[3] == -1);

    double x_orig[] = {2, 2, 1, 1};
    double y_orig[] = {1, 0, -1, 1, -1, 0};
    double x_red[2], y_red[4], z_red[2];
    map_original_sol_to_reduced(presolver, x_orig, y_orig, x_red, y_red, z_red);
    mu_assert("mapped point", x_red[0] == 2 && x_red[1] == 1 && y_red[0] == 1 &&
                                  y_red[1] == -1 && y_red[2] == 1 && y_red[3] == -1);
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));

    postsolve(presolver, x_red, y_red, z_red);
    const double *x = presolver->sol->x, *y = presolver->sol->y;
    mu_assert("original KKT point",
              is_kkt_point(&prob, x, y, presolver->sol->z, 1e-9));
    mu_assert("postsolved point", x[0] == 2 && x[2] == 1 && y[1] == 0 && y[5] == 0 &&
                                      y[2] == -1 && y[4] == -1);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*  min x0 + 2 x1  s.t.  x0 + x1 >= 1,  x0 - x1 >= -3,  (x0, x1, 0) in ROTATED,
    x free. The block is x0 >= 0, x1 >= 0 and goes. Reduced: the two linear
    rows with lbs = 0, optimum x = (1, 0), y = (1, 0), z = (0, 1). Postsolve
    moves z1 to the second head: y = (1, 0, 0, 1, 0). */
static char *test_cones_block_to_bounds()
{
    double Ax[] = {1, 1, 1, -1, 1, 1};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 5, 6, 6};
    double lhs[] = {1, -3, -INF, -INF, -INF};
    double rhs[] = {INF, INF, 0, 0, 0};
    double lbs[] = {-INF, -INF};
    double ubs[] = {INF, INF};
    double c[] = {1, 2};
    ConeBlock cones[] = {{PSLP_CONE_ROTATED, 2, 3}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 5, 2, 6, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("reduced dims",
              reduced->m == 2 && reduced->n == 2 && reduced->n_cones == 0);
    mu_assert("heads became bounds", reduced->lbs[0] == 0 && reduced->lbs[1] == 0);

    double x_orig[] = {1, 0};
    double y_orig[] = {1, 0, 0, 1, 0};
    double x_red[2], y_red[2], z_red[2];
    map_original_sol_to_reduced(presolver, x_orig, y_orig, x_red, y_red, z_red);
    mu_assert("mapped point", x_red[0] == 1 && x_red[1] == 0 && y_red[0] == 1 &&
                                  y_red[1] == 0 && z_red[0] == 0 && z_red[1] == 1);
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));

    postsolve(presolver, x_red, y_red, z_red);
    const double *y = presolver->sol->y, *z = presolver->sol->z;
    mu_assert("original KKT point",
              is_kkt_point(&prob, presolver->sol->x, y, z, 1e-9));
    mu_assert("bound dual moved to the head",
              y[2] == 0 && y[3] == 1 && y[4] == 0 && z[0] == 0 && z[1] == 0);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*  A block of empty rows: (-r) in SOC decides feasibility. With r = (1, 0, 0)
    the problem is infeasible, with r = (-1, 0, 0) the block goes with y = 0.
    The one column has x0 >= 1 and cost 1, so the reduced problem is empty. */
static char *test_cones_constant_block()
{
    double Ax[] = {1};
    int Ai[] = {0};
    int Ap[] = {0, 1, 1, 1, 1};
    double lhs[] = {2, -INF, -INF, -INF};
    double rhs[] = {INF, 1, 0, 0};
    double lbs[] = {-INF};
    double ubs[] = {INF};
    double c[] = {1};
    ConeBlock cones[] = {{PSLP_CONE_SOC, 1, 3}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 4, 1, 1, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("infeasible", run_presolver(presolver) == INFEASIBLE);
    free_presolver(presolver);

    rhs[1] = -1;
    presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init feasible", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    mu_assert("solved by presolve", presolver->reduced_prob->m == 0 &&
                                        presolver->reduced_prob->n == 0 &&
                                        presolver->reduced_prob->n_cones == 0);

    double x_red[1], y_red[1], z_red[1]; /* the reduced problem is empty */
    postsolve(presolver, x_red, y_red, z_red);
    const double *y = presolver->sol->y;
    mu_assert("original KKT point",
              is_kkt_point(&prob, presolver->sol->x, y, presolver->sol->z, 1e-9));
    mu_assert("block multipliers", y[1] == 0 && y[2] == 0 && y[3] == 0);

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

/*  test_cones_rotated_to_P with a zero tail appended to the block. One visit
    deletes the tail and converts the block, and postsolve gives the deleted
    tail the multiplier 0. */
static char *test_cones_chain()
{
    double Ax[] = {1, 1, 1, 3, 4};
    int Ai[] = {1, 2, 0, 1, 2};
    int Ap[] = {0, 2, 2, 3, 4, 5, 5};
    double lhs[] = {25, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {INF, -1, 0, 0, 0, 0};
    double lbs[] = {-INF, -INF, -INF};
    double ubs[] = {INF, INF, INF};
    double c[] = {1, 0, 0};
    ConeBlock cones[] = {{PSLP_CONE_ROTATED, 1, 5}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 6, 3, 5, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("reduced dims", reduced->m == 1 && reduced->n == 2 &&
                                  reduced->n_cones == 0 && reduced->Pnnz == 2);
    mu_assert("reduced P", reduced->Px[0] == 9 && reduced->Px[1] == 16);

    double x_red[] = {16, 9};
    double y_red[] = {144};
    double z_red[] = {0, 0};
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));

    postsolve(presolver, x_red, y_red, z_red);
    const double *x = presolver->sol->x, *y = presolver->sol->y;
    mu_assert("original KKT point",
              is_kkt_point(&prob, x, y, presolver->sol->z, 1e-9));
    mu_assert("postsolved point", x[0] == 1800 && y[0] == 144 && y[1] == 1800 &&
                                      y[2] == 1 && y[3] == -48 && y[4] == -36 &&
                                      y[5] == 0);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*  The block of test_cones_chain visited by hand with the conversion off: the
    zero tail goes and stays a hole (the block keeps its dim until the
    cleanup), a second visit changes nothing, and a visit with the conversion
    on converts across the hole. The loop alone, no presolver run after it. */
static char *test_cones_hole()
{
    double Ax[] = {1, 1, 1, 3, 4};
    int Ai[] = {1, 2, 0, 1, 2};
    int Ap[] = {0, 2, 2, 3, 4, 5, 5};
    double lhs[] = {25, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {INF, -1, 0, 0, 0, 0};
    double lbs[] = {-INF, -INF, -INF};
    double ubs[] = {INF, INF, INF};
    double c[] = {1, 0, 0};
    ConeBlock cones[] = {{PSLP_CONE_ROTATED, 1, 5}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 6, 3, 5, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->rotated_to_P = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    Problem *p = presolver->prob;
    const Constraints *constraints = p->constraints;
    const int *row_sizes = constraints->state->row_sizes;

    mu_assert("first visit", process_cones(p, stgs, presolver->stats) == UNCHANGED);
    mu_assert("hole", row_sizes[5] == SIZE_INACTIVE_ROW && row_sizes[1] == 0 &&
                          row_sizes[2] == 1 && row_sizes[3] == 1 &&
                          row_sizes[4] == 1 && constraints->cones[0].dim == 5);
    mu_assert("second visit", process_cones(p, stgs, presolver->stats) == UNCHANGED);
    mu_assert("nothing more",
              row_sizes[1] == 0 && row_sizes[4] == 1 && p->obj->P == NULL);

    stgs->rotated_to_P = true;
    mu_assert("third visit", process_cones(p, stgs, presolver->stats) == UNCHANGED);
    mu_assert("converted across the hole", row_sizes[1] == SIZE_INACTIVE_ROW &&
                                               row_sizes[2] == SIZE_INACTIVE_ROW &&
                                               row_sizes[3] == SIZE_INACTIVE_ROW &&
                                               row_sizes[4] == SIZE_INACTIVE_ROW &&
                                               p->obj->P != NULL);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*  min x0 - x3  s.t.  x1 - 2 x3 = 0,  (x0, x2, x1) in EXP,  x2 = 1 (lb = ub),
    0.5 <= x3 <= 1. The LP reductions work on the block's rows: the fixed x2
    leaves the second cone row empty with r = -1, and the doubleton row
    eliminates x1 or x3 into the third. The block itself is passed through
    with its three rows. Optimum x = (e, 1, 1, 1/2): s = (e, 1, 1) is on the
    boundary with the supporting dual y_B = (1, 0, -e), y_0 = e,
    z = (0, 0, 0, 2e - 1), objective e - 1/2. */
static char *test_cones_exp_flow()
{
    double Ax[] = {1, -2, 1, 1, 1};
    int Ai[] = {1, 3, 0, 2, 1};
    int Ap[] = {0, 2, 3, 4, 5};
    double lhs[] = {0, -INF, -INF, -INF};
    double rhs[] = {0, 0, 0, 0};
    double lbs[] = {-INF, -INF, 1, 0.5};
    double ubs[] = {INF, INF, 1, 1};
    double c[] = {1, 0, 0, -1};
    ConeBlock cones[] = {{PSLP_CONE_EXP, 1, 3}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 4, 4, 5, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    double e = exp(1.0);
    double x[] = {e, 1, 1, 0.5};
    double y[] = {e, 1, 0, -e};
    double z[] = {0, 0, 0, 2 * e - 1};
    mu_assert("original KKT point", is_kkt_point(&prob, x, y, z, 1e-9));
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("reduced dims", reduced->m == 3 && reduced->n == 2);
    mu_assert("block passed through",
              reduced->n_cones == 1 && reduced->cones[0].type == PSLP_CONE_EXP &&
                  reduced->cones[0].start == 0 && reduced->cones[0].dim == 3);
    mu_assert("constant row kept", reduced->Ap[2] == reduced->Ap[1] &&
                                       reduced->lhs[1] == -1 &&
                                       reduced->rhs[1] == -1);
    mu_assert("cone stats", presolver->stats->n_cones_reduced == 1 &&
                                presolver->stats->n_cone_rows_reduced == 3);

    double x_red[2], y_red[3], z_red[2];
    map_original_sol_to_reduced(presolver, x, y, x_red, y_red, z_red);
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));

    postsolve(presolver, x_red, y_red, z_red);
    const double *xs = presolver->sol->x, *ys = presolver->sol->y;
    mu_assert("original KKT point after postsolve",
              is_kkt_point(&prob, xs, ys, presolver->sol->z, 1e-9));
    mu_assert("postsolved point", ABS(xs[0] - e) <= 1e-12 && xs[1] == 1 &&
                                      xs[2] == 1 && xs[3] == 0.5 &&
                                      ABS(ys[0] - e) <= 1e-12 && ys[1] == 1 &&
                                      ys[2] == 0 && ABS(ys[3] + e) <= 1e-12);

    free_presolver(presolver);
    free_settings(stgs);
    return 0;
}

/*  min x0  s.t.  x0 >= 1,  (3, 1, 1) in EXP as a block of three empty rows
    (r = (-3, -1, -1)). A constant SOC block goes (test_cones_constant_block);
    a constant EXP block stays, its multiplier 0, and the column is solved by
    presolve. With r = (-1, -1, -1) the block is infeasible (1 < e), which
    the presolver does not detect: the solver sees the constant block. */
static char *test_cones_exp_constant_block()
{
    double Ax[] = {1};
    int Ai[] = {0};
    int Ap[] = {0, 1, 1, 1, 1};
    double lhs[] = {1, -INF, -INF, -INF};
    double rhs[] = {INF, -3, -1, -1};
    double lbs[] = {-INF};
    double ubs[] = {INF};
    double c[] = {1};
    ConeBlock cones[] = {{PSLP_CONE_EXP, 1, 3}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 4, 1, 1, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("reduced", run_presolver(presolver) == REDUCED);
    const ProblemData *reduced = presolver->reduced_prob;
    mu_assert("block kept",
              reduced->m == 3 && reduced->n == 0 && reduced->n_cones == 1 &&
                  reduced->cones[0].type == PSLP_CONE_EXP &&
                  reduced->cones[0].start == 0 && reduced->cones[0].dim == 3);
    mu_assert("constant rows", reduced->nnz == 0 && reduced->rhs[0] == -3 &&
                                   reduced->rhs[1] == -1 && reduced->rhs[2] == -1);

    double x_red[1];
    double y_red[] = {0, 0, 0};
    double z_red[1];
    mu_assert("reduced KKT point", is_kkt_point(reduced, x_red, y_red, z_red, 1e-9));
    postsolve(presolver, x_red, y_red, z_red);
    const double *y = presolver->sol->y;
    mu_assert("original KKT point",
              is_kkt_point(&prob, presolver->sol->x, y, presolver->sol->z, 1e-9));
    mu_assert("block multipliers", y[1] == 0 && y[2] == 0 && y[3] == 0);
    free_presolver(presolver);

    rhs[1] = -1;
    presolver = new_presolver_from_problem(&prob, stgs);
    mu_assert("init infeasible block", presolver != NULL);
    mu_assert("not detected", run_presolver(presolver) == REDUCED);
    mu_assert("infeasible block kept", presolver->stats->n_cones_reduced == 1);
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
    mu_run_test(test_cones_rotated_to_P, counter_cones);
    mu_run_test(test_cones_rotated_to_P_signs, counter_cones);
    mu_run_test(test_cones_zero_tails, counter_cones);
    mu_run_test(test_cones_block_to_bounds, counter_cones);
    mu_run_test(test_cones_constant_block, counter_cones);
    mu_run_test(test_cones_chain, counter_cones);
    mu_run_test(test_cones_hole, counter_cones);
    mu_run_test(test_cones_exp_flow, counter_cones);
    mu_run_test(test_cones_exp_constant_block, counter_cones);
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
