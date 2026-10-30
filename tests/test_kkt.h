#ifndef TEST_KKT_H
#define TEST_KKT_H

#include "Numerics.h"
#include "glbopts.h"
#include "kkt.h"
#include "minunit.h"
#include <stdio.h>

static int counter_kkt = 0;

/*  min. [-3  -1  4] x + 2.5
    s.t.        [1  1  0] x <= 4
        -4 <=   [0  1  1] x
        0 <= x1 <= 3,  0 <= x2,  -5 <= x3 <= 2

    Optimum x = (3, 1, -5) with row 0 active at its rhs, row 1 at its lhs,
    x1 at its upper bound, x2 interior and x3 at its lower bound:
    y = (-2, 1), z = (-1, 0, 3). */
#define KKT_TEST_LP                                                                 \
    double Ax[] = {1, 1, 1, 1};                                                     \
    int Ai[] = {0, 1, 1, 2};                                                        \
    int Ap[] = {0, 2, 4};                                                           \
    int nnz = 4;                                                                    \
    int n_rows = 2;                                                                 \
    int n_cols = 3;                                                                 \
    double lhs[] = {-INF, -4};                                                      \
    double rhs[] = {4, INF};                                                        \
    double lbs[] = {0, 0, -5};                                                      \
    double ubs[] = {3, INF, 2};                                                     \
    double c[] = {-3, -1, 4};                                                       \
    ProblemData prob =                                                              \
        problem_from_csr(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c);   \
    double x[] = {3, 1, -5};                                                        \
    double y[] = {-2, 1};                                                           \
    double z[] = {-1, 0, 3};

static char *test_kkt_optimal_point()
{
    KKT_TEST_LP
    prob.obj_offset = 2.5;

    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("primal objective", ABS(r.primal_obj - (-27.5)) <= 1e-12);
    mu_assert("dual residual", r.dual_res_abs <= 1e-12 && r.dual_res_rel <= 1e-12);
    mu_assert("primal residual",
              r.primal_res_abs <= 1e-12 && r.primal_res_rel <= 1e-12);
    mu_assert("gap", r.gap_abs <= 1e-12 && r.gap_rel <= 1e-12);
    mu_assert("comp slack", r.y_comp_slack <= 1e-12 && r.z_comp_slack <= 1e-12);
    mu_assert("bound violation", r.bound_viol <= 1e-12);
    mu_assert("norms", ABS(r.x_norm * r.x_norm - 35) <= 1e-12 &&
                           ABS(r.y_norm * r.y_norm - 5) <= 1e-12 &&
                           ABS(r.z_norm * r.z_norm - 10) <= 1e-12);
    mu_assert("is_kkt_point", is_kkt_point(&prob, x, y, z, 1e-9));
    return 0;
}

static char *test_kkt_bound_violation()
{
    KKT_TEST_LP
    x[0] = 3.5;
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("bound violation not detected", ABS(r.bound_viol - 0.5) <= 1e-12);
    mu_assert("is_kkt_point accepted a bound violation",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

static char *test_kkt_primal_residual()
{
    KKT_TEST_LP
    x[1] = 2; // row 0 activity becomes 5 > 4, x stays within its bounds
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("primal residual not detected", ABS(r.primal_res_abs - 1) <= 1e-12);
    mu_assert("bound violation reported", r.bound_viol <= 1e-12);
    mu_assert("is_kkt_point accepted an infeasible x",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

static char *test_kkt_dual_residual()
{
    KKT_TEST_LP
    z[1] = 0.5; // x2 is interior, so this also breaks complementary slackness
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("dual residual not detected", ABS(r.dual_res_abs - 0.5) <= 1e-12);
    mu_assert("z comp slack not detected", ABS(r.z_comp_slack - 0.5) <= 1e-12);
    mu_assert("primal residual reported", r.primal_res_abs <= 1e-12);
    mu_assert("is_kkt_point accepted a wrong z",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

static char *test_kkt_y_comp_slack()
{
    KKT_TEST_LP
    y[1] = -1; // wrong sign on the lhs-active row
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("y comp slack not detected", ABS(r.y_comp_slack - 1) <= 1e-12);
    mu_assert("is_kkt_point accepted a wrong-signed y",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

static char *test_kkt_gap()
{
    KKT_TEST_LP
    y[0] = -4; // scaled multiplier: c'x = -30, p(-y) + p(-z) = 42
    y[1] = 2;
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("gap not detected", ABS(r.gap_abs - 12) <= 1e-12);
    mu_assert("comp slack reported", r.y_comp_slack <= 1e-12);
    mu_assert("is_kkt_point accepted a wrong gap",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

/* A problem without rows: min x1 - x2 over the unit box. */
static char *test_kkt_no_rows()
{
    int Ap[] = {0};
    double lbs[] = {0, 0};
    double ubs[] = {1, 1};
    double c[] = {1, -1};
    ProblemData prob =
        problem_from_csr(NULL, NULL, Ap, 0, 2, 0, NULL, NULL, lbs, ubs, c);
    double x[] = {0, 1};
    double z[] = {1, -1};
    double y[1] = {0};
    mu_assert("no rows optimal point", is_kkt_point(&prob, x, y, z, 1e-9));
    x[0] = 1;
    mu_assert("no rows wrong point", !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

/*  Infeasible: x1 + x2 <= 5, 2x1 + 2x2 <= 4, x1 + 2x2 >= 3.5, 2x1 + x2 >= 3.5,
    x >= 0 (the problem of test_infeasible_rhs_change_primal_ray_postsolve). */
static char *test_kkt_primal_ray_certificate()
{
    double Ax[] = {1, 1, 2, 2, 1, 2, 2, 1};
    int Ai[] = {0, 1, 0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6, 8};
    int nnz = 8;
    int n_rows = 4;
    int n_cols = 2;
    double lhs[] = {-INF, -INF, 3.5, 3.5};
    double rhs[] = {5, 4, INF, INF};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};
    ProblemData prob =
        problem_from_csr(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c);

    // A'y = 0 and support = 4 * 1.5 - 3.5 - 3.5 = -1
    double y[] = {0, 1.5, -1, -1};
    mu_assert("valid certificate rejected",
              is_primal_ray_certificate(&prob, y, FEAS_TOL));

    double y_neg[] = {0, -1.5, 1, 1};
    mu_assert("negated certificate accepted",
              !is_primal_ray_certificate(&prob, y_neg, FEAS_TOL));

    // z = (-2, -1) pairs with the finite lower bounds, but support = 2.5
    double y_pos[] = {0, 1.5, -1, 0};
    mu_assert("positive support accepted",
              !is_primal_ray_certificate(&prob, y_pos, FEAS_TOL));

    // positive multiplier on a row without rhs
    double y_inf[] = {0, 1.5, -1, 1};
    mu_assert("multiplier on infinite side accepted",
              !is_primal_ray_certificate(&prob, y_inf, FEAS_TOL));
    return 0;
}

/*  Unbounded: x1 - x2 <= 1, -x1 + x2 <= 1, 2x1 - x3 <= 0, x1, x2 >= 0,
    min -x1 - x2 (the problem of test_unbounded_fix_col_inf_dual_ray_postsolve). */
static char *test_kkt_dual_ray_certificate()
{
    double Ax[] = {1, -1, -1, 1, 2, -1};
    int Ai[] = {0, 1, 0, 1, 0, 2};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 3;
    double lhs[] = {-INF, -INF, -INF};
    double rhs[] = {1, 1, 0};
    double lbs[] = {0, 0, -INF};
    double ubs[] = {INF, INF, INF};
    double c[] = {-1, -1, 0};
    ProblemData prob =
        problem_from_csr(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c);

    double d[] = {1, 1, 2};
    mu_assert("valid certificate rejected",
              is_dual_ray_certificate(&prob, d, FEAS_TOL));

    double d_neg[] = {-1, -1, -2};
    mu_assert("negated certificate accepted",
              !is_dual_ray_certificate(&prob, d_neg, FEAS_TOL));

    double d_row[] = {1, 1, 1}; // (A d)_2 = 1 > 0 on an rhs-only row
    mu_assert("non-recession direction accepted",
              !is_dual_ray_certificate(&prob, d_row, FEAS_TOL));

    double d_zero[] = {0, 0, 1}; // recession direction with c'd = 0
    mu_assert("direction with zero objective accepted",
              !is_dual_ray_certificate(&prob, d_zero, FEAS_TOL));
    return 0;
}

/* FEAS_TOL scaled by the point, max(1, ||s||_inf). */
static double scaled_tol(const double *s, int dim)
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
    return cone_violation(type, s, dim) <= scaled_tol(s, dim);
}

static bool in_dual_cone(ConeType type, const double *y, int dim)
{
    return cone_dual_violation(type, y, dim) <= scaled_tol(y, dim);
}

static char *test_kkt_soc_violation()
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
    double two_off[] = {3, 3, 4};

    mu_assert("inside", in_cone(PSLP_CONE_SOC, inside, 3));
    mu_assert("on", in_cone(PSLP_CONE_SOC, on, 3));
    mu_assert("outside", !in_cone(PSLP_CONE_SOC, outside, 3));
    mu_assert("negative head", !in_cone(PSLP_CONE_SOC, negative_head, 3));
    mu_assert("scaled on", in_cone(PSLP_CONE_SOC, scaled_on, 3));
    mu_assert("scaled outside", !in_cone(PSLP_CONE_SOC, scaled_outside, 3));
    mu_assert("dim 1 zero", in_cone(PSLP_CONE_SOC, &zero, 1));
    mu_assert("dim 1 small negative", in_cone(PSLP_CONE_SOC, &small_negative, 1));
    mu_assert("dim 1 negative", !in_cone(PSLP_CONE_SOC, &negative, 1));
    mu_assert("violation value",
              ABS(cone_violation(PSLP_CONE_SOC, two_off, 3) - 2) <= 1e-12);
    mu_assert("unknown type has no formula",
              isinf(cone_violation((ConeType) 99, inside, 3)));

    mu_assert("dual inside", in_dual_cone(PSLP_CONE_SOC, inside, 3));
    mu_assert("dual on", in_dual_cone(PSLP_CONE_SOC, on, 3));
    mu_assert("dual outside", !in_dual_cone(PSLP_CONE_SOC, outside, 3));
    mu_assert("dual negative head", !in_dual_cone(PSLP_CONE_SOC, negative_head, 3));
    return 0;
}

static char *test_kkt_rotated_violation()
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
    double three_off[] = {1, 2, 2, 3};

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
    mu_assert("violation value",
              ABS(cone_violation(PSLP_CONE_ROTATED, three_off, 4) - 3) <= 1e-12);

    mu_assert("dual inside", in_dual_cone(PSLP_CONE_ROTATED, inside, 4));
    mu_assert("dual on", in_dual_cone(PSLP_CONE_ROTATED, on, 4));
    mu_assert("dual outside", !in_dual_cone(PSLP_CONE_ROTATED, outside, 4));
    mu_assert("dual negative head",
              !in_dual_cone(PSLP_CONE_ROTATED, negative_first_head, 4));
    return 0;
}

/*  EXP is (t, s, r) with t >= s exp(r / s); its dual (u, v, w) has
    u >= -w exp(v / w - 1) for w < 0. The point (e, 1, 1) is on the boundary
    and (1, 0, -e) is its supporting dual, so the two are complementary. */
static char *test_kkt_exp_violation()
{
    double inside[] = {3, 1, 1};
    double on[] = {exp(1.0), 1, 1};
    double outside[] = {2.7, 1, 1};
    double negative_t[] = {-1, 1, -5};
    double closure_inside[] = {1, 0, -1};
    double closure_positive_r[] = {1, 0, 1};
    double closure_negative_t[] = {-1, 0, -1};
    double negative_s[] = {1, -1e-3, -1};
    double zero[] = {0, 0, 0};
    double scaled_inside[] = {3e6, 1e6, 1e6};
    double scaled_outside[] = {2.7e6, 1e6, 1e6};
    double one_off[] = {1, 1, 1};

    mu_assert("inside", in_cone(PSLP_CONE_EXP, inside, 3));
    mu_assert("on", in_cone(PSLP_CONE_EXP, on, 3));
    mu_assert("outside", !in_cone(PSLP_CONE_EXP, outside, 3));
    mu_assert("negative t", !in_cone(PSLP_CONE_EXP, negative_t, 3));
    mu_assert("closure inside", in_cone(PSLP_CONE_EXP, closure_inside, 3));
    mu_assert("closure positive r", !in_cone(PSLP_CONE_EXP, closure_positive_r, 3));
    mu_assert("closure negative t", !in_cone(PSLP_CONE_EXP, closure_negative_t, 3));
    mu_assert("negative s", !in_cone(PSLP_CONE_EXP, negative_s, 3));
    mu_assert("zero", in_cone(PSLP_CONE_EXP, zero, 3));
    mu_assert("scaled inside", in_cone(PSLP_CONE_EXP, scaled_inside, 3));
    mu_assert("scaled outside", !in_cone(PSLP_CONE_EXP, scaled_outside, 3));
    mu_assert("violation value", ABS(cone_violation(PSLP_CONE_EXP, one_off, 3) -
                                     (exp(1.0) - 1)) <= 1e-12);

    double dual_inside[] = {2, 0, -1};
    double dual_on[] = {exp(-1.0), 0, -1};
    double dual_outside[] = {0.3, 0, -1};
    double dual_of_on[] = {1, 0, -exp(1.0)};
    double dual_closure_inside[] = {1, 1, 0};
    double dual_closure_negative_u[] = {-1, 1, 0};
    double dual_closure_negative_v[] = {1, -1, 0};
    double dual_positive_w[] = {1, 1, 1e-3};
    double dual_one_off[] = {0, 0, -exp(1.0)};
    double dual_negative_v[] = {1, -1, -1};

    mu_assert("dual inside", in_dual_cone(PSLP_CONE_EXP, dual_inside, 3));
    mu_assert("dual on", in_dual_cone(PSLP_CONE_EXP, dual_on, 3));
    mu_assert("dual outside", !in_dual_cone(PSLP_CONE_EXP, dual_outside, 3));
    mu_assert("dual of the boundary point",
              in_dual_cone(PSLP_CONE_EXP, dual_of_on, 3));
    mu_assert("dual closure inside",
              in_dual_cone(PSLP_CONE_EXP, dual_closure_inside, 3));
    mu_assert("dual closure negative u",
              !in_dual_cone(PSLP_CONE_EXP, dual_closure_negative_u, 3));
    mu_assert("dual closure negative v",
              !in_dual_cone(PSLP_CONE_EXP, dual_closure_negative_v, 3));
    mu_assert("dual positive w", !in_dual_cone(PSLP_CONE_EXP, dual_positive_w, 3));
    mu_assert("dual violation value",
              ABS(cone_dual_violation(PSLP_CONE_EXP, dual_one_off, 3) - 1) <= 1e-12);
    mu_assert("dual negative v", in_dual_cone(PSLP_CONE_EXP, dual_negative_v, 3));
    mu_assert("the cones differ", !in_dual_cone(PSLP_CONE_EXP, on, 3) &&
                                      !in_cone(PSLP_CONE_EXP, dual_negative_v, 3));
    mu_assert("complementary pair",
              ABS(on[0] * dual_of_on[0] + on[1] * dual_of_on[1] +
                  on[2] * dual_of_on[2]) <= 1e-12);
    return 0;
}

/*  Columns x0 .. x4, all free except x4 with lb = ub = 1. One linear row
    x0 - x3 = 1 and the SOC rows s = (x0, x1 + x3, x2 + x4 - 1). The lhs of the
    cone rows is -INF to show that it is ignored.

    Primal x = (sqrt2, 2 - sqrt2, 1, sqrt2 - 1, 1), so s = (sqrt2, 1, 1) on the
    boundary. Dual y = (1; sqrt2, -1, -1) with y_B = (s0, -s_u) in K*,
    z = (0, 0, 0, 0, 1), so c = A'y + z = (1 + sqrt2, -1, -1, -2, 0) and the
    objective is 1. */
#define KKT_TEST_SOC                                                                \
    double Ax[] = {1, -1, 1, 1, 1, 1, 1};                                           \
    int Ai[] = {0, 3, 0, 1, 3, 2, 4};                                               \
    int Ap[] = {0, 2, 3, 5, 7};                                                     \
    double lhs[] = {1, -INF, -INF, -INF};                                           \
    double rhs[] = {1, 0, 0, 1};                                                    \
    double lbs[] = {-INF, -INF, -INF, -INF, 1};                                     \
    double ubs[] = {INF, INF, INF, INF, 1};                                         \
    double c[] = {1 + sqrt(2.0), -1, -1, -2, 0};                                    \
    ConeBlock cones[] = {{PSLP_CONE_SOC, 1, 3}};                                    \
    ProblemData prob =                                                              \
        problem_from_csr(Ax, Ai, Ap, 4, 5, 7, lhs, rhs, lbs, ubs, c);               \
    prob.cones = cones;                                                             \
    prob.n_cones = 1;                                                               \
    double x[] = {sqrt(2.0), 2 - sqrt(2.0), 1, sqrt(2.0) - 1, 1};                   \
    double y[] = {1, sqrt(2.0), -1, -1};                                            \
    double z[] = {0, 0, 0, 0, 1};

static char *test_kkt_soc_block()
{
    KKT_TEST_SOC
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("primal objective", ABS(r.primal_obj - 1) <= 1e-12);
    mu_assert("dual residual", r.dual_res_abs <= 1e-12);
    mu_assert("primal residual", r.primal_res_abs <= 1e-12);
    mu_assert("gap", r.gap_abs <= 1e-12);
    mu_assert("comp slack", r.y_comp_slack <= 1e-12 && r.z_comp_slack <= 1e-12);
    mu_assert("cone quantities", r.cone_primal_viol <= 1e-12 &&
                                     r.cone_dual_viol <= 1e-12 &&
                                     r.cone_comp <= 1e-12);
    mu_assert("is_kkt_point", is_kkt_point(&prob, x, y, z, 1e-9));
    return 0;
}

static char *test_kkt_cone_primal_violation()
{
    KKT_TEST_SOC
    x[0] = 1; // s = (1, 1, 1), violation sqrt2 - 1 over 1 + ||s||_inf
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("cone primal violation not detected",
              ABS(r.cone_primal_viol - (sqrt(2.0) - 1) / 2) <= 1e-12);
    mu_assert("is_kkt_point accepted a slack outside the cone",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

static char *test_kkt_cone_dual_violation()
{
    KKT_TEST_SOC
    y[1] = 1; // y_B = (1, -1, -1) is outside K*, violation sqrt2 - 1 over 2
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("cone dual violation not detected",
              ABS(r.cone_dual_viol - (sqrt(2.0) - 1) / 2) <= 1e-12);
    mu_assert("is_kkt_point accepted a multiplier outside the dual cone",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

static char *test_kkt_cone_complementarity()
{
    KKT_TEST_SOC
    y[1] = 2 * sqrt(2.0); // y_B stays in K*, s'y = 2, ||s|| = 2, ||y|| = sqrt10
    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("dual violation reported", r.cone_dual_viol <= 1e-12);
    mu_assert("cone complementarity not detected",
              ABS(r.cone_comp - 2 / (1 + 2 * sqrt(10.0))) <= 1e-12);
    mu_assert("is_kkt_point accepted a non-complementary pair",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

/*  Four free columns, A = I, one ROTATED block of dimension 4 with r = 0.
    x = s = (1, 2, 2, 0) on the boundary, y = (s1, s0, -s_u) = (2, 1, -2, 0)
    in K* with s'y = 0, z = 0, c = y, objective 0. */
static char *test_kkt_rotated_block()
{
    double Ax[] = {1, 1, 1, 1};
    int Ai[] = {0, 1, 2, 3};
    int Ap[] = {0, 1, 2, 3, 4};
    double lhs[] = {-INF, -INF, -INF, -INF};
    double rhs[] = {0, 0, 0, 0};
    double lbs[] = {-INF, -INF, -INF, -INF};
    double ubs[] = {INF, INF, INF, INF};
    double c[] = {2, 1, -2, 0};
    ConeBlock cones[] = {{PSLP_CONE_ROTATED, 0, 4}};
    ProblemData prob = problem_from_csr(Ax, Ai, Ap, 4, 4, 4, lhs, rhs, lbs, ubs, c);
    prob.cones = cones;
    prob.n_cones = 1;
    double x[] = {1, 2, 2, 0};
    double y[] = {2, 1, -2, 0};
    double z[] = {0, 0, 0, 0};

    KKTResiduals r;
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("primal objective", ABS(r.primal_obj) <= 1e-12);
    mu_assert("dual residual", r.dual_res_abs <= 1e-12);
    mu_assert("gap", r.gap_abs <= 1e-12);
    mu_assert("cone quantities", r.cone_primal_viol <= 1e-12 &&
                                     r.cone_dual_viol <= 1e-12 &&
                                     r.cone_comp <= 1e-12);
    mu_assert("is_kkt_point", is_kkt_point(&prob, x, y, z, 1e-9));

    x[3] = 1; // s = (1, 2, 2, 1): defect sqrt(5 - 4) = 1 over 1 + ||s||_inf = 3
    kkt_residuals(&prob, x, y, z, &r);
    mu_assert("cone primal violation not detected",
              ABS(r.cone_primal_viol - 1.0 / 3) <= 1e-12);
    mu_assert("complementarity kept", r.cone_comp <= 1e-12);
    mu_assert("is_kkt_point accepted a slack outside the cone",
              !is_kkt_point(&prob, x, y, z, 1e-6));
    return 0;
}

static const char *all_tests_kkt()
{
    mu_run_test(test_kkt_optimal_point, counter_kkt);
    mu_run_test(test_kkt_bound_violation, counter_kkt);
    mu_run_test(test_kkt_primal_residual, counter_kkt);
    mu_run_test(test_kkt_dual_residual, counter_kkt);
    mu_run_test(test_kkt_y_comp_slack, counter_kkt);
    mu_run_test(test_kkt_gap, counter_kkt);
    mu_run_test(test_kkt_no_rows, counter_kkt);
    mu_run_test(test_kkt_primal_ray_certificate, counter_kkt);
    mu_run_test(test_kkt_dual_ray_certificate, counter_kkt);
    mu_run_test(test_kkt_soc_violation, counter_kkt);
    mu_run_test(test_kkt_rotated_violation, counter_kkt);
    mu_run_test(test_kkt_exp_violation, counter_kkt);
    mu_run_test(test_kkt_soc_block, counter_kkt);
    mu_run_test(test_kkt_cone_primal_violation, counter_kkt);
    mu_run_test(test_kkt_cone_dual_violation, counter_kkt);
    mu_run_test(test_kkt_cone_complementarity, counter_kkt);
    mu_run_test(test_kkt_rotated_block, counter_kkt);
    return 0;
}

int test_kkt()
{
    const char *result = all_tests_kkt();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("kkt: TEST FAILED!\n");
    }
    else
    {
        printf("kkt: ALL TESTS PASSED\n");
    }
    printf("kkt: Tests run: %d\n", counter_kkt);
    return result == 0;
}

#endif
