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
    PresolvedProblem prob =                                                         \
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
    PresolvedProblem prob =
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
    PresolvedProblem prob =
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
    PresolvedProblem prob =
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
