#ifndef TEST_RAY_POSTSOLVE_H
#define TEST_RAY_POSTSOLVE_H

#include "Numerics.h"
#include "PSLP_API.h"
#include "Postsolver.h"
#include "minunit.h"

static int counter_ray_postsolve = 0;

static bool is_ray_correct(const double *ray, const double *correct_ray, int len)
{
    bool correct = true;
    for (int i = 0; i < len; ++i)
    {
        if (!IS_EQUAL_FEAS_TOL(ray[i], correct_ray[i]))
        {
            correct = false;
        }
    }

    return correct;
}

static bool is_primal_ray_z_correct(const double *Ax, const int *Ai, const int *Ap,
                                    const double *y, const double *z, int n_rows,
                                    int n_cols)
{
    for (int j = 0; j < n_cols; ++j)
    {
        double z_j = 0.0;
        for (int i = 0; i < n_rows; ++i)
        {
            for (int p = Ap[i]; p < Ap[i + 1]; ++p)
            {
                if (Ai[p] == j)
                {
                    z_j -= Ax[p] * y[i];
                }
            }
        }

        if (!IS_EQUAL_FEAS_TOL(z_j, z[j]))
        {
            return false;
        }
    }

    return true;
}

static bool is_primal_ray_valid_certificate(const double *Ax, const int *Ai,
                                            const int *Ap, const double *lhs,
                                            const double *rhs, const double *lbs,
                                            const double *ubs, const double *y,
                                            int n_rows, int n_cols)
{
    double support = 0.0;
    for (int i = 0; i < n_rows; ++i)
    {
        if (y[i] > FEAS_TOL)
        {
            if (IS_POS_INF(rhs[i]))
            {
                return false;
            }
            support += rhs[i] * y[i];
        }
        else if (y[i] < -FEAS_TOL)
        {
            if (IS_NEG_INF(lhs[i]))
            {
                return false;
            }
            support += lhs[i] * y[i];
        }
    }

    for (int j = 0; j < n_cols; ++j)
    {
        double z_j = 0.0;
        for (int i = 0; i < n_rows; ++i)
        {
            for (int p = Ap[i]; p < Ap[i + 1]; ++p)
            {
                if (Ai[p] == j)
                {
                    z_j -= Ax[p] * y[i];
                }
            }
        }

        if (z_j > FEAS_TOL)
        {
            if (IS_POS_INF(ubs[j]))
            {
                return false;
            }
            support += ubs[j] * z_j;
        }
        else if (z_j < -FEAS_TOL)
        {
            if (IS_NEG_INF(lbs[j]))
            {
                return false;
            }
            support += lbs[j] * z_j;
        }
    }

    return support < -FEAS_TOL;
}

static bool is_dual_ray_valid_certificate(const double *Ax, const int *Ai,
                                          const int *Ap, const double *lhs,
                                          const double *rhs, const double *lbs,
                                          const double *ubs, const double *c,
                                          const double *x, int n_rows, int n_cols)
{
    for (int i = 0; i < n_rows; ++i)
    {
        double ax_i = 0.0;
        for (int p = Ap[i]; p < Ap[i + 1]; ++p)
        {
            ax_i += Ax[p] * x[Ai[p]];
        }

        if (!IS_NEG_INF(lhs[i]) && !IS_POS_INF(rhs[i]))
        {
            if (!IS_ZERO_FEAS_TOL(ax_i))
            {
                return false;
            }
        }
        else if (!IS_POS_INF(rhs[i]))
        {
            if (ax_i > FEAS_TOL)
            {
                return false;
            }
        }
        else if (!IS_NEG_INF(lhs[i]))
        {
            if (ax_i < -FEAS_TOL)
            {
                return false;
            }
        }
    }

    for (int j = 0; j < n_cols; ++j)
    {
        if (!IS_NEG_INF(lbs[j]) && !IS_POS_INF(ubs[j]))
        {
            if (!IS_ZERO_FEAS_TOL(x[j]))
            {
                return false;
            }
        }
        else if (!IS_POS_INF(ubs[j]))
        {
            if (x[j] > FEAS_TOL)
            {
                return false;
            }
        }
        else if (!IS_NEG_INF(lbs[j]))
        {
            if (x[j] < -FEAS_TOL)
            {
                return false;
            }
        }
    }

    double cx = 0.0;
    for (int j = 0; j < n_cols; ++j)
    {
        cx += c[j] * x[j];
    }

    return cx < -FEAS_TOL;
}

static char *test_0_primal_ray_postsolve()
{
    double Ax[] = {1.0, 2.0, -1.0, 1.0};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 2;

    double lhs[] = {-1.0, -2.0};
    double rhs[] = {3.0, 4.0};
    double lbs[] = {-10.0, -10.0};
    double ubs[] = {10.0, 10.0};
    double c[] = {1.0, -1.0};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("primal ray postsolve presolve status", status == UNCHANGED);

    double y[] = {1.0, -2.0};
    double y_orig[] = {0.0, 0.0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("primal ray postsolve y0", IS_EQUAL_FEAS_TOL(y_orig[0], y[0]));
    mu_assert("primal ray postsolve y1", IS_EQUAL_FEAS_TOL(y_orig[1], y[1]));
    mu_assert("primal ray postsolve z",
              is_primal_ray_z_correct(Ax, Ai, Ap, y_orig, presolver->sol->z, n_rows,
                                      n_cols));

    free_settings(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_0_dual_ray_postsolve()
{
    double Ax[] = {1.0, 2.0, -1.0, 1.0};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int nnz = 4;
    int n_rows = 2;
    int n_cols = 2;

    double lhs[] = {-1.0, -2.0};
    double rhs[] = {3.0, 4.0};
    double lbs[] = {-10.0, -10.0};
    double ubs[] = {10.0, 10.0};
    double c[] = {1.0, -1.0};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("dual ray postsolve presolve status", status == UNCHANGED);

    double x[] = {0.0, 0.0};
    double x_orig[] = {0.0, 0.0};
    postsolve_dual_infeas_ray(presolver, x, x_orig);

    mu_assert("dual ray postsolve x0", IS_EQUAL_FEAS_TOL(x_orig[0], x[0]));
    mu_assert("dual ray postsolve x1", IS_EQUAL_FEAS_TOL(x_orig[1], x[1]));

    free_settings(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_singleton_eq_primal_ray_postsolve()
{
    double Ax[] = {0.1, 0.2, -0.3, -0.1, 0.4,  -0.4, 0.1,  0.3, 0,
                   1,   0.3, 0.5,  -0.6, -0.3, 0.2,  -0.2, 0.4, 0.1};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 0, 1, 2, 3, 4, 5, 6, 7};
    int Ap[] = {0, 8, 10, 18};
    int nnz = 18;
    int n_rows = 3;
    int n_cols = 8;

    double lhs[] = {0.5, 3, -2};
    double rhs[] = {0.5, 3, 2};
    double lbs[] = {-10, -10, -10, -10, -10, -10, -10, -10};
    double ubs[] = {10, 10, 10, 10, 10, 10, 10, 10};
    double c[] = {-1, 1, -1, 1, -1, 2, -0.2, 0.7};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray
    double y[] = {-2.57142857, 0.14285714};
    double y_orig[] = {0.0, 0.0, 0.0};
    double correct_y[] = {-2.57142857, -0.68571429, 0.14285714};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("singleton eq primal ray postsolve error",
              is_ray_correct(y_orig, correct_y, n_rows));
    mu_assert("singleton eq primal ray z",
              is_primal_ray_z_correct(Ax, Ai, Ap, y_orig, presolver->sol->z, n_rows,
                                      n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_dton_primal_ray_postsolve()
{
    double Ax[] = {0.1, 0.2, -0.3, -0.1, 0.4,  -0.4, 0.1,  0.3, 1,  1,
                   1,   0.3, 0.5,  -0.6, -0.3, 0.2,  -0.2, 0.4, 0.1};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 0, 1, 2, 3, 4, 5, 6, 7};
    int Ap[] = {0, 8, 9, 11, 19};
    int nnz = 19;
    int n_rows = 4;
    int n_cols = 8;

    double lhs[] = {0.5, 3, 4, -2};
    double rhs[] = {0.5, 3, 4, 2};
    double lbs[] = {-10, -10, -10, -10, -10, -10, -10, -10};
    double ubs[] = {10, 10, 10, 10, 10, 10, 10, 10};
    double c[] = {-1, 1, -1, 1, -1, 2, -0.2, 0.7};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray
    double y[] = {-4.66666667, 0.66666667};
    double y_orig[] = {0.0, 0.0, 0.0, 0.0};
    double correct_y[] = {-4.66666667, -1., -0.26666667, 0.66666667};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("dton primal ray postsolve error",
              is_ray_correct(y_orig, correct_y, n_rows));
    mu_assert("dton primal ray z",
              is_primal_ray_z_correct(Ax, Ai, Ap, y_orig, presolver->sol->z, n_rows,
                                      n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_bound_change_primal_ray_postsolve()
{
    double Ax[] = {0.1, 0.2, -0.3, -0.1, 0.4,  -0.4, 0.1,  0.3, 1,  1,
                   1,   0.3, 0.5,  -0.6, -0.3, 0.2,  -0.2, 0.4, 0.1};
    int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 0, 1, 2, 3, 4, 5, 6, 7};
    int Ap[] = {0, 8, 9, 11, 19};
    int nnz = 19;
    int n_rows = 4;
    int n_cols = 8;

    double lhs[] = {0.5, 3, 4, -2};
    double rhs[] = {0.5, 3, 4, 2};
    double lbs[] = {-10, -10, -10, -10, -10, -10, -10, -10};
    double ubs[] = {10, 10, 10, 10, 4.0 - 8.0 / 5.0, 10, 10, 10};
    double c[] = {-1, 1, -1, 1, -1, 2, -0.2, 0.7};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray
    double y[] = {-2.0, 0.0};
    double y_orig[] = {0.0, 0.0, 0.0, 0.0};
    double correct_y[] = {-2., -0.6, -0.2, 0.0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("bound change primal ray postsolve error",
              is_ray_correct(y_orig, correct_y, n_rows));
    mu_assert("bound change primal ray z",
              is_primal_ray_z_correct(Ax, Ai, Ap, y_orig, presolver->sol->z, n_rows,
                                      n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_bound_change_dual_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(1, 2);

    save_retrieval_bound_change_no_row(info, 1, 3.0, 0.0, true);

    int rows[] = {0};
    double vals[] = {2.0};
    save_retrieval_fixed_col(info, 1, 3.0, 0.0, vals, rows, 1);

    int cols[] = {0, 1};
    double row_vals[] = {1.0, 2.0};
    save_retrieval_bound_change_the_row(info, 0, cols, row_vals, 2, 1);

    int col_map[] = {0, -1};
    int row_map[] = {0};
    postsolver_update(info, 1, 1, col_map, row_map);

    double x[] = {4.0};
    double sol_x[2];
    Solution sol = {sol_x, NULL, NULL, 2, 0};

    postsolver_run_dual_infeas_ray(info, &sol, x);

    mu_assert("bound change dual ray x0", IS_EQUAL_FEAS_TOL(sol.x[0], 4.0));
    mu_assert("bound change dual ray x1", IS_EQUAL_FEAS_TOL(sol.x[1], 0.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_fix_col_inf_primal_ray_postsolve()
{
    double Ax[] = {1.0, 1, 1, 2, 2, 1, 1, 1, 2, 1};
    int Ai[] = {0, 1, 2, 0, 1, 2, 0, 2, 0, 2};
    int Ap[] = {0, 3, 6, 8, 10};
    int nnz = 10;
    int n_rows = 4;
    int n_cols = 3;

    double lhs[] = {1, 3, 1, 2};
    double rhs[] = {INF, INF, INF, INF};
    double lbs[] = {-1.0, -INF, 0.5};
    double ubs[] = {INF, INF, INF};
    double c[] = {1.0, 0, 2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray
    double y[] = {0.0, 0.5};
    double y_orig[] = {0.0, 0.0, 0.0, 0.0};
    double correct_y[] = {0.0, 0.0, 0.0, 0.5};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("fix col inf primal ray postsolve error",
              is_ray_correct(y_orig, correct_y, n_rows));
    mu_assert("fix col inf primal ray z",
              is_primal_ray_z_correct(Ax, Ai, Ap, y_orig, presolver->sol->z, n_rows,
                                      n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_fix_col_inf_dual_ray_postsolve()
{
    double Ax[] = {1.0, 1, 1, 2, 2, 1, 1, 1, 2, 1};
    int Ai[] = {0, 1, 2, 0, 1, 2, 0, 2, 0, 2};
    int Ap[] = {0, 3, 6, 8, 10};
    int nnz = 10;
    int n_rows = 4;
    int n_cols = 3;

    double lhs[] = {1, 3, 1, 2};
    double rhs[] = {INF, INF, INF, INF};
    double lbs[] = {-1.0, -INF, 0.5};
    double ubs[] = {INF, INF, INF};
    double c[] = {1.0, 0, 2};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced dual infeasibility ray
    double x[] = {-5.0, 0.0};
    double x_orig[] = {0.0, 0.0, 0.0};
    double correct_x[] = {-5.0, 5.0, 0.0};
    postsolve_dual_infeas_ray(presolver, x, x_orig);

    mu_assert("fix col inf dual ray postsolve error",
              is_ray_correct(x_orig, correct_x, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_fix_col_primal_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(1, 2);

    int rows[] = {0};
    double vals[] = {2.0};
    save_retrieval_fixed_col(info, 0, 0.0, 0.0, vals, rows, 1);

    int col_map[] = {-1, 0};
    int row_map[] = {0};
    postsolver_update(info, 1, 1, col_map, row_map);

    double y[] = {3.0};
    double z[] = {-1.0};
    double sol_x[2];
    double sol_y[1];
    double sol_z[2];
    Solution sol = {sol_x, sol_y, sol_z, 2, 1};

    postsolver_run_primal_infeas_ray(info, &sol, y, z);

    mu_assert("fix col primal ray y0", IS_EQUAL_FEAS_TOL(sol.y[0], 3.0));
    mu_assert("fix col primal ray z0", IS_EQUAL_FEAS_TOL(sol.z[0], -6.0));
    mu_assert("fix col primal ray z1", IS_EQUAL_FEAS_TOL(sol.z[1], -1.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_fix_col_dual_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(1, 2);

    int rows[] = {0};
    double vals[] = {2.0};
    save_retrieval_fixed_col(info, 0, 0.0, 0.0, vals, rows, 1);

    int col_map[] = {-1, 0};
    int row_map[] = {0};
    postsolver_update(info, 1, 1, col_map, row_map);

    double x[] = {2.0};
    double sol_x[2];
    Solution sol = {sol_x, NULL, NULL, 2, 0};

    postsolver_run_dual_infeas_ray(info, &sol, x);

    mu_assert("fix col dual ray x0", IS_EQUAL_FEAS_TOL(sol.x[0], 0.0));
    mu_assert("fix col dual ray x1", IS_EQUAL_FEAS_TOL(sol.x[1], 2.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_sub_col_primal_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(1, 2);

    int cols[] = {0, 1};
    double vals[] = {2.0, 5.0};
    save_retrieval_sub_col(info, 0, cols, vals, 2, 0.0, 0, 0.0);

    int col_map[] = {-1, 0};
    int row_map[] = {0};
    postsolver_update(info, 1, 1, col_map, row_map);

    double y[] = {4.0};
    double z[] = {7.0};
    double sol_x[2];
    double sol_y[1];
    double sol_z[2];
    Solution sol = {sol_x, sol_y, sol_z, 2, 1};

    postsolver_run_primal_infeas_ray(info, &sol, y, z);

    mu_assert("sub col primal ray y0", IS_EQUAL_FEAS_TOL(sol.y[0], 4.0));
    mu_assert("sub col primal ray z0", IS_EQUAL_FEAS_TOL(sol.z[0], -8.0));
    mu_assert("sub col primal ray z1", IS_EQUAL_FEAS_TOL(sol.z[1], 7.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_sub_col_dual_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(1, 2);

    int cols[] = {0, 1};
    double vals[] = {2.0, 5.0};
    save_retrieval_sub_col(info, 0, cols, vals, 2, 0.0, 0, 0.0);

    int col_map[] = {-1, 0};
    int row_map[] = {0};
    postsolver_update(info, 1, 1, col_map, row_map);

    double x[] = {4.0};
    double sol_x[2];
    Solution sol = {sol_x, NULL, NULL, 2, 0};

    postsolver_run_dual_infeas_ray(info, &sol, x);

    mu_assert("sub col dual ray x0", IS_EQUAL_FEAS_TOL(sol.x[0], -10.0));
    mu_assert("sub col dual ray x1", IS_EQUAL_FEAS_TOL(sol.x[1], 4.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_parallel_col_primal_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(0, 2);

    save_retrieval_parallel_col(info, INF, -INF, -INF, INF, 2.0, 0, 1, C_TAG_NONE,
                                C_TAG_NONE);

    int col_map[] = {0, -1};
    postsolver_update(info, 1, 0, col_map, NULL);

    double z[] = {3.0};
    double sol_x[2];
    double sol_z[2];
    Solution sol = {sol_x, NULL, sol_z, 2, 0};

    postsolver_run_primal_infeas_ray(info, &sol, NULL, z);

    mu_assert("parallel col primal ray z0", IS_EQUAL_FEAS_TOL(sol.z[0], 3.0));
    mu_assert("parallel col primal ray z1", IS_EQUAL_FEAS_TOL(sol.z[1], 6.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_parallel_col_dual_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(0, 2);

    save_retrieval_parallel_col(info, 1.0, 0.0, -INF, INF, 2.0, 0, 1, C_TAG_UB_INF,
                                C_TAG_LB_INF | C_TAG_UB_INF);

    int col_map[] = {0, -1};
    postsolver_update(info, 1, 0, col_map, NULL);

    double x[] = {4.0};
    double sol_x[2];
    Solution sol = {sol_x, NULL, NULL, 2, 0};

    postsolver_run_dual_infeas_ray(info, &sol, x);

    mu_assert("parallel col dual ray x0", IS_EQUAL_FEAS_TOL(sol.x[0], 0.0));
    mu_assert("parallel col dual ray x1", IS_EQUAL_FEAS_TOL(sol.x[1], 2.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_lhs_change_primal_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(2, 2);

    int cols[] = {0, 1};
    double vals[] = {1.0, 2.0};
    save_retrieval_rhs_or_lhs_change(info, 0, vals, cols, 2, 1.0, 1, 3.0, true);

    int col_map[] = {0, 1};
    int row_map[] = {0, -1};
    postsolver_update(info, 2, 1, col_map, row_map);

    double y[] = {-2.0};
    double z[] = {0.0, 0.0};
    double sol_x[2];
    double sol_y[2];
    double sol_z[2];
    Solution sol = {sol_x, sol_y, sol_z, 2, 2};

    postsolver_run_primal_infeas_ray(info, &sol, y, z);

    mu_assert("lhs change primal ray y0", IS_EQUAL_FEAS_TOL(sol.y[0], 0.0));
    mu_assert("lhs change primal ray y1", IS_EQUAL_FEAS_TOL(sol.y[1], -6.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_lhs_change_primal_ray_postsolve_inactive()
{
    PostsolveInfo *info = postsolve_info_new(2, 2);

    int cols[] = {0, 1};
    double vals[] = {1.0, 2.0};
    save_retrieval_rhs_or_lhs_change(info, 0, vals, cols, 2, 1.0, 1, 3.0, true);

    int col_map[] = {0, 1};
    int row_map[] = {0, -1};
    postsolver_update(info, 2, 1, col_map, row_map);

    double y[] = {2.0};
    double z[] = {0.0, 0.0};
    double sol_x[2];
    double sol_y[2];
    double sol_z[2];
    Solution sol = {sol_x, sol_y, sol_z, 2, 2};

    postsolver_run_primal_infeas_ray(info, &sol, y, z);

    mu_assert("lhs change inactive primal ray y0", IS_EQUAL_FEAS_TOL(sol.y[0], 2.0));
    mu_assert("lhs change inactive primal ray y1", IS_EQUAL_FEAS_TOL(sol.y[1], 0.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_rhs_change_primal_ray_postsolve()
{
    PostsolveInfo *info = postsolve_info_new(2, 2);

    int cols[] = {0, 1};
    double vals[] = {1.0, 2.0};
    save_retrieval_rhs_or_lhs_change(info, 0, vals, cols, 2, 1.0, 1, 3.0, false);

    int col_map[] = {0, 1};
    int row_map[] = {0, -1};
    postsolver_update(info, 2, 1, col_map, row_map);

    double y[] = {2.0};
    double z[] = {0.0, 0.0};
    double sol_x[2];
    double sol_y[2];
    double sol_z[2];
    Solution sol = {sol_x, sol_y, sol_z, 2, 2};

    postsolver_run_primal_infeas_ray(info, &sol, y, z);

    mu_assert("rhs change primal ray y0", IS_EQUAL_FEAS_TOL(sol.y[0], 0.0));
    mu_assert("rhs change primal ray y1", IS_EQUAL_FEAS_TOL(sol.y[1], 6.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_rhs_change_primal_ray_postsolve_inactive()
{
    PostsolveInfo *info = postsolve_info_new(2, 2);

    int cols[] = {0, 1};
    double vals[] = {1.0, 2.0};
    save_retrieval_rhs_or_lhs_change(info, 0, vals, cols, 2, 1.0, 1, 3.0, false);

    int col_map[] = {0, 1};
    int row_map[] = {0, -1};
    postsolver_update(info, 2, 1, col_map, row_map);

    double y[] = {-2.0};
    double z[] = {0.0, 0.0};
    double sol_x[2];
    double sol_y[2];
    double sol_z[2];
    Solution sol = {sol_x, sol_y, sol_z, 2, 2};

    postsolver_run_primal_infeas_ray(info, &sol, y, z);

    mu_assert("rhs change inactive primal ray y0",
              IS_EQUAL_FEAS_TOL(sol.y[0], -2.0));
    mu_assert("rhs change inactive primal ray y1", IS_EQUAL_FEAS_TOL(sol.y[1], 0.0));

    postsolve_info_free(info);
    return 0;
}

static char *test_infeasible_rhs_change_primal_ray_postsolve()
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

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_cols = false;
    stgs->dual_fix = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray
    double y[] = {3.0, -1.0, -1.0};
    double y_orig[] = {0.0, 0.0, 0.0, 0.0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("infeasible rhs change primal ray certificate",
              is_primal_ray_valid_certificate(Ax, Ai, Ap, lhs, rhs, lbs, ubs, y_orig,
                                              n_rows, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_infeasible_lhs_change_primal_ray_postsolve()
{
    double Ax[] = {1, 1, 2, 2, 1, 2, 2, 1};
    int Ai[] = {0, 1, 0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6, 8};
    int nnz = 8;
    int n_rows = 4;
    int n_cols = 2;

    double lhs[] = {0.5, 4, -INF, -INF};
    double rhs[] = {INF, INF, 1.5, 1.5};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_cols = false;
    stgs->dual_fix = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray
    double y[] = {-3.0, 1.0, 1.0};
    double y_orig[] = {0.0, 0.0, 0.0, 0.0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("infeasible lhs change primal ray certificate",
              is_primal_ray_valid_certificate(Ax, Ai, Ap, lhs, rhs, lbs, ubs, y_orig,
                                              n_rows, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_infeasible_both_sides_changed_primal_ray_postsolve()
{
    // row0 survives the parallel-row merge and BOTH of its sides are
    // tightened from row1, which is then deleted. A reduced ray that uses
    // the merged lhs must end up with its multiplier on row1.
    double Ax[] = {1, 1, 1, 1, 1, 2};
    int Ai[] = {0, 1, 0, 1, 0, 1};
    int Ap[] = {0, 2, 4, 6};
    int nnz = 6;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {0.0, 3.5, -INF};
    double rhs[] = {10.0, 5.0, 1.0};
    double lbs[] = {0, 0};
    double ubs[] = {INF, INF};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_cols = false;
    stgs->dual_fix = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray using the lhs of the merged row
    double y[] = {-1.0, 1.0};
    double y_orig[] = {0.0, 0.0, 0.0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("both sides changed primal ray certificate",
              is_primal_ray_valid_certificate(Ax, Ai, Ap, lhs, rhs, lbs, ubs, y_orig,
                                              n_rows, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_infeasible_bound_change_primal_ray_postsolve()
{
    double Ax[] = {1, -2, 1, 1, 1};
    int Ai[] = {0, 1, 0, 1, 0};
    int Ap[] = {0, 2, 4, 5};
    int nnz = 5;
    int n_rows = 3;
    int n_cols = 2;

    double lhs[] = {-INF, -INF, 2};
    double rhs[] = {0, 2, INF};
    double lbs[] = {0, 0};
    double ubs[] = {10, 10};
    double c[] = {0, 0};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->dual_fix = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray
    double y[] = {0.5, 1.0};
    double y_orig[] = {0.0, 0.0, 0.0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("infeasible bound change primal ray certificate",
              is_primal_ray_valid_certificate(Ax, Ai, Ap, lhs, rhs, lbs, ubs, y_orig,
                                              n_rows, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_infeasible_parallel_col_primal_ray_postsolve()
{
    double Ax[] = {1, 2, 1, 1, 2, -1, 1, 2};
    int Ai[] = {0, 1, 2, 0, 1, 2, 0, 1};
    int Ap[] = {0, 3, 6, 8};
    int nnz = 8;
    int n_rows = 3;
    int n_cols = 3;

    double lhs[] = {-INF, -INF, 2};
    double rhs[] = {1, 1, INF};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {0, 0, 0};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->dual_fix = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray
    double y[] = {1.0, 0.0};
    double y_orig[] = {0.0, 0.0, 0.0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("infeasible parallel col primal ray certificate",
              is_primal_ray_valid_certificate(Ax, Ai, Ap, lhs, rhs, lbs, ubs, y_orig,
                                              n_rows, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_infeasible_fix_col_inf_primal_ray_postsolve()
{
    double Ax[] = {1, 2, 2, 1, 1, 1, 1, -1, 1, -1};
    int Ai[] = {0, 1, 0, 1, 0, 1, 0, 2, 1, 2};
    int Ap[] = {0, 2, 4, 6, 8, 10};
    int nnz = 10;
    int n_rows = 5;
    int n_cols = 3;

    double lhs[] = {-INF, -INF, 2.5, -INF, -INF};
    double rhs[] = {3, 3, INF, 0, 0};
    double lbs[] = {0, 0, -INF};
    double ubs[] = {INF, INF, INF};
    double c[] = {0, 0, 0};

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced primal infeasibility ray, zero-cost free col x2 fixed to +inf
    double y[] = {1.0, 1.0, -3.0};
    double y_orig[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);

    mu_assert("infeasible fix col inf primal ray certificate",
              is_primal_ray_valid_certificate(Ax, Ai, Ap, lhs, rhs, lbs, ubs, y_orig,
                                              n_rows, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_unbounded_fix_col_inf_dual_ray_postsolve()
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

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_rows = false;
    stgs->parallel_cols = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced dual infeasibility ray
    double x[] = {1.0, 1.0};
    double x_orig[] = {0.0, 0.0, 0.0};
    postsolve_dual_infeas_ray(presolver, x, x_orig);

    mu_assert("unbounded fix col inf dual ray certificate",
              is_dual_ray_valid_certificate(Ax, Ai, Ap, lhs, rhs, lbs, ubs, c,
                                            x_orig, n_rows, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static char *test_unbounded_negated_fix_col_inf_dual_ray_postsolve()
{
    double Ax[] = {1, -1, -1, 1, -2, 1};
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

    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->primal_propagation = false;
    stgs->parallel_rows = false;
    stgs->parallel_cols = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    run_presolver(presolver);

    // reduced dual infeasibility ray
    double x[] = {1.0, 1.0};
    double x_orig[] = {0.0, 0.0, 0.0};
    postsolve_dual_infeas_ray(presolver, x, x_orig);

    mu_assert("unbounded negated fix col inf dual ray certificate",
              is_dual_ray_valid_certificate(Ax, Ai, Ap, lhs, rhs, lbs, ubs, c,
                                            x_orig, n_rows, n_cols));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_ray_postsolve()
{
    mu_run_test(test_0_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_0_dual_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_singleton_eq_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_dton_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_bound_change_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_bound_change_dual_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_fix_col_inf_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_fix_col_inf_dual_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_fix_col_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_fix_col_dual_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_sub_col_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_sub_col_dual_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_parallel_col_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_parallel_col_dual_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_lhs_change_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_lhs_change_primal_ray_postsolve_inactive,
                counter_ray_postsolve);
    mu_run_test(test_rhs_change_primal_ray_postsolve, counter_ray_postsolve);
    mu_run_test(test_rhs_change_primal_ray_postsolve_inactive,
                counter_ray_postsolve);
    mu_run_test(test_infeasible_rhs_change_primal_ray_postsolve,
                counter_ray_postsolve);
    mu_run_test(test_infeasible_lhs_change_primal_ray_postsolve,
                counter_ray_postsolve);
    mu_run_test(test_infeasible_both_sides_changed_primal_ray_postsolve,
                counter_ray_postsolve);
    mu_run_test(test_infeasible_bound_change_primal_ray_postsolve,
                counter_ray_postsolve);
    mu_run_test(test_infeasible_parallel_col_primal_ray_postsolve,
                counter_ray_postsolve);
    mu_run_test(test_infeasible_fix_col_inf_primal_ray_postsolve,
                counter_ray_postsolve);
    mu_run_test(test_unbounded_fix_col_inf_dual_ray_postsolve,
                counter_ray_postsolve);
    mu_run_test(test_unbounded_negated_fix_col_inf_dual_ray_postsolve,
                counter_ray_postsolve);
    return NULL;
}

int test_ray_postsolve()
{
    const char *result = all_tests_ray_postsolve();
    if (result != NULL)
    {
        printf("ray_postsolve: TEST FAILED!\n");
        printf("%s\n", result);
    }
    else
    {
        printf("ray_postsolve: ALL TESTS PASSED\n");
    }
    printf("ray_postsolve: Tests run: %d\n", counter_ray_postsolve);
    return result == 0;
}

#endif // TEST_RAY_POSTSOLVE_H
