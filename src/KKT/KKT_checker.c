/*
 * Copyright 2025 Daniel Cederberg
 *
 * This file is part of the PSLP project (LP Presolver).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#include "KKT_checker.h"
#include "KKT_utils.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "glbopts.h"
#include <stdio.h>
#include <string.h>

KKT_checker *new_KKT_checker(const double *c, const double *lhs, const double *rhs,
                             const double *Ax, const int *Ai, const int *Ap,
                             const double *lbs, const double *ubs, size_t m,
                             size_t n, size_t nnz)
{
    KKT_checker *checker = (KKT_checker *) ps_malloc(1, sizeof(KKT_checker));
    RETURN_PTR_IF_NULL(checker, NULL);

    checker->c = c;
    checker->lhs = lhs;
    checker->rhs = rhs;
    checker->lbs = lbs;
    checker->ubs = ubs;
    checker->m = m;
    checker->n = n;

    checker->A = matrix_new_no_extra_space(Ax, Ai, Ap, m, n, nnz);
    if (checker->A == NULL)
    {
        PS_FREE(checker);
        return NULL;
    }

    int *work = (int *) ps_calloc(n, sizeof(int));
    if (work == NULL)
    {
        free_matrix(checker->A);
        PS_FREE(checker);
        return NULL;
    }

    checker->AT = transpose(checker->A, work);
    PS_FREE(work);

    if (checker->AT == NULL)
    {
        free_matrix(checker->A);
        PS_FREE(checker);
        return NULL;
    }

    return checker;
}

void free_KKT_checker(KKT_checker *checker)
{
    if (checker == NULL) return;
    free_matrix(checker->A);
    free_matrix(checker->AT);
    PS_FREE(checker);
}

static double project_scalar(double val, double low, double high)
{
    if (!IS_NEG_INF(low) && val < low)
    {
        return low;
    }

    if (!IS_POS_INF(high) && val > high)
    {
        return high;
    }
    return val;
}

static void project_z_comp_slack(const double *z, const double *x, const double *lbs,
                                 const double *ubs, size_t n, double tol,
                                 double *z_proj)
{
    for (size_t j = 0; j < n; ++j)
    {
        bool is_lb_finite = !IS_ABS_INF(lbs[j]);
        bool is_ubs_finite = !IS_ABS_INF(ubs[j]);

        if (is_lb_finite && is_ubs_finite && fabs(lbs[j] - ubs[j]) <= tol)
        {
            z_proj[j] = z[j];
        }
        else if (is_lb_finite && fabs(x[j] - lbs[j]) <= tol)
        {
            z_proj[j] = MAX(z[j], 0.0);
        }
        else if (is_ubs_finite && fabs(x[j] - ubs[j]) <= tol)
        {
            z_proj[j] = MIN(z[j], 0.0);
        }
        else
        {
            z_proj[j] = 0.0;
        }
    }
}

static void project_y_comp_slack(const double *y, const double *Ax_val,
                                 const double *lhs, const double *rhs, size_t m,
                                 double tol, double *y_proj)
{
    for (size_t i = 0; i < m; ++i)
    {
        bool is_lhs_finite = !IS_ABS_INF(lhs[i]);
        bool is_rhs_finite = !IS_ABS_INF(rhs[i]);

        if (is_lhs_finite && is_rhs_finite && fabs(lhs[i] - rhs[i]) <= tol)
        {
            y_proj[i] = y[i];
        }
        else if (is_lhs_finite && fabs(Ax_val[i] - lhs[i]) <= tol)
        {
            y_proj[i] = MAX(y[i], 0.0);
        }
        else if (is_rhs_finite && fabs(Ax_val[i] - rhs[i]) <= tol)
        {
            y_proj[i] = MIN(y[i], 0.0);
        }
        else
        {
            y_proj[i] = 0.0;
        }
    }
}

void KKT_checker_compute_residuals(KKT_checker *checker, const double *x,
                                   const double *y, const double *z)
{
    size_t m = checker->m;
    size_t n = checker->n;

    double *Ax_val = (double *) ps_calloc(m, sizeof(double));
    double *ATy = (double *) ps_calloc(n, sizeof(double));
    double *work_n = (double *) ps_calloc(n, sizeof(double));
    double *work_m = (double *) ps_calloc(m, sizeof(double));

    if (!Ax_val || !ATy || !work_n || !work_m)
    {
        free(Ax_val);
        free(ATy);
        free(work_n);
        free(work_m);
        checker->dual_res = INFINITY;
        checker->primal_res = INFINITY;
        checker->z_comp_slack = INFINITY;
        checker->y_comp_slack = INFINITY;
        return;
    }

    csr_matvec(checker->A, x, Ax_val);
    csr_matvec(checker->AT, y, ATy);

    /* dual residual: ATy + z - c */
    for (size_t j = 0; j < n; ++j)
    {
        work_n[j] = ATy[j] + z[j] - checker->c[j];
    }
    checker->dual_res = norm2(work_n, n);

    /* primal residual: Ax - proj_{[bl,bu]}(Ax) */
    for (size_t i = 0; i < m; ++i)
    {
        double proj = project_scalar(Ax_val[i], checker->lhs[i], checker->rhs[i]);
        work_m[i] = Ax_val[i] - proj;
    }
    checker->primal_res = norm2(work_m, m);

    /* z complementary slackness */
    project_z_comp_slack(z, x, checker->lbs, checker->ubs, n, FEAS_TOL, work_n);
    for (size_t j = 0; j < n; ++j)
    {
        work_n[j] = z[j] - work_n[j];
    }
    checker->z_comp_slack = norm2(work_n, n);

    /* y complementary slackness */
    project_y_comp_slack(y, Ax_val, checker->lhs, checker->rhs, m, FEAS_TOL, work_m);
    for (size_t i = 0; i < m; ++i)
    {
        work_m[i] = y[i] - work_m[i];
    }
    checker->y_comp_slack = norm2(work_m, m);

    free(Ax_val);
    free(ATy);
    free(work_n);
    free(work_m);
}

bool KKT_checker_abs(KKT_checker *checker, const double *x, const double *y,
                     const double *z, double eps_dual_res_abs,
                     double eps_primal_res_abs, double eps_z_comp_slack_abs,
                     double eps_y_comp_slack_abs)
{
    KKT_checker_compute_residuals(checker, x, y, z);

    bool ok = true;
    if (checker->dual_res > eps_dual_res_abs)
    {
        printf("KKT FAIL: dual_res = %.2e "
               "(tol = %.2e)\n",
               checker->dual_res, eps_dual_res_abs);
        ok = false;
    }
    if (checker->primal_res > eps_primal_res_abs)
    {
        printf("KKT FAIL: primal_res = %.2e "
               "(tol = %.2e)\n",
               checker->primal_res, eps_primal_res_abs);
        ok = false;
    }
    if (checker->z_comp_slack > eps_z_comp_slack_abs)
    {
        printf("KKT FAIL: z_comp_slack = %.2e "
               "(tol = %.2e)\n",
               checker->z_comp_slack, eps_z_comp_slack_abs);
        ok = false;
    }
    if (checker->y_comp_slack > eps_y_comp_slack_abs)
    {
        printf("KKT FAIL: y_comp_slack = %.2e "
               "(tol = %.2e)\n",
               checker->y_comp_slack, eps_y_comp_slack_abs);
        ok = false;
    }

    return ok;
}
