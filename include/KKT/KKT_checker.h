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

#ifndef KKT_CHECKER_H
#define KKT_CHECKER_H

#include "Matrix.h"
#include <stdbool.h>
#include <stddef.h>

/* We consider an LP of the form
 *   min c^T x
 *   s.t. lhs <= Ax <= rhs
 *        lbs <= x <= ubs
 *
 * The absolute residuals are defined like this:
 * 1. dual_res = ||A^T y + z - c||_2
 * 2. primal_res = ||Ax - proj_{[lhs,rhs]}(Ax)||_2
 * 3. gap = |c^T x + p(-y; lhs, rhs) + p(-z; lbs, ubs)|, where p is defined as
 *    p(s; l, u) = sum_i p(s_i; l_i, u_i) and p(s_i; l_i, u_i) = s_i * l_i if
 *    s_i < 0, s_i * u_i if s_i > 0, and 0 if s_i = 0.
 * 4. y_comp_slack = ||y - y(x)||_2, where y_i(x) is defined as
 *    (equalities checked with FEAS_TOL tolerance):
 *    - y_i          if lhs_i = rhs_i (both finite)
 *    - max(y_i, 0)  if (Ax)_i = lhs_i (finite)
 *    - min(y_i, 0)  if (Ax)_i = rhs_i (finite)
 *    - 0            otherwise
 * 5. z_comp_slack = ||z - z(x)||_2, where z_j(x) is defined as
 *    (equalities checked with FEAS_TOL tolerance):
 *    - z_j          if lbs_j = ubs_j (both finite)
 *    - max(z_j, 0)  if x_j = lbs_j (finite)
 *    - min(z_j, 0)  if x_j = ubs_j (finite)
 *    - 0            otherwise
 */

typedef struct KKT_checker
{
    const double *c;
    const double *lhs;
    const double *rhs;
    const double *lbs;
    const double *ubs;
    size_t m;
    size_t n;
    Matrix *A;
    Matrix *AT;
    double dual_res;
    double primal_res;
    double gap;
    double z_comp_slack;
    double y_comp_slack;
} KKT_checker;

KKT_checker *new_KKT_checker(const double *c, const double *lhs, const double *rhs,
                             const double *Ax, const int *Ai, const int *Ap,
                             const double *lbs, const double *ubs, size_t m,
                             size_t n, size_t nnz);

void free_KKT_checker(KKT_checker *checker);

void KKT_checker_compute_residuals(KKT_checker *checker, const double *x,
                                   const double *y, const double *z);

bool KKT_checker_abs(KKT_checker *checker, const double *x, const double *y,
                     const double *z, double eps_dual_res_abs,
                     double eps_primal_res_abs, double eps_gap_abs,
                     double eps_z_comp_slack_abs, double eps_y_comp_slack_abs);

#endif // KKT_CHECKER_H
