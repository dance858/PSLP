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
 * Absolute residuals:
 * 1. dual_res_abs = ||A^T y + z - c||_2
 * 2. primal_res_abs = ||Ax - proj_{[lhs,rhs]}(Ax)||_2
 * 3. gap_abs = |c^T x + p(-y; lhs, rhs) + p(-z; lbs, ubs)|,
 *    where p(s; l, u) = sum_i p_i and p_i = s_i * l_i if
 *    s_i < 0, s_i * u_i if s_i > 0, 0 otherwise.
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
 *
 * Relative residuals:
 * 1. dual_res_rel = dual_res_abs / (1 + ||c||_2)
 * 2. primal_res_rel = primal_res_abs /
 *        (1 + ||(lhs_fin, rhs_fin)||_2)
 *    where lhs_fin, rhs_fin are the finite components.
 * 3. gap_rel = gap_abs /
 *        (1 + |p(-y; lhs, rhs) + p(-z; lbs, ubs)| + |c^T x|)
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
    double dual_res_abs;
    double primal_res_abs;
    double gap_abs;
    double z_comp_slack;
    double y_comp_slack;
    double dual_res_rel;
    double primal_res_rel;
    double gap_rel;
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
