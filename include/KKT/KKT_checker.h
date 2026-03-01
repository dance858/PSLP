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
    double z_comp_slack;
    double y_comp_slack;
} KKT_checker;

KKT_checker *new_KKT_checker(
    const double *c, const double *lhs,
    const double *rhs, const double *Ax,
    const int *Ai, const int *Ap,
    const double *lbs, const double *ubs, size_t m,
    size_t n, size_t nnz);

void free_KKT_checker(KKT_checker *checker);

void KKT_checker_compute_residuals(
    KKT_checker *checker, const double *x,
    const double *y, const double *z);

bool KKT_checker_abs(KKT_checker *checker,
    const double *x, const double *y,
    const double *z, double eps_dual_res_abs,
    double eps_primal_res_abs,
    double eps_z_comp_slack_abs,
    double eps_y_comp_slack_abs);

#endif // KKT_CHECKER_H
