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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* Public header containing the outward facing API. It includes all the
   input/output data structs and the API functions.  Make sure this file is
   somewhere appropriate and then use `#include <PSLP/"API.h"
>` to access the
   public API. */
#ifndef PRESOLVER_H
#define PRESOLVER_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "PresolveStatus.h"
#include "Sol.h"
#include <stdbool.h>

    /* forward declaration */
    struct Problem;
    struct PresolveStats;

    typedef struct
    {
        bool ston_cols;
        bool dton_eq;
        bool parallel_rows;
        bool parallel_cols;
        bool primal_propagation;
        bool clean_small_coeff;
        bool finite_bound_tightening;
        bool dual_fix;
        bool relax_bounds;
        int max_shift;
        double max_time;
        bool verbose;
    } Settings;

    /* struct corresponding to the presolved problem */
    typedef struct
    {
        // CSR format of a matrix A of size m x n with nnz non-zeros
        double *Ax;
        int *Ai;
        int *Ap;
        int m;
        int n;
        int nnz;

        // lhs and rhs in the form lhs <= Ax <= rhs
        double *lhs;
        double *rhs;
        double *c;

        // variable bounds bounds[k].lb <= x_k <= bounds[k].ub
        // Bound *bounds;
        double *lbs;
        double *ubs;
    } PresolvedProblem;

    /* struct corresponding to the presolver*/
    typedef struct
    {
        struct PresolveStats *stats;
        const Settings *stgs;
        struct Problem *prob;
        PresolvedProblem *reduced_prob;
        Solution *sol;
    } Presolver;

    /* The user is responsible for freeing the settings struct using standard free.
     */
    Settings *default_settings();
    void free_settings(Settings *stgs);
    void set_settings_true(Settings *stgs);
    void set_settings_false(Settings *stgs);

    /* Initialize presolver, allocate memory, and build internal data structures.
       The presolver maintains internal deep copies of Ax, Ai, Ap, lhs, rhs, lbs,
       ubs, and c. The user is responsible for freeing this memory using
       'free_presolver'. If the allocation fails, the function returns NULL. */
    Presolver *new_presolver(const double *Ax, const int *Ai, const int *Ap, int m,
                             int n, int nnz, const double *lhs, const double *rhs,
                             const double *lbs, const double *ubs, const double *c,
                             const Settings *stgs, bool CSR);

    /* Free the memory allocated for the presolver. */
    void free_presolver(Presolver *presolver);

    /* Runs the presolver. At completion, the 'problem' field of the presolver
       contains the presolved problem. */
    PresolveStatus run_presolver(Presolver *presolver);

    /* Postsolve the problem given the primal-dual solution (x, y, z) of the
       reduced problem. The optimal value of the reduced problem is 'obj'.
       The function populates presolver->sol. */
    void postsolve(Presolver *presolver, const double *x, const double *y,
                   const double *z, double obj);

#ifdef __cplusplus
}
#endif

#endif // PRESOLVER_H
