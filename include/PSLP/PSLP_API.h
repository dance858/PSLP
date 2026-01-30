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

/* Public header containing the outward facing API. Together with the other
   files in the folder "PSLP", it includes the input/output data structs
   and the API functions. */
#ifndef PRESOLVER_H
#define PRESOLVER_H

#ifdef __cplusplus
#include <cstdbool>
#include <cstddef> // size_t
extern "C"
{
#else
#include <stdbool.h>
#include <stddef.h> // size_t
#endif

#include "PSLP_status.h"

    /* forward declarations */
    struct Solution;
    struct Problem;
    struct PresolveStats;

    typedef struct
    {
        bool ston_cols;
        bool dton_eq;
        bool parallel_rows;
        bool parallel_cols;
        bool primal_propagation;
        // bool clean_small_coeff;
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
        size_t m;
        size_t n;
        size_t nnz;

        // lhs and rhs in the form lhs <= Ax <= rhs
        double *lhs;
        double *rhs;
        double *c;

        // variable bounds lbs <= x <= ubs
        double *lbs;
        double *ubs;

        /* offset to be added to the objective value. (When the presolver fixes
           variables, it adds an offset to the objective. This offset should
           possibly be taken into account on the solver side when evaluating the
           relative optimality gap of the reduced problem.) */
        double obj_offset;
    } PresolvedProblem;

    /* struct corresponding to the presolver:
         - 'stats' contains statistics about the presolving process
         - 'stgs' contains the settings used for presolving
         - 'prob' contains the internal problem representation used during
            presolving (never needs to be accessed by the user)
         - 'reduced_prob' contains the presolved problem after running
            'run_presolver'
         - 'sol' contains the solution to the original problem after running
            'postsolve'
    */
    typedef struct
    {
        struct PresolveStats *stats;
        const Settings *stgs;
        struct Problem *prob;
        PresolvedProblem *reduced_prob;
        struct Solution *sol;
    } Presolver;

    /* The user is responsible for freeing the settings using 'free_settings'. */
    Settings *default_settings();
    void free_settings(Settings *stgs);
    void set_settings_true(Settings *stgs);
    void set_settings_false(Settings *stgs);

    /* Initialize presolver, allocate memory, and build internal data structures.
       The presolver maintains internal deep copies of Ax, Ai, Ap, lhs, rhs, lbs,
       ubs, and c. The user is responsible for freeing the presolver using
       'free_presolver'. If the allocation fails, the function returns NULL.
       The matrix should be given in CSR form.*/
    Presolver *new_presolver(const double *Ax, const int *Ai, const int *Ap,
                             size_t m, size_t n, size_t nnz, const double *lhs,
                             const double *rhs, const double *lbs, const double *ubs,
                             const double *c, const Settings *stgs);

    /* Free the memory allocated for the presolver. */
    void free_presolver(Presolver *presolver);

    /* Runs the presolver. At completion, the 'reduced_prob' field of the
       presolver contains the presolved problem. */
    PresolveStatus run_presolver(Presolver *presolver);

    /* Postsolve the problem given the primal-dual solution (x, y, z) of the
       reduced problem. The function populates presolver->sol, so if you're
       looking for the solution to the original problem, you want to look there.
       If the solver has added the offset to the objective when solving the reduced
       problem, the optimal value of the original problem is the same as that of
       the reduced problem. */
    void postsolve(Presolver *presolver, const double *x, const double *y,
                   const double *z);

#ifdef __cplusplus
}
#endif

#endif // PRESOLVER_H
