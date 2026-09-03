/*
 * Copyright 2025-2026 Daniel Cederberg
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
    Settings *default_settings(void);
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

    /* Free the presolved problem data owned by 'presolver' after the caller has
       copied all arrays in 'reduced_prob'. This does not free 'presolver'; the
       caller must still call 'free_presolver' after postsolve. Normal postsolve
       and dual-infeasibility-ray postsolve remain available. Primal-
       infeasibility-ray postsolve is no longer available because it requires
       the reduced constraint matrix. */
    void free_presolver_reduced_problem(Presolver *presolver);

    /* Postsolve the problem given the primal-dual solution (x, y, z) of the
       reduced problem. The function populates presolver->sol, so if you're
       looking for the solution to the original problem, you want to look there.
       If the solver has added the offset to the objective when solving the reduced
       problem, the optimal value of the original problem is the same as that of
       the reduced problem. */
    void postsolve(Presolver *presolver, const double *x, const double *y,
                   const double *z);

    /* Map a primal-dual point (x, y) of the original problem to the reduced
       problem, e.g. to warm start a solver on the reduced problem from a
       solution of the original one. Must be called after 'run_presolver' has
       returned UNCHANGED or REDUCED.

       x has length n and y has length m (original dimensions). On return,
       x_red (length reduced_prob->n) and y_red (length reduced_prob->m) hold
       the mapped point, and z_red (length reduced_prob->n) holds the reduced
       dual slack z_red = c_red - A_red^T y_red. Any of the output pointers may
       be NULL to skip that part; x_red requires x, y_red requires y, and z_red
       requires y_red (a violated requirement prints a warning and leaves that
       output untouched). The function does not modify presolver->sol.

       The map replays the recorded reductions forward: removed rows/columns
       are dropped, merged parallel columns are aggregated, and multipliers of
       merged/eliminated rows are transferred to the rows that replaced them.
       x_red is projected onto the bounds of the reduced problem. Note that the
       result is a sensible starting point, not necessarily a feasible or
       optimal point of the reduced problem. In particular, if (x, y) is optimal
       for the original problem, (x_red, y_red) is optimal for the reduced
       problem up to the multipliers of rows that presolve removed as redundant,
       which are discarded.

       After 'free_presolver_reduced_problem', x_red and y_red are still
       computed, but x_red is not projected onto the reduced bounds and z_red is
       skipped (with a warning) since the reduced problem data is gone. */
    void map_solution_to_reduced(Presolver *presolver, const double *x,
                                 const double *y, double *x_red, double *y_red,
                                 double *z_red);

    /* Postsolve a primal infeasibility ray y of the reduced problem.
       The function writes the corresponding ray for the original problem
       to y_orig. It does not check whether y is a valid ray.

       y uses the Farkas sign convention, yi >= 0 when row i is active at its
       rhs, the opposite of postsolve(). It matches Gurobi's FarkasDual. */
    void postsolve_primal_infeas_ray(Presolver *presolver, const double *y,
                                     double *y_orig);

    /* Postsolve a dual infeasibility ray x of the reduced problem.
       The function writes the corresponding ray for the original problem
       to x_orig. It does not check whether x is a valid ray.

       x is an unbounded direction with c'x < 0; no sign convention involved. */
    void postsolve_dual_infeas_ray(Presolver *presolver, const double *x,
                                   double *x_orig);

#ifdef __cplusplus
}
#endif

#endif // PRESOLVER_H
