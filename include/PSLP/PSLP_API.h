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
   and the API functions.

   The presolver handles problems of the form

       minimize    (1/2) x'Px + c'x
       subject to  lhs <= Ax <= rhs
                   lbs <= x <= ubs

   with P symmetric positive semidefinite. An LP is the special case P = 0
   (see new_presolver). */
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
        bool lp_fold;
        double max_time;
        bool verbose;
        /* Threads PSLP may use. 0 is automatic: the smaller of 8 and the
           processors the system reports. Clamped to [1, 64]. Ignored on
           Windows. */
        int n_threads;
    } Settings;

    /* struct corresponding to the problem to presolve,
           min (1/2) x'Px + c'x  s.t.  lhs <= Ax <= rhs,  lbs <= x <= ubs.
       Optional parts are left NULL or 0, so an LP leaves P unset. */
    typedef struct
    {
        // CSR format of a matrix A of size m x n with nnz non-zeros
        const double *Ax;
        const int *Ai;
        const int *Ap;
        size_t m;
        size_t n;
        size_t nnz;

        // lhs and rhs in the form lhs <= Ax <= rhs
        const double *lhs;
        const double *rhs;
        const double *c;

        // variable bounds lbs <= x <= ubs
        const double *lbs;
        const double *ubs;

        /* CSR format of the quadratic objective matrix P (n x n, full
           symmetric matrix with sorted column indices, Pnnz non-zeros).
           Pp == NULL or Pnnz == 0 means P = 0, i.e. an LP. */
        const double *Px;
        const int *Pi;
        const int *Pp;
        size_t Pnnz;
    } OriginalProblem;

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

        /* CSR format of the reduced quadratic objective matrix P (n x n,
           full symmetric matrix, Pnnz non-zeros). For an LP the three
           arrays are NULL and Pnnz is 0. */
        double *Px;
        int *Pi;
        int *Pp;
        size_t Pnnz;

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
            'postsolve' (it is NULL until the first postsolve call)
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
       The presolver maintains internal deep copies of all arrays of 'problem'
       and drops explicit zeros of P. The user is responsible for freeing the
       presolver using 'free_presolver'. Since the row pointers are ints, nnz and
       Pnnz must be at most INT_MAX. The function returns NULL if the allocation
       fails, nnz or Pnnz is too large, Pnnz > 0 with a NULL array of P, the CSR
       structure of P is invalid (see matrix_valid_csr_input) or P is not exactly
       symmetric. Positive semidefiniteness of P is not checked. lp_fold is
       skipped when P has entries. */
    Presolver *new_presolver_from_problem(const OriginalProblem *problem,
                                          const Settings *stgs);

    /* Same as new_presolver_from_problem for an LP, with the problem data
       given as arguments. */
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
       reduced problem, where z = Px + c - A'y are the multipliers of the
       variable bounds (z = c - A'y for an LP): z >= 0 for a variable at its
       lower bound, z <= 0 at its upper bound, y >= 0 for a row active at its
       lhs and y <= 0 at its rhs. The function populates presolver->sol, so if
       you're looking for the solution to the original problem, you want to
       look there. If the solver has added the offset to the objective when
       solving the reduced problem, the optimal value of the original problem
       is the same as that of the reduced problem. */
    void postsolve(Presolver *presolver, const double *x, const double *y,
                   const double *z);

    /* Map a primal-dual point (x, y) of the original problem to the reduced
       problem. This can be used to warm start a solver on the reduced problem.
       Must be called after 'run_presolver' has returned UNCHANGED or REDUCED.

       x has length n and y has length m (original dimensions). On return, x_red
       (length reduced_prob->n) and y_red (length reduced_prob->m) hold the
       mapped point and z_red (length reduced_prob->n) holds z_red = c_red -
       A_red^T y_red. Any output pointer may be NULL to skip that part; x_red
       requires x, y_red requires y, and z_red requires y_red. The function does
       not modify presolver->sol.

       The semantics are as follows. If (x, y) is an optimal primal-dual point
       of the original problem and 'relax_bounds' is false, (x_red, y_red,
       z_red) is an optimal primal-dual point of the reduced problem. */
    void map_original_sol_to_reduced(Presolver *presolver, const double *x,
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

       x is an unbounded direction: a recession direction of the constraint
       set with Px = 0 and c'x < 0 (for an LP simply c'x < 0); no sign
       convention involved. The recovered ray satisfies the same conditions
       for the original problem: fixed variables get a zero component and a
       substitution is a linear change of variables, so d'Pd = 0 carries
       over; for a PSD P that forces Pd = 0, and then c'd carries over too. */
    void postsolve_dual_infeas_ray(Presolver *presolver, const double *x,
                                   double *x_orig);

#ifdef __cplusplus
}
#endif

#endif // PRESOLVER_H
