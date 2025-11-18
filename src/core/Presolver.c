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

#include "Activity.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "DTonsEq.h"
#include "Debugger.h"
#include "Locks.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "PSLP_API.h"
#include "Parallel_cols.h"
#include "Parallel_rows.h"
#include "Postsolver.h"
#include "PresolveStats.h"
#include "Primal_propagation.h"
#include "Problem.h"
#include "RowColViews.h"
#include "SimpleReductions.h"
#include "Simple_dual_fix.h"
#include "State.h"
#include "StonCols.h"
#include "Tags.h"
#include "Timer.h"
#include "Workspace.h"
#include "dVec.h"
#include "glbopts.h"
#include "iVec.h"
#include <stdio.h>
#include <stdlib.h>

#define NNZ_REMOVED_FAST 0.95
#define NNZ_REMOVED_CYCLE 0.95

typedef uint8_t Complexity;

enum Complexity
{
    FAST = 0,
    MEDIUM = 1 << 0
};

PresolveStats *init_stats(int n_rows, int n_cols, int nnz)
{
    PresolveStats *stats = (PresolveStats *) ps_calloc(1, sizeof(PresolveStats));
    stats->n_rows_original = n_rows;
    stats->n_cols_original = n_cols;
    stats->nnz_original = nnz;
    return stats;
}

void free_settings(Settings *stgs)
{
    PS_FREE(stgs);
}

Settings *default_settings()
{
    Settings *stgs = (Settings *) ps_malloc(1, sizeof(Settings));
    RETURN_PTR_IF_NULL(stgs, NULL);
    stgs->ston_cols = true;
    stgs->dton_eq = true;
    stgs->parallel_rows = true;
    stgs->parallel_cols = true;
    stgs->primal_propagation = true;
    stgs->dual_fix = true;
    stgs->clean_small_coeff = false;
    stgs->finite_bound_tightening = true;
    stgs->relax_bounds = true;
    stgs->max_shift = 10;
    stgs->max_time = 60.0;
    stgs->verbose = true;
    return stgs;
}

void set_settings_true(Settings *stgs)
{
    stgs->ston_cols = true;
    stgs->dton_eq = true;
    stgs->parallel_rows = true;
    stgs->parallel_cols = true;
    stgs->primal_propagation = true;
    stgs->dual_fix = true;
    stgs->clean_small_coeff = true;
    stgs->finite_bound_tightening = true;
    stgs->relax_bounds = true;
    stgs->verbose = true;
}

void set_settings_false(Settings *stgs)
{
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    stgs->parallel_rows = false;
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    stgs->dual_fix = false;
    stgs->clean_small_coeff = false;
    stgs->finite_bound_tightening = false;
    stgs->relax_bounds = false;
    stgs->verbose = false;
}

typedef struct clean_up_scope
{
    Matrix *A, *AT;
    double *lhs_copy, *rhs_copy, *c_copy;
    int *row_sizes, *col_sizes;
    RowTag *row_tags;
    Lock *locks;
    Activity *activities;
    State *data;
    Constraints *constraints;
    Presolver *presolver;
    Objective *obj;
    ColTag *col_tags;
    Bound *bounds;
    Work *work;
} clean_up_scope;

void presolver_clean_up(clean_up_scope scope)
{
    free_matrix(scope.A);
    free_matrix(scope.AT);
    PS_FREE(scope.lhs_copy);
    PS_FREE(scope.rhs_copy);
    PS_FREE(scope.c_copy);
    PS_FREE(scope.row_sizes);
    PS_FREE(scope.col_sizes);
    PS_FREE(scope.row_tags);
    PS_FREE(scope.col_tags);
    PS_FREE(scope.bounds);
    free_activities(scope.activities);
    free_locks(scope.locks);
    free_work(scope.work);
    free_state(scope.data);
    free_constraints(scope.constraints);
    free_objective(scope.obj);
    free_presolver(scope.presolver);
}

Presolver *new_presolver(const double *Ax, const int *Ai, const int *Ap, int m,
                         int n, int nnz, const double *lhs, const double *rhs,
                         const double *lbs, const double *ubs, const double *c,
                         const Settings *stgs, bool CSR)
{
    Timer timer;
    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    Matrix *A = NULL, *AT = NULL;
    double *lhs_copy = NULL, *rhs_copy = NULL, *c_copy = NULL;
    int *row_sizes = NULL, *col_sizes = NULL;
    RowTag *row_tags = NULL;
    Lock *locks = NULL;
    Activity *activities = NULL;
    State *data = NULL;
    Constraints *constraints = NULL;
    Presolver *presolver = NULL;
    Objective *obj = NULL;
    ColTag *col_tags = NULL;
    Bound *bounds = NULL;
    Work *work = NULL;

    //  ---------------------------------------------------------------------------
    //   Copy data and allocate memory. For an exact specification of the
    //   workspace, see the workspace class. The problem object owns all the
    //   memory that is allocated in this block of code (and therefore frees it).
    //  ---------------------------------------------------------------------------
    int n_rows = m;
    int n_cols = n;
    lhs_copy = (double *) ps_malloc(n_rows, sizeof(double));
    rhs_copy = (double *) ps_malloc(n_rows, sizeof(double));
    c_copy = (double *) ps_malloc(n_cols, sizeof(double));
    col_tags = (ColTag *) ps_calloc(n_cols, sizeof(ColTag));
    bounds = (Bound *) ps_malloc(n_cols, sizeof(Bound));
    work = new_work(n_rows, n_cols);
    row_sizes = (int *) ps_malloc(n_rows, sizeof(int));
    col_sizes = (int *) ps_malloc(n_cols, sizeof(int));

    if (!lhs_copy || !rhs_copy || !c_copy || !col_tags || !bounds || !work ||
        !row_sizes || !col_sizes)
    {
        goto cleanup;
    }

    memcpy(lhs_copy, lhs, n_rows * sizeof(double));
    memcpy(rhs_copy, rhs, n_rows * sizeof(double));
    memcpy(c_copy, c, n_cols * sizeof(double));

    // ---------------------------------------------------------------------------
    //                  Build bounds, row tags, A and AT.
    // ---------------------------------------------------------------------------
    row_tags = new_rowtags(lhs_copy, rhs_copy, n_rows);
    if (!row_tags) goto cleanup;

    for (int i = 0; i < n_cols; i++)
    {
        bounds[i].lb = lbs[i];
        bounds[i].ub = ubs[i];

        if (IS_NEG_INF(lbs[i]))
        {
            UPDATE_TAG(col_tags[i], C_TAG_LB_INF);
        }

        if (IS_POS_INF(ubs[i]))
        {
            UPDATE_TAG(col_tags[i], C_TAG_UB_INF);
        }
    }

    if (CSR)
    {
        A = matrix_new_no_extra_space(Ax, Ai, Ap, n_rows, n_cols, nnz);
        if (!A) goto cleanup;

        if (stgs->clean_small_coeff)
        {
            clean_small_coeff_A(A, bounds, row_tags, col_tags, rhs_copy, lhs_copy);
        }

        AT = transpose(A, work->iwork_n_cols);
        if (!AT) goto cleanup;
    }
    else
    {
        AT = matrix_new(Ax, Ai, Ap, n_cols, n_rows, nnz);
        if (!AT) goto cleanup;
        A = transpose(AT, work->iwork_n_rows);
        if (!A) goto cleanup;
    }

    // ---------------------------------------------------------------------------
    //           locks, activities, row sizes and column sizes
    // ---------------------------------------------------------------------------
    locks = new_locks(A, row_tags);
    activities = new_activities(A, col_tags, bounds);
    if (!locks || !activities) goto cleanup;
    count_rows(A, row_sizes);
    count_rows(AT, col_sizes);

    // ---------------------------------------------------------------------------
    //  Initialize internal data and constraints
    // ---------------------------------------------------------------------------
    data = new_state(row_sizes, col_sizes, locks, n_rows, n_cols, activities, work,
                     row_tags);

    if (!data) goto cleanup;
    constraints =
        constraints_new(A, AT, lhs_copy, rhs_copy, bounds, data, row_tags, col_tags);
    if (!constraints) goto cleanup;

    // ---------------------------------------------------------------------------
    //             Allocate the actual presolver
    // ---------------------------------------------------------------------------
    presolver = (Presolver *) ps_malloc(1, sizeof(Presolver));
    obj = objective_new(c_copy);
    if (!presolver || !obj) goto cleanup;
    presolver->stgs = stgs;
    presolver->prob = new_problem(constraints, obj);
    presolver->stats = init_stats(A->m, A->n, nnz);
    presolver->stats->nnz_removed_trivial = nnz - A->nnz; // due to clean_small_coeff
    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    presolver->stats->ps_time_init = GET_ELAPSED_SECONDS(timer);
    presolver->reduced_prob =
        (PresolvedProblem *) ps_calloc(1, sizeof(PresolvedProblem));
    DEBUG(run_debugger(constraints, false));

    // ---------------------------------------------------------------------------
    //           Allocate space for returning the solution
    // ---------------------------------------------------------------------------
    presolver->sol = (Solution *) ps_malloc(1, sizeof(Solution));
    if (!presolver->sol) goto cleanup;
    presolver->sol->x = (double *) ps_malloc(n_cols, sizeof(double));
    presolver->sol->y = (double *) ps_malloc(n_rows, sizeof(double));
    presolver->sol->z = (double *) ps_malloc(n_cols, sizeof(double));
    presolver->sol->dim_x = n_cols;
    presolver->sol->dim_y = n_rows;
    if (!presolver->sol->x || !presolver->sol->y || !presolver->sol->z)
    {
        goto cleanup;
    }

    return presolver;

cleanup:
{
    struct clean_up_scope scope = {
        A,         AT,       lhs_copy, rhs_copy,   c_copy, row_sizes,
        col_sizes, row_tags, locks,    activities, data,   constraints,
        presolver, obj,      col_tags, bounds,     work};
    presolver_clean_up(scope);
}
    return NULL;
}

/* This function updates the termination criterion for the presolver.
   We terminate if one of the following conditions is satisfied:
   1. The problem is infeasible or unbounded.
   2. A cycle reduced the number of non-zeros by less than
     (1 - NNZ_REMOVED_CYCLE)%.
   3. A maximum time limit is reached.
*/
static inline bool update_termination(int nnz_after_cycle, int nnz_before_cycle,
                                      Complexity complexity, PresolveStatus status,
                                      double max_time, Timer outer_timer)
{
    if (HAS_STATUS(status, UNBNDORINFEAS))
    {
        return true;
    }

    if (complexity == MEDIUM &&
        nnz_after_cycle >= NNZ_REMOVED_CYCLE * nnz_before_cycle)
    {
        return true;
    }

    clock_gettime(CLOCK_MONOTONIC, &outer_timer.end);

    if (GET_ELAPSED_SECONDS(outer_timer) >= max_time)
    {
        printf("Maximum time limit of %.2f seconds reached.\n", max_time);
        return true;
    }

    return false;
}

/* This function updates the complexity class according to
   the following rules:
    * If the current complexity is fast, the next complexity is
      fast if sufficiently many non-zeros were removed, otherwise
      it is medium.
    * If the current complexity is medium, the next complexity is
      fast.
*/
static inline Complexity update_complexity(Complexity curr_complexity,
                                           int nnz_before_phase, int nnz_after_phase)
{
    if (curr_complexity == FAST)
    {
        if (nnz_after_phase < NNZ_REMOVED_FAST * nnz_before_phase)
        {
            return FAST;
        }

        return MEDIUM;
    }
    else if (curr_complexity == MEDIUM)
    {
        return FAST;
    }
    else
    {
        assert(false);
    }

    return FAST; // to suppress compiler warning
}

/* This function exhaustively removes singleton rows from the problem. It
   also eliminates empty rows and empty columns. If the problem is feasible
   and bounded, it returns UNCHANGED (even if the problem was reduced).
*/
static inline PresolveStatus run_trivial_explorers(Problem *prob,
                                                   const Settings *stgs)
{
    PresolveStatus status = UNCHANGED;

    remove_variables_with_close_bounds(prob);

    if (stgs->dual_fix)
    {
        //   after simple dual fix there can be new empty rows and singleton rows,
        //   and empty columns (if rows are marked as inactive due to a variable
        //   being set to inf)
        status |= remove_empty_cols(prob);
        status |= simple_dual_fix(prob);
        RETURN_IF_UNBNDORINFEAS(status);
    }

    do
    {
        status = remove_ston_rows(prob);
        RETURN_IF_INFEASIBLE(status);
    } while (status == REDUCED);

    status = remove_empty_rows(prob->constraints);
    RETURN_IF_INFEASIBLE(status);

    status = remove_empty_cols(prob);
    RETURN_IF_UNBNDORINFEAS(status);

    assert(prob->constraints->state->ston_rows->len == 0);
    assert(prob->constraints->state->empty_rows->len == 0);
    assert(prob->constraints->state->empty_cols->len == 0);
    return UNCHANGED;
}

/* This function applies the fastest reductions, which are doubleton equality
   rows, singleton columns, and simple dual fix. It returns UNCHANGED even if
   the problem was reduced (provided that the problem does not seem infeasible or
   bounded). */
static inline PresolveStatus run_fast_explorers(Problem *prob, const Settings *stgs)
{
    assert(prob->constraints->state->ston_rows->len == 0);
    assert(prob->constraints->state->empty_rows->len == 0);
    assert(prob->constraints->state->empty_cols->len == 0);
    PresolveStatus status = UNCHANGED;

    if (stgs->ston_cols)
    {
        status |= remove_ston_cols(prob);

        // after removing singleton columns, there can be new empty columns
        // and singleton rows, but no empty rows
        assert(prob->constraints->state->empty_rows->len == 0);
        status |= run_trivial_explorers(prob, stgs);
    }

    if (stgs->dton_eq)
    {
        status |= remove_dton_eq_rows(prob, stgs->max_shift);
        // after removing doubleton equality rows, there can be new empty rows,
        // new singleton rows, and new empty columns
        status |= run_trivial_explorers(prob, stgs);
    }
    return status;
}

/* This function applies the medium complexity presolvers, which are domain
   propagation and parallel rows. It returns UNCHANGED even if the problem
   was reduced (provided that the problem is infeasible and bounded). */
static inline PresolveStatus
run_medium_explorers(Problem *prob, const Settings *stgs, PresolveStats *stats)
{
    assert(prob->constraints->state->ston_rows->len == 0);
    assert(prob->constraints->state->empty_rows->len == 0);
    assert(prob->constraints->state->empty_cols->len == 0);
    DEBUG(verify_no_duplicates_sort(prob->constraints->state->updated_activities));
    int nnz_before;
    int *nnz = &prob->constraints->A->nnz;
    PresolveStatus status = UNCHANGED;
    Timer timer;

    if (stgs->primal_propagation)
    {
        // stats
        clock_gettime(CLOCK_MONOTONIC, &timer.start);
        nnz_before = *nnz;

        // the actual reduction
        status |= check_activities(prob);
        status |= propagate_primal(prob, stgs->finite_bound_tightening);

        // after dom prop propagation there can be new empty and ston rows rows
        status |= run_trivial_explorers(prob, stgs);

        // stats
        stats->nnz_removed_primal_propagation += nnz_before - *nnz;
        clock_gettime(CLOCK_MONOTONIC, &timer.end);
        stats->ps_time_primal_propagation += GET_ELAPSED_SECONDS(timer);
    }

    if (stgs->parallel_rows)
    {
        // stats
        nnz_before = *nnz;
        clock_gettime(CLOCK_MONOTONIC, &timer.start);

        // the actual reduction (after removing parallel rows there can
        // be no new empty or ston rows so no need to run trivial presolvers)
        status |= remove_parallel_rows(prob->constraints);
        assert(prob->constraints->state->empty_rows->len == 0);
        assert(prob->constraints->state->ston_rows->len == 0);

        // stats
        stats->nnz_removed_parallel_rows += nnz_before - *nnz;
        clock_gettime(CLOCK_MONOTONIC, &timer.end);
        stats->ps_time_parallel_rows += GET_ELAPSED_SECONDS(timer);
    }

    if (stgs->parallel_cols)
    {
        // stats
        nnz_before = *nnz;
        clock_gettime(CLOCK_MONOTONIC, &timer.start);

        // the actual reduction (after removing parallel columns there can
        // be no new empty rows or empty cols but might be new stonrows, probably
        // although that's rare but it does happen
        status |= remove_parallel_cols(prob);
        assert(prob->constraints->state->empty_rows->len == 0);
        assert(prob->constraints->state->empty_cols->len == 0);
        status |= run_trivial_explorers(prob, stgs);
        assert(prob->constraints->state->ston_rows->len == 0);

        // stats
        stats->nnz_removed_parallel_cols += nnz_before - *nnz;
        clock_gettime(CLOCK_MONOTONIC, &timer.end);
        stats->ps_time_parallel_cols += GET_ELAPSED_SECONDS(timer);
    }

    return status;
}

void populate_presolved_problem(Presolver *presolver)
{
    PresolvedProblem *ps_prob = presolver->reduced_prob;
    Constraints *constraints = presolver->prob->constraints;
    Matrix *A = constraints->A;
    ps_prob->m = A->m;
    ps_prob->n = A->n;
    ps_prob->nnz = A->nnz;
    ps_prob->Ax = A->x;
    ps_prob->Ai = A->i;
    ps_prob->rhs = constraints->rhs;
    ps_prob->lhs = constraints->lhs;
    ps_prob->c = presolver->prob->obj->c;

    // create bounds arrays
    ps_prob->lbs = (double *) malloc(A->n * sizeof(double));
    ps_prob->ubs = (double *) malloc(A->n * sizeof(double));
    for (int i = 0; i < A->n; i++)
    {
        ps_prob->lbs[i] = constraints->bounds[i].lb;
        ps_prob->ubs[i] = constraints->bounds[i].ub;
    }

    // create row pointers
    ps_prob->Ap = (int *) malloc((A->m + 1) * sizeof(int));
    for (int i = 0; i < A->m + 1; i++)
    {
        ps_prob->Ap[i] = A->p[i].start;
    }
}

static inline void print_start_message(const PresolveStats *stats)
{
    printf("\n\t       PSLP v%s - LP presolver \n\t(c) Daniel "
           "Cederberg, Stanford University, 2025\n",
           PSLP_presolve_VERSION);
    printf("Original problem:  %d rows, %d columns, %d nnz\n",
           stats->n_rows_original, stats->n_cols_original, stats->nnz_original);
}

static inline void print_end_message(const Matrix *A, const PresolveStats *stats)
{
    printf("Presolved problem: %d rows, %d columns, %d nnz\n", A->m, A->n, A->nnz);

    printf("Presolver init & run time : %.3f seconds, %.3f \n", stats->ps_time_init,
           stats->presolve_total_time);
}

PresolveStatus run_presolver(Presolver *presolver)
{
    Timer inner_timer, outer_timer;
    int nnz_before_cycle, nnz_after_cycle;
    int nnz_before_phase, nnz_after_phase;
    int nnz_before_reduction;
    Problem *prob = presolver->prob;
    PresolveStats *stats = presolver->stats;
    Matrix *A = presolver->prob->constraints->A;
    const Settings *stgs = presolver->stgs;
    Complexity curr_complexity = FAST;
    bool terminate = false;
    PresolveStatus status = UNCHANGED;
    clock_gettime(CLOCK_MONOTONIC, &outer_timer.start);

    if (stgs->verbose)
    {
        print_start_message(stats);
    }

    DEBUG(run_debugger(prob->constraints, false));

    // ------------------------------------------------------------------------
    // Main loop for the presolver. The logic is organized into phases and
    // cycles. A phase is a sequence of presolvers belonging to a certain
    // complexity class. The cycle resets after the medium presolvers.
    // ------------------------------------------------------------------------
    nnz_before_cycle = A->nnz;
    while (!terminate)
    {
        // before each phase we run the trivial presolvers
        nnz_before_reduction = A->nnz;
        status = run_trivial_explorers(prob, stgs);
        stats->nnz_removed_trivial += nnz_before_reduction - A->nnz;
        RETURN_IF_NOT_UNCHANGED(status); // TODO: this name is misleading!
        DEBUG(run_debugger(prob->constraints, false));
        nnz_before_phase = A->nnz;

        // apply presolvers belonging to a certain complexity class
        if (curr_complexity == FAST)
        {
            nnz_before_reduction = A->nnz;
            RUN_AND_TIME(run_fast_explorers, inner_timer, stats->ps_time_fast,
                         status, prob, stgs);
            stats->nnz_removed_fast += nnz_before_reduction - A->nnz;
        }
        else if (curr_complexity == MEDIUM)
        {
            RUN_AND_TIME(run_medium_explorers, inner_timer, stats->ps_time_medium,
                         status, prob, stgs, stats);
            nnz_after_cycle = A->nnz;
        }

        nnz_after_phase = A->nnz;

        terminate =
            update_termination(nnz_after_cycle, nnz_before_cycle, curr_complexity,
                               status, stgs->max_time, outer_timer);

        if (curr_complexity == MEDIUM)
        {
            // a new cycle starts after the medium presolvers
            nnz_before_cycle = nnz_after_cycle;
        }

        curr_complexity =
            update_complexity(curr_complexity, nnz_before_phase, nnz_after_phase);
    }

    if (stgs->relax_bounds)
    {
        remove_redundant_bounds(prob->constraints);
    }

    problem_clean(prob, true);
    DEBUG(run_debugger(prob->constraints, true));

    clock_gettime(CLOCK_MONOTONIC, &outer_timer.end);
    stats->n_rows_reduced = A->m;
    stats->n_cols_reduced = A->n;
    stats->nnz_reduced = A->nnz;
    stats->presolve_total_time = GET_ELAPSED_SECONDS(outer_timer);
    DEBUG(run_debugger_stats_consistency_check(stats));
    populate_presolved_problem(presolver);

    if (stgs->verbose)
    {
        print_end_message(A, stats);
    }

    return status;
}

void postsolve(Presolver *presolver, const double *x, const double *y,
               const double *z, double obj)
{
    Timer timer;
    Solution *sol = presolver->sol;
    PresolveStats *stats = presolver->stats;
    State *data = presolver->prob->constraints->state;
    PostsolveInfo *postsolve_info = data->postsolve_info;
    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    postsolver_update(postsolve_info, stats->n_cols_reduced, stats->n_rows_reduced,
                      data->work->mappings->cols, data->work->mappings->rows);
    postsolver_run(postsolve_info, sol, x, y, z);
    sol->obj = obj + presolver->prob->obj->offset;
    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    stats->ps_time_post_solve = GET_ELAPSED_SECONDS(timer);

    if (presolver->stgs->verbose)
    {
        printf("Postsolve time: %.4f seconds\n", stats->ps_time_post_solve);
    }
}

void free_presolver(Presolver *presolver)
{
    if (presolver == NULL)
    {
        return;
    }

    free_problem(presolver->prob);
    PS_FREE(presolver->stats);

    if (presolver->reduced_prob)
    {
        PS_FREE(presolver->reduced_prob->Ap);
        PS_FREE(presolver->reduced_prob->lbs);
        PS_FREE(presolver->reduced_prob->ubs);
        PS_FREE(presolver->reduced_prob);
    }

    if (presolver->sol)
    {
        PS_FREE(presolver->sol->x);
        PS_FREE(presolver->sol->y);
        PS_FREE(presolver->sol->z);
        PS_FREE(presolver->sol);
    }

    PS_FREE(presolver);
}
