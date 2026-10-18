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

#include "Activity.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "DtonsEq.h"
#include "Locks.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "PSLP_API.h"
#include "PSLP_sol.h"
#include "PSLP_stats.h"
#include "PSLP_status.h"
#include "Parallel_cols.h"
#include "Parallel_rows.h"
#include "Postsolver.h"
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
#include "fold.h"
#include "fold_apply.h"
#include "fold_postsolve.h"
#include "glbopts.h"
#include "iVec.h"
#include "pslp_thread.h"
#include "radix_sort.h"
#include "thread_pool.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define NNZ_REMOVED_FAST 0.95
#define NNZ_REMOVED_CYCLE 0.95
#define LP_FOLD_MIN_REDUCTION 0.01 /* removable nnz required, as a share of nnz */
#define LP_FOLD_WORK_BUDGET 100.0
#define LP_FOLD_MAX_ROUNDS 8192

typedef uint8_t Complexity;

enum Complexity
{
    FAST = 0,
    MEDIUM = 1 << 0
};

PresolveStats *init_stats(size_t n_rows, size_t n_cols, size_t nnz)
{
    PresolveStats *stats = (PresolveStats *) ps_calloc(1, sizeof(PresolveStats));
    stats->n_rows_original = n_rows;
    stats->n_cols_original = n_cols;
    stats->nnz_original = nnz;
    stats->n_rows_folded = n_rows;
    stats->n_cols_folded = n_cols;
    return stats;
}

void free_settings(Settings *stgs)
{
    PS_FREE(stgs);
}

Settings *default_settings(void)
{
    Settings *stgs = (Settings *) ps_malloc(1, sizeof(Settings));
    RETURN_PTR_IF_NULL(stgs, NULL);
    stgs->ston_cols = true;
    stgs->dton_eq = true;
    stgs->parallel_rows = true;
    stgs->parallel_cols = true;
    stgs->primal_propagation = true;
    stgs->dual_fix = true;
    stgs->finite_bound_tightening = true;
    stgs->relax_bounds = false;
    stgs->lp_fold = false;
    stgs->max_time = 60.0;
    stgs->verbose = true;
    stgs->n_threads = 0;
    return stgs;
}

// Helper struct for parallel initialization work
typedef struct
{
    Matrix *A;
    Work *work;
    size_t n_cols;
    size_t n_rows;
    const double *lbs;
    const double *ubs;
    double *lhs_copy;
    double *rhs_copy;
    Bound *bounds;
    ColTag *col_tags;
    RowTag *row_tags;
    Lock *locks;
    Activity *activities;
    int *row_sizes;
    bool derived_data; /* locks, activities and row sizes; false with the fold */
} ParallelInitData;

// thread function for initializing data
static void *init_thread_func(void *arg)
{
    ParallelInitData *data = (ParallelInitData *) arg;

    data->row_tags = new_rowtags(data->lhs_copy, data->rhs_copy, data->n_rows);

    for (int i = 0; i < data->n_cols; i++)
    {
        data->bounds[i].lb = data->lbs[i];
        data->bounds[i].ub = data->ubs[i];

        if (IS_NEG_INF(data->lbs[i]))
        {
            UPDATE_TAG(data->col_tags[i], C_TAG_LB_INF);
        }

        if (IS_POS_INF(data->ubs[i]))
        {

            UPDATE_TAG(data->col_tags[i], C_TAG_UB_INF);
        }
    }

    /* With the fold the matrix may still change; the caller builds these
       afterwards. */
    if (data->derived_data)
    {
        data->locks = new_locks(data->A, data->row_tags);
        data->activities = new_activities(data->A, data->col_tags, data->bounds);
        count_rows(data->A, data->row_sizes);
    }

    return NULL;
}

void set_settings_true(Settings *stgs)
{
    stgs->ston_cols = true;
    stgs->dton_eq = true;
    stgs->parallel_rows = true;
    stgs->parallel_cols = true;
    stgs->primal_propagation = true;
    stgs->dual_fix = true;
    stgs->finite_bound_tightening = true;
    stgs->relax_bounds = true;
    stgs->lp_fold = false;
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
    stgs->finite_bound_tightening = false;
    stgs->relax_bounds = false;
    stgs->lp_fold = false;
    stgs->verbose = false;
}

typedef struct clean_up_scope
{
    Matrix *A, *AT, *P;
    double *lhs_copy, *rhs_copy, *c_copy;
    int *row_sizes, *col_sizes, *P_row_sizes;
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
    free_matrix(scope.P);
    PS_FREE(scope.lhs_copy);
    PS_FREE(scope.rhs_copy);
    PS_FREE(scope.c_copy);
    PS_FREE(scope.row_sizes);
    PS_FREE(scope.col_sizes);
    PS_FREE(scope.P_row_sizes);
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

Presolver *new_presolver(const double *Ax, const int *Ai, const int *Ap, size_t m,
                         size_t n, size_t nnz, const double *lhs, const double *rhs,
                         const double *lbs, const double *ubs, const double *c,
                         const Settings *stgs)
{
    OriginalProblem problem = {.Ax = Ax,
                               .Ai = Ai,
                               .Ap = Ap,
                               .m = m,
                               .n = n,
                               .nnz = nnz,
                               .lhs = lhs,
                               .rhs = rhs,
                               .c = c,
                               .lbs = lbs,
                               .ubs = ubs};
    return new_presolver_from_problem(&problem, stgs);
}

Presolver *new_presolver_from_problem(const OriginalProblem *problem,
                                      const Settings *stgs)
{
    const double *Ax = problem->Ax;
    const int *Ai = problem->Ai;
    const int *Ap = problem->Ap;
    size_t m = problem->m;
    size_t n = problem->n;
    size_t nnz = problem->nnz;
    const double *lhs = problem->lhs;
    const double *rhs = problem->rhs;
    const double *c = problem->c;
    const double *lbs = problem->lbs;
    const double *ubs = problem->ubs;
    const double *Px = problem->Px;
    const int *Pi = problem->Pi;
    const int *Pp = problem->Pp;
    size_t Pnnz = problem->Pnnz;

    /* the row pointers of the internal matrices are ints, so a valid CSR
       input has at most INT_MAX non-zeros */
    if (nnz > (size_t) INT_MAX || Pnnz > (size_t) INT_MAX)
    {
        return NULL;
    }

    if (Pnnz > 0 && (Px == NULL || Pi == NULL || Pp == NULL ||
                     !matrix_valid_csr_input(Pi, Pp, n, n, Pnnz)))
    {
        return NULL;
    }

    Timer timer;
    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    Matrix *A = NULL, *AT = NULL, *P = NULL;
    double *lhs_copy = NULL, *rhs_copy = NULL, *c_copy = NULL;
    int *row_sizes = NULL, *col_sizes = NULL, *P_row_sizes = NULL;
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
    Timer fold_timer;
    double fold_time = 0.0;
    size_t nnz_stripped;

    //  ---------------------------------------------------------------------------
    //   Copy data and allocate memory. For an exact specification of the
    //   workspace, see the workspace class. The problem object owns all the
    //   memory that is allocated in this block of code (and therefore frees it).
    //  ---------------------------------------------------------------------------
    size_t n_rows = m;
    size_t n_cols = n;
    lhs_copy = (double *) ps_malloc(n_rows, sizeof(double));
    rhs_copy = (double *) ps_malloc(n_rows, sizeof(double));
    c_copy = (double *) ps_malloc(n_cols, sizeof(double));
    col_tags = (ColTag *) ps_calloc(n_cols, sizeof(ColTag));
    bounds = (Bound *) ps_malloc(n_cols, sizeof(Bound));
    work = new_work(n_rows, n_cols, stgs->dton_eq);
    row_sizes = (int *) ps_malloc(n_rows, sizeof(int));
    col_sizes = (int *) ps_malloc(n_cols, sizeof(int));

    if (!lhs_copy || !rhs_copy || !c_copy || !col_tags || !bounds || !work ||
        !row_sizes || !col_sizes)
    {
        goto cleanup;
    }

    /* The pool of the fold and of the parallel sort. Without the fold, a
       problem below the sort's threshold never uses workers, so it gets none. */
    bool needs_workers = stgs->lp_fold || n_rows >= PARALLEL_SORT_THRESHOLD ||
                         n_cols >= PARALLEL_SORT_THRESHOLD;
    work->pool = thread_pool_new(
        needs_workers ? thread_pool_resolve_n_threads(stgs->n_threads) - 1 : 0);
    if (!work->pool) goto cleanup;

    memcpy(lhs_copy, lhs, n_rows * sizeof(double));
    memcpy(rhs_copy, rhs, n_rows * sizeof(double));
    memcpy(c_copy, c, n_cols * sizeof(double));

    // ---------------------------------------------------------------------------
    //  Build and validate P (NULL for an LP)
    // ---------------------------------------------------------------------------
    if (Pnnz > 0)
    {
        P = matrix_new_no_extra_space(Px, Pi, Pp, n_cols, n_cols, Pnnz);
        if (!P) goto cleanup;
        if (!matrix_is_symmetric(P)) goto cleanup;
        if (P->nnz == 0)
        {
            /* All entries were explicit zeros. */
            free_matrix(P);
            P = NULL;
        }
    }
    if (P)
    {
        P_row_sizes = (int *) ps_malloc(n_cols, sizeof(int));
        if (!P_row_sizes) goto cleanup;
        count_rows(P, P_row_sizes);
    }
    size_t Pnnz_stripped =
        P ? P->nnz : 0; /* after the explicit zeros were dropped */
    bool lp_fold = stgs->lp_fold && P == NULL; /* no fold of a QP */

    // ---------------------------------------------------------------------------
    //  Build bounds, row tags, A and AT. We first build A, then in parallel
    //  we build AT (main thread) and some other things (second thread).
    // ---------------------------------------------------------------------------
    A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);
    if (!A) goto cleanup;

    ps_thread_t thread_id;
    ParallelInitData parallel_data = {A,    work,     n_cols,   n_rows,    lbs,
                                      ubs,  lhs_copy, rhs_copy, bounds,    col_tags,
                                      NULL, NULL,     NULL,     row_sizes, !lp_fold};

    ps_thread_create(&thread_id, NULL, init_thread_func, &parallel_data);

    // Main thread: Transpose A
    AT = transpose(A, work->iwork_n_cols,
                   stgs->dton_eq ? dton_extra_memory(A->nnz) : 0);
    if (!AT)
    {
        ps_thread_join(&thread_id, NULL);
        goto cleanup;
    }

    // sync threads
    ps_thread_join(&thread_id, NULL);

    row_tags = parallel_data.row_tags;
    if (!row_tags) goto cleanup;
    locks = parallel_data.locks;
    activities = parallel_data.activities;

    nnz_stripped = A->nnz; /* after the explicit zeros were dropped */
    if (lp_fold)
    {
        FoldOptions fold_options;
        fold_options.min_reduction = LP_FOLD_MIN_REDUCTION;
        fold_options.work_budget = LP_FOLD_WORK_BUDGET;
        fold_options.max_rounds = LP_FOLD_MAX_ROUNDS;
        clock_gettime(CLOCK_MONOTONIC, &fold_timer.start);
        if (fold_problem(&A, AT, lhs_copy, rhs_copy, row_tags, bounds, col_tags,
                         c_copy, &fold_options, work->pool, &work->fold_maps, NULL))
        {
            free_matrix(AT);
            AT = transpose(A, work->iwork_n_cols,
                           stgs->dton_eq ? dton_extra_memory(A->nnz) : 0);
            if (!AT) goto cleanup;
            /* postsolve replays the records into the folded dimensions */
            work->folded_sol = new_solution(A->m, A->n);
            if (!work->folded_sol) goto cleanup;
        }
        clock_gettime(CLOCK_MONOTONIC, &fold_timer.end);
        fold_time = GET_ELAPSED_SECONDS(fold_timer);

        /* The init thread left these to us. */
        locks = new_locks(A, row_tags);
        activities = new_activities(A, col_tags, bounds);
        count_rows(A, row_sizes);
    }
    if (!locks || !activities) goto cleanup;
    count_rows(AT, col_sizes);

    // ---------------------------------------------------------------------------
    //  Initialize internal data and constraints
    // ---------------------------------------------------------------------------
    data = new_state(row_sizes, col_sizes, P_row_sizes, locks, A->m, A->n,
                     activities, work, row_tags);

    if (!data) goto cleanup;
    constraints =
        constraints_new(A, AT, lhs_copy, rhs_copy, bounds, data, row_tags, col_tags);
    if (!constraints) goto cleanup;

    // ---------------------------------------------------------------------------
    //             Allocate the actual presolver
    // ---------------------------------------------------------------------------
    presolver = (Presolver *) ps_malloc(1, sizeof(Presolver));
    obj = objective_new(c_copy, P);
    if (!presolver || !obj) goto cleanup;
    presolver->sol = NULL; /* allocated by the first postsolve call */
    presolver->stgs = stgs;
    presolver->prob = new_problem(constraints, obj);
    presolver->stats = init_stats(n_rows, n_cols, nnz);
    presolver->stats->nnz_removed_trivial = nnz - nnz_stripped; /* explicit 0's */
    presolver->stats->nnz_removed_lp_fold = nnz_stripped - A->nnz;
    presolver->stats->n_rows_folded = A->m;
    presolver->stats->n_cols_folded = A->n;
    presolver->stats->Pnnz_original = Pnnz_stripped;
    presolver->reduced_prob =
        (PresolvedProblem *) ps_calloc(1, sizeof(PresolvedProblem));
    DEBUG(run_debugger_problem(presolver->prob, false));

    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    presolver->stats->time_lp_fold = fold_time;
    presolver->stats->time_init = GET_ELAPSED_SECONDS(timer) - fold_time;

    return presolver;

cleanup:
{
    struct clean_up_scope scope = {.A = A,
                                   .AT = AT,
                                   .P = P,
                                   .lhs_copy = lhs_copy,
                                   .rhs_copy = rhs_copy,
                                   .c_copy = c_copy,
                                   .row_sizes = row_sizes,
                                   .col_sizes = col_sizes,
                                   .P_row_sizes = P_row_sizes,
                                   .row_tags = row_tags,
                                   .locks = locks,
                                   .activities = activities,
                                   .data = data,
                                   .constraints = constraints,
                                   .presolver = presolver,
                                   .obj = obj,
                                   .col_tags = col_tags,
                                   .bounds = bounds,
                                   .work = work};
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
static inline bool update_termination(size_t nnz_after_cycle,
                                      size_t nnz_before_cycle, Complexity complexity,
                                      PresolveStatus status, double max_time,
                                      Timer outer_timer)
{
    if (HAS_STATUS(status, UNBNDORINFEAS))
    {
        return true;
    }

    if (complexity == MEDIUM &&
        (double) nnz_after_cycle >= NNZ_REMOVED_CYCLE * (double) nnz_before_cycle)
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
                                           size_t nnz_before_phase,
                                           size_t nnz_after_phase)
{
    if (curr_complexity == FAST)
    {
        if ((double) nnz_after_phase < NNZ_REMOVED_FAST * (double) nnz_before_phase)
        {
            return FAST;
        }

        return MEDIUM;
    }
    else if (curr_complexity == MEDIUM)
    {
        return FAST;
    }

    // should never reach here
    assert(false);
    return FAST;
}

/* This function exhaustively removes singleton rows from the problem. It
   also eliminates empty rows and empty columns. If the problem is feasible
   and bounded, it returns UNCHANGED (even if the problem was reduced).
*/
static inline PresolveStatus run_trivial_explorers(Problem *prob,
                                                   const Settings *stgs)
{
    PresolveStatus status = UNCHANGED;

    // TODO: should we just run this once in the beginning of presolve?
    // Very important question.
    status = remove_variables_with_close_bounds(prob);
    RETURN_IF_INFEASIBLE(status);

    if (stgs->dual_fix)
    {
        //   after simple dual fix there can be new empty rows and singleton rows,
        //   and empty columns (if rows are marked as inactive due to a variable
        //   being set to inf)
        status |= remove_empty_cols(prob);
        RETURN_IF_UNBNDORINFEAS(status);
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
    assert(!has_removable_empty_cols(prob->constraints->state));
    return UNCHANGED;
}

/* This function applies the fastest reductions, which are doubleton equality
   rows, singleton columns, and simple dual fix. It returns UNCHANGED even if
   the problem was reduced (provided that the problem does not seem infeasible or
   bounded). */
static inline PresolveStatus run_fast_explorers(Problem *prob, const Settings *stgs,
                                                PresolveStats *stats)
{
    assert(prob->constraints->state->ston_rows->len == 0);
    assert(prob->constraints->state->empty_rows->len == 0);
    assert(!has_removable_empty_cols(prob->constraints->state));
    PresolveStatus status = UNCHANGED;
    Timer timer;

    if (stgs->ston_cols)
    {
        clock_gettime(CLOCK_MONOTONIC, &timer.start);
        status |= remove_ston_cols(prob);
        clock_gettime(CLOCK_MONOTONIC, &timer.end);
        stats->time_ston_cols += GET_ELAPSED_SECONDS(timer);

        // after removing singleton columns, there can be new empty columns
        // and singleton rows, but no empty rows
        assert(prob->constraints->state->empty_rows->len == 0);
        status |= run_trivial_explorers(prob, stgs);
    }

    if (stgs->dton_eq)
    {
        clock_gettime(CLOCK_MONOTONIC, &timer.start);
        status |= remove_dton_eq_rows(prob);
        clock_gettime(CLOCK_MONOTONIC, &timer.end);
        stats->time_dton_rows += GET_ELAPSED_SECONDS(timer);
        RETURN_IF_INFEASIBLE(status);
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
    assert(!has_removable_empty_cols(prob->constraints->state));
    DEBUG(verify_no_duplicates_sort(prob->constraints->state->updated_activities));
    size_t nnz_before;
    size_t *nnz = &prob->constraints->A->nnz;
    PresolveStatus status = UNCHANGED;
    Timer timer;

    if (stgs->primal_propagation)
    {
        nnz_before = *nnz;

        // the actual reduction
        clock_gettime(CLOCK_MONOTONIC, &timer.start);
        status |= check_activities(prob);
        status |= propagate_primal(prob, stgs->finite_bound_tightening);
        clock_gettime(CLOCK_MONOTONIC, &timer.end);
        stats->time_primal_propagation += GET_ELAPSED_SECONDS(timer);

        // after dom prop propagation there can be new empty and ston rows
        status |= run_trivial_explorers(prob, stgs);
        assert(!(status & REDUCED));
        RETURN_IF_NOT_UNCHANGED(status);
        stats->nnz_removed_primal_propagation += nnz_before - *nnz;
    }

    if (stgs->parallel_rows)
    {
        nnz_before = *nnz;

        // the actual reduction (after removing parallel rows there can
        // be no new empty or ston rows so no need to run trivial presolvers)
        clock_gettime(CLOCK_MONOTONIC, &timer.start);
        status |= remove_parallel_rows(prob->constraints);
        clock_gettime(CLOCK_MONOTONIC, &timer.end);
        stats->time_parallel_rows += GET_ELAPSED_SECONDS(timer);
        assert(prob->constraints->state->empty_rows->len == 0);
        assert(prob->constraints->state->ston_rows->len == 0);
        assert(!(status & REDUCED));
        RETURN_IF_NOT_UNCHANGED(status);
        stats->nnz_removed_parallel_rows += nnz_before - *nnz;
    }

    if (stgs->parallel_cols)
    {
        nnz_before = *nnz;

        // the actual reduction (after removing parallel columns there can
        // be no new empty rows or empty cols but might be new stonrows, probably
        // although that's rare but it does happen)
        clock_gettime(CLOCK_MONOTONIC, &timer.start);
        status |= remove_parallel_cols(prob);
        clock_gettime(CLOCK_MONOTONIC, &timer.end);
        stats->time_parallel_cols += GET_ELAPSED_SECONDS(timer);
        assert(prob->constraints->state->empty_rows->len == 0);
        assert(!has_removable_empty_cols(prob->constraints->state));
        status |= run_trivial_explorers(prob, stgs);
        assert(prob->constraints->state->ston_rows->len == 0);
        assert(!(status & REDUCED));
        RETURN_IF_NOT_UNCHANGED(status);
        stats->nnz_removed_parallel_cols += nnz_before - *nnz;
    }

    return status;
}

void populate_presolved_problem(Presolver *presolver)
{
    PresolvedProblem *reduced_prob = presolver->reduced_prob;
    Constraints *constraints = presolver->prob->constraints;
    Matrix *A = constraints->A;
    reduced_prob->m = A->m;
    reduced_prob->n = A->n;
    reduced_prob->nnz = A->nnz;
    reduced_prob->Ax = A->x;
    reduced_prob->Ai = A->i;
    reduced_prob->rhs = constraints->rhs;
    reduced_prob->lhs = constraints->lhs;
    reduced_prob->c = presolver->prob->obj->c;
    reduced_prob->obj_offset = presolver->prob->obj->offset;

    // create bounds arrays
    reduced_prob->lbs = (double *) malloc(A->n * sizeof(double));
    reduced_prob->ubs = (double *) malloc(A->n * sizeof(double));
    for (int i = 0; i < A->n; i++)
    {
        reduced_prob->lbs[i] = constraints->bounds[i].lb;
        reduced_prob->ubs[i] = constraints->bounds[i].ub;
    }

    // create row pointers
    reduced_prob->Ap = (int *) malloc((A->m + 1) * sizeof(int));
    for (int i = 0; i < A->m + 1; i++)
    {
        reduced_prob->Ap[i] = A->p[i].start;
    }

    /* For an LP the arrays of P stay NULL and Pnnz stays 0. */
    const Matrix *P = presolver->prob->obj->P;
    if (P == NULL)
    {
        return;
    }
    assert(P->m == A->n && P->n == A->n);
    reduced_prob->Px = P->x;
    reduced_prob->Pi = P->i;
    reduced_prob->Pnnz = P->nnz;
    reduced_prob->Pp = (int *) malloc((P->m + 1) * sizeof(int));
    for (int i = 0; i < P->m + 1; i++)
    {
        reduced_prob->Pp[i] = P->p[i].start;
    }
}

static inline void print_start_message(const PresolveStats *stats)
{
    printf("\n\t       PSLP v%s - LP/QP presolver \n\t(c) Daniel "
           "Cederberg, Stanford University, 2025\n",
           PSLP_VERSION);
    printf("Original problem:  %zu rows, %zu columns, %zu nnz, %zu nnz in P\n",
           stats->n_rows_original, stats->n_cols_original, stats->nnz_original,
           stats->Pnnz_original);
}

static inline void print_end_message(const Matrix *A, const PresolveStats *stats)
{
    printf("Presolved problem: %zu rows, %zu columns, %zu nnz, %zu nnz in P\n", A->m,
           A->n, A->nnz, stats->Pnnz_reduced);
    printf("PSLP init, fold & run time : %.3f, %.3f, %.3f seconds\n",
           stats->time_init, stats->time_lp_fold, stats->time_presolve);
}

static inline void print_infeas_or_unbnd_message(PresolveStatus status)
{
    if (status == INFEASIBLE)
    {
        printf("PSLP declares problem as infeasible.\n");
        return;
    }
    else if (status == UNBNDORINFEAS)
    {
        printf("PSLP declares problem as infeasible or unbounded.\n");
        return;
    }

    // should never reach here
    assert(false);
}

static inline PresolveStatus final_status(PresolveStats *stats)
{
    if (stats->n_cols_original == stats->n_cols_reduced &&
        stats->n_rows_original == stats->n_rows_reduced &&
        stats->nnz_original == stats->nnz_reduced)
    {
        return UNCHANGED;
    }
    else
    {
        return REDUCED;
    }
}

/* Frees the memory that only presolve uses. Postsolve keeps the recorded
   reductions and the row/column mappings. */
static void free_presolve_scratch(Problem *prob)
{
    free_matrix(prob->constraints->AT);
    prob->constraints->AT = NULL;
    free_work_scratch(prob->constraints->state->work);
}

PresolveStatus run_presolver(Presolver *presolver)
{
    Timer inner_timer, outer_timer;
    size_t nnz_before_cycle, nnz_after_cycle;
    size_t nnz_before_phase, nnz_after_phase;
    size_t nnz_before_reduction;
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

    DEBUG(run_debugger_problem(prob, false));

    // ------------------------------------------------------------------------
    // Main loop for the presolver. The logic is organized into phases and
    // cycles. A phase is a sequence of presolvers belonging to a certain
    // complexity class. The cycle resets after the medium presolvers.
    // ------------------------------------------------------------------------
    nnz_before_cycle = A->nnz;
    nnz_after_cycle = A->nnz; /* to avoid unitialized variable warning */
    while (!terminate)
    {
        // before each phase we run the trivial presolvers
        nnz_before_reduction = A->nnz;
        status = run_trivial_explorers(prob, stgs);

        // run_trivial_explorers returns something else than UNCHANGED only
        // if the problem is detected to be infeasible or unbounded
        if (status != UNCHANGED)
        {
            break;
        }

        stats->nnz_removed_trivial += nnz_before_reduction - A->nnz;
        DEBUG(run_debugger_problem(prob, false));
        nnz_before_phase = A->nnz;

        // apply presolvers belonging to a certain complexity class
        if (curr_complexity == FAST)
        {
            nnz_before_reduction = A->nnz;
            RUN_AND_TIME(run_fast_explorers, inner_timer,
                         stats->time_fast_reductions, status, prob, stgs, stats);
            stats->nnz_removed_fast += nnz_before_reduction - A->nnz;
        }
        else if (curr_complexity == MEDIUM)
        {
            RUN_AND_TIME(run_medium_explorers, inner_timer,
                         stats->time_medium_reductions, status, prob, stgs, stats);
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

    if (status != UNCHANGED)
    {
        // problem detected to be infeasible or unbounded
        if (stgs->verbose)
        {
            print_infeas_or_unbnd_message(status);
        }
        return status;
    }

    if (stgs->relax_bounds)
    {
        remove_redundant_bounds(prob->constraints);
    }

    problem_clean(prob);
    DEBUG(run_debugger_problem(prob, true));
    stats->n_rows_reduced = A->m;
    stats->n_cols_reduced = A->n;
    stats->nnz_reduced = A->nnz;
    stats->Pnnz_reduced = prob->obj->P ? prob->obj->P->nnz : 0;
    DEBUG(run_debugger_stats_consistency_check(stats));
    populate_presolved_problem(presolver);

    free_presolve_scratch(prob);

    clock_gettime(CLOCK_MONOTONIC, &outer_timer.end);
    stats->time_presolve = GET_ELAPSED_SECONDS(outer_timer);

    if (stgs->verbose)
    {
        print_end_message(A, stats);
    }

    return final_status(stats);
}

void free_presolver_reduced_problem(Presolver *presolver)
{

    // 'Ap' is non-null when 'run_presolver' has populated the reduced problem.
    // Bailing out here keeps this function safe to call before a successful
    // presolve.
    if (!presolver || !presolver->reduced_prob || !presolver->reduced_prob->Ap)
    {
        return;
    }

    PresolvedProblem *reduced_prob = presolver->reduced_prob;
    Problem *prob = presolver->prob;
    if (prob)
    {
        Constraints *constraints = prob->constraints;
        if (constraints)
        {
            free_matrix(constraints->A);
            constraints->A = NULL;
            free_matrix(constraints->AT);
            constraints->AT = NULL;
            PS_FREE(constraints->lhs);
            PS_FREE(constraints->rhs);
            PS_FREE(constraints->bounds);
            PS_FREE(constraints->row_tags);
            PS_FREE(constraints->col_tags);
        }
        if (prob->obj)
        {
            PS_FREE(prob->obj->c);
            free_matrix(prob->obj->P);
            prob->obj->P = NULL;
        }
    }

    reduced_prob->Ax = NULL;
    reduced_prob->Ai = NULL;
    reduced_prob->lhs = NULL;
    reduced_prob->rhs = NULL;
    reduced_prob->c = NULL;
    reduced_prob->Px = NULL;
    reduced_prob->Pi = NULL;
    PS_FREE(reduced_prob->Ap);
    PS_FREE(reduced_prob->Pp);
    PS_FREE(reduced_prob->lbs);
    PS_FREE(reduced_prob->ubs);
}

/* Allocates presolver->sol on the first postsolve call. Returns false (and
   prints an error) if the allocation fails. */
static bool ensure_solution(Presolver *presolver)
{
    if (presolver->sol)
    {
        return true;
    }

    presolver->sol = new_solution(presolver->stats->n_rows_original,
                                  presolver->stats->n_cols_original);
    if (!presolver->sol)
    {
        fprintf(stderr, "PSLP error: out of memory when allocating the "
                        "postsolved solution. Outputs are left untouched. \n");
        return false;
    }
    return true;
}

void postsolve(Presolver *presolver, const double *x, const double *y,
               const double *z)
{
    Timer timer;
    PresolveStats *stats = presolver->stats;
    if (!ensure_solution(presolver))
    {
        return;
    }
    Solution *sol = presolver->sol;
    State *data = presolver->prob->constraints->state;
    PostsolveInfo *postsolve_info = data->postsolve_info;
    const FoldMaps *fold_maps = data->work->fold_maps;
    /* With a fold the records live in the folded problem: replay them there,
       then lift to the original. */
    Solution *target = fold_maps ? data->work->folded_sol : sol;
    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    postsolver_update(postsolve_info, stats->n_cols_reduced, stats->n_rows_reduced,
                      data->work->mappings->cols, data->work->mappings->rows);
    postsolver_run(postsolve_info, target, x, y, z);
    if (fold_maps)
    {
        fold_postsolve(fold_maps, (int) sol->dim_y, (int) sol->dim_x, target->x,
                       target->y, target->z, sol->x, sol->y, sol->z);
    }
    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    stats->time_postsolve = GET_ELAPSED_SECONDS(timer);

    if (presolver->stgs->verbose)
    {
        printf("PSLP postsolve time: %.4f seconds\n", stats->time_postsolve);
    }
}

static void compute_primal_infeas_ray_z(double *z, const double *Ax, const int *Ai,
                                        const int *Ap, const double *y,
                                        size_t n_rows, size_t n_cols)
{
    for (size_t j = 0; j < n_cols; ++j)
    {
        z[j] = 0.0;
    }

    for (size_t i = 0; i < n_rows; ++i)
    {
        for (int p = Ap[i]; p < Ap[i + 1]; ++p)
        {
            z[Ai[p]] -= Ax[p] * y[i];
        }
    }
}

void postsolve_primal_infeas_ray(Presolver *presolver, const double *y,
                                 double *y_orig)
{
    Timer timer;
    PresolveStats *stats = presolver->stats;
    State *data = presolver->prob->constraints->state;
    PostsolveInfo *postsolve_info = data->postsolve_info;
    PresolvedProblem *reduced_prob = presolver->reduced_prob;
    size_t len_z = MAX((size_t) 1, stats->n_cols_reduced);
    double *z;
    assert(reduced_prob != NULL);
    if (data->work->fold_maps)
    {
        fprintf(stderr, "PSLP error: postsolve_primal_infeas_ray is not supported "
                        "with lp_fold yet! y_orig is left untouched. \n");
        return;
    }
    if (!reduced_prob || !reduced_prob->Ax)
    {
        // The reduced constraint matrix was released by
        // 'free_presolver_reduced_problem', so the ray cannot be recovered.
        fprintf(stderr, "PSLP warning: postsolve_primal_infeas_ray called after "
                        "free_presolver_reduced_problem! Don't do this! \n");
        return;
    }
    if (!ensure_solution(presolver))
    {
        return;
    }
    Solution *sol = presolver->sol;
    z = (double *) ps_malloc(len_z, sizeof(double));
    assert(z != NULL);

    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    compute_primal_infeas_ray_z(z, reduced_prob->Ax, reduced_prob->Ai,
                                reduced_prob->Ap, y, stats->n_rows_reduced,
                                stats->n_cols_reduced);
    postsolver_update(postsolve_info, stats->n_cols_reduced, stats->n_rows_reduced,
                      data->work->mappings->cols, data->work->mappings->rows);
    postsolver_run_primal_infeas_ray(postsolve_info, sol, y, z);
    for (int i = 0; i < sol->dim_y; ++i)
    {
        y_orig[i] = sol->y[i];
    }
    PS_FREE(z);
    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    stats->time_postsolve = GET_ELAPSED_SECONDS(timer);

    if (presolver->stgs->verbose)
    {
        printf("PSLP postsolve time: %.4f seconds\n", stats->time_postsolve);
    }
}

void postsolve_dual_infeas_ray(Presolver *presolver, const double *x, double *x_orig)
{
    Timer timer;
    PresolveStats *stats = presolver->stats;
    State *data = presolver->prob->constraints->state;
    PostsolveInfo *postsolve_info = data->postsolve_info;
    assert(presolver->reduced_prob != NULL);
    if (data->work->fold_maps)
    {
        fprintf(stderr, "PSLP error: postsolve_dual_infeas_ray is not supported "
                        "with lp_fold yet! x_orig is left untouched. \n");
        return;
    }
    if (!ensure_solution(presolver))
    {
        return;
    }
    Solution *sol = presolver->sol;

    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    postsolver_update(postsolve_info, stats->n_cols_reduced, stats->n_rows_reduced,
                      data->work->mappings->cols, data->work->mappings->rows);
    postsolver_run_dual_infeas_ray(postsolve_info, sol, x);
    for (int i = 0; i < sol->dim_x; ++i)
    {
        x_orig[i] = sol->x[i];
    }
    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    stats->time_postsolve = GET_ELAPSED_SECONDS(timer);

    if (presolver->stgs->verbose)
    {
        printf("PSLP postsolve time: %.4f seconds\n", stats->time_postsolve);
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
        PS_FREE(presolver->reduced_prob->Pp);
        PS_FREE(presolver->reduced_prob->lbs);
        PS_FREE(presolver->reduced_prob->ubs);
        PS_FREE(presolver->reduced_prob);
    }

    free_solution(presolver->sol);

    PS_FREE(presolver);
}
