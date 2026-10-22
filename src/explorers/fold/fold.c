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

/* Entry point of the fold: colour the rows and columns by their line data,
   refine, extract the maps. */

#include "fold.h"
#include "Bounds.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Timer.h"
#include "fold_internal.h"
#include "hash.h"
#include <assert.h>
#include <string.h>

const char *fold_reason_string(FoldReason reason)
{
    switch (reason)
    {
        case FOLD_FOLDABLE:
            return "foldable";
        case FOLD_NOTHING:
            return "nothing to fold";
        case FOLD_LIMIT:
            return "limit reached";
        case FOLD_FAILED:
            return "failed";
        case FOLD_COLLISION:
            return "hash collision";
    }
    return "unknown";
}

FoldWorkspace *fold_workspace_new(const Matrix *A, const Matrix *AT,
                                  ThreadPool *pool)
{
    FoldWorkspace *work = (FoldWorkspace *) ps_calloc(1, sizeof(FoldWorkspace));
    if (!work)
    {
        return NULL;
    }
    size_t n_lines = A->m > A->n ? A->m : A->n;
    work->pool = pool;
    work->refine = colour_refinement_new(A, AT, pool);
    work->row_colours = (uint64_t *) ps_malloc(A->m + 1, sizeof(uint64_t));
    work->col_colours = (uint64_t *) ps_malloc(A->n + 1, sizeof(uint64_t));
    work->row_classes = (int *) ps_malloc(A->m + 1, sizeof(int));
    work->col_classes = (int *) ps_malloc(A->n + 1, sizeof(int));
    work->scratch = (int *) ps_malloc(n_lines + 1, sizeof(int));
    if (!work->refine || !work->row_colours || !work->col_colours ||
        !work->row_classes || !work->col_classes || !work->scratch)
    {
        fold_workspace_free(work);
        return NULL;
    }
    return work;
}

void fold_workspace_free(FoldWorkspace *work)
{
    if (!work)
    {
        return;
    }
    colour_refinement_free(work->refine);
    PS_FREE(work->row_colours);
    PS_FREE(work->col_colours);
    PS_FREE(work->row_classes);
    PS_FREE(work->col_classes);
    PS_FREE(work->scratch);
    PS_FREE(work);
}

static FoldMaps *maps_new(int n_rows, int n_cols)
{
    FoldMaps *maps = (FoldMaps *) ps_calloc(1, sizeof(FoldMaps));
    if (!maps)
    {
        return NULL;
    }
    maps->row_rep = (int *) ps_malloc((size_t) n_rows + 1, sizeof(int));
    maps->col_rep = (int *) ps_malloc((size_t) n_cols + 1, sizeof(int));
    maps->folded_row = (int *) ps_malloc((size_t) n_rows + 1, sizeof(int));
    maps->folded_col = (int *) ps_malloc((size_t) n_cols + 1, sizeof(int));
    maps->row_class_size = (int *) ps_malloc((size_t) n_rows + 1, sizeof(int));
    maps->col_class_size = (int *) ps_malloc((size_t) n_cols + 1, sizeof(int));
    if (!maps->row_rep || !maps->col_rep || !maps->folded_row || !maps->folded_col ||
        !maps->row_class_size || !maps->col_class_size)
    {
        fold_maps_free(maps);
        return NULL;
    }
    return maps;
}

void fold_maps_free(FoldMaps *maps)
{
    if (!maps)
    {
        return;
    }
    PS_FREE(maps->row_rep);
    PS_FREE(maps->col_rep);
    PS_FREE(maps->folded_row);
    PS_FREE(maps->folded_col);
    PS_FREE(maps->row_class_size);
    PS_FREE(maps->col_class_size);
    PS_FREE(maps);
}

/* Initialise the start colours for the fold. */
void fold_start_colours(FoldWorkspace *work, int n_rows, int n_cols,
                        const double *lhs, const double *rhs, const Bound *bounds,
                        const double *c)
{
    for (int row = 0; row < n_rows; row++)
    {
        work->row_colours[row] = hash_two_values(lhs[row], rhs[row]);
    }

    for (int col = 0; col < n_cols; col++)
    {
        work->col_colours[col] =
            hash_three_values(bounds[col].lb, bounds[col].ub, c[col]);
    }
}

/* Numbers the classes 0, 1, 2, ... in the order of their representatives,
   writes the representative and folded line of every line, counts the members
   of each class and returns the number of classes. The lines are visited in
   order, so the first line seen in a class is its representative (the lowest
   line of the class). The class ids are in [0, n_lines), so 'first' has
   n_lines entries. */
static int number_classes(int n_lines, const int *class_id, int *first, int *rep,
                          int *folded, int *class_size)
{
    for (int k = 0; k < n_lines; k++)
    {
        first[k] = -1;
    }

    int n_folded = 0;
    for (int line = 0; line < n_lines; line++)
    {
        int c = class_id[line];
        if (first[c] < 0)
        {
            /* first line of a new class */
            first[c] = line;
            rep[line] = line;
            folded[line] = n_folded;
            class_size[n_folded] = 1;
            n_folded++;
        }
        else
        {
            /* join the class of its representative */
            rep[line] = first[c];
            folded[line] = folded[first[c]];
            class_size[folded[line]]++;
        }
    }
    return n_folded;
}

/* The maps of the refined partition. */
static void extract_maps(const FoldWorkspace *work, int n_rows, int n_cols,
                         FoldMaps *maps)
{
    maps->n_rows_folded =
        number_classes(n_rows, work->row_classes, work->scratch, maps->row_rep,
                       maps->folded_row, maps->row_class_size);
    maps->n_cols_folded =
        number_classes(n_cols, work->col_classes, work->scratch, maps->col_rep,
                       maps->folded_col, maps->col_class_size);
}

bool fold_check(const FoldMaps *maps, int n_rows, int n_cols, const double *lhs,
                const double *rhs, const Bound *bounds, const double *c)
{
    for (int row = 0; row < n_rows; row++)
    {
        int rep = maps->row_rep[row];
        if (lhs[row] != lhs[rep] || rhs[row] != rhs[rep])
        {
            return false;
        }
    }

    for (int col = 0; col < n_cols; col++)
    {
        int rep = maps->col_rep[col];
        if (bounds[col].lb != bounds[rep].lb || bounds[col].ub != bounds[rep].ub ||
            c[col] != c[rep])
        {
            return false;
        }
    }
    return true;
}

/* Colours, refines, gates, extracts the maps and checks them. */
static FoldReason find_fold(FoldWorkspace *work, const Matrix *A, const double *lhs,
                            const double *rhs, const Bound *bounds, const double *c,
                            const FoldOptions *options, FoldMaps *maps,
                            FoldStats *stats)
{
    int n_rows = (int) A->m, n_cols = (int) A->n;

    /* assign initial colours */
    fold_start_colours(work, n_rows, n_cols, lhs, rhs, bounds, c);

    // --------------------------------------------------------------------
    // /        find coarsest equitable partition
    // --------------------------------------------------------------------
    ColourRefinementOptions refine_options = {options->work_budget,
                                              options->max_rounds};
    ColourRefinementStats refine_stats;
    ColourRefinementResult result = colour_refinement_run(
        work->refine, work->row_colours, work->col_colours, &refine_options,
        work->row_classes, work->col_classes, &refine_stats);

    // --------------------------------------------------------------------
    //      process the result to see if folding seems worthwhile
    // --------------------------------------------------------------------
    stats->rounds = refine_stats.rounds;
    stats->work = refine_stats.work;
    if (result == COLOUR_REFINEMENT_ALLOCATION_FAILED)
    {
        return FOLD_FAILED;
    }
    if (result == COLOUR_REFINEMENT_NOT_CONVERGED)
    {
        return FOLD_LIMIT;
    }
    /* Debug builds check the partition exactly. */
    assert(colour_refinement_verify(work->refine, work->row_classes,
                                    work->col_classes) !=
           COLOUR_REFINEMENT_NOT_EQUITABLE);
    bool removes_nothing = refine_stats.duplicate_entries == 0;
    bool too_small =
        refine_stats.duplicate_entries < options->min_reduction * (double) A->nnz;
    if (removes_nothing || too_small)
    {
        return FOLD_NOTHING;
    }

    // --------------------------------------------------------------------
    //         extract the maps of the refined partition
    // --------------------------------------------------------------------
    extract_maps(work, n_rows, n_cols, maps);
    if (!fold_check(maps, n_rows, n_cols, lhs, rhs, bounds, c))
    {
        return FOLD_COLLISION;
    }
    return FOLD_FOLDABLE;
}

FoldReason fold_detect(const Matrix *A, const Matrix *AT, const double *lhs,
                       const double *rhs, const Bound *bounds, const double *c,
                       const FoldOptions *options, ThreadPool *pool, FoldMaps **maps,
                       FoldStats *stats)
{
    memset(stats, 0, sizeof(*stats));
    *maps = NULL;
    if (A->nnz == 0)
    {
        stats->reason = FOLD_NOTHING;
        return FOLD_NOTHING;
    }

    Timer timer;
    clock_gettime(CLOCK_MONOTONIC, &timer.start);

    FoldReason reason;
    FoldWorkspace *work = fold_workspace_new(A, AT, pool);
    FoldMaps *found = maps_new((int) A->m, (int) A->n);
    if (!work || !found)
    {
        reason = FOLD_FAILED;
    }
    else
    {
        reason = find_fold(work, A, lhs, rhs, bounds, c, options, found, stats);
    }
    fold_workspace_free(work);
    if (reason != FOLD_FOLDABLE)
    {
        fold_maps_free(found);
        found = NULL;
    }
    *maps = found;

    stats->reason = reason;
    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    stats->time = GET_ELAPSED_SECONDS(timer);
    return reason;
}
