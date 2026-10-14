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
#include <stdlib.h>
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
    }
    return "unknown";
}

FoldWorkspace *fold_workspace_new(const Matrix *A, const Matrix *AT)
{
    FoldWorkspace *work = (FoldWorkspace *) ps_calloc(1, sizeof(FoldWorkspace));
    if (!work)
    {
        return NULL;
    }
    size_t n_lines = A->m > A->n ? A->m : A->n;
    work->refine = colour_refinement_new(A, AT);
    work->rows = (RowData *) ps_malloc(A->m + 1, sizeof(RowData));
    work->cols = (ColData *) ps_malloc(A->n + 1, sizeof(ColData));
    work->row_colours = (int *) ps_malloc(A->m + 1, sizeof(int));
    work->col_colours = (int *) ps_malloc(A->n + 1, sizeof(int));
    work->row_classes = (int *) ps_malloc(A->m + 1, sizeof(int));
    work->col_classes = (int *) ps_malloc(A->n + 1, sizeof(int));
    work->scratch = (int *) ps_malloc(n_lines + 1, sizeof(int));
    if (!work->refine || !work->rows || !work->cols || !work->row_colours ||
        !work->col_colours || !work->row_classes || !work->col_classes ||
        !work->scratch)
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
    PS_FREE(work->rows);
    PS_FREE(work->cols);
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

/* -1, 0 or 1 as a is below, equal to or above b; -0.0 equals 0.0. */
static int compare_doubles(double a, double b)
{
    return (a > b) - (a < b);
}

/* Orders rows by lhs, then by rhs. Returns 0 only when both are equal. */
static int compare_row_data(const void *a, const void *b)
{
    const RowData *row_a = (const RowData *) a;
    const RowData *row_b = (const RowData *) b;

    int cmp_lhs = compare_doubles(row_a->lhs, row_b->lhs);
    if (cmp_lhs != 0)
    {
        return cmp_lhs;
    }

    return compare_doubles(row_a->rhs, row_b->rhs);
}

/* Orders columns by lb, then ub, then c. Returns 0 only when all are equal. */
static int compare_col_data(const void *a, const void *b)
{
    const ColData *col_a = (const ColData *) a;
    const ColData *col_b = (const ColData *) b;

    int cmp_lb = compare_doubles(col_a->lb, col_b->lb);
    if (cmp_lb != 0)
    {
        return cmp_lb;
    }

    int cmp_ub = compare_doubles(col_a->ub, col_b->ub);
    if (cmp_ub != 0)
    {
        return cmp_ub;
    }

    return compare_doubles(col_a->c, col_b->c);
}

/* Sort rows so that rows with equal (lhs, rhs) get the same colour. */
static void colour_rows(RowData *rows, int n_rows, int *colours)
{
    qsort(rows, (size_t) n_rows, sizeof(RowData), compare_row_data);
    int colour = 0;
    for (int ii = 0; ii < n_rows; ii++)
    {
        if (ii > 0 && compare_row_data(&rows[ii - 1], &rows[ii]) != 0)
        {
            colour++;
        }
        colours[rows[ii].row] = colour;
    }
}

/* Sort columns so that columns with equal (lb, ub, c) get the same colour. */
static void colour_cols(ColData *cols, int n_cols, int *colours)
{
    qsort(cols, (size_t) n_cols, sizeof(ColData), compare_col_data);
    int colour = 0;
    for (int ii = 0; ii < n_cols; ii++)
    {
        if (ii > 0 && compare_col_data(&cols[ii - 1], &cols[ii]) != 0)
        {
            colour++;
        }
        colours[cols[ii].col] = colour;
    }
}

/* Initialise the start colours for the fold. */
void fold_start_colours(FoldWorkspace *work, int n_rows, int n_cols,
                        const double *lhs, const double *rhs, const Bound *bounds,
                        const double *c)
{

    /* assign initial row colours */
    for (int row = 0; row < n_rows; row++)
    {
        RowData data = {lhs[row], rhs[row], row};
        work->rows[row] = data;
    }
    colour_rows(work->rows, n_rows, work->row_colours);

    /* assign initial column colours */
    for (int col = 0; col < n_cols; col++)
    {
        ColData data = {bounds[col].lb, bounds[col].ub, c[col], col};
        work->cols[col] = data;
    }
    colour_cols(work->cols, n_cols, work->col_colours);
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

/* Colours, refines, gates and extracts the maps. */
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
    if (result == COLOUR_REFINEMENT_NOT_CONVERGED)
    {
        return FOLD_LIMIT;
    }
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
    return FOLD_FOLDABLE;
}

FoldReason fold_detect(const Matrix *A, const Matrix *AT, const double *lhs,
                       const double *rhs, const Bound *bounds, const double *c,
                       const FoldOptions *options, FoldMaps **maps, FoldStats *stats)
{
    Timer timer;
    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    memset(stats, 0, sizeof(*stats));

    FoldReason reason;
    FoldWorkspace *work = fold_workspace_new(A, AT);
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
