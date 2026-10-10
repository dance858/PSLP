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
   refine, verify, extract the maps. */

#include "fold.h"
#include "Bounds.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Timer.h"
#include "fold_internal.h"
#include <stdlib.h>
#include <string.h>

double fold_memory_estimate(size_t n_rows, size_t n_cols, size_t nnz)
{
    double lines = (double) n_rows + (double) n_cols;
    double larger = (double) (n_rows > n_cols ? n_rows : n_cols);
    double bytes = colour_refinement_memory_estimate(n_rows, n_cols, nnz);
    bytes += (sizeof(LineData) + sizeof(int)) * larger; /* Line records, scratch. */
    bytes += sizeof(int) * lines;                       /* The maps. */
    return bytes;
}

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

/* ------------------------------------------------------------------------- */
/*  Workspace and maps                                                        */
/* ------------------------------------------------------------------------- */

bool fold_workspace_new(FoldWorkspace *work, const Matrix *A, const Matrix *AT)
{
    size_t n_lines = A->m > A->n ? A->m : A->n;
    work->refine = colour_refinement_new(A, AT);
    work->lines = (LineData *) ps_malloc(n_lines + 1, sizeof(LineData));
    work->scratch = (int *) ps_malloc(n_lines + 1, sizeof(int));
    return work->refine && work->lines && work->scratch;
}

void fold_workspace_free(FoldWorkspace *work)
{
    colour_refinement_free(work->refine);
    PS_FREE(work->lines);
    PS_FREE(work->scratch);
}

static bool maps_new(FoldMaps *maps, int n_rows, int n_cols)
{
    maps->row_class = (int *) ps_malloc((size_t) n_rows + 1, sizeof(int));
    maps->col_class = (int *) ps_malloc((size_t) n_cols + 1, sizeof(int));
    return maps->row_class && maps->col_class;
}

void fold_maps_free(FoldMaps *maps)
{
    PS_FREE(maps->row_class);
    PS_FREE(maps->col_class);
}

/* ------------------------------------------------------------------------- */
/*  Line data                                                                 */
/* ------------------------------------------------------------------------- */

/* Bit pattern of a value with -0.0 read as 0.0. */
static uint64_t value_bits(double value)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits == 0x8000000000000000ULL ? 0 : bits;
}

static int compare_line_data(const void *a, const void *b)
{
    return memcmp(((const LineData *) a)->data, ((const LineData *) b)->data,
                  sizeof(((const LineData *) a)->data));
}

/* Sorts the lines by their data and colours every run of equal data. A
   comparison sort; look into a radix sort on the three keys when the
   refinement is fast enough for this to show. */
static void colour_lines(LineData *lines, int n_lines, uint64_t *colours)
{
    qsort(lines, (size_t) n_lines, sizeof(LineData), compare_line_data);
    uint64_t colour = 0;
    for (int ii = 0; ii < n_lines; ii++)
    {
        if (ii > 0 && compare_line_data(&lines[ii - 1], &lines[ii]) != 0)
        {
            colour++;
        }
        colours[lines[ii].line] = colour;
    }
}

void fold_start_colours(FoldWorkspace *work, int n_rows, int n_cols,
                        const double *lhs, const double *rhs, const Bound *bounds,
                        const double *c)
{
    for (int row = 0; row < n_rows; row++)
    {
        LineData line = {{value_bits(lhs[row]), value_bits(rhs[row]), 0}, row};
        work->lines[row] = line;
    }
    colour_lines(work->lines, n_rows, colour_refinement_row_colours(work->refine));
    for (int col = 0; col < n_cols; col++)
    {
        LineData line = {{value_bits(bounds[col].lb), value_bits(bounds[col].ub),
                          value_bits(c[col])},
                         col};
        work->lines[col] = line;
    }
    colour_lines(work->lines, n_cols, colour_refinement_col_colours(work->refine));
}

/* ------------------------------------------------------------------------- */
/*  The fold                                                                  */
/* ------------------------------------------------------------------------- */

/* The maps of the refined partition. */
static void extract_maps(const FoldWorkspace *work, FoldMaps *maps)
{
    const ColourRefinement *refine = work->refine;
    colour_refinement_lowest_member_ids(refine, true, work->scratch,
                                        maps->row_class);
    colour_refinement_lowest_member_ids(refine, false, work->scratch,
                                        maps->col_class);
    maps->n_rows_folded = colour_refinement_n_row_classes(refine);
    maps->n_cols_folded = colour_refinement_n_col_classes(refine);
    maps->removable_entries = colour_refinement_duplicate_entries(refine);
}

/* Colours, refines, gates, verifies and extracts the maps. */
static FoldReason detect(FoldWorkspace *work, const Matrix *A, const double *lhs,
                         const double *rhs, const Bound *bounds, const double *c,
                         const FoldOptions *options, FoldMaps *maps,
                         FoldStats *stats)
{
    ColourRefinement *refine = work->refine;
    int n_rows = (int) A->m, n_cols = (int) A->n;
    fold_start_colours(work, n_rows, n_cols, lhs, rhs, bounds, c);

    ColourRefinementOptions refine_options = {options->work_budget,
                                              options->max_rounds};
    ColourRefinementStats refine_stats;
    ColourRefinementResult result =
        colour_refinement_run(refine, &refine_options, &refine_stats);
    stats->rounds = refine_stats.rounds;
    stats->work = refine_stats.work;
    if (result == COLOUR_REFINEMENT_NOT_CONVERGED)
    {
        return FOLD_LIMIT;
    }
    if ((colour_refinement_n_row_classes(refine) == n_rows &&
         colour_refinement_n_col_classes(refine) == n_cols) ||
        colour_refinement_duplicate_entries(refine) <
            options->min_reduction * (double) A->nnz)
    {
        return FOLD_NOTHING;
    }

    /* One more exact round. It cannot split after an exact refinement; it is
       the check a refinement that compares hashes needs. */
    if (!colour_refinement_verify(refine))
    {
        return FOLD_FAILED;
    }
    extract_maps(work, maps);
    return FOLD_FOLDABLE;
}

FoldReason fold_problem(const Matrix *A, const Matrix *AT, const double *lhs,
                        const double *rhs, const Bound *bounds, const double *c,
                        const FoldOptions *options, FoldMaps *maps, FoldStats *stats)
{
    Timer timer;
    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    memset(stats, 0, sizeof(*stats));
    memset(maps, 0, sizeof(*maps));
    stats->memory_estimate = fold_memory_estimate(A->m, A->n, A->nnz);

    FoldReason reason;
    FoldWorkspace work = {NULL, NULL, NULL};
    if (A->nnz == 0)
    {
        reason = FOLD_NOTHING;
    }
    else if (options->max_memory > 0 && stats->memory_estimate > options->max_memory)
    {
        reason = FOLD_LIMIT;
    }
    else if (!fold_workspace_new(&work, A, AT) ||
             !maps_new(maps, (int) A->m, (int) A->n))
    {
        reason = FOLD_FAILED;
    }
    else
    {
        reason = detect(&work, A, lhs, rhs, bounds, c, options, maps, stats);
    }
    fold_workspace_free(&work);
    if (reason != FOLD_FOLDABLE)
    {
        fold_maps_free(maps);
    }

    stats->reason = reason;
    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    stats->time = GET_ELAPSED_SECONDS(timer);
    return reason;
}
