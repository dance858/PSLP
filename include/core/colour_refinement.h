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

#ifndef COLOUR_REFINEMENT_H
#define COLOUR_REFINEMENT_H

/* Colour refinement (1-dimensional Weisfeiler-Leman) of the bipartite graph of
   a sparse matrix: a node per line (a row or a column) and an edge per entry,
   labelled with the bit pattern of its value (-0.0 read as 0.0). Starting from
   colours given by the caller, the rows and the columns are split into the
   coarsest equitable partition: two rows of a class have, for every column
   class, the same multiset of labels on their edges into it, and likewise for
   columns.

   The multisets are compared exactly, so the partition is equitable whenever
   the run converges. Every round recomputes and sorts the signature of every
   line, which costs about nnz log nnz.

   The numbering of the classes is unspecified; only the lowest member ids
   identify a class across implementations. */

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

struct Matrix;

typedef struct ColourRefinement ColourRefinement;

typedef enum ColourRefinementResult
{
    COLOUR_REFINEMENT_CONVERGED,     /* Stable partition. */
    COLOUR_REFINEMENT_NOT_CONVERGED, /* Work budget or round limit reached. */
    COLOUR_REFINEMENT_EQUITABLE,     /* From verify. */
    COLOUR_REFINEMENT_NOT_EQUITABLE, /* From verify: two lines of a class differ. */
    COLOUR_REFINEMENT_ALLOCATION_FAILED
} ColourRefinementResult;

typedef struct ColourRefinementOptions
{
    double work_budget; /* Entries touched before giving up, in multiples of nnz. */
    int max_rounds;     /* A round refines the rows, then the columns. */
} ColourRefinementOptions;

typedef struct ColourRefinementStats
{
    int rounds;
    double work; /* Entries touched, in multiples of nnz. */
} ColourRefinementStats;

/* State for A and its transpose. Keeps pointers to both. NULL on failure. */
ColourRefinement *colour_refinement_new(const struct Matrix *A,
                                        const struct Matrix *AT);
void colour_refinement_free(ColourRefinement *refine);

/* Start colours, one per row and one per column, written by the caller before
   the run. Lines with equal colours and equally many entries start in one
   class. */
uint64_t *colour_refinement_row_colours(ColourRefinement *refine);
uint64_t *colour_refinement_col_colours(ColourRefinement *refine);

/* Refines to a fixpoint or a limit, once per state. The class ids are valid
   either way. */
ColourRefinementResult colour_refinement_run(ColourRefinement *refine,
                                             const ColourRefinementOptions *options,
                                             ColourRefinementStats *stats);

/* Class ids of the rows and of the columns, in [0, n_classes). */
const uint32_t *colour_refinement_row_classes(const ColourRefinement *refine);
const uint32_t *colour_refinement_col_classes(const ColourRefinement *refine);
int colour_refinement_n_row_classes(const ColourRefinement *refine);
int colour_refinement_n_col_classes(const ColourRefinement *refine);

/* Sum over the classes of both sides of (size - 1) times the entries of a
   member line: the entries that duplicate those of another line of the class. */
int colour_refinement_duplicate_entries(const ColourRefinement *refine);

/* Class ids of the rows (or the columns) as lowest member indices. 'scratch'
   has n_classes entries. */
void colour_refinement_lowest_member_ids(const ColourRefinement *refine, bool rows,
                                         int *scratch, int *ids);

/* Exact check that the partition is equitable: one more round splits nothing.
   Returns EQUITABLE, NOT_EQUITABLE or ALLOCATION_FAILED. */
ColourRefinementResult colour_refinement_verify(const ColourRefinement *refine);

#endif // COLOUR_REFINEMENT_H
