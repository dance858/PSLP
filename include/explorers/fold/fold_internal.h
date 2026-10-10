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

#ifndef FOLD_INTERNAL_H
#define FOLD_INTERNAL_H

/* The workspace and the first step of fold_problem, exposed for the unit
   tests. */

#include "Bounds.h"
#include "colour_refinement.h"
#include <stdbool.h>
#include <stdint.h>

struct Matrix;

/* The line data of a line as bit patterns, with the line for the sort. */
typedef struct LineData
{
    uint64_t data[3]; /* (lhs, rhs, 0) of a row or (lb, ub, c) of a column. */
    int line;
} LineData;

/* Everything the fold allocates, apart from the maps. */
typedef struct FoldWorkspace
{
    ColourRefinement *refine;
    LineData *lines; /* Of the larger side, for the colouring. */
    int *scratch;    /* Of the larger side, for the lowest member ids. */
} FoldWorkspace;

/* False when an allocation fails. Freeing is safe either way. */
bool fold_workspace_new(FoldWorkspace *work, const struct Matrix *A,
                        const struct Matrix *AT);
void fold_workspace_free(FoldWorkspace *work);

/* Start colours from the line data: rows with equal (lhs, rhs) get one colour,
   columns with equal (lb, ub, c) likewise. */
void fold_start_colours(FoldWorkspace *work, int n_rows, int n_cols,
                        const double *lhs, const double *rhs, const Bound *bounds,
                        const double *c);

#endif // FOLD_INTERNAL_H
