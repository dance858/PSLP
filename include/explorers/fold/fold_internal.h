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

/* The workspace, the first step and the check of fold_detect, exposed for
   the unit tests. */

#include "Bounds.h"
#include "colour_refinement.h"
#include "fold.h"
#include "thread_pool.h"
#include <stdbool.h>
#include <stdint.h>

struct Matrix;

/* Everything the fold allocates, apart from the maps. */
typedef struct FoldWorkspace
{
    ThreadPool *pool; /* The caller's; it runs the refinement. */
    ColourRefinement *refine;
    uint64_t *row_colours; /* size: n_rows, the start colours of the rows. */
    uint64_t *col_colours; /* size: n_cols, the start colours of the columns. */
    int *row_classes;      /* size: n_rows, the refined classes of the rows. */
    int *col_classes;      /* size: n_cols, the refined classes of the columns. */
    int *scratch; /* size: max(n_rows, n_cols), the first line of each class. */
} FoldWorkspace;

/* NULL when an allocation fails. Freeing NULL is safe; the pool stays. */
FoldWorkspace *fold_workspace_new(const struct Matrix *A, const struct Matrix *AT,
                                  ThreadPool *pool);
void fold_workspace_free(FoldWorkspace *work);

/* Start colours from the line data: a hash of (lhs, rhs) for a row and of
   (lb, ub, c) for a column. Lines with equal data get one colour. Lines with
   different data may share one. */
void fold_start_colours(FoldWorkspace *work, int n_rows, int n_cols,
                        const double *lhs, const double *rhs, const Bound *bounds,
                        const double *c);

/* True when every line has the line data of its representative. */
bool fold_check(const FoldMaps *maps, int n_rows, int n_cols, const double *lhs,
                const double *rhs, const Bound *bounds, const double *c);

#endif // FOLD_INTERNAL_H
