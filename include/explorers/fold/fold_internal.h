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

/* The workspace and the first step of fold_detect, exposed for the unit
   tests. */

#include "Bounds.h"
#include "colour_refinement.h"
#include <stdbool.h>

struct Matrix;

/* The line data of a row, with the row for the sort. */
typedef struct RowData
{
    double lhs;
    double rhs;
    int row;
} RowData;

/* The line data of a column, with the column for the sort. */
typedef struct ColData
{
    double lb;
    double ub;
    double c;
    int col;
} ColData;

/* Everything the fold allocates, apart from the maps. */
typedef struct FoldWorkspace
{
    ColourRefinement *refine;
    RowData *rows;    /* size: n_rows, for the start colours. */
    ColData *cols;    /* size: n_cols, for the start colours. */
    int *row_colours; /* size: n_rows, the start colours of the rows. */
    int *col_colours; /* size: n_cols, the start colours of the columns. */
    int *row_classes; /* size: n_rows, the refined classes of the rows. */
    int *col_classes; /* size: n_cols, the refined classes of the columns. */
    int *scratch;     /* size: max(n_rows, n_cols), the first line of each class. */
} FoldWorkspace;

/* NULL when an allocation fails. Freeing NULL is safe. */
FoldWorkspace *fold_workspace_new(const struct Matrix *A, const struct Matrix *AT);
void fold_workspace_free(FoldWorkspace *work);

/* Start colours from the line data: rows with equal (lhs, rhs) get one colour,
   columns with equal (lb, ub, c) likewise. */
void fold_start_colours(FoldWorkspace *work, int n_rows, int n_cols,
                        const double *lhs, const double *rhs, const Bound *bounds,
                        const double *c);

#endif // FOLD_INTERNAL_H
