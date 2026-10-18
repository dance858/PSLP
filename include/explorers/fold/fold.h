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

#ifndef FOLD_H
#define FOLD_H

/* LP folding. Rows whose entries agree, class by class of the columns, and
   columns whose entries agree, class by class of the rows, form the coarsest
   equitable partition that refines the partition by line data (lhs, rhs) and
   (lb, ub, c). Values are equal when they compare equal as doubles; matrix
   entries are compared by their bits, with -0.0 read as 0.0. Indices and line
   counts are int, like the matrix indices. */

struct Matrix;
struct ThreadPool;
struct Bound;

typedef enum FoldReason
{
    FOLD_FOLDABLE, /* Partition found. */
    FOLD_NOTHING,  /* No entries to remove, or fewer than min_reduction asks for. */
    FOLD_LIMIT,    /* Work budget or round limit reached. */
    FOLD_FAILED,   /* Out of memory. */
    FOLD_COLLISION /* Lines with different data share a start colour. */
} FoldReason;

typedef struct FoldOptions
{
    double min_reduction; /* Removable-entry bound required, as a share of nnz. */
    double work_budget;   /* Work allowed, in multiples of nnz. */
    int max_rounds;
} FoldOptions;

/*  FoldMaps is used to keep track of the row and column classes.

    * row_rep[i]:        idx of the row that represents the class of row i.
    * col_rep[j]:        idx of the col that represents the class of col j.
    * folded_row[i]:     idx of the folded row for the class of row i.
    * folded_col[j]:     idx of the folded col for the class of col j.
    * row_class_size[k]: number of rows in the class of folded row k.
    * col_class_size[l]: number of cols in the class of folded col l.

    The representative of a class is its lowest-indexed member, so row i is a
    representative exactly when row_rep[i] == i (columns likewise).
*/
typedef struct FoldMaps
{
    int *row_rep;        /* size: n_rows */
    int *col_rep;        /* size: n_cols */
    int *folded_row;     /* size: n_rows */
    int *folded_col;     /* size: n_cols */
    int *row_class_size; /* size: n_rows_folded */
    int *col_class_size; /* size: n_cols_folded */
    int n_rows_folded;
    int n_cols_folded;
} FoldMaps;

typedef struct FoldStats
{
    FoldReason reason;
    int rounds;
    double work; /* In multiples of nnz. */
    double time;
} FoldStats;

/* Finds the partition of (A, lhs, rhs, bounds, c). AT is the transpose of A
   and 'pool' runs the refinement. *maps is set to the maps when the result is
   FOLD_FOLDABLE and to NULL otherwise. */
FoldReason fold_detect(const struct Matrix *A, const struct Matrix *AT,
                       const double *lhs, const double *rhs,
                       const struct Bound *bounds, const double *c,
                       const FoldOptions *options, struct ThreadPool *pool,
                       FoldMaps **maps, FoldStats *stats);

/* Freeing NULL is safe. */
void fold_maps_free(FoldMaps *maps);

const char *fold_reason_string(FoldReason reason);

#endif // FOLD_H
