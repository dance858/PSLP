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
   (lb, ub, c). Equality is bitwise. */

#include <stddef.h>

struct Matrix;
struct Bound;

typedef enum FoldReason
{
    FOLD_FOLDABLE, /* Partition found and verified. */
    FOLD_NOTHING,  /* Every line alone in its class, or too few removable entries. */
    FOLD_LIMIT,    /* Work budget, round limit or memory limit reached. */
    FOLD_FAILED /* The partition did not verify, or no memory for the workspace. */
} FoldReason;

typedef struct FoldOptions
{
    double min_reduction; /* Removable-entry bound required, as a share of nnz. */
    double work_budget;   /* Work allowed, in multiples of nnz. */
    int max_rounds;
    double max_memory; /* Bytes allowed, 0 for no limit. */
} FoldOptions;

/* Partition of the rows and columns. A class id is the lowest member index, so
   a line whose class id equals its own index is the representative of its row
   class or the keeper of its column class. */
typedef struct FoldMaps
{
    int *row_class; /* [n_rows] */
    int *col_class; /* [n_cols] */
    int n_rows_folded;
    int n_cols_folded;
    int removable_entries; /* Upper bound, rows and columns counted separately. */
} FoldMaps;

typedef struct FoldStats
{
    FoldReason reason;
    int rounds;
    double work; /* In multiples of nnz. */
    double time;
    double memory_estimate; /* Bytes allocated by the fold. */
} FoldStats;

/* Finds the partition of (A, lhs, rhs, bounds, c) and verifies it exactly. AT
   is the transpose of A. 'maps' is filled when the result is FOLD_FOLDABLE and
   left with NULL arrays otherwise. */
FoldReason fold_problem(const struct Matrix *A, const struct Matrix *AT,
                        const double *lhs, const double *rhs,
                        const struct Bound *bounds, const double *c,
                        const FoldOptions *options, FoldMaps *maps,
                        FoldStats *stats);

void fold_maps_free(FoldMaps *maps);

/* Bytes the fold allocates. */
double fold_memory_estimate(size_t n_rows, size_t n_cols, size_t nnz);

const char *fold_reason_string(FoldReason reason);

#endif // FOLD_H
