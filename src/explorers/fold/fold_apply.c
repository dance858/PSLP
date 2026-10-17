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

/* Builds the folded problem into a new matrix. The entries of a representative
   row are summed by folded column in an accumulator. The line data is compacted
   only once the matrix is complete. */

#include "fold_apply.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include <stdlib.h>

/* Sums the entries of one row by folded column. stamp[col] is the folded row
   that last touched the column, so nothing is cleared between rows. */
typedef struct RowAccumulator
{
    int *stamp;    /* size: n_cols_folded, -1 at the start. */
    double *value; /* size: n_cols_folded. */
    int *touched;  /* size: longest representative row, in first-touch order. */
} RowAccumulator;

static void row_accumulator_free(RowAccumulator *acc)
{
    if (!acc)
    {
        return;
    }
    PS_FREE(acc->stamp);
    PS_FREE(acc->value);
    PS_FREE(acc->touched);
    PS_FREE(acc);
}

static RowAccumulator *row_accumulator_new(int n_cols_folded, int max_len)
{
    RowAccumulator *acc = (RowAccumulator *) ps_calloc(1, sizeof(RowAccumulator));
    if (!acc)
    {
        return NULL;
    }
    acc->stamp = (int *) ps_malloc((size_t) n_cols_folded + 1, sizeof(int));
    acc->value = (double *) ps_malloc((size_t) n_cols_folded + 1, sizeof(double));
    acc->touched = (int *) ps_malloc((size_t) max_len + 1, sizeof(int));
    if (!acc->stamp || !acc->value || !acc->touched)
    {
        row_accumulator_free(acc);
        return NULL;
    }
    for (int col = 0; col < n_cols_folded; col++)
    {
        acc->stamp[col] = -1;
    }
    return acc;
}

static int compare_ints(const void *a, const void *b)
{
    int first = *(const int *) a;
    int second = *(const int *) b;
    return (first > second) - (first < second);
}

/* Writes the folded row 'folded' of the representative 'row': its entries
   summed over each column class it touches, sorted by folded column, without
   exact-zero sums. Returns the count. */
static int fold_row(const Matrix *A, int row, int folded, const FoldMaps *maps,
                    RowAccumulator *acc, int *out_index, double *out_value)
{
    int n_touched = 0;
    for (int idx = A->p[row].start; idx < A->p[row].end; idx++)
    {
        int col = maps->folded_col[A->i[idx]];
        if (acc->stamp[col] != folded)
        {
            acc->stamp[col] = folded;
            acc->value[col] = A->x[idx];
            acc->touched[n_touched++] = col;
        }
        else
        {
            acc->value[col] += A->x[idx];
        }
    }

    /* The touched columns are nearly always ascending already. */
    for (int ii = 1; ii < n_touched; ii++)
    {
        if (acc->touched[ii] < acc->touched[ii - 1])
        {
            qsort(acc->touched, (size_t) n_touched, sizeof(int), compare_ints);
            break;
        }
    }

    int count = 0;
    for (int ii = 0; ii < n_touched; ii++)
    {
        int col = acc->touched[ii];
        if (acc->value[col] != 0.0)
        {
            out_index[count] = col;
            out_value[count] = acc->value[col];
            count++;
        }
    }
    return count;
}

/* lhs, rhs and row tags of the row representing each row class. */
static void fold_row_sides(const FoldMaps *maps, int n_rows, double *lhs,
                           double *rhs, RowTag *row_tags)
{
    int idx = 0;
    for (int row = 0; row < n_rows; row++)
    {
        if (maps->row_rep[row] == row)
        {
            lhs[idx] = lhs[row];
            rhs[idx] = rhs[row];
            row_tags[idx] = row_tags[row];
            idx++;
        }
    }
}

/* Bounds and column tags of columns representing each column class. */
static void fold_bounds(const FoldMaps *maps, int n_cols, Bound *bounds,
                        ColTag *col_tags)
{
    int idx = 0;
    for (int col = 0; col < n_cols; col++)
    {
        if (maps->col_rep[col] == col)
        {
            bounds[idx] = bounds[col];
            col_tags[idx] = col_tags[col];
            idx++;
        }
    }
}

/* Cost of the column representing each column class. */
static void fold_objective(const FoldMaps *maps, int n_cols, double *c)
{
    int idx = 0;
    for (int col = 0; col < n_cols; col++)
    {
        if (maps->col_rep[col] == col)
        {
            c[idx] = maps->col_class_size[idx] * c[col];
            idx++;
        }
    }
}

FoldReason fold_apply(const Matrix *A, double *lhs, double *rhs, RowTag *row_tags,
                      Bound *bounds, ColTag *col_tags, double *c,
                      const FoldMaps *maps, Matrix **folded)
{
    int n_rows = (int) A->m, n_cols = (int) A->n;
    size_t nnz_max = 0;
    int max_len = 0, len, new_idx = 0, nnz = 0;
    Matrix *A_new;
    RowAccumulator *acc;
    *folded = NULL;

    // ----------------------------------------------------------------------
    //   compute nnz bound on folded matrix and longest representative row
    // ----------------------------------------------------------------------
    for (int row = 0; row < n_rows; row++)
    {
        if (maps->row_rep[row] == row)
        {
            len = A->p[row].end - A->p[row].start;
            nnz_max += (size_t) len;
            max_len = MAX(max_len, len);
        }
    }

    // ----------------------------------------------------------------------
    //             allocate folded matrix and temporary storage
    // ----------------------------------------------------------------------
    A_new = matrix_alloc((size_t) maps->n_rows_folded, (size_t) maps->n_cols_folded,
                         nnz_max);
    acc = row_accumulator_new(maps->n_cols_folded, max_len);
    if (!A_new || !acc)
    {
        free_matrix(A_new);
        row_accumulator_free(acc);
        return FOLD_FAILED;
    }

    // ----------------------------------------------------------------------
    //     form the constraint matrix of the folded problem, row by row
    // ----------------------------------------------------------------------
    for (int row = 0; row < n_rows; row++)
    {
        if (maps->row_rep[row] == row)
        {
            A_new->p[new_idx].start = nnz;
            nnz +=
                fold_row(A, row, new_idx, maps, acc, A_new->i + nnz, A_new->x + nnz);
            A_new->p[new_idx].end = nnz;
            new_idx++;
        }
    }
    A_new->p[new_idx].start = nnz;
    A_new->p[new_idx].end = nnz;
    A_new->nnz = (size_t) nnz;
    row_accumulator_free(acc);

    // ----------------------------------------------------------------------
    //              fold the row sides, bounds, and objective
    // ----------------------------------------------------------------------
    fold_row_sides(maps, n_rows, lhs, rhs, row_tags);
    fold_bounds(maps, n_cols, bounds, col_tags);
    fold_objective(maps, n_cols, c);
    *folded = A_new;
    return FOLD_FOLDABLE;
}

bool fold_problem(Matrix **A, const Matrix *AT, double *lhs, double *rhs,
                  RowTag *row_tags, Bound *bounds, ColTag *col_tags, double *c,
                  const FoldOptions *options, FoldMaps **maps, FoldStats *stats)
{
    FoldStats local_stats;
    if (!stats)
    {
        stats = &local_stats;
    }
    Matrix *folded = NULL;
    FoldReason reason =
        fold_detect(*A, AT, lhs, rhs, bounds, c, options, maps, stats);
    if (reason == FOLD_FOLDABLE)
    {
        reason =
            fold_apply(*A, lhs, rhs, row_tags, bounds, col_tags, c, *maps, &folded);
    }
    bool applied = reason == FOLD_FOLDABLE;
    if (applied)
    {
        free_matrix(*A);
        *A = folded;
    }
    else
    {
        fold_maps_free(*maps);
        *maps = NULL;
    }
    return applied;
}
