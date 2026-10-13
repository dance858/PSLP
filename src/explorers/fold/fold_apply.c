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
   row are sorted by folded column and each run of one folded column is summed.
   The line data is compacted only once the matrix is complete. */

#include "fold_apply.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include <stdlib.h>

typedef struct FoldEntry
{
    int col; /* Folded column. */
    double val;
} FoldEntry;

static int compare_entries(const void *a, const void *b)
{
    const FoldEntry *first = a;
    const FoldEntry *second = b;
    return (first->col > second->col) - (first->col < second->col);
}

/* Writes the folded row of the representative 'row', its entries summed over
   each column class it touches, without exact-zero sums, and returns the count. */
static int fold_row(const Matrix *A, int row, const FoldMaps *maps,
                    FoldEntry *entries, int *out_index, double *out_value)
{
    int start = A->p[row].start, n_entries = A->p[row].end - start;
    for (int k = 0; k < n_entries; k++)
    {
        entries[k].col = maps->col_fold_index[A->i[start + k]];
        entries[k].val = A->x[start + k];
    }
    qsort(entries, (size_t) n_entries, sizeof(FoldEntry), compare_entries);

    int count = 0, idx = 0;
    while (idx < n_entries)
    {
        int col = entries[idx].col;
        double sum = 0.0;
        while (idx < n_entries && entries[idx].col == col)
        {
            sum += entries[idx].val;
            idx++;
        }
        if (sum != 0.0)
        {
            out_index[count] = col;
            out_value[count] = sum;
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
        if (maps->row_class[row] == row)
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
        if (maps->col_class[col] == col)
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
        if (maps->col_class[col] == col)
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
    FoldEntry *entries;
    *folded = NULL;

    /* compute nnz bound on folded matrix and longest representative row */
    for (int row = 0; row < n_rows; row++)
    {
        if (maps->row_class[row] == row)
        {
            len = A->p[row].end - A->p[row].start;
            nnz_max += (size_t) len;
            max_len = MAX(max_len, len);
        }
    }

    /* allocate folded matrix and temporary storage */
    A_new = matrix_alloc((size_t) maps->n_rows_folded, (size_t) maps->n_cols_folded,
                         nnz_max);
    entries = ps_malloc((size_t) max_len + 1, sizeof(FoldEntry));
    if (!A_new || !entries)
    {
        free_matrix(A_new);
        PS_FREE(entries);
        return FOLD_FAILED;
    }

    /* fold the matrix */
    for (int row = 0; row < n_rows; row++)
    {
        if (maps->row_class[row] == row)
        {
            A_new->p[new_idx].start = nnz;
            nnz += fold_row(A, row, maps, entries, A_new->i + nnz, A_new->x + nnz);
            A_new->p[new_idx].end = nnz;
            new_idx++;
        }
    }
    A_new->p[new_idx].start = nnz;
    A_new->p[new_idx].end = nnz;
    A_new->nnz = (size_t) nnz;
    PS_FREE(entries);

    /* fold the row sides, bounds, and objective */
    fold_row_sides(maps, n_rows, lhs, rhs, row_tags);
    fold_bounds(maps, n_cols, bounds, col_tags);
    fold_objective(maps, n_cols, c);
    *folded = A_new;
    return FOLD_FOLDABLE;
}
