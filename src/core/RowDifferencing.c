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
#include "RowDifferencing.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "Memory_wrapper.h"
#include "Numerics.h"
#include "glbopts.h"

#define ROW_DIFFERENCING_MIN_LEN 16

/* The share of the entries that must go. */
#define ROW_DIFFERENCING_MIN_GAIN 0.25

/* initial capacity of the lists of candidate pairs */
#define ROW_DIFFERENCING_INIT_PAIRS 1024

typedef struct DiffPass
{
    /* input matrix */
    const double *Ax;
    const int *Ai;
    const int *Ap;

    /* output matrix */
    Matrix *A;

    /* The candidate pairs, with rows ascending. Pair k replaces row rows[k] by
       the difference rows[k] - alphas[k] * parents[k]. After the count only
       the accepted pairs are left, and they become the result. */
    RowDifferencing *rd;

    /* diff_lens[k] is the number of entries of the difference of pair k, or -1
       if the pair was rejected. */
    int *diff_lens;
} DiffPass;

/* Checks whether row 1 can be replaced by the difference row1 - alpha * row2,
   and counts the entries the difference would have. Nothing is written.

   Each row is given by its column indices (cols1, cols2), its values (vals1,
   vals2) and its length (len1, len2). The two rows are walked in step, in
   increasing column order:
     - a column that appears in only one row survives in the difference and
       is counted (unless the stored value is an explicit zero);
     - a column that appears in both rows must cancel exactly, that is
       vals1 == alpha * vals2. It then contributes nothing.

   Returns the number of entries of the difference, or -1 if the pair must be
   rejected:
     - some shared column does not cancel exactly;
     - the difference has more than max_len entries (checked during the walk,
       so a bad pair is given up early).

   The columns of each row must be strictly increasing, as the API requires of
   the input. This is not checked. */
static int diff_count(const int *cols1, const double *vals1, int len1,
                      const int *cols2, const double *vals2, int len2, double alpha,
                      int max_len)
{
    int idx1 = 0;
    int idx2 = 0;
    int diff_len = 0;

    // ------------------------------------------------------------------------
    //          Two-pointer merge to find number of elements that would cancel
    // ------------------------------------------------------------------------
    while (idx1 < len1 && idx2 < len2)
    {
        int col1 = cols1[idx1];
        int col2 = cols2[idx2];

        if (col1 < col2)
        {
            // the column is only in row 1
            if (vals1[idx1] != 0.0)
            {
                ++diff_len;
            }
            ++idx1;
        }
        else if (col2 < col1)
        {
            // the column is only in row 2
            if (vals2[idx2] != 0.0)
            {
                ++diff_len;
            }
            ++idx2;
        }
        else
        {
            // the column is in both rows and must cancel
            if (vals1[idx1] != alpha * vals2[idx2])
            {
                return -1;
            }
            ++idx1;
            ++idx2;
        }

        // check for excessive fill-in
        if (diff_len > max_len)
        {
            return -1;
        }
    }

    // ------------------------------------------------------------------------
    //               Process the rest of the longer row.
    // ------------------------------------------------------------------------
    if (idx1 < len1)
    {
        for (; idx1 < len1; ++idx1)
        {
            if (vals1[idx1] != 0.0)
            {
                ++diff_len;
            }
        }
    }
    else
    {
        for (; idx2 < len2; ++idx2)
        {
            if (vals2[idx2] != 0.0)
            {
                ++diff_len;
            }
        }
    }

    return diff_len <= max_len ? diff_len : -1;
}

/* Writes the difference row1 - alpha * row2 to (cols_out, vals_out), sorted
   by column and without explicit zeros. The pair must have passed diff_count,
   so shared columns cancel and are skipped. The output has exactly as many
   entries as diff_count returned. */
static void diff_write(const int *cols1, const double *vals1, int len1,
                       const int *cols2, const double *vals2, int len2, double alpha,
                       int *cols_out, double *vals_out)
{
    int idx1 = 0;
    int idx2 = 0;
    int n_out = 0;

    // ------------------------------------------------------------------------
    //          Two-pointer merge that writes the entries that do not cancel
    // ------------------------------------------------------------------------
    while (idx1 < len1 && idx2 < len2)
    {
        int col1 = cols1[idx1];
        int col2 = cols2[idx2];

        if (col1 < col2)
        {
            // the column is only in row 1
            if (vals1[idx1] != 0.0)
            {
                cols_out[n_out] = col1;
                vals_out[n_out++] = vals1[idx1];
            }
            ++idx1;
        }
        else if (col2 < col1)
        {
            // the column is only in row 2
            if (vals2[idx2] != 0.0)
            {
                cols_out[n_out] = col2;
                vals_out[n_out++] = -alpha * vals2[idx2];
            }
            ++idx2;
        }
        else
        {
            // the column is in both rows and cancels
            ++idx1;
            ++idx2;
        }
    }

    // ------------------------------------------------------------------------
    //               Process the rest of the longer row.
    // ------------------------------------------------------------------------
    if (idx1 < len1)
    {
        for (; idx1 < len1; ++idx1)
        {
            if (vals1[idx1] != 0.0)
            {
                cols_out[n_out] = cols1[idx1];
                vals_out[n_out++] = vals1[idx1];
            }
        }
    }
    else
    {
        for (; idx2 < len2; ++idx2)
        {
            if (vals2[idx2] != 0.0)
            {
                cols_out[n_out] = cols2[idx2];
                vals_out[n_out++] = -alpha * vals2[idx2];
            }
        }
    }
}

/* Copies the unchanged input rows first_row, ..., end_row - 1 to the output
   matrix. No row in the range is differenced, so all of them are shifted by
   the same amount and one block copy is enough. */
static void copy_rows(const DiffPass *pass, int first_row, int end_row)
{
    if (end_row <= first_row)
    {
        return;
    }

    const int *Ap = pass->Ap;
    const Matrix *A = pass->A;
    size_t count = (size_t) (Ap[end_row] - Ap[first_row]);
    int dest = A->p[first_row].start;
    memcpy(A->i + dest, pass->Ai + Ap[first_row], count * sizeof(int));
    memcpy(A->x + dest, pass->Ax + Ap[first_row], count * sizeof(double));
}

/* A row is eligible for differencing, as a row or as a parent, if it is an
   equality with at least ROW_DIFFERENCING_MIN_LEN entries. */
static bool is_eligible_row(int len, double lhs, double rhs)
{
    return len >= ROW_DIFFERENCING_MIN_LEN && lhs == rhs && !IS_ABS_INF(rhs);
}

/* The alpha of a pair, from the leading coefficients of the two rows: +1 if
   they are equal, -1 if they are opposite, and 0 if the pair is no candidate. */
static int alpha_of_pair(double lead_row, double lead_parent)
{
    if (lead_row == lead_parent)
    {
        return 1;
    }
    if (lead_row == -lead_parent)
    {
        return -1;
    }
    return 0;
}

/* True if the new side rhs_row - alpha * rhs_parent is computed without
   rounding: either it is zero or one of the two terms is zero. */
static bool side_is_exact(double rhs_row, double rhs_parent, int alpha)
{
    return rhs_row == alpha * rhs_parent || rhs_row == 0.0 || rhs_parent == 0.0;
}

/* The screen: one pass over the linear rows that collects the candidate
   pairs. Only the first entry and the sides of each row are read.

   The parent of an eligible row is the latest earlier eligible row whose first
   entry is in the same column. The pair is a candidate if the two leading
   coefficients are equal up to sign, which gives alpha, and the new side is
   exact.

   Returns false when out of memory. */
static bool screen_candidates(DiffPass *pass, size_t n_lin_rows, size_t n_cols,
                              const double *lhs, const double *rhs)
{
    const double *Ax = pass->Ax;
    const int *Ai = pass->Ai;
    const int *Ap = pass->Ap;

    // last_row_of_col[col] - 1 is the latest eligible row whose first entry is
    // in col, or -1 if there is none. The shift by one lets the array start as
    // zeros, which calloc gives without touching the memory.
    int *last_row_of_col = (int *) ps_calloc(MAX(n_cols, 1), sizeof(int));
    if (!last_row_of_col)
    {
        return false;
    }

    for (int row = 0; row < (int) n_lin_rows; ++row)
    {
        int start = Ap[row];
        if (!is_eligible_row(Ap[row + 1] - start, lhs[row], rhs[row]))
        {
            continue;
        }

        // the row replaces its parent as the latest row of its first column
        int first_col = Ai[start];
        int parent = last_row_of_col[first_col] - 1;
        last_row_of_col[first_col] = row + 1;
        if (parent < 0)
        {
            continue;
        }

        int alpha = alpha_of_pair(Ax[start], Ax[Ap[parent]]);
        if (alpha == 0 || !side_is_exact(rhs[row], rhs[parent], alpha))
        {
            continue;
        }

        iVec_append(pass->rd->rows, row);
        iVec_append(pass->rd->parents, parent);
        iVec_append(pass->rd->alphas, alpha);
    }

    PS_FREE(last_row_of_col);
    return true;
}

/* The count: tests every candidate pair exactly and keeps the accepted pairs (rows,
   parents, alphas and diff_lens are compacted in place). A difference is accepted
   only if it is at most half as long as the row, and not empty. A row is read at
   most twice, as a row and as a parent. nnz_removed is set to the number of entries
   the accepted pairs remove.

   Returns false when out of memory. */
static bool count_and_filter(DiffPass *pass, size_t *nnz_removed)
{
    const int *Ap = pass->Ap;
    int *rows = pass->rd->rows->data;
    int *parents = pass->rd->parents->data;
    int *alphas = pass->rd->alphas->data;
    int n_pairs = (int) pass->rd->rows->len;

    pass->diff_lens = (int *) ps_malloc((size_t) n_pairs, sizeof(int));
    if (!pass->diff_lens)
    {
        return false;
    }

    int n_accepted = 0;
    *nnz_removed = 0;
    for (int k = 0; k < n_pairs; ++k)
    {
        int row = rows[k];
        int parent = parents[k];
        int row_len = Ap[row + 1] - Ap[row];
        int parent_len = Ap[parent + 1] - Ap[parent];

        int diff_len = diff_count(pass->Ai + Ap[row], pass->Ax + Ap[row], row_len,
                                  pass->Ai + Ap[parent], pass->Ax + Ap[parent],
                                  parent_len, alphas[k], row_len / 2);
        if (diff_len <= 0)
        {
            continue;
        }

        *nnz_removed += (size_t) (row_len - diff_len);
        rows[n_accepted] = row;
        parents[n_accepted] = parent;
        alphas[n_accepted] = alphas[k];
        pass->diff_lens[n_accepted] = diff_len;
        n_accepted++;
    }
    pass->rd->rows->len = (size_t) n_accepted;
    pass->rd->parents->len = (size_t) n_accepted;
    pass->rd->alphas->len = (size_t) n_accepted;
    return true;
}

/* The build: allocates the output matrix with exactly nnz - nnz_removed
   entries, sets the new row starts, and writes the rows: for each accepted
   pair the unchanged rows since the previous pair, then the difference, and
   at the end the unchanged rows after the last pair. Returns NULL when out of
   memory. */
static Matrix *build_matrix(DiffPass *pass, size_t n_rows, size_t n_cols, size_t nnz,
                            size_t nnz_removed)
{
    const int *Ap = pass->Ap;
    const int *rows = pass->rd->rows->data;
    const int *parents = pass->rd->parents->data;
    const int *alphas = pass->rd->alphas->data;
    int n_pairs = (int) pass->rd->rows->len;

    Matrix *A = matrix_alloc(n_rows, n_cols, nnz - nnz_removed);
    if (!A)
    {
        return NULL;
    }

    // the new row starts: every row moves down by the number of entries
    // removed from the rows before it
    int shift = 0;
    int k = 0;
    for (int row = 0; row < (int) n_rows; ++row)
    {
        A->p[row].start = Ap[row] - shift;
        if (k < n_pairs && rows[k] == row)
        {
            shift += Ap[row + 1] - Ap[row] - pass->diff_lens[k];
            ++k;
        }
        A->p[row].end = Ap[row + 1] - shift;
    }
    A->p[n_rows].start = A->p[n_rows].end = (int) (nnz - nnz_removed);

    pass->A = A;
    int first_unchanged = 0;
    for (k = 0; k < n_pairs; ++k)
    {
        int row = rows[k];
        int parent = parents[k];
        copy_rows(pass, first_unchanged, row);

        diff_write(pass->Ai + Ap[row], pass->Ax + Ap[row], Ap[row + 1] - Ap[row],
                   pass->Ai + Ap[parent], pass->Ax + Ap[parent],
                   Ap[parent + 1] - Ap[parent], alphas[k], A->i + A->p[row].start,
                   A->x + A->p[row].start);
        first_unchanged = row + 1;
    }
    copy_rows(pass, first_unchanged, (int) n_rows);

    return A;
}

/* Sets the sides of the differenced rows to rhs[row] - alpha * rhs[parent].
   The subtraction is exact: the screen only admits a pair if the result is
   zero or one of the two terms is zero. The pairs are visited with rows
   descending: a parent may be differenced too, and its input side must be
   read before it is replaced. */
static void update_sides(const DiffPass *pass, double *lhs, double *rhs)
{
    for (int k = (int) pass->rd->rows->len - 1; k >= 0; --k)
    {
        int row = pass->rd->rows->data[k];
        int parent = pass->rd->parents->data[k];
        double alpha = pass->rd->alphas->data[k];

        rhs[row] -= alpha * rhs[parent];
        lhs[row] = rhs[row];
    }
}

/* An empty list of pairs. Returns NULL when out of memory. */
static RowDifferencing *row_differencing_new(void)
{
    RowDifferencing *rd = (RowDifferencing *) ps_calloc(1, sizeof(RowDifferencing));
    if (!rd)
    {
        return NULL;
    }

    rd->rows = iVec_new(ROW_DIFFERENCING_INIT_PAIRS);
    rd->parents = iVec_new(ROW_DIFFERENCING_INIT_PAIRS);
    rd->alphas = iVec_new(ROW_DIFFERENCING_INIT_PAIRS);
    if (!rd->rows || !rd->parents || !rd->alphas)
    {
        row_differencing_free(rd);
        return NULL;
    }
    return rd;
}

Matrix *row_differencing_matrix_new(const double *Ax, const int *Ai, const int *Ap,
                                    size_t n_rows, size_t n_lin_rows, size_t n_cols,
                                    size_t nnz, double *lhs, double *rhs,
                                    RowDifferencing **out)
{
    *out = NULL;
    Matrix *A = NULL;
    size_t nnz_removed = 0;

    RowDifferencing *rd = row_differencing_new();
    if (!rd)
    {
        return NULL;
    }
    DiffPass pass = {.Ax = Ax, .Ai = Ai, .Ap = Ap, .rd = rd};

    if (!screen_candidates(&pass, n_lin_rows, n_cols, lhs, rhs) ||
        rd->rows->len == 0)
    {
        goto done;
    }

    if (!count_and_filter(&pass, &nnz_removed))
    {
        goto done;
    }

    // not worth it unless a large share of the entries goes (this also rejects
    // the case that no pair was accepted)
    if ((double) nnz_removed < ROW_DIFFERENCING_MIN_GAIN * (double) nnz)
    {
        goto done;
    }

    A = build_matrix(&pass, n_rows, n_cols, nnz, nnz_removed);
    if (!A)
    {
        goto done;
    }

    update_sides(&pass, lhs, rhs);

    // the unchanged rows may hold explicit zeros
    matrix_remove_explicit_zeros(A);

    rd->nnz_removed = nnz_removed;

done:
    PS_FREE(pass.diff_lens);
    if (!A)
    {
        row_differencing_free(rd);
        return NULL;
    }
    *out = rd;
    return A;
}

void row_differencing_postsolve_duals(const RowDifferencing *rd, double *y)
{
    // Pairs ascending: y[row] is still the dual of the differenced row when it
    // is read, since only a later row can have it as its parent.
    for (size_t k = 0; k < rd->rows->len; ++k)
    {
        y[rd->parents->data[k]] -= rd->alphas->data[k] * y[rd->rows->data[k]];
    }
}

void row_differencing_map_duals(const RowDifferencing *rd, double *y)
{
    // the steps of row_differencing_postsolve_duals undone in reverse order
    for (int k = (int) rd->rows->len - 1; k >= 0; --k)
    {
        y[rd->parents->data[k]] += rd->alphas->data[k] * y[rd->rows->data[k]];
    }
}

void row_differencing_free(RowDifferencing *rd)
{
    RETURN_IF_NULL(rd);
    iVec_free(rd->rows);
    iVec_free(rd->parents);
    iVec_free(rd->alphas);
    PS_FREE(rd);
}
