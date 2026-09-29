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

#ifndef SPARSE_LA_MATRIX_H
#define SPARSE_LA_MATRIX_H

#include <stdbool.h>

#include "debug_macros.h"

struct RowView;

typedef struct
{
    int start;
    int end;
} RowRange;

// Sparse matrix in CSR format with an explicit range per row; the rows need
// not be contiguous or in index order.
typedef struct Matrix
{
    size_t m;
    size_t n;
    size_t nnz;
    size_t n_alloc;

    int *i;
    RowRange *p;
    double *x;
} Matrix;

/* Constructs a matrix from CSR input (Ax, Ai, Ap); explicit zeros are
   dropped. */
Matrix *matrix_new(const double *Ax, const int *Ai, const int *Ap, size_t n_rows,
                   size_t n_cols, size_t nnz);

/* True if (Ai, Ap) is a valid n_rows x n_cols CSR structure with nnz entries:
   Ap[0] == 0, Ap non-decreasing, Ap[n_rows] == nnz, and every row has strictly
   increasing column indices in [0, n_cols). Works on the raw input arrays, so
   it can run before a matrix is built. */
bool matrix_valid_csr_input(const int *Ai, const int *Ap, size_t n_rows,
                            size_t n_cols, size_t nnz);
// Allocates a matrix with the given dimensions and room for nnz entries.
Matrix *matrix_alloc(size_t n_rows, size_t n_cols, size_t nnz);

/* Returns the transpose of A with room for 'tail' spare entries after its
   rows, or NULL on failure. */
Matrix *transpose(const Matrix *A, int *work_n_cols, size_t tail);

/* Like transpose(), but writes into AT's existing allocation, which must have
   the transposed dimensions and room for A->nnz entries. */
void transpose_into(const Matrix *A, Matrix *AT, int *work_n_cols);

// frees all allocated memory
void free_matrix(Matrix *A);

/* Removes 'col' from a row. Updates the length of the row, but not column
   sizes. It is assumed  that col exists in the row. */
void remove_coeff(struct RowView *row, int col);

void count_rows(const Matrix *A, int *row_sizes);

/*
When the presolving is finished we want to remove all redundant space, i.e.
the inactive rows (marked SIZE_INACTIVE_ROW in row_sizes).

It may also be beneficial to remove the space associated with inactive
rows/inactive variables in the middle of the presolving process. In this
case we must keep track of inactive rows/columns that are deleted, since
the corresponding rows must also be removed from other data structures
(eg. rowSize, stonRows, rowActivities etc).

The column indices are renumbered through col_idxs_map at the end of the
function (columns[j] = colsmap[columns[j]]). Assumes the rows are stored in
index order.
*/
void remove_extra_space(Matrix *A, const int *row_sizes, const int *col_idxs_map,
                        size_t new_n_cols);

void print_row_starts(const RowRange *row_ranges, size_t len);

#ifdef TESTING
Matrix *random_matrix_new(size_t n_rows, size_t n_cols, double density);

// replace_row_A assumes the matrix has sufficient with space to shift
// rows; otherwise it throws an assertion
void replace_row_A(Matrix *A, int row, double ratio, double *new_vals, int *cols_new,
                   int new_len);

#endif

#endif // SPARSE_LA_MATRIX_H
