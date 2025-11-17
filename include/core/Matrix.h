/*
 * Copyright 2025 Daniel Cederberg
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

// Represents a sparse matrix stored in a modified CSR format.
// The modification is that every row includes some extra space.
typedef struct Matrix
{
    int m;
    int n;
    int nnz;
    int n_alloc;

    int *i;
    RowRange *p;
    double *x;
} Matrix;

// Constructs a new matrix with the given values but extra space.
// Assumes that (Ax, Ai, Ap) are CSR.
Matrix *matrix_new(const double *Ax, const int *Ai, const int *Ap, int n_rows,
                   int n_cols, int nnz);

Matrix *matrix_new_no_extra_space(const double *Ax, const int *Ai, const int *Ap,
                                  int n_rows, int n_cols, int nnz);

// Allocate a new matrix with the given dimensions and nnz.
// write matrix_alloc
Matrix *matrix_alloc(int n_rows, int n_cols, int nnz);

// Returns the transpose of the given matrix
Matrix *transpose(const Matrix *A, int *work_n_cols);

// Computes the number of entries allocated for a row with 'size' nnzs
int calc_memory_row(int size, int extra_row_space, double memory_ratio);

// Computes the total number of entries allocated for A
int calc_memory(int nnz, int n_rows, int extra_row_space, double memory_ratio);

// frees all allocated memory
void free_matrix(Matrix *A);

// Shifts row obtain extra space. Returns true if the shift was
// successful, and false otherwise.
bool shift_row(Matrix *A, int row, int extra_space, int max_shift);

/* Updates a coefficient of the matrix. If the coefficient does not exists,
   it inserts it. This function assumes there is space. Returns the old
   value. The function updates the row size but not any other internal
   data structures.

   This function can handle the case when val is zero. In this case,
   if the coefficient of 'col' is nonzero in the row, it will be removed.

   This function expects that if 'val' is zero, then 'col' exists in the
   row.
*/
double insert_or_update_coeff(Matrix *A, int row, int col, double val,
                              int *row_size);

/* Removes 'col' from a row. Updates the length of the row, but not column
   sizes. It is assumed  that col exists in the row. */
void remove_coeff(struct RowView *row, int col);

void count_rows(const Matrix *A, int *row_sizes);

/*
When the presolving is finished we want to remove all redundant space,
regardless of its nature.

It may also be beneficial to remove redundant space associated with
inactive rows/inactive variables in the middle of the presolving process.
In this case we must keep track of inactive rows/columns that are deleted,
since the corresponding rows must also be removed from other data
structures (eg. rowSize, stonRows, rowActivities etc).

This function will be called on both A and AT. An empty row of AT, ie. an
empty column of A, will be removed when this function is called on AT.
To have consistency between A and AT we must update the column indices of A.
This is done in the end of the function (columns[j] = colsmap[columns[j]]).
*/
void remove_extra_space(Matrix *A, const int *row_sizes, const int *col_sizes,
                        bool remove_all, int *col_idxs_map);

void print_row_starts(const RowRange *row_ranges, size_t len);

#ifdef TESTING
Matrix *random_matrix_new(int n_rows, int n_cols, double density);

// replace_row_A assumes the matrix has sufficient with space to shift
// rows; otherwise it throws an assertion
void replace_row_A(Matrix *A, int row, double ratio, double *new_vals, int *cols_new,
                   int new_len);
#endif

#endif // SPARSE_LA_MATRIX_H
