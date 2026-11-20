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
#ifndef CORE_ROW_DIFFERENCING_H
#define CORE_ROW_DIFFERENCING_H

#include <stddef.h>

#include "Matrix.h"
#include "iVec.h"

/* Row differencing. An equality row r that repeats an earlier equality row p
   and adds a few entries is replaced by r - alpha * p, with alpha = +1 or -1.
   The shared entries must cancel exactly. The primal solution is unchanged;
   the duals of the original rows are recovered by y[p] -= alpha * y[r]. */
typedef struct RowDifferencing
{
    iVec *rows;    /* the rewritten rows, ascending */
    iVec *parents; /* parents[k] < rows[k] */
    iVec *alphas;  /* +1 or -1 */
    size_t nnz_removed;
} RowDifferencing;

/* Builds the matrix of the presolver from CSR input with the rows differenced
   and their sides (lhs = rhs) updated. Only the first n_lin_rows rows are
   considered. Returns NULL and sets *out to NULL when nothing is rewritten
   (less than a quarter of the entries would go) or on an allocation
   failure: the caller then builds the matrix with matrix_new. */
Matrix *row_differencing_matrix_new(const double *Ax, const int *Ai, const int *Ap,
                                    size_t n_rows, size_t n_lin_rows, size_t n_cols,
                                    size_t nnz, double *lhs, double *rhs,
                                    RowDifferencing **out);

/* y of the differenced rows -> y of the original rows, in place. */
void row_differencing_postsolve_duals(const RowDifferencing *rd, double *y);

/* The inverse of row_differencing_postsolve_duals: maps the duals of the input rows
   to duals of the differenced rows, in place. */
void row_differencing_map_duals(const RowDifferencing *rd, double *y);

void row_differencing_free(RowDifferencing *rd);

#endif // CORE_ROW_DIFFERENCING_H
