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

#ifndef PARALLEL_ROWS_H
#define PARALLEL_ROWS_H

#include "PresolveStatus.h"
#include "Tags.h"

// forward declarations
struct Constraints;
struct Matrix;
struct iVec;

// this is static and inlined in .c file, but we expose it for testing
#ifdef TESTING
void compute_supp_and_coeff_hash(const struct Matrix *A, const RowTag *rowtags,
                                 int *sparsity_IDs, int *coeff_hashes,
                                 RowTag INACTIVE_TAG);
#endif

/* This function finds all parallel rows among the active rows of a matrix A.
   If you want to find all parallel columns among the active columns, call
   it on AT. OBS! When calling this function, it is assumed that coefficients
   corresponding to inactive columns have been removed from A.

   The function works by first hashing the sparsity patterns. It then hashes
   the coefficients of rows with the same sparsity pattern? Finally, ....

   After completiton:
   * 'group_starts' is a vector of size 'n_groups', where 'parallel_rows[i]' is
   the index of the first row in group 'i'.
   * The indices of the rows in group 'i' are stored in 'parallel_rows' starting
   at 'group_starts[i]' and ending at 'group_starts[i+1] - 1'.

    Memory requirements:
    * 'parallel_rows' must be of size A->m. (It is used as internal workspace)
    * 'sparsity_IDs' must be of size A->m.
    * 'coeff_hashes' must be of size A->m.
*/
void find_parallel_rows(const struct Matrix *A, const RowTag *r_Tags,
                        struct iVec *group_starts, int *parallel_rows,
                        int *sparsity_IDs, int *coeff_hashes, RowTag INACTIVE_TAG);

/* This function finds parallel rows and deletes the redundant ones, or declares
   the problem as infeasible if the parallel rows are not consistent. It returns
   INFEASIBLE or UNCHANGED. Even if the problem is reduced, UNCHANGED is
   returned.

   When this function finishes, the problem is up to date, and no extra
   synchronization is needed.
*/
PresolveStatus remove_parallel_rows(struct Constraints *constraints);

#endif // PARALLEL_ROWS_H