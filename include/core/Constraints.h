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

#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <stdbool.h>

#include "Tags.h"
#include <stddef.h>

// forward declarations
struct State;
struct Mapping;
struct Matrix;
struct Bound;

typedef struct Constraints
{
    // left and right hand sides for each row
    double *lhs;
    double *rhs;
    size_t m; // number of constraints
    size_t n; // number of variables
    RowTag *row_tags;
    ColTag *col_tags;

    // CSR representations of A and A transpose. CSR representation of A
    // transpose is basically CSC representation of A.
    struct Matrix *A;
    struct Matrix *AT;

    // bounds for each variable
    struct Bound *bounds;

    // the constraints must update internal data structures
    struct State *state;
} Constraints;

/* Constructor and destructor */
Constraints *constraints_new(struct Matrix *A, struct Matrix *AT, double *lhs,
                             double *rhs, struct Bound *bounds, struct State *state,
                             RowTag *row_tags, ColTag *col_tags);
void free_constraints(Constraints *constraints);

/* Cleans the problem and removes free space from eg. rows or variables that
   have been deleted. For example, after this function, the pointers lhs and rhs
   only contain active rows, and the pointer bounds contain bounds only for
   active variables.

   After this call, 'mappings' contains the mapping from the old indices to the
   new indices for both rows and columns. For example, if mapping->rows[3] = 2,
   then row 3 in the old problem is now row 2 in the new problem. If
   mapping->rows[1] = -1, then row 1 in the old problem has been removed. */
void constraints_clean(Constraints *constraints, struct Mapping *mappings,
                       bool remove_all);

/* This function processes the list of deleted rows.
   1. It updates col sizes and the list of empty and singleton columns.
   2. The row size of a deleted row is set to SIZEINACTIVEROW.
   3. It updates locks. */
void delete_inactive_rows(Constraints *constraints);

/* This function processes the lists of fixed/substituted columns.
    1. It updates A and AT.
    1. It updates row sizes and the lists of empty/ston rows.
    2. The column size of a deleted column is set to SIZEINACTIVECOL.

    Note that it only updates interal data structures that are
    immediately linked to A and AT. In other words, it does not update
    the lhs or rhs, or activities etc.

    This function clears the lists of fixed and substituted columns to
    delete.
*/
void delete_inactive_cols_from_A_and_AT(Constraints *constraints);

#endif // CONSTRAINTS_H
