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

#ifndef CORE_PRESOLVER_DEBUGGER_H
#define CORE_PRESOLVER_DEBUGGER_H

// forward declarations
struct Constraints;
struct State;
struct Matrix;
struct Bound;
struct Lock;
struct iVec;
struct PresolveStats;

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>

#include "Activity.h"
#include "Tags.h"

int ASSERT_NO_ZEROS_D(const double *x, size_t len);
int ASSERT_INCREASING_I(const int *x, size_t len);
void ASSERT_NO_ACTIVE_STON_ROWS(const struct Matrix *A, const RowTag *row_tags);
void print_int_array(const int *arr, size_t len);
void print_double_array(const double *arr, size_t len);

// runs all functions below
void run_debugger(const struct Constraints *constraints, bool finished);

void verify_problem_up_to_date(const struct Constraints *constraints);

// write a verify activity function. The max and min should be zero unless it is
// valid. void verify_activity(const Constraints *constraints, int row);
// TODO: should maybe take large bounds into account here?
void verify_activities(const struct Constraints *constraints);
void verify_activity(const ColTag *col_tags, const struct Bound *bounds,
                     Activity activity, RowTag row_tag, const int *cols,
                     const double *vals, int len);

// Verifies that the states of the activities in updated_activities are ADDED
// and that they have n_inf_min = 0 or n_inf_max = 0.
void verify_row_states(const Activity *acts, const iVec *updated_activities);

/* Coefficients in the problem data that are too large may cause
   numerical issues. This function logs a warning if any entry in
   A, lb, ub, lhs, rhs, c has magnitude larger than 10^7.
   void verifyMagnitudeOfEntries(); */

/* verifies that all empty rows have been appended to the list of
   empty rows, and that all rows in the list are empty */
void verify_empty_rows(const struct Constraints *constraints);

/* verifies that all singleton rows have been appended to the list of
   singleton rows, and that all rows in the list are singleton */
void verify_ston_rows(const struct Constraints *constraints);

/* verifies that all doubleton equality rows have been appended
   to the list of doubleton rows, and that all rows in the list
   are doubleton equality rows. */
void verify_doubleton_rows(const struct Constraints *constraints);

/* verify that all empty columns have been appended to the list of
   empty columns, and that all columns in the list are empty */
void verify_empty_cols(const struct Constraints *constraints);

/* verify that all singleton columns have been appended to the list of
   singleton columns, and that all columns in the list are singleton */
void verify_ston_cols(const struct Constraints *constraints);

// verify row and column sizes
void verify_row_and_col_sizes(const struct Constraints *constraints);

/* verifies that there are no duplicates in the list of updated activities,
   list of dton_rows and ston_rows */
// void verify_no_duplicates(const iVec *vec);
void verify_no_duplicates_lists(const struct State *data);
// void verify_no_duplicates_ptr(const int *data, int len);
void verify_no_duplicates_sort(const struct iVec *vec);
void verify_no_duplicates_sort_ptr(const int *data, int len);

/* verifies that the input is a valid CSR matrix, no explicit zeros are
   allowed, nnz are counted (compressed equal to true means that there
   is no extra space between rows)
   Currently this is static inside debugger.c */
// void verify_CSR_matrix(const Matrix *A, bool compressed);

/* verifies that A and AT are consistent */
bool verify_A_and_AT_consistency(const struct Matrix *A, const struct Matrix *AT);

/* verifies that
    1. A and AT are consistent
    2. both are valid CSR matrices
    3. A does not store the coefficients of inactive columns
    4. AT does not store the coefficients of inactive rows
*/
void verify_A_and_AT(const struct Constraints *constraints, bool compressed);

// verifies that the locks are correct
void verify_locks(const struct Matrix *AT, const struct Lock *locks,
                  const ColTag *col_tags, const RowTag *row_tags);

// verifies that the row_tags are correct
void verify_row_tags(const struct Constraints *constraints);

/* verify that there are no empty rows, empty columns, inactive rows, or
   inactive columns when finished */
void verify_empty_when_finished(const struct Constraints *constraints);

void run_debugger_stats_consistency_check(const struct PresolveStats *stats);

void print_matrix(const struct Matrix *A);

#endif // CORE_PRESOLVER_DEBUGGER_H
