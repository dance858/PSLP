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

#ifndef SIMPLEREDUCTIONS_H
#define SIMPLEREDUCTIONS_H

#include "PSLP_status.h"
#include "Tags.h"

// forward declarations
struct Problem;
struct Constraints;
struct Bound;
struct Matrix;

/* Closely related to 'delete_inactive_cols_from_A_and_AT'. While that
   function updates A and AT, this function updates the lhs, rhs and
   activities when deleting fixed variables. This function does not
   update A or AT.

   This function does not clear the list of fixed columns to delete.
*/
void delete_fixed_cols_from_problem(struct Problem *prob);

/* Eliminates singleton rows. A variable is fixed if the row is an equality.
   If the row is an inequality the bounds are updated and the row is marked
   as inactive. Note that we only call 'set_row_to_inactive' for inequalities
   and not for equalities, since for equalities fixing the variable will
   append the column to fixedColsToDelete, and this ensures that the
   coefficient of the column is removed from A and AT. This is a very important
   implementation detail.

   When this function finishes, the problem is up to date, and no
   extra synchronization is needed. New singleton rows may have been created.
*/
PresolveStatus remove_ston_rows(struct Problem *prob);

/* Loops through the empty rows and checks feasibility with respect to the
   feasibility tolerance. An empty row is marked as inactive. Note, however,
   that 'set_row_to_inactive' isn't called, since there is no need to append
   an empty row to the list of rows that should be deleted (since there are
   no coefficients of an empty row that must be deleted).
*/
PresolveStatus remove_empty_rows(const struct Constraints *constraints);

/* Sets small coefficients of A to zero according to the rules
   in "Presolve Reductions in Mixed Integer Programming Section 3.1".

   Note that this includes fixing variables with very close lower and upper
   bounds.

   If this function is called, clean_small_coefficients_AT should never be
   called.

   Note that the remaining part of the code only calls this code if
   A is in CSR format.
*/
void clean_small_coeff_A(struct Matrix *A, const struct Bound *bounds,
                         const RowTag *row_tags, const ColTag *col_tags, double *rhs,
                         double *lhs);

/* Checks for redundant and infeasible constraints. Note that this
   function does *not* check for forcing constraints. We use
   similar identical numerics as in Presolve Reductions
   in Mixed Integer Programming Section 3.1.

   This function only processes activities inside
   internaldata->updated_activities.

   TODO (don't remove this todo, it's not really a todo):
         As Gurobi, we are very aggresive with declaring constraints
         as redundant here. If you ever get weird bugs, check this
         function.

   TODO: (very important, don't remove!) How should we treat equalities?
         We probably want to try different strategies when we are done.
         One simple strategy is to just ignore equalities. I think this
         makes sense from a numerical perspective. Try this!

   When this function finishes, the problem is up to date, and no
   extra synchronization is needed.

   Note that this function does not clear updated_activities!
*/
PresolveStatus check_activities(struct Problem *prob);

/* Loops through the empty columns and fixes them. An empty column is marked
   as fixed and the column size is set to SIZEINACTIVECOL. Note, however,
   that 'set_col_to_fixed' isn't called, since there is no need to append an
   empty column to the list of fixed columns that should be deleted (since there
   are no coefficients of an empty column that must be deleted). We DO NOT even
   have to call remove_fixed_cols().
*/
PresolveStatus remove_empty_cols(struct Problem *prob);

PresolveStatus remove_variables_with_close_bounds(struct Problem *prob);

// only used for debugging purposes
// void fix_cols_with_equal_bounds(struct Constraints *constraints);

#endif // SIMPLEREDUCTIONS_H
