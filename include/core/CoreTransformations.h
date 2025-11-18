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

#ifndef CORE_TRANSFORMATIONS_H
#define CORE_TRANSFORMATIONS_H

#include "PSLP_status.h"
#include "Postsolver.h"
#include "Tags.h"
#include "debug_macros.h"
#include "glbopts.h"
#include "iVec.h"
#include <assert.h>

struct Constraints;
struct Lock;
struct Activity;
struct RowView;
struct ConstColView;
struct ConstRowView;
struct Matrix;
struct Bound;

/* Checks if fixing a column to 'val' is feasible with respect to its bounds.
   If yes, it fixes it and updates the activities. When a variable is fixed we
   set both its lb and ub to its value. */
PresolveStatus fix_col(struct Constraints *constraints, int col, double val,
                       double ck);

/* Fixes a variable to infinity and marks the constraints it appears in as
   inactive. */
void fix_col_to_negative_inf(struct Constraints *constraints, int col);
void fix_col_to_positive_inf(struct Constraints *constraints, int col);

// These should not be called with large values of new_lb or new_ub.
// 'row' is the row that causes the bound change.
PresolveStatus update_lb(struct Constraints *constraints, int col, double new_lb,
                         int row HUGE_BOUND_PARAM);
PresolveStatus update_ub(struct Constraints *constraints, int col, double new_ub,
                         int row HUGE_BOUND_PARAM);

/* Updates the rhs/lhs of a row. These functions check for infeasibility
   and also if the row becomes an equality constraint after the bound change.
   The locks are updated. It also appends to doubleton equality rows
   if a doubleton inequality row becomes an equality row.
*/
PresolveStatus change_rhs_of_ineq(struct RowView *row, double new_rhs,
                                  struct Lock *locks, iVec *dton_rows,
                                  struct PostsolveInfo *postsolve_info,
                                  int other_row_idx, double ratio);

PresolveStatus change_lhs_of_ineq(struct RowView *row, double new_lhs,
                                  struct Lock *locks, iVec *dton_rows,
                                  struct PostsolveInfo *postsolve_info,
                                  int other_row_idx, double ratio);

/* Uses activities of a row to deduce if a variable is implied free. */
bool implied_free_from_above_by_row(double Aik, const struct ConstRowView *row,
                                    double lb, double ub, ColTag col_tag,
                                    const struct Activity *act,
                                    const ColTag *col_tags,
                                    const struct Bound *bounds);
bool implied_free_from_below_by_row(double Aik, const struct ConstRowView *row,
                                    double lb, double ub, ColTag col_tag,
                                    const struct Activity *act,
                                    const ColTag *col_tags,
                                    const struct Bound *bounds);

/* Checks if the upper or lower bound of a column is implied free by any of the
   constraints it appears in. */
bool is_ub_implied_free(const struct ConstColView *col, struct Activity *activities,
                        const double *lhs, const double *rhs, const RowTag *row_tags,
                        const struct Matrix *A, const ColTag *col_tags,
                        const struct Bound *bounds);
bool is_lb_implied_free(const struct ConstColView *col, struct Activity *activities,
                        const double *lhs, const double *rhs, const RowTag *row_tags,
                        const struct Matrix *A, const ColTag *col_tags,
                        const struct Bound *bounds);

static inline void set_col_to_fixed(int col, ColTag *col_tag,
                                    iVec *fixed_cols_to_delete)
{
    assert(!HAS_TAG(*col_tag, C_TAG_INACTIVE));
    UPDATE_TAG(*col_tag, C_TAG_FIXED);
    iVec_append(fixed_cols_to_delete, col);
}

static inline void set_row_to_inactive(int row, RowTag *row_tag,
                                       iVec *rows_to_delete, PostsolveInfo *info,
                                       double val)
{
    assert(!HAS_TAG(*row_tag, R_TAG_INACTIVE));
    iVec_append(rows_to_delete, row);

    // very important that we use UPDATE_TAG here instead of RESET_TAG,
    // since the old infinity flags are needed for updating the locks
    // correctly when the row is deleted.
    UPDATE_TAG(*row_tag, R_TAG_INACTIVE);

    save_retrieval_deleted_row(info, row, val);
}

static inline void set_col_to_substituted(int col, ColTag *col_tag,
                                          iVec *substituted_cols_to_delete)
{
    assert(!HAS_TAG(*col_tag, C_TAG_INACTIVE));
    UPDATE_TAG(*col_tag, C_TAG_SUBSTITUTED);
    iVec_append(substituted_cols_to_delete, col);
}

#endif // CORE_TRANSFORMATIONS_H
