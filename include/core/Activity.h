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

#ifndef ACTIVITY_H
#define ACTIVITY_H

#include <stdbool.h>

#include "Tags.h"
#include "debug_macros.h"
#include "iVec.h"

// Forward declarations.
struct Bound;
struct Matrix;
struct ConstColView;

#define INVALID_ACT_DEBUG 2 * INF // only used in debug mode

typedef enum
{
    NOT_ADDED = 0,
    ADDED = 1,
    PROPAGATED_THIS_ROUND = 2,
    PROPAGATE_NEXT_ROUND = 3,
} Status;

typedef struct Activity
{
    // minimal and maximal activities
    double min;
    double max;

    // number of variables contributing with an infinite bound to the lower
    // activity of this row
    int n_inf_min;

    // number of variables contributing with an infinite bound to the upper
    // activity of this row
    int n_inf_max;

    // if the constraint has already been propagated in the current round
    uint8_t status;

} Activity;

typedef uint8_t Altered_Activity;

enum Altered_Activity
{
    NO_RECOMPUTE = 1 << 0,
    MAX_ALTERED = 1 << 1,
    MIN_ALTERED = 1 << 2,
    // define like this for easy bit operations
    MAX_ALTERED_RECOMPUTE = (1 << 3) | MAX_ALTERED,
    MIN_ALTERED_RECOMPUTE = (1 << 4) | MIN_ALTERED
};

/* Constructor. Returns an activity struct. If the activity is valid and
   can be used for deduction of something useful, then min and max are
   computed. Otherwise, only n_inf_min and n_inf_max are updated. Indices
   of rows with potentially useful activities are appended to
   updated_activities. */
Activity *new_activities(const struct Matrix *A, const ColTag *col_tags,
                         const struct Bound *bounds);

void free_activities(Activity *activities);

/* Other type of constructor used for computing dual activities. Unlike
   new_activities, this function does not allocate any memory and it does not
   compute activities for inactive rows */
void Activities_init(const struct Matrix *A, const ColTag *col_tags,
                     const RowTag *row_tags, const struct Bound *bounds,
                     Activity *activities);

// Initializes a single activity. Used for dual propagation.
void Activity_init(Activity *act, const double *vals, const int *cols, int len,
                   const struct Bound *bounds, const ColTag *col_tags);

/* Recomputes n_inf_min and n_inf_max of a single activity. */
void recompute_n_infs(Activity *act, const double *vals, const int *cols, int len,
                      const ColTag *coltags);

/* Computes minimal and maximal activities of a row, taking column tags into
 * account*/
double compute_min_act_tags(const double *vals, const int *cols, int len,
                            const struct Bound *bounds, const ColTag *col_tags);
double compute_max_act_tags(const double *vals, const int *cols, int len,
                            const struct Bound *bounds, const ColTag *col_tags);

/* Computes minimal and maximal activities, not taking column tags into account.
   These functions should only be called when you know that the n_inf_min = 0
   or n_inf_max = 0. */
double compute_min_act_no_tags(const double *vals, const int *cols, int len,
                               const struct Bound *bounds);
double compute_max_act_no_tags(const double *vals, const int *cols, int len,
                               const struct Bound *bounds);

/* Computes minimal and maximal activities, taking ONE infinite bound into
   account. These functions should only be called when you know that
   n_inf_min = 1 or n_inf_max = 1. */
double compute_min_act_one_tag(const double *vals, const int *cols, int len,
                               const struct Bound *bounds, int col_to_avoid);
double compute_max_act_one_tag(const double *vals, const int *cols, int len,
                               const struct Bound *bounds, int col_to_avoid);

/* (vals, rows, len) represent the column with a bound change. finite_bound is
   true if the old bound was finite. lower is true if a lower bound has been
   updated. The new bound must have been updated in the bounds struct before
   calling this function.

   This function allows the new_bound to be huge. It is the responsibility of
   the caller to make sure that this only happens when it is supposed to happen
   (eg. when we update a bound because we have fixed a variable to a huge value,
    or because of a bound update due to a singleton inequality row).
*/
void update_activities_bound_change(Activity *activities, const struct Matrix *A,
                                    const struct Bound *bounds, const double *vals,
                                    const int *rows, int len, double old_bound,
                                    double new_bound, double finite_bound,
                                    bool lower,
                                    iVec *updated_activities HUGE_BOUND_PARAM);

void update_activities_fixed_col(Activity *activities, const struct Matrix *A,
                                 const struct Bound *bounds, const double *vals,
                                 const int *rows, int len, double old_ub,
                                 double old_lb, double val, bool ub_update,
                                 bool lb_update, bool is_ub_inf, bool is_lb_inf,
                                 iVec *updated_activities);

/* Updates the activity struct when a coefficient in the matrix has changed.
   Depending on the return value, we might have to recompute a max/min
   activity value from scratch outside this function */
Altered_Activity update_activity_coeff_change(Activity *act, double lb, double ub,
                                              double old_coeff, double new_coeff,
                                              ColTag cTag);

/* Updates the activity struct when a bound changes. */
Altered_Activity update_act_bound_change(Activity *act, double coeff,
                                         double old_bound, double new_bound,
                                         bool finite_bound, bool lower,
                                         const double *vals, const int *cols,
                                         int len, const struct Bound *bounds);

void remove_finite_lb_from_activities(const struct ConstColView *col, Activity *acts,
                                      double bound);

void remove_finite_ub_from_activities(const struct ConstColView *col, Activity *acts,
                                      double bound);
#endif // ACTIVITY_H