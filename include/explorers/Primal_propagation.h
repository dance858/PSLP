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

#ifndef PRIMAL_PROPAGATION_H
#define PRIMAL_PROPAGATION_H

#include <stdbool.h>

#include "PresolveStatus.h"
#include "Tags.h"

struct Constraints;
struct Bound;
struct Problem;
struct Activity;
struct ConstRowView;

typedef PresolveStatus (*BoundUpdater)(double, double *, double, int, ColTag *,
                                       struct Problem *, void *);

/* This function tries to tighten bounds of variables appearing in rows
   listed in data->updated_activities. We always tight infinite bounds, and
   if finite_bound_tightening is true we also tight finite bounds. We reject too
   small bound changes. When an implied bound is not an integer, we loosen
   the bound by the derived bound by 'BOUND_MARGINAL'. Note that this
   function implicitly takes care of forced rows.

   The function returns INFEASIBLE if infeasibility was discovered and
   otherwise UNCHANGED, even if the problem was reduced. This is a design
   choice.

   Like Gurobi, we only perform at most one round of domain propagation on
   each constraint in each presolver round.

   After this function, data->updated_activities is empty. Furthermore,
   when this function finishes, the problem is up to date, and no extra
   synchronization is needed.
*/
PresolveStatus propagate_primal(struct Problem *prob, bool finite_bound_tightening);

/* This function uses the row defined by 'row' to tighten bounds of variables
   appearing in the row. The arguments arg1 and arg2 are passed to the bound
   updaters.

    We expose this function because it is used in primal
    propagation and will be used dual propagation. In primal propagation we use
    it as

    status = bound_tightening_single_row(
                row, act, problem, bounds, col_tags,
                &finite_bound_tightening, &bound_marginal,
                (BoundUpdater)update_lb_within_propagation,
                (BoundUpdater)update_ub_within_propagation);

    and in dual propagation we'll use it differently.
: */
PresolveStatus bound_tightening_single_row(const struct ConstRowView *row,
                                           const struct Activity *act,
                                           struct Problem *problem,
                                           struct Bound *bounds, ColTag *col_tags,
                                           void *arg1, BoundUpdater updater_lb,
                                           BoundUpdater updater_ub);

void remove_redundant_bounds(struct Constraints *constraints);
#endif // PRIMAL_PROPAGATION_H