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
#include "Primal_propagation.h"
#include "Activity.h"
#include "Bounds.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "Matrix.h"
#include "Numerics.h"
#include "Problem.h"
#include "RowColViews.h"
#include "SimpleReductions.h"
#include "State.h"
#include "Workspace.h"
#include "utils.h"

static PresolveStatus update_lb_within_propagation(double new_lb, double *lb,
                                                   double ub, int col, ColTag *cTag,
                                                   Problem *prob,
                                                   bool *finite_bound_tightening)
{
    assert(ABS(new_lb) != INF && !HAS_TAG(*cTag, C_TAG_INACTIVE));
    assert((HAS_TAG(*cTag, C_TAG_LB_INF) || new_lb >= *lb) && !IS_HUGE(new_lb));

    bool is_ub_inf = HAS_TAG(*cTag, C_TAG_UB_INF);
    bool is_lb_inf = HAS_TAG(*cTag, C_TAG_LB_INF);

    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT;
    const double *vals = AT->x + AT->p[col].start;
    const int *rows = AT->i + AT->p[col].start;
    int len = AT->p[col].end - AT->p[col].start;

    // -----------------------------------------------------------------------
    // If ub is finite we compare it to the new lower bound (we check for
    // infeasibility and if the lower bound is so close to the ub that we
    // can fix the variable).
    // -----------------------------------------------------------------------
    if (!is_ub_inf)
    {
        assert(ABS(ub) != -INF);

        // check for infeasibility
        if (new_lb >= ub + FEAS_TOL)
        {
            return INFEASIBLE;
        }

        // check if variable can be fixed
        if (new_lb >= ub || (ub - new_lb) * get_max_abs(vals, len) <= FEAS_TOL)
        {
            save_retrieval_bound_change_no_row(constraints->state->postsolve_info,
                                               col, new_lb, ub, false);

            // doing this later would make the postsolve more elegant
            fix_col(constraints, col, ub, prob->obj->c[col]);
            return REDUCED;
        }
    }

    // -----------------------------------------------------------------------
    //  We always tighten infinite bounds, and like Gurobi we reject too small
    //  bound changes.
    // -----------------------------------------------------------------------
    if (is_lb_inf || (*finite_bound_tightening && (new_lb - *lb > FEAS_TOL * 1e4) &&
                      (new_lb - *lb > 1e-2 * ABS(*lb))))
    {
        // very important to scale bound marginal with ABS
        if (!IS_INTEGRAL(new_lb))
        {
            new_lb -= BOUND_MARGINAL * ABS(new_lb);
        }

        REMOVE_TAG(*cTag, C_TAG_LB_INF);

        double old_bound = *lb;
        *lb = new_lb;

        update_activities_bound_change(
            constraints->state->activities, constraints->A, constraints->bounds,
            vals, rows, len, old_bound, new_lb, !is_lb_inf, true,
            constraints->state->updated_activities HUGE_BOUND_IS_NOT_OK);

        save_retrieval_bound_change_no_row(constraints->state->postsolve_info, col,
                                           new_lb, ub, false);

        return REDUCED;
    }

    return UNCHANGED;
}

static PresolveStatus update_ub_within_propagation(double new_ub, double *ub,
                                                   double lb, int col, ColTag *cTag,
                                                   Problem *prob,
                                                   bool *finite_bound_tightening)
{
    assert(new_ub != -INF && new_ub != INF && !HAS_TAG(*cTag, C_TAG_INACTIVE));
    assert((HAS_TAG(*cTag, C_TAG_UB_INF) || new_ub <= *ub) && !IS_HUGE(new_ub));

    bool is_lb_inf = HAS_TAG(*cTag, C_TAG_LB_INF);
    bool is_ub_inf = HAS_TAG(*cTag, C_TAG_UB_INF);

    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT;
    const double *vals = AT->x + AT->p[col].start;
    const int *rows = AT->i + AT->p[col].start;
    int len = AT->p[col].end - AT->p[col].start;

    // -----------------------------------------------------------------------
    // If lb is finite we compare it to the new upper bound (we check for
    // infeasibility and if the upper bound is so close to the lb that we
    // can fix the variable).
    // -----------------------------------------------------------------------
    if (!is_lb_inf)
    {
        assert(ABS(lb) != -INF);

        if (new_ub <= lb - FEAS_TOL)
        {
            return INFEASIBLE;
        }

        // check if variable can be fixed
        if (new_ub <= lb || (new_ub - lb) * get_max_abs(vals, len) <= FEAS_TOL)
        {
            save_retrieval_bound_change_no_row(constraints->state->postsolve_info,
                                               col, new_ub, lb, true);

            /// doing this later would make the postsolve more elegant
            fix_col(constraints, col, lb, prob->obj->c[col]);
            return REDUCED;
        }
    }

    // -----------------------------------------------------------------------
    //  We always tighten infinite bounds, and like Gurobi we reject too small
    //  bound changes.
    // -----------------------------------------------------------------------
    if (is_ub_inf || (*finite_bound_tightening && (*ub - new_ub > FEAS_TOL * 1e4) &&
                      (*ub - new_ub > 1e-2 * ABS(*ub))))
    {
        // very important to scale bound marginal with ABS
        if (!IS_INTEGRAL(new_ub))
        {
            new_ub += BOUND_MARGINAL * ABS(new_ub);
        }

        REMOVE_TAG(*cTag, C_TAG_UB_INF);

        double old_bound = *ub;
        *ub = new_ub;

        update_activities_bound_change(
            constraints->state->activities, constraints->A, constraints->bounds,
            vals, rows, len, old_bound, new_ub, !is_ub_inf, false,
            constraints->state->updated_activities HUGE_BOUND_IS_NOT_OK);

        save_retrieval_bound_change_no_row(constraints->state->postsolve_info, col,
                                           new_ub, lb, true);

        return REDUCED;
    }

    return UNCHANGED;
}

/* Bound tightening based on the right-hand side of the
   constraint.

    Bound tightening can occur in several scenarios:
 1) the RHS is finite and the min activity is valid (ie. when
    nInfMin == 0). In this case every variable in the constraint
    is a candidate for bound tightening.
 2) the RHS is finite and nInfMin == 1. In this case only the variable
    contributing with the infinite bound is a candidate for bound
    tightening.
 3) If the RHS is infinite, and the max activity is valid, and
    there is one variable contributing with an infinite bound to the
    min activity, then we may be able to tight the finite bound of
    the variable contributing the infinite bound.

    This function does not directly update act, but it may be updated
    by other functions called by this function.
*/
static inline PresolveStatus
bound_tightening_single_row_rhs(const ConstRowView *row, const Activity *act,
                                Problem *problem, Bound *bounds, ColTag *col_tags,
                                void *arg1, BoundUpdater updater_lb,
                                BoundUpdater updater_ub, int *num_of_bounds_changed)
{
    int j = 1;
    int k = -1;
    double Aik = INF;
    double implied_ub, implied_lb, implied_bound, min_act;
    PresolveStatus status = UNCHANGED;
    PresolveStatus temp_status = UNCHANGED;
    bool is_rhs_inf = HAS_TAG(*row->tag, R_TAG_RHS_INF);
    bool is_lb_inf = false;
    bool is_ub_inf = false;
    bool is_col_inactive;

    // ------------------------------------------------------------------
    // Case 1: assume RHS is finite and the min activity is valid. Then
    //         the following bounds are valid for every xk in row i:
    //         * if aik > 0: xk <= lb[k] + (rhs[i] - minAct) / aik
    //         * if aik < 0: xk >= ub[k] + (rhs[i] - minAct) / aik.
    // ------------------------------------------------------------------
    if (!is_rhs_inf && act->n_inf_min == 0)
    {
        assert(ABS(*row->rhs) != INF && ABS(act->min) != INF);

        for (j = 0; j < *row->len; ++j)
        {
            k = row->cols[j];

            // In primal propagation this might happen: suppose row 1 is forcing
            // and fixes x1. Assume that x1 appears in row 2. When we loop
            // through row 2
            //  we want to skip x1.
            if (HAS_TAG(col_tags[k], C_TAG_INACTIVE))
            {
                continue;
            }

            Aik = row->vals[j];

            // --------------------------------------------------------------------
            // if Aik > 0 we might be able to tighten the upper bound if
            // implied bound is tighter (and not too large)
            // --------------------------------------------------------------------
            if (Aik > 0)
            {
                assert(!IS_ABS_INF(bounds[k].lb) &&
                       !HAS_TAG(col_tags[k], C_TAG_LB_INF));
                implied_ub = bounds[k].lb + (*row->rhs - act->min) / Aik;
                is_ub_inf = HAS_TAG(col_tags[k], C_TAG_UB_INF);

                if ((is_ub_inf || implied_ub < bounds[k].ub) && !IS_HUGE(implied_ub))
                {
                    temp_status =
                        updater_ub(implied_ub, &(bounds[k].ub), bounds[k].lb, k,
                                   col_tags + k, problem, arg1);

                    if (temp_status == REDUCED)
                    {
                        (*num_of_bounds_changed)++;
                    }

                    status |= temp_status;
                }
            }
            // ---------------------------------------------------------------------
            // if Aik < 0 we might be able to tighten the lower bound if
            // implied bound is tighter (and not too large)
            // ---------------------------------------------------------------------
            else
            {
                assert(Aik < 0 && ABS(bounds[k].ub) != INF);
                assert(!HAS_TAG(col_tags[k], C_TAG_UB_INF));

                implied_lb = bounds[k].ub + (*row->rhs - act->min) / Aik;
                is_lb_inf = HAS_TAG(col_tags[k], C_TAG_LB_INF);

                if ((is_lb_inf || implied_lb > bounds[k].lb) && !IS_HUGE(implied_lb))
                {
                    temp_status =
                        updater_lb(implied_lb, &(bounds[k].lb), bounds[k].ub, k,
                                   col_tags + k, problem, arg1);

                    if (temp_status == REDUCED)
                    {
                        (*num_of_bounds_changed)++;
                    }

                    status |= temp_status;
                }
            }
        }
    }
    // -----------------------------------------------------------------------------
    // Case 2: assume RHS is finite and nInfMin = 1. Let xk be the variable
    //         contributing with the infinite bound. Then the following bounds
    //         are valid for xk:
    // * if aik > 0: xk <= (rhs[i] - \sum_{j != k} aij xj) / aik
    //                  <= (rhs[i] - \sum_{j != k : aij > 0} aij lb[j]
    //                             - \sum_{j != k : aij < 0} aij ub[j]) / aik
    //
    // * if aik < 0: xk >= (rhs[i] - \sum_{j \neq k} a_ij xj) / aik
    //                  >= (rhs[i] - \sum_{j != k : aij > 0} aij lb[j]
    //                             - \sum_{j != k : aij < 0} aij ub[j]) / aik
    //
    // Case 3: assume RHS is infinite, that the max activity is valid and
    //         nInfMin = 1. Let xk be the variable contributing with the
    //         infinite bound.
    //     * if aik > 0, then lb[k] = -inf and ub[k] < inf. The following bound
    //       is valid: xk <= (maxAct - minAct) / aik
    //     * if aik < 0, then ub[k] = inf and lb[k] > -inf. The following bound
    //       is valid: xk >= (maxAct - minAct) / aik.
    // -----------------------------------------------------------------------------
    else if ((!is_rhs_inf || act->n_inf_max == 0) && act->n_inf_min == 1)
    {
        // find the variable contributing with the infinite bound
        for (j = 0; j < *row->len; ++j)
        {
            Aik = row->vals[j];
            k = row->cols[j];
            is_lb_inf = HAS_TAG(col_tags[k], C_TAG_LB_INF);
            is_ub_inf = HAS_TAG(col_tags[k], C_TAG_UB_INF);
            is_col_inactive = HAS_TAG(col_tags[k], C_TAG_INACTIVE);

            if (((Aik > 0 && is_lb_inf) || (Aik < 0 && is_ub_inf)) &&
                !is_col_inactive)
            {
                break;
            }
        }

        assert(j != *row->len && Aik != 0 && !is_col_inactive);

        // ------------------------------------------------------------------------------
        // compute an implied bound and use it if it is tighter
        // ------------------------------------------------------------------------------
        min_act =
            compute_min_act_one_tag(row->vals, row->cols, *row->len, bounds, k);
        assert((is_rhs_inf && act->n_inf_max == 0 && ABS(act->max) != INF) ||
               (!is_rhs_inf && ABS(*row->rhs) != INF));
        assert(ABS(min_act) != INF && Aik != INF && Aik != 0);
        implied_bound =
            (is_rhs_inf) ? (act->max - min_act) / Aik : (*row->rhs - min_act) / Aik;

        if (Aik > 0)
        {
            assert(HAS_TAG(col_tags[k], C_TAG_LB_INF));
            if ((is_ub_inf || implied_bound < bounds[k].ub) &&
                !IS_HUGE(implied_bound))
            {
                temp_status =
                    updater_ub(implied_bound, &(bounds[k].ub), bounds[k].lb, k,
                               col_tags + k, problem, arg1);

                if (temp_status == REDUCED)
                {
                    (*num_of_bounds_changed)++;
                }
                status |= temp_status;
            }
        }
        else
        {
            assert(Aik < 0 && HAS_TAG(col_tags[k], C_TAG_UB_INF));
            if ((is_lb_inf || implied_bound > bounds[k].lb) &&
                !IS_HUGE(implied_bound))
            {
                temp_status =
                    updater_lb(implied_bound, &(bounds[k].lb), bounds[k].ub, k,
                               col_tags + k, problem, arg1);

                if (temp_status == REDUCED)
                {
                    (*num_of_bounds_changed)++;
                }

                status |= temp_status;
            }
        }
    }

    return status;
}

/* Bound tightening based on the left-hand side of the
   constraint. Bound tightening can occur in several scenarios:
   1) the LHS is finite and the max activity is valid (ie. when
      nInfMax == 0). In this case every variable in the constraint
      is a candidate for bound tightening.
   2) the LHS is finite and nInfMax == 1. In this case only the variable
      contributing with the infinite bound is a candidate for bound
      tightening.
   3) If the LHS is infinite, and the min activity is valid, and
      there is one variable contributing with an infinite bound to the
      max activity, then we may be able to tight the finite bound of
      the variable contributing with the infinite bound.
*/
static inline PresolveStatus
bound_tightening_single_row_lhs(const ConstRowView *row, const Activity *act,
                                Problem *problem, Bound *bounds, ColTag *col_tags,
                                void *arg1, BoundUpdater updater_lb,
                                BoundUpdater updater_ub, int *num_of_bounds_changed)
{
    int j = -1;
    int k = -1;
    double Aik = INF;
    double implied_ub, implied_lb, implied_bound, max_act;
    PresolveStatus status = UNCHANGED;
    PresolveStatus temp_status = UNCHANGED;
    bool is_lhs_inf = HAS_TAG(*row->tag, R_TAG_LHS_INF);
    bool is_lb_inf = false;
    bool is_ub_inf = false;
    bool is_col_inactive;

    // ------------------------------------------------------------------
    // Case 1: assume LHS is finite and the max activity is valid. Then
    //         the following bounds are valid for every xk in row i:
    //         * if aik > 0: xk >= ub[k] + (lhs[i] - maxAct) / aik
    //         * if aik < 0: xk <= lb[k] + (lhs[i] - maxAct) / aik.
    // ------------------------------------------------------------------
    if (!is_lhs_inf && act->n_inf_max == 0)
    {
        assert(ABS(*row->lhs) != INF && ABS(act->max) != INF);

        for (j = 0; j < *row->len; ++j)
        {
            k = row->cols[j];

            // Consider the following scenario. Suppose row 1 is forcing and
            // fixes x1. Assume that x1 appears in row 2. When we loop
            // through row 2 we want to skip x1.
            if (HAS_TAG(col_tags[k], C_TAG_INACTIVE))
            {
                continue;
            }

            Aik = row->vals[j];

            // --------------------------------------------------------------------
            // if Aik > 0 we might be able to tighten the upper bound if
            // implied bound is tighter (and not too large)
            // --------------------------------------------------------------------
            if (Aik > 0)
            {
                assert(ABS(bounds[k].ub) != INF &&
                       !HAS_TAG(col_tags[k], C_TAG_UB_INF));

                implied_lb = bounds[k].ub + (*row->lhs - act->max) / Aik;
                is_lb_inf = HAS_TAG(col_tags[k], C_TAG_LB_INF);

                if ((is_lb_inf || implied_lb > bounds[k].lb) && !IS_HUGE(implied_lb))
                {
                    temp_status =
                        updater_lb(implied_lb, &(bounds[k].lb), bounds[k].ub, k,
                                   col_tags + k, problem, arg1);

                    if (temp_status == REDUCED)
                    {
                        (*num_of_bounds_changed)++;
                    }

                    status |= temp_status;
                }
            }
            // ---------------------------------------------------------------------
            // if Aik < 0 we might be able to tighten the lower bound if
            // implied bound is tighter (and not too large)
            // ---------------------------------------------------------------------
            else
            {
                assert(Aik < 0 && ABS(bounds[k].lb) != -INF);
                assert(!HAS_TAG(col_tags[k], C_TAG_LB_INF));

                implied_ub = bounds[k].lb + (*row->lhs - act->max) / Aik;
                is_ub_inf = HAS_TAG(col_tags[k], C_TAG_UB_INF);

                if ((is_ub_inf || implied_ub < bounds[k].ub) && !IS_HUGE(implied_ub))
                {
                    temp_status =
                        updater_ub(implied_ub, &(bounds[k].ub), bounds[k].lb, k,
                                   col_tags + k, problem, arg1);

                    if (temp_status == REDUCED)
                    {
                        (*num_of_bounds_changed)++;
                    }

                    status |= temp_status;
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // Case 2: assume LHS is finite and nInfMax = 1. Let xk be the variable
    //         contributing with the infinite bound. Then the following bounds
    //         are valid for xk:
    // * if aik > 0: xk >= (lhs[i] - \sum_{j != k} aij xj) / aik
    //                  >= (lhs[i] - \sum_{j != k : aij > 0} aij ub[j]
    //                             - \sum_{j != k : aij < 0} aij lb[j]) / aik
    //
    // * if aik < 0: xk <= (lhs[i] - \sum_{j \neq k} a_ij xj) / aik
    //                  <= (lhs[i] - \sum_{j != k : aij > 0} aij ub[j]
    //                             - \sum_{j != k : aij < 0} aij lb[j]) / aik
    //
    // Case 3: assume LHS is infinite, that the min activity is valid and
    //         nInfMax = 1. Let xk be the variable contributing with the
    //         infinite bound.
    //     * if aik > 0, then ub[k] = inf and lb[k] > -inf. The following bound
    //       is valid: xk >= (minAct - maxAct) / aik
    //     * if aik < 0, then lb[k] = -inf and ub[k] < inf. The following bound
    //       is valid: xk <= (minAct - maxAct) / aik
    // -----------------------------------------------------------------------
    else if ((!is_lhs_inf || act->n_inf_min == 0) && act->n_inf_max == 1)
    {
        // find the variable contributing with the infinite bound
        for (j = 0; j < *row->len; ++j)
        {
            Aik = row->vals[j];
            k = row->cols[j];
            is_lb_inf = HAS_TAG(col_tags[k], C_TAG_LB_INF);
            is_ub_inf = HAS_TAG(col_tags[k], C_TAG_UB_INF);
            is_col_inactive = HAS_TAG(col_tags[k], C_TAG_INACTIVE);

            if (((Aik > 0 && is_ub_inf) || (Aik < 0 && is_lb_inf)) &&
                !is_col_inactive)
            {
                break;
            }
        }

        assert(j != *row->len && Aik != 0 && !is_col_inactive);

        // ------------------------------------------------------------------------------
        // compute an implied bound and use it if it is tighter
        // ------------------------------------------------------------------------------
        max_act =
            compute_max_act_one_tag(row->vals, row->cols, *row->len, bounds, k);
        assert((act->n_inf_min == 0 && ABS(act->min) != INF && is_lhs_inf) ||
               (!is_lhs_inf && ABS(*row->lhs) != INF));
        assert(max_act != INF);
        assert(Aik != INF && Aik != 0);
        implied_bound =
            (is_lhs_inf) ? (act->min - max_act) / Aik : (*row->lhs - max_act) / Aik;

        if (Aik > 0)
        {
            assert(HAS_TAG(col_tags[k], C_TAG_UB_INF));

            if ((is_lb_inf || implied_bound > bounds[k].lb) &&
                !IS_HUGE(implied_bound))
            {
                temp_status =
                    updater_lb(implied_bound, &(bounds[k].lb), bounds[k].ub, k,
                               col_tags + k, problem, arg1);

                if (temp_status == REDUCED)
                {
                    (*num_of_bounds_changed)++;
                }

                status |= temp_status;
            }
        }
        else
        {
            assert(Aik < 0 && HAS_TAG(col_tags[k], C_TAG_LB_INF));

            if ((is_ub_inf || implied_bound < bounds[k].ub) &&
                !IS_HUGE(implied_bound))
            {
                temp_status =
                    updater_ub(implied_bound, &(bounds[k].ub), bounds[k].lb, k,
                               col_tags + k, problem, arg1);

                if (temp_status == REDUCED)
                {
                    (*num_of_bounds_changed)++;
                }

                status |= temp_status;
            }
        }
    }

    return status;
}

/* This function is guaranteed to return infeasible if the problem is
   infeasible, but not necessarily reduced if the problem is reduced. This
   is a design choice. */
PresolveStatus bound_tightening_single_row(const ConstRowView *row,
                                           const Activity *act, Problem *problem,
                                           Bound *bounds, ColTag *col_tags,
                                           void *arg1, BoundUpdater updater_lb,
                                           BoundUpdater updater_ub)
{
    assert(!HAS_TAG(*row->tag, R_TAG_INACTIVE) &&
           (act->n_inf_min <= 1 || act->n_inf_max <= 1));

    PresolveStatus status = UNCHANGED;

    int num_of_bound_changes = 0;

    status |= bound_tightening_single_row_rhs(row, act, problem, bounds, col_tags,
                                              arg1, updater_lb, updater_ub,
                                              &num_of_bound_changes);

    status |= bound_tightening_single_row_lhs(row, act, problem, bounds, col_tags,
                                              arg1, updater_lb, updater_ub,
                                              &num_of_bound_changes);

    // dual postsolve if a bound change occured (this assumes primal propagation. We
    // must make sure this is not called during dual propagation)
    if (num_of_bound_changes > 0)
    {
        save_retrieval_bound_change_the_row(
            problem->constraints->state->postsolve_info, row->i, row->cols,
            row->vals, *row->len, num_of_bound_changes);
    }
    return status;
}

PresolveStatus propagate_primal(Problem *prob, bool finite_bound_tightening)
{
    assert(prob->constraints->state->ston_rows->len == 0);
    assert(prob->constraints->state->empty_rows->len == 0);
    DEBUG(verify_problem_up_to_date(prob->constraints));

    // empty columns may occur because of check_activities
    // assert(prob->constraints->state->empty_cols->len == 0);

    Constraints *constraints = prob->constraints;
    const iVec *updated_activities = constraints->state->updated_activities;
    Activity *acts = constraints->state->activities;
    const int *row_sizes = constraints->state->row_sizes;
    const RowTag *row_tags = constraints->row_tags;
    const Matrix *A = constraints->A;
    Bound *bounds = constraints->bounds;
    ColTag *col_tags = constraints->col_tags;
    const double *rhs = constraints->rhs;
    const double *lhs = constraints->lhs;
    int i, ii, current_len;
    PresolveStatus status = UNCHANGED;

    // check that updated_activities has no duplicates, and that all
    // activities in updated_activities have status ADDED, and that these
    // have n_inf_min = 0 or n_inf_max = 0.
    DEBUG(verify_no_duplicates_sort(updated_activities));
    DEBUG(verify_row_states(acts, updated_activities));

#ifndef NDEBUG
    bool *HAVE_ROWS_BEEN_PROP = (bool *) ps_calloc(A->m, sizeof(bool));
#endif

    // -------------------------------------------------------------------------
    // Loop through and propagate the rows in order to tighten variable bounds.
    // We skip rows that are that are inactive, or will be handled by
    // empty rows / stonrows, or do not have the potential to offer a reduction.
    //
    // Note that bound_tightening_single_row might append to updated_activities,
    // so we must be careful when we loop through it. This is why we use the
    // double-loop below. The debug code below is used to verify that we
    // perform at most one round of domain propagation on each constraint.
    // -------------------------------------------------------------------------
    ii = 0;
    while (ii < updated_activities->len)
    {
        current_len = updated_activities->len;
        for (; ii < current_len; ++ii)
        {
            i = updated_activities->data[ii];

            // can it ever happen that n_inf_min and n_inf_max are both not
            // zero?
            assert(acts[i].n_inf_min == 0 || acts[i].n_inf_max == 0);
            assert(acts[i].status != NOT_ADDED);

            if (row_sizes[i] <= 1 || acts[i].status != ADDED ||
                !(acts[i].n_inf_min == 0 || acts[i].n_inf_max == 0 ||
                  (acts[i].n_inf_min == 1 && !HAS_TAG(row_tags[i], R_TAG_RHS_INF)) ||
                  (acts[i].n_inf_max == 1 && !HAS_TAG(row_tags[i], R_TAG_LHS_INF))))
            {
                continue;
            }

            assert(row_sizes[i] > 1 && !HAVE_ROWS_BEEN_PROP[i]);
            DEBUG(HAVE_ROWS_BEEN_PROP[i] = true;);

            ConstRowView row = new_const_rowview(
                A->x + A->p[i].start, A->i + A->p[i].start, row_sizes + i, A->p + i,
                lhs + i, rhs + i, row_tags + i, i);

            acts[i].status = PROPAGATED_THIS_ROUND;

            status = bound_tightening_single_row(
                &row, acts + i, prob, bounds, col_tags, &finite_bound_tightening,
                (BoundUpdater) update_lb_within_propagation,
                (BoundUpdater) update_ub_within_propagation);

            if (HAS_STATUS(status, INFEASIBLE))
            {
                return INFEASIBLE;
            }
        }
    }

#ifndef NDEBUG
    for (ii = 0; ii < A->m; ++ii)
    {
        if (HAS_TAG(row_tags[ii], R_TAG_INACTIVE))
        {
            continue;
        }

        if (acts[ii].status != NOT_ADDED)
        {
            assert(HAVE_ROWS_BEEN_PROP[ii]);
        }
        else
        {
            assert(!HAVE_ROWS_BEEN_PROP[ii]);
        }
    }
#endif

    DEBUG(PS_FREE(HAVE_ROWS_BEEN_PROP););

    // -----------------------------------------------------------------------------
    // Loop through updated activities and add some of them to the next round.
    // Also reset the status of all activities.
    // -----------------------------------------------------------------------------
    int *iwork_n_rows = constraints->state->work->iwork_n_rows;
    int new_len = 0;
    for (ii = 0; ii < updated_activities->len; ++ii)
    {
        if (acts[updated_activities->data[ii]].status == PROPAGATE_NEXT_ROUND)
        {
            iwork_n_rows[new_len++] = updated_activities->data[ii];
            acts[updated_activities->data[ii]].status = ADDED;
        }
    }
    DEBUG(verify_no_duplicates_sort_ptr(iwork_n_rows, new_len));

    for (ii = 0; ii < A->m; ++ii)
    {
        acts[ii].status = NOT_ADDED;
    }

    for (ii = 0; ii < new_len; ++ii)
    {
        acts[iwork_n_rows[ii]].status = ADDED;
    }

    iVec_clear_no_resize(constraints->state->updated_activities);
    iVec_append_array(constraints->state->updated_activities, iwork_n_rows, new_len);
    //----------------------------------------------------------------------------

    // delete contribution of fixed columns from rhs and lhs etc.
    delete_fixed_cols_from_problem(prob);

    // The domain propagation might fix columns, so there may be inactive cols.
    // When we fix these columns we might get zero rows etc.
    delete_inactive_cols_from_A_and_AT(constraints);

    // We have not yet deleted the zero rows so there should be no rows
    // to delete right now.
    assert(constraints->state->rows_to_delete->len == 0);

    // the rows that were forcing are empty now so we process them
    if (remove_empty_rows(constraints) == INFEASIBLE)
    {
        return INFEASIBLE;
    }

    DEBUG(verify_problem_up_to_date(constraints));
    return UNCHANGED;
}

// Static variable to hold col_sizes (for qsort)
static const int *global_col_sizes;

int compare_col_len(const void *a, const void *b)
{
    int idx_a = *(const int *) a;
    int idx_b = *(const int *) b;

    if (global_col_sizes[idx_a] < global_col_sizes[idx_b])
    {
        return -1;
    }
    if (global_col_sizes[idx_a] > global_col_sizes[idx_b])
    {
        return 1;
    }

    // if equal len, sort by index
    return idx_a - idx_b;
}

void remove_redundant_bounds(Constraints *constraints)
{
    const Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT;
    double *lhs = constraints->lhs;
    double *rhs = constraints->rhs;
    RowTag *row_tags = constraints->row_tags;
    ColTag *col_tags = constraints->col_tags;
    Bound *bounds = constraints->bounds;
    int n_cols = constraints->n;
    Activity *acts = constraints->state->activities;
    const int *col_sizes = constraints->state->col_sizes;

    bool is_ub_inf, is_lb_inf;
    int ii, jj;
    int removed_bounds = 0;

    // ------------------------------------------------------------------------
    // Sort columns by column sizes. We want to process the shortest columns
    // first, since when we remove a bound from a variable, the activities of
    // all rows containing that variable become invalid.
    // ------------------------------------------------------------------------
    int *col_order = constraints->state->work->iwork_n_cols;
    for (ii = 0; ii < n_cols; ++ii)
    {
        col_order[ii] = ii;
    }

    global_col_sizes = col_sizes;
    qsort(col_order, n_cols, sizeof(int), compare_col_len);
    global_col_sizes = NULL;

    for (jj = 0; jj < n_cols; jj++)
    {
        ii = col_order[jj];
        if (HAS_TAG(col_tags[ii], C_TAG_INACTIVE))
        {
            continue;
        }

        is_ub_inf = HAS_TAG(col_tags[ii], C_TAG_UB_INF);
        is_lb_inf = HAS_TAG(col_tags[ii], C_TAG_LB_INF);

        if (is_ub_inf && is_lb_inf)
        {
            continue;
        }

        ConstColView col = new_const_colview(
            AT->x + AT->p[ii].start, AT->i + AT->p[ii].start, col_sizes + ii,
            AT->p + ii, &bounds[ii].lb, &bounds[ii].ub, col_tags + ii, ii);

        // TODO: can we check if we can remove both bounds?
        if (is_ub_inf)
        {
            assert(!is_lb_inf);

            if (is_lb_implied_free(&col, acts, lhs, rhs, row_tags, A, col_tags,
                                   bounds))
            {
                remove_finite_lb_from_activities(&col, acts, bounds[ii].lb);
                removed_bounds++;
                DEBUG(bounds[ii].lb = -INF);
                UPDATE_TAG(col_tags[ii], C_TAG_LB_INF);
            }
        }
        else
        {
            assert(!is_ub_inf);

            if (is_ub_implied_free(&col, acts, lhs, rhs, row_tags, A, col_tags,
                                   bounds))
            {
                remove_finite_ub_from_activities(&col, acts, bounds[ii].ub);
                removed_bounds++;
                DEBUG(bounds[ii].ub = INF);
                UPDATE_TAG(col_tags[ii], C_TAG_UB_INF);
            }
        }
    }

    DEBUG(verify_activities(constraints));
}
