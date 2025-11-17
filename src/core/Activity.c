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

#include "Activity.h"
#include "Bounds.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "RowColViews.h"
#include "glbopts.h"

Activity *new_activities(const Matrix *A, const ColTag *col_tags,
                         const Bound *bounds)
{
    Activity *activities = (Activity *) ps_malloc(A->m, sizeof(Activity));
    RETURN_PTR_IF_NULL(activities, NULL);
    int i, j, start, end;
    int *cols = A->i;
    double *vals = A->x;

    for (i = 0; i < A->m; ++i)
    {
        Activity *act = activities + i;
        act->n_inf_min = 0;
        act->n_inf_max = 0;
        act->max = 0.0;
        act->min = 0.0;
        act->status = NOT_ADDED;

        start = A->p[i].start;
        end = A->p[i].end;

        // record infinite contributions
        for (j = start; j < end; ++j)
        {
            if (vals[j] > 0)
            {
                if (HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF))
                {
                    act->n_inf_max++;
                }

                if (HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF))
                {
                    act->n_inf_min++;
                }
            }
            else
            {
                if (HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF))
                {
                    act->n_inf_max++;
                }

                if (HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF))
                {
                    act->n_inf_min++;
                }
            }
        }

        // compute max act if it is useful
        if (act->n_inf_max == 0)
        {
            for (j = start; j < end; ++j)
            {
                if (vals[j] > 0)
                {
                    assert(!HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF));
                    assert(bounds[cols[j]].ub != INF);
                    act->max += vals[j] * bounds[cols[j]].ub;
                }
                else
                {
                    assert(!HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF));
                    assert(bounds[cols[j]].lb != -INF);
                    act->max += vals[j] * bounds[cols[j]].lb;
                }
            }
        }
#ifndef NDEBUG
        else
        {
            act->max = INVALID_ACT_DEBUG;
        }
#endif

        // compute min act if it is useful
        if (act->n_inf_min == 0)
        {
            for (j = start; j < end; ++j)
            {
                if (vals[j] > 0)
                {
                    assert(!HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF));
                    assert(bounds[cols[j]].lb != -INF);
                    act->min += vals[j] * bounds[cols[j]].lb;
                }
                else
                {
                    assert(!HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF));
                    assert(bounds[cols[j]].ub != INF);
                    act->min += vals[j] * bounds[cols[j]].ub;
                }
            }
        }
#ifndef NDEBUG
        else
        {
            act->min = INVALID_ACT_DEBUG;
        }
#endif
    }

    return activities;
}

void free_activities(Activity *activities)
{
    PS_FREE(activities);
}

void Activities_init(const Matrix *A, const ColTag *col_tags, const RowTag *row_tags,
                     const Bound *bounds, Activity *acts)
{

    for (int i = 0; i < A->m; ++i)
    {
        if (HAS_TAG(row_tags[i], R_TAG_INACTIVE))
        {
            DEBUG(acts[i].status = NOT_ADDED;);
            continue;
        }

        Activity_init(acts + i, A->x + A->p[i].start, A->i + A->p[i].start,
                      A->p[i].end - A->p[i].start, bounds, col_tags);
    }
}

void Activity_init(Activity *act, const double *vals, const int *cols, int len,
                   const Bound *bounds, const ColTag *col_tags)
{
    act->n_inf_min = 0;
    act->n_inf_max = 0;
    act->max = 0.0;
    act->min = 0.0;
    act->status = NOT_ADDED;
    int j;

    // record infinite contributions
    for (j = 0; j < len; ++j)
    {
        if (vals[j] > 0)
        {
            if (HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF))
            {
                act->n_inf_max++;
            }

            if (HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF))
            {
                act->n_inf_min++;
            }
        }
        else
        {
            if (HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF))
            {
                act->n_inf_max++;
            }

            if (HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF))
            {
                act->n_inf_min++;
            }
        }
    }

    // compute max act if it is useful
    if (act->n_inf_max == 0)
    {
        for (j = 0; j < len; ++j)
        {
            if (vals[j] > 0)
            {
                assert(!HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF));
                assert(bounds[cols[j]].ub != INF);
                act->max += vals[j] * bounds[cols[j]].ub;
            }
            else
            {
                assert(!HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF));
                assert(bounds[cols[j]].lb != -INF);
                act->max += vals[j] * bounds[cols[j]].lb;
            }
        }
    }
#ifndef NDEBUG
    else
    {
        act->max = INVALID_ACT_DEBUG;
    }
#endif

    // compute min act if it is useful
    if (act->n_inf_min == 0)
    {
        for (j = 0; j < len; ++j)
        {
            if (vals[j] > 0)
            {
                assert(!HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF));
                assert(bounds[cols[j]].lb != -INF);
                act->min += vals[j] * bounds[cols[j]].lb;
            }
            else
            {
                assert(!HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF));
                assert(bounds[cols[j]].ub != INF);
                act->min += vals[j] * bounds[cols[j]].ub;
            }
        }
    }
#ifndef NDEBUG
    else
    {
        act->min = INVALID_ACT_DEBUG;
    }
#endif
}

double compute_min_act_tags(const double *vals, const int *cols, int len,
                            const Bound *bounds, const ColTag *col_tags)
{
    double min_act = 0.0;
    for (int j = 0; j < len; ++j)
    {
        if (vals[j] > 0)
        {
            if (!HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF))
            {
                assert(bounds[cols[j]].lb != -INF);
                min_act += vals[j] * bounds[cols[j]].lb;
            }
        }
        else
        {
            if (!HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF))
            {
                assert(bounds[cols[j]].ub != INF);
                min_act += vals[j] * bounds[cols[j]].ub;
            }
        }
    }
    return min_act;
}

double compute_max_act_tags(const double *vals, const int *cols, int len,
                            const Bound *bounds, const ColTag *col_tags)
{
    double max_act = 0.0;
    for (int j = 0; j < len; ++j)
    {
        if (vals[j] > 0)
        {
            if (!HAS_TAG(col_tags[cols[j]], C_TAG_UB_INF))
            {
                assert(bounds[cols[j]].ub != INF);
                max_act += vals[j] * bounds[cols[j]].ub;
            }
        }
        else
        {
            if (!HAS_TAG(col_tags[cols[j]], C_TAG_LB_INF))
            {
                assert(bounds[cols[j]].lb != -INF);
                max_act += vals[j] * bounds[cols[j]].lb;
            }
        }
    }
    return max_act;
}

double compute_min_act_no_tags(const double *vals, const int *cols, int len,
                               const Bound *bounds)
{
    double min_act = 0.0;
    for (int j = 0; j < len; ++j)
    {
        if (vals[j] > 0)
        {
            assert(bounds[cols[j]].lb != -INF);
            min_act += vals[j] * bounds[cols[j]].lb;
        }
        else
        {
            assert(bounds[cols[j]].ub != INF);
            min_act += vals[j] * bounds[cols[j]].ub;
        }
    }
    return min_act;
}

double compute_max_act_no_tags(const double *vals, const int *cols, int len,
                               const Bound *bounds)
{
    double max_act = 0.0;
    for (int j = 0; j < len; ++j)
    {
        if (vals[j] > 0)
        {
            assert(bounds[cols[j]].ub != INF);
            max_act += vals[j] * bounds[cols[j]].ub;
        }
        else
        {
            assert(bounds[cols[j]].lb != -INF);
            max_act += vals[j] * bounds[cols[j]].lb;
        }
    }
    return max_act;
}

double compute_min_act_one_tag(const double *vals, const int *cols, int len,
                               const Bound *bounds, int col_to_avoid)
{
    assert(col_to_avoid >= 0);
    double min_act = 0.0;
    for (int j = 0; j < len; ++j)
    {
        if (cols[j] == col_to_avoid)
        {
            continue;
        }

        if (vals[j] > 0)
        {
            assert(bounds[cols[j]].lb != -INF);
            min_act += vals[j] * bounds[cols[j]].lb;
        }
        else
        {
            assert(bounds[cols[j]].ub != INF);
            min_act += vals[j] * bounds[cols[j]].ub;
        }
    }
    return min_act;
}

double compute_max_act_one_tag(const double *vals, const int *cols, int len,
                               const Bound *bounds, int col_to_avoid)
{
    assert(col_to_avoid >= 0);
    double max_act = 0.0;
    for (int j = 0; j < len; ++j)
    {
        if (cols[j] == col_to_avoid)
        {
            continue;
        }

        if (vals[j] > 0)
        {
            assert(bounds[cols[j]].ub != INF);
            max_act += vals[j] * bounds[cols[j]].ub;
        }
        else
        {
            assert(bounds[cols[j]].lb != -INF);
            max_act += vals[j] * bounds[cols[j]].lb;
        }
    }
    return max_act;
}

Altered_Activity update_act_bound_change(Activity *act, double coeff,
                                         double old_bound, double new_bound,
                                         bool finite_bound, bool lower,
                                         const double *vals, const int *cols,
                                         int len, const Bound *bounds)
{
    // ---------------------------------------------------------------------
    //              If a lower bound has been updated
    // ---------------------------------------------------------------------
    if (lower)
    {
        if (coeff < 0)
        {
            if (act->n_inf_max == 0)
            {
                assert(finite_bound);
                act->max += (new_bound - old_bound) * coeff;
            }
            else if (!finite_bound)
            {
                act->n_inf_max--;

                if (act->n_inf_max == 0)
                {
                    act->max = compute_max_act_no_tags(vals, cols, len, bounds);
                }
            }

            return MAX_ALTERED;
        }
        else
        {
            if (act->n_inf_min == 0)
            {
                assert(finite_bound);
                act->min += (new_bound - old_bound) * coeff;
            }
            else if (!finite_bound)
            {
                act->n_inf_min--;

                if (act->n_inf_min == 0)
                {
                    act->min = compute_min_act_no_tags(vals, cols, len, bounds);
                }
            }
            return MIN_ALTERED;
        }
    }
    // ------------------------------------------------------------------------
    //                 If an upper bound has been updated
    // ------------------------------------------------------------------------
    else
    {
        if (coeff < 0)
        {
            if (act->n_inf_min == 0)
            {
                assert(finite_bound);
                act->min += (new_bound - old_bound) * coeff;
            }
            else if (!finite_bound)
            {
                act->n_inf_min--;

                if (act->n_inf_min == 0)
                {
                    act->min = compute_min_act_no_tags(vals, cols, len, bounds);
                }
            }
            return MIN_ALTERED;
        }
        else
        {
            if (act->n_inf_max == 0)
            {
                assert(finite_bound);
                act->max += (new_bound - old_bound) * coeff;
            }
            else if (!finite_bound)
            {
                act->n_inf_max--;

                if (act->n_inf_max == 0)
                {
                    act->max = compute_max_act_no_tags(vals, cols, len, bounds);
                }
            }
            return MAX_ALTERED;
        }
    }
}

void update_activities_bound_change(Activity *activities, const Matrix *A,
                                    const Bound *bounds, const double *vals,
                                    const int *rows, int len, double old_bound,
                                    double new_bound, double finite_bound,
                                    bool lower,
                                    iVec *updated_activities HUGE_BOUND_PARAM)
{
    assert(!IS_HUGE(new_bound) || huge_bound_ok);
    assert(!finite_bound || (lower && new_bound > old_bound) ||
           (!lower && new_bound < old_bound));

    Altered_Activity altered;
    for (int ii = 0; ii < len; ++ii)
    {
        int i = rows[ii];
        Activity *act = activities + i;
        altered = update_act_bound_change(act, vals[ii], old_bound, new_bound,
                                          finite_bound, lower, A->x + A->p[i].start,
                                          A->i + A->p[i].start,
                                          A->p[i].end - A->p[i].start, bounds);

        if (((act->n_inf_max == 0 && (altered & MAX_ALTERED)) ||
             (act->n_inf_min == 0 && (altered & MIN_ALTERED))))
        {
            if (act->status == NOT_ADDED)
            {
                act->status = ADDED;
                iVec_append(updated_activities, rows[ii]);
            }
            else if (act->status == PROPAGATED_THIS_ROUND)
            {
                act->status = PROPAGATE_NEXT_ROUND;
                iVec_append(updated_activities, rows[ii]);
            }
        }
    }
}

void update_activities_fixed_col(Activity *activities, const Matrix *A,
                                 const Bound *bounds, const double *vals,
                                 const int *rows, int len, double old_ub,
                                 double old_lb, double val, bool ub_update,
                                 bool lb_update, bool is_ub_inf, bool is_lb_inf,
                                 iVec *updated_activities)
{
    Altered_Activity altered1;
    Altered_Activity altered2;
    Altered_Activity altered;

    for (int ii = 0; ii < len; ++ii)
    {
        altered1 = NO_RECOMPUTE;
        altered2 = NO_RECOMPUTE;
        int i = rows[ii];
        Activity *act = activities + i;

        if (ub_update)
        {
            altered1 = update_act_bound_change(
                act, vals[ii], old_ub, val, !is_ub_inf, false, A->x + A->p[i].start,
                A->i + A->p[i].start, A->p[i].end - A->p[i].start, bounds);
        }

        if (lb_update)
        {
            altered2 = update_act_bound_change(
                act, vals[ii], old_lb, val, !is_lb_inf, true, A->x + A->p[i].start,
                A->i + A->p[i].start, A->p[i].end - A->p[i].start, bounds);
        }

        altered = altered1 | altered2;
        if (((act->n_inf_max == 0 && (altered & MAX_ALTERED)) ||
             (act->n_inf_min == 0 && (altered & MIN_ALTERED))))
        {
            if (act->status == NOT_ADDED)
            {
                act->status = ADDED;
                iVec_append(updated_activities, i);
            }
            else if (act->status == PROPAGATED_THIS_ROUND)
            {
                act->status = PROPAGATE_NEXT_ROUND;
                iVec_append(updated_activities, i);
            }
        }
    }
}

Altered_Activity update_activity_coeff_change(Activity *act, double lb, double ub,
                                              double old_coeff, double new_coeff,
                                              ColTag cTag)
{
    assert(!HAS_TAG(cTag, C_TAG_INACTIVE));
    bool is_lb_finite = !HAS_TAG(cTag, C_TAG_LB_INF);
    bool is_ub_finite = !HAS_TAG(cTag, C_TAG_UB_INF);

    int n_inf_min_before = act->n_inf_min;
    int n_inf_max_before = act->n_inf_max;

    // ----------------------------------------
    // remove contribution from old coefficient
    // ----------------------------------------
    if (old_coeff > 0)
    {
        if (is_lb_finite)
        {
            if (act->n_inf_min == 0)
            {
                act->min -= old_coeff * lb;
            }
        }
        else
        {
            act->n_inf_min--;
        }

        if (is_ub_finite)
        {
            if (act->n_inf_max == 0)
            {
                act->max -= old_coeff * ub;
            }
        }
        else
        {
            act->n_inf_max--;
        }
    }
    else if (old_coeff < 0)
    {
        if (is_lb_finite)
        {
            if (act->n_inf_max == 0)
            {
                act->max -= old_coeff * lb;
            }
        }
        else
        {
            act->n_inf_max--;
        }

        if (is_ub_finite)
        {
            if (act->n_inf_min == 0)
            {
                act->min -= old_coeff * ub;
            }
        }
        else
        {
            act->n_inf_min--;
        }
    }

    // -----------------------------------
    // add contribution of new coefficient
    // -----------------------------------
    if (new_coeff > 0)
    {
        if (is_lb_finite)
        {
            if (act->n_inf_min == 0)
            {
                act->min += new_coeff * lb;
            }
        }
        else
        {
            act->n_inf_min++;
        }

        if (is_ub_finite)
        {
            if (act->n_inf_max == 0)
            {
                act->max += new_coeff * ub;
            }
        }
        else
        {
            act->n_inf_max++;
        }
    }
    else if (new_coeff < 0)
    {
        if (is_lb_finite)
        {
            if (act->n_inf_max == 0)
            {
                act->max += new_coeff * lb;
            }
        }
        else
        {
            act->n_inf_max++;
        }

        if (is_ub_finite)
        {
            if (act->n_inf_min == 0)
            {
                act->min += new_coeff * ub;
            }
        }
        else
        {
            act->n_inf_min++;
        }
    }

#ifndef NDEBUG
    if (n_inf_min_before == 0 && act->n_inf_min == 1)
    {
        act->min = INVALID_ACT_DEBUG;
    }

    if (n_inf_max_before == 0 && act->n_inf_max == 1)
    {
        act->max = INVALID_ACT_DEBUG;
    }
#endif

    // recompute from scratch if necessary
    if (act->n_inf_min == 0 && n_inf_min_before == 1)
    {
        return MIN_ALTERED_RECOMPUTE;
    }

    if (act->n_inf_max == 0 && n_inf_max_before == 1)
    {
        return MAX_ALTERED_RECOMPUTE;
    }

    return NO_RECOMPUTE;
}

void remove_finite_lb_from_activities(const ConstColView *col, Activity *acts,
                                      double bound)
{

    for (int ii = 0; ii < *col->len; ii++)
    {
        int i = col->rows[ii];
        Activity *act = acts + i;

        if (col->vals[ii] > 0)
        {
            if (act->n_inf_min == 0)
            {
                act->min -= col->vals[ii] * bound;
            }

            act->n_inf_min++;
        }
        else if (col->vals[ii] < 0)
        {
            if (act->n_inf_max == 0)
            {
                act->max -= col->vals[ii] * bound;
            }
            act->n_inf_max++;
        }

        DEBUG(acts[i].min =
                  (acts[i].n_inf_min != 0) ? INVALID_ACT_DEBUG : acts[i].min);
        DEBUG(acts[i].max =
                  (acts[i].n_inf_max != 0) ? INVALID_ACT_DEBUG : acts[i].max);
    }

    return;
}

void remove_finite_ub_from_activities(const ConstColView *col, Activity *acts,
                                      double bound)
{

    for (int ii = 0; ii < *col->len; ii++)
    {
        int i = col->rows[ii];
        Activity *act = acts + i;

        if (col->vals[ii] > 0)
        {
            if (act->n_inf_max == 0)
            {
                act->max -= col->vals[ii] * bound;
            }

            act->n_inf_max++;
        }
        else if (col->vals[ii] < 0)
        {
            if (act->n_inf_min == 0)
            {
                act->min -= col->vals[ii] * bound;
            }
            act->n_inf_min++;
        }

        DEBUG(acts[i].min =
                  (acts[i].n_inf_min != 0) ? INVALID_ACT_DEBUG : acts[i].min);
        DEBUG(acts[i].max =
                  (acts[i].n_inf_max != 0) ? INVALID_ACT_DEBUG : acts[i].max);
    }

    return;
}

void recompute_n_infs(Activity *act, const double *vals, const int *cols, int len,
                      const ColTag *coltags)
{
    int n_inf_max = 0;
    int n_inf_min = 0;

    for (int j = 0; j < len; ++j)
    {
        if (HAS_TAG(coltags[cols[j]], C_TAG_INACTIVE))
        {
            continue;
        }

        if (vals[j] > 0)
        {
            if (HAS_TAG(coltags[cols[j]], C_TAG_UB_INF))
            {
                n_inf_max++;
            }

            if (HAS_TAG(coltags[cols[j]], C_TAG_LB_INF))
            {
                n_inf_min++;
            }
        }
        else
        {
            if (HAS_TAG(coltags[cols[j]], C_TAG_LB_INF))
            {
                n_inf_max++;
            }

            if (HAS_TAG(coltags[cols[j]], C_TAG_UB_INF))
            {
                n_inf_min++;
            }
        }
    }

    act->n_inf_max = n_inf_max;
    act->n_inf_min = n_inf_min;
}
