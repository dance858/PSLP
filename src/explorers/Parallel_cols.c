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

#include "Parallel_cols.h"
#include "Activity.h"
#include "Bounds.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "Matrix.h"
#include "Numerics.h"
#include "Parallel_rows.h"
#include "Problem.h"
#include "RowColViews.h"
#include "SimpleReductions.h"
#include "State.h"
#include "Tags.h"
#include "Workspace.h"

static inline PresolveStatus process_single_bin(const Problem *prob, const int *bin,
                                                int bin_size)
{
    assert(bin_size > 1);
    Constraints *constraints = prob->constraints;
    ColTag *col_tags = constraints->col_tags;
    Bound *bounds = constraints->bounds;
    const Matrix *AT = constraints->AT;
    const RowRange *col_ranges = AT->p;
    const double *c = prob->obj->c;
    Activity *acts = constraints->state->activities;
    iVec *sub_cols_to_delete = constraints->state->sub_cols_to_delete;
    PostsolveInfo *postsolve_info = constraints->state->postsolve_info;

    bool recount_ninfs = false;

    // column j stays, column k goes
    int ii, jj, j, k, len_j;
    const int *aj_cols;
    double ratio, cj, ck, ak0, lb_j_old, ub_j_old;
    const double *aj_vals;

    for (ii = 0; ii < bin_size - 1; ++ii)
    {
        j = bin[ii];

        if (HAS_TAG(col_tags[j], C_TAG_INACTIVE))
        {
            continue;
        }

        cj = c[j];
        aj_vals = AT->x + col_ranges[j].start;
        aj_cols = AT->i + col_ranges[j].start;
        len_j = col_ranges[j].end - col_ranges[j].start;
        recount_ninfs = false;

        for (jj = ii + 1; jj < bin_size; ++jj)
        {
            k = bin[jj];

            if (HAS_TAG(col_tags[k], C_TAG_INACTIVE))
            {
                continue;
            }

            ck = c[k];
            ak0 = AT->x[col_ranges[k].start];
            ratio = ak0 / aj_vals[0];
            assert(ratio != 0);

            // if column j and column k are parallel, remove column k
            // if you run into issues with parallel cols activated, try to
            // change this
            if (IS_EQUAL_FEAS_TOL(ck * aj_vals[0], cj * ak0))
            {
                // mark that we must recount n_infs for rows in which column j
                // appears in
                recount_ninfs = true;

                // ---------------------------------------------------------------------
                //                    Update the bounds of column j
                // ---------------------------------------------------------------------
                lb_j_old = bounds[j].lb;
                ub_j_old = bounds[j].ub;
                if (ratio > 0)
                {
                    // update lb
                    if (HAS_TAG((col_tags[j] | col_tags[k]), C_TAG_LB_INF))
                    {
                        bounds[j].lb = -INF;
                        UPDATE_TAG(col_tags[j], C_TAG_LB_INF);
                    }
                    else
                    {
                        assert(!IS_ABS_INF(bounds[j].lb) &&
                               !IS_ABS_INF(bounds[k].lb));
                        bounds[j].lb += bounds[k].lb * ratio;
                    }

                    // update ub
                    if (HAS_TAG((col_tags[j] | col_tags[k]), C_TAG_UB_INF))
                    {
                        bounds[j].ub = INF;
                        UPDATE_TAG(col_tags[j], C_TAG_UB_INF);
                    }
                    else
                    {
                        assert(!IS_ABS_INF(bounds[j].ub) &&
                               !IS_ABS_INF(bounds[k].ub));
                        bounds[j].ub += bounds[k].ub * ratio;
                    }
                }
                else
                {
                    // update lb
                    if (HAS_TAG(col_tags[j], C_TAG_LB_INF) ||
                        HAS_TAG(col_tags[k], C_TAG_UB_INF))
                    {
                        bounds[j].lb = -INF;
                        UPDATE_TAG(col_tags[j], C_TAG_LB_INF);
                    }
                    else
                    {
                        assert(!IS_ABS_INF(bounds[j].lb) &&
                               !IS_ABS_INF(bounds[k].ub));
                        bounds[j].lb += bounds[k].ub * ratio;
                    }

                    // update ub
                    if (HAS_TAG(col_tags[j], C_TAG_UB_INF) ||
                        HAS_TAG(col_tags[k], C_TAG_LB_INF))
                    {
                        bounds[j].ub = INF;
                        UPDATE_TAG(col_tags[j], C_TAG_UB_INF);
                    }
                    else
                    {
                        assert(!IS_ABS_INF(bounds[j].ub) &&
                               !IS_ABS_INF(bounds[k].lb));
                        bounds[j].ub += bounds[k].lb * ratio;
                    }
                }

                // ---------------------------------------------------------------------
                //                  mark col as substituted
                // ---------------------------------------------------------------------
                set_col_to_substituted(k, col_tags + k, sub_cols_to_delete);
                save_retrieval_parallel_col(postsolve_info, ub_j_old, lb_j_old,
                                            bounds[k].lb, bounds[k].ub, ratio, j, k,
                                            col_tags[j], col_tags[k]);
            }
            else
            {
                // we could add an implied check here but it doesn't seem to
                // make a difference on miplib. If you do it, don't forget that
                // the implied free check implicitly corresponds to a bound
                // change that must be dual postsolved

                bool fix_xk_to_lower = false;
                bool fix_xk_to_upper = false;
                bool fix_xj_to_upper = false;
                bool fix_xj_to_lower = false;

                assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
                assert(!HAS_TAG(col_tags[j], C_TAG_INACTIVE));

                if (ck > ratio * cj)
                {
                    if (ratio > 0)
                    {
                        fix_xk_to_lower = HAS_TAG(col_tags[j], C_TAG_UB_INF);
                        fix_xj_to_upper = HAS_TAG(col_tags[k], C_TAG_LB_INF);
                    }
                    else
                    {
                        fix_xk_to_lower = HAS_TAG(col_tags[j], C_TAG_LB_INF);
                        fix_xj_to_lower = HAS_TAG(col_tags[k], C_TAG_LB_INF);
                    }
                }
                else
                {
                    assert(ck < ratio * cj);
                    if (ratio > 0)
                    {
                        fix_xk_to_upper = HAS_TAG(col_tags[j], C_TAG_LB_INF);
                        fix_xj_to_lower = HAS_TAG(col_tags[k], C_TAG_UB_INF);
                    }
                    else
                    {
                        fix_xk_to_upper = HAS_TAG(col_tags[j], C_TAG_UB_INF);
                        fix_xj_to_upper = HAS_TAG(col_tags[k], C_TAG_UB_INF);
                    }
                }

                // check if we can fix xk
                if (fix_xk_to_lower)
                {
                    if (HAS_TAG(col_tags[k], C_TAG_LB_INF))
                    {
                        return UNBNDORINFEAS;
                    }

                    assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
                    assert(!IS_ABS_INF(bounds[k].lb));
                    fix_col(constraints, k, bounds[k].lb, ck);
                    continue;
                }
                else if (fix_xk_to_upper)
                {
                    if (HAS_TAG(col_tags[k], C_TAG_UB_INF))
                    {
                        return UNBNDORINFEAS;
                    }

                    assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
                    assert(!IS_ABS_INF(bounds[k].ub));
                    fix_col(constraints, k, bounds[k].ub, ck);
                    continue;
                }

                // check if we can fix xj
                if (fix_xj_to_lower)
                {
                    if (HAS_TAG(col_tags[j], C_TAG_LB_INF))
                    {
                        return UNBNDORINFEAS;
                    }

                    assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
                    assert(!IS_ABS_INF(bounds[j].lb));
                    fix_col(constraints, j, bounds[j].lb, cj);
                    // break from checking if xj is parallel to other cols since
                    // we fix it
                    break;
                }
                else if (fix_xj_to_upper)
                {
                    if (HAS_TAG(col_tags[j], C_TAG_UB_INF))
                    {
                        return UNBNDORINFEAS;
                    }

                    assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
                    assert(!IS_ABS_INF(bounds[j].ub));
                    fix_col(constraints, j, bounds[j].ub, cj);
                    // break from checking if xj is parallel to other cols since
                    // we fix it
                    break;
                }
            }
        }

        // we might have to recount n_inf_min and n_inf_max for the rows column
        // j appears in due to bound change
        if (recount_ninfs)
        {
            const Matrix *A = constraints->A;

            for (jj = 0; jj < len_j; jj++)
            {
                int row = aj_cols[jj];

                recompute_n_infs(acts + row, A->x + A->p[row].start,
                                 A->i + A->p[row].start,
                                 A->p[row].end - A->p[row].start, col_tags);
            }
        }
    }

    return UNCHANGED;
}

static PresolveStatus process_all_bins(const Problem *prob, const int *parallel_cols,
                                       const iVec *groups)
{
    int n_groups = groups->len - 1;
    PresolveStatus status = UNCHANGED;

    for (int i = 0; i < n_groups; ++i)
    {
        int n_cols_this_group = groups->data[i + 1] - groups->data[i];
        status |= process_single_bin(prob, parallel_cols + groups->data[i],
                                     n_cols_this_group);
    }

    return status;
}

PresolveStatus remove_parallel_cols(Problem *prob)
{
    assert(prob->constraints->state->ston_rows->len == 0);
    assert(prob->constraints->state->empty_rows->len == 0);
    assert(prob->constraints->state->empty_cols->len == 0);
    DEBUG(verify_problem_up_to_date(prob->constraints));

    Constraints *constraints = prob->constraints;
    int *parallel_cols = constraints->state->work->iwork_n_cols;
    int *sparsity_IDs = constraints->state->work->iwork1_max_nrows_ncols;
    int *coeff_hashes = constraints->state->work->iwork2_max_nrows_ncols;
    iVec *group_starts = constraints->state->work->int_vec;

    // finding parallel rows of AT is the same as finding parallel cols of A
    find_parallel_rows(constraints->AT, constraints->col_tags, group_starts,
                       parallel_cols, sparsity_IDs, coeff_hashes, C_TAG_INACTIVE);
#ifndef NDEBUG
    // there should be no inactive cols in group_starts here
    for (int i = 0; i < group_starts->len - 1; ++i)
    {
        for (int j = group_starts->data[i]; j < group_starts->data[i + 1]; ++j)
        {
            assert(
                !HAS_TAG(constraints->col_tags[parallel_cols[j]], C_TAG_INACTIVE));
            assert(constraints->state->col_sizes[parallel_cols[j]] !=
                   SIZE_INACTIVE_COL);
        }
    }
#endif

    PresolveStatus status = process_all_bins(prob, parallel_cols, group_starts);
    assert(constraints->state->empty_rows->len == 0);

    delete_fixed_cols_from_problem(prob);
    delete_inactive_cols_from_A_and_AT(prob->constraints);

    DEBUG(verify_problem_up_to_date(constraints));

    return status;
}
