/*
 * Copyright 2025-2026 Daniel Cederberg
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

#include "Bounds.h"
#include "Chains.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "DtonsEq_internal.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "Postsolver.h"
#include "Problem.h"
#include "State.h"
#include "Tags.h"
#include <limits.h>

/* Decides which of the two entries of a doubleton equality row to
   substitute. Returns the index (0 or 1). */
int dton_choose_subst(const double *row_vals, int col_size0, int col_size1)
{
    double a_abs = ABS(row_vals[0]);
    double b_abs = ABS(row_vals[1]);

    /* choose column with larger absolute value */
    if (a_abs != b_abs)
    {
        return (a_abs > b_abs) ? 0 : 1;
    }

    /* equal magnitude: the sparser column */
    return (col_size0 < col_size1) ? 0 : 1;
}

/* Reserves the round's per-record arrays for 'need' records. False if an
   allocation fails. */
static inline bool dton_reserve_records(DtonWorkspace *ws, size_t need)
{
    if (need <= (size_t) ws->substs.cap)
    {
        return true;
    }
    assert(need <= (size_t) INT_MAX);
    bool ok = true;
    ok = ps_grow(&ws->substs.recs, need, sizeof(DtonSubst)) && ok;
    ok = ps_grow(&ws->substs.order, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->substs.succ, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->substs.drop_priority, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->substs.depth, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->substs.stamp, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->targets.list, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->targets.old_size, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->log.start, need + 1, sizeof(int)) && ok;
    if (!ok)
    {
        return false;
    }
    ws->substs.cap = (int) need;
    return true;
}

/* Phase 1, no mutation: walks state->dton_rows and claims the substituted
   column of every eliminable row into ws. A row whose chosen column is
   already claimed is appended to 'deferred' and retried next round. Stale
   worklist entries are dropped. False if the round's records cannot be
   reserved (nothing is claimed then). */
bool dton_claim(Problem *prob, DtonWorkspace *ws, int *deferred, int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const double *rhs = constraints->rhs;
    const int *row_sizes = constraints->state->row_sizes;
    const int *col_sizes = constraints->state->col_sizes;
    const iVec *dton_rows = constraints->state->dton_rows;

    ws->substs.n = 0;

    int i, ii, col0, col1, subst, stay, k;

    // every worklist row claims at most one record, and the targets and log
    // segments are bounded by the records
    if (!dton_reserve_records(ws, (size_t) dton_rows->len))
    {
        return false;
    }

    for (ii = 0; ii < dton_rows->len; ++ii)
    {
        i = dton_rows->data[ii];

        /* a row that used to be dton might have been modified */
        if (row_sizes[i] < 2)
        {
            continue;
        }

        assert(row_sizes[i] == 2 && A->p[i].end - A->p[i].start == 2 &&
               constraints->lhs[i] == rhs[i] && !IS_ABS_INF(rhs[i]) &&
               HAS_TAG(constraints->row_tags[i], R_TAG_EQ) &&
               !HAS_TAG(constraints->row_tags[i], R_TAG_INACTIVE));

        const int *cols = A->i + A->p[i].start;
        const double *vals = A->x + A->p[i].start;
        col0 = cols[0];
        col1 = cols[1];

        assert(!HAS_TAG(constraints->col_tags[col0] | constraints->col_tags[col1],
                        C_TAG_INACTIVE));

        subst = dton_choose_subst(vals, col_sizes[col0], col_sizes[col1]);
        stay = 1 - subst;

        k = cols[subst];

        /* the column this dton row wants to substitute has already been claimed
         * by another dtonrow in this round */
        if (ws->substs.col_subst[k] >= 0)
        {
            deferred[(*n_deferred)++] = i;
            continue;
        }

        assert(ws->substs.n < ws->substs.cap);
        DtonSubst *rec = ws->substs.recs + ws->substs.n;
        rec->k = k;
        rec->owner = i;
        rec->j = cols[stay];
        rec->aik = vals[subst];
        rec->aij = vals[stay];
        rec->dir_mult = -rec->aij / rec->aik;
        rec->dir_shift = rhs[i] / rec->aik;
        ws->substs.col_subst[k] = ws->substs.n++;
    }
    return true;
}

/* Phase 2, no mutation: composes the per-link maps into maps onto the round's
   survivors. A cycle is broken by un-eliminating the record with the largest
   owner row on it, and that owner row is appended to 'deferred'. */
void dton_compose(DtonWorkspace *ws, int *deferred, int *n_deferred)
{
    DtonSubst *recs = ws->substs.recs;
    int *col_subst = ws->substs.col_subst;
    int *succ = ws->substs.succ;
    int *depth = ws->substs.depth;
    int *order = ws->substs.order;
    int n = ws->substs.n;

    /* set drop priority for cycle breaking to the owner row of each record */
    for (int idx = 0; idx < n; ++idx)
    {
        succ[idx] = col_subst[recs[idx].j]; // -1 when the stay column survives
        ws->substs.drop_priority[idx] = recs[idx].owner;
    }

    /* break cycles */
    int *dropped = deferred + *n_deferred;
    int n_dropped =
        compute_chain_depths(n, succ, ws->substs.drop_priority, depth, order,
                             dropped, ws->substs.stamp, ws->acc.touched);

    /* update the deferred list with the dropped records */
    for (int ii = 0; ii < n_dropped; ++ii)
    {
        int idx = dropped[ii];
        col_subst[recs[idx].k] = -1;
        dropped[ii] = recs[idx].owner;
    }
    *n_deferred += n_dropped;

    /* map every eliminated col onto the surviving column at the end of its chain */
    int n_kept = n - n_dropped;
    for (int ii = n_kept - 1; ii >= 0; --ii)
    {
        int idx = order[ii];
        DtonSubst *rec = recs + idx;
        if (depth[idx] == 0)
        {
            rec->target = rec->j;
            rec->mult = rec->dir_mult;
            rec->shift = rec->dir_shift;
        }
        else
        {
            const DtonSubst *child = recs + succ[idx];
            rec->target = child->target;
            rec->mult = rec->dir_mult * child->mult;
            rec->shift = rec->dir_mult * child->shift + rec->dir_shift;
        }
    }

    /* remove dropped records, keep claim order, and renumber col_subst and order */
    int *map = ws->substs.stamp;
    int out = 0;
    for (int idx = 0; idx < n; ++idx)
    {
        if (depth[idx] == CHAINS_DROPPED)
        {
            continue;
        }
        recs[out] = recs[idx];
        depth[out] = depth[idx];
        col_subst[recs[out].k] = out;
        map[idx] = out++;
    }
    ws->substs.n = out;
    for (int ii = 0; ii < n_kept; ++ii)
    {
        order[ii] = map[order[ii]];
    }
}

/* Transfers the bounds of x_subst onto x_stay through the doubleton equality
   row i: aij * x_stay + aik * x_subst = rhs. Returns INFEASIBLE when a
   transferred bound contradicts the bounds of x_stay, UNCHANGED otherwise. */
static PresolveStatus dton_transfer_bounds_link(Constraints *constraints, int i,
                                                double aij, double aik, double rhs,
                                                double lb_subst, double ub_subst,
                                                int stay, ColTag col_tag_subst)
{
    assert(aij != 0 && aik != 0);
    bool same_sign = (aik * aij > 0.0);
    PresolveStatus status = UNCHANGED;

    if (same_sign)
    {
        if (!HAS_TAG(col_tag_subst, C_TAG_LB_INF))
        {
            double new_ub_cand = (rhs - aik * lb_subst) / aij;
            status =
                update_ub(constraints, stay, new_ub_cand, i HUGE_BOUND_IS_NOT_OK);
            RETURN_IF_INFEASIBLE(status);
        }

        if (!HAS_TAG(col_tag_subst, C_TAG_UB_INF))
        {
            double new_lb_cand = (rhs - aik * ub_subst) / aij;
            status =
                update_lb(constraints, stay, new_lb_cand, i HUGE_BOUND_IS_NOT_OK);
        }
    }
    else
    {
        if (!HAS_TAG(col_tag_subst, C_TAG_LB_INF))
        {
            double new_lb_cand = (rhs - aik * lb_subst) / aij;
            status =
                update_lb(constraints, stay, new_lb_cand, i HUGE_BOUND_IS_NOT_OK);
            RETURN_IF_INFEASIBLE(status);
        }

        if (!HAS_TAG(col_tag_subst, C_TAG_UB_INF))
        {
            double new_ub_cand = (rhs - aik * ub_subst) / aij;
            status =
                update_ub(constraints, stay, new_ub_cand, i HUGE_BOUND_IS_NOT_OK);
        }
    }
    return status;
}

/* Phase 3: transfers the bounds of each eliminated column onto its direct stay
   column, one link at a time in descending chain depth, so each transfer sees
   the bounds tightened by the deeper links. A composed transfer would break
   dual complementarity on chains. Returns INFEASIBLE at the first
   contradicting transfer, UNCHANGED otherwise. */
PresolveStatus dton_transfer_bounds(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    const double *rhs = constraints->rhs;
    const Bound *bounds = constraints->bounds;
    const ColTag *col_tags = constraints->col_tags;
    const DtonSubst *recs = ws->substs.recs;
    const int *order = ws->substs.order;

    for (int ii = 0; ii < ws->substs.n; ++ii)
    {
        const DtonSubst *rec = recs + order[ii];
        int i = rec->owner;
        int k = rec->k;
        int j = rec->j;
        assert(ii == 0 ||
               ws->substs.depth[order[ii - 1]] >= ws->substs.depth[order[ii]]);

        PresolveStatus status =
            dton_transfer_bounds_link(constraints, i, rec->aij, rec->aik, rhs[i],
                                      bounds[k].lb, bounds[k].ub, j, col_tags[k]);
        RETURN_IF_INFEASIBLE(status);

        // update_lb/update_ub cannot fix or deactivate columns
        assert(!HAS_TAG(col_tags[j], C_TAG_INACTIVE));
        assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
    }
    return UNCHANGED;
}

/* Phase 3b: postsolve records, three per eliminated column k with owner row i
   and survivor s: ADDED_ROWS(i, column k, aik), SUB_COL_DTON(k, s, mult,
   shift) and DELETED_ROW(i, c_k / aik). Replayed in reverse they give x_k and
       y_i = (c_k - sum_{r != i} a_rk y_r) / a_ik,
   where the sum includes the owner rows whose stay column is k. Those records
   are deeper, so the columns are emitted by ascending depth and replayed
   deepest first. */
void dton_record(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT; // pre-round transpose
    const double *c = prob->obj->c;     // pre-round objective
    PostsolveInfo *info = constraints->state->postsolve_info;

    // ascending depth (the order lists the records by descending depth)
    for (int ii = ws->substs.n - 1; ii >= 0; --ii)
    {
        const DtonSubst *rec = ws->substs.recs + ws->substs.order[ii];
        int k = rec->k;
        int i = rec->owner;
        double aik = rec->aik;

        const int *rows = AT->i + AT->p[k].start;
        const double *col_vals = AT->x + AT->p[k].start;
        size_t len = (size_t) (AT->p[k].end - AT->p[k].start);

#ifndef NDEBUG
        // every row of the pre-round column is active
        for (size_t jj = 0; jj < len; ++jj)
        {
            assert(!HAS_TAG(constraints->row_tags[rows[jj]], R_TAG_INACTIVE));
        }
#endif

        save_retrieval_added_rows(info, i, rows, col_vals, len, aik);
        save_retrieval_sub_col_dton(info, k, rec->target, rec->mult, rec->shift);
        save_retrieval_deleted_row(info, i, c[k] / aik);
    }
}
