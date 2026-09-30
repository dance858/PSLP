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

/* Decides which of the two entries of a doubleton equality row to
   substitute. Returns the index (0 or 1). */
int dton_choose_subst(const double *row_vals, int col_size0, int col_size1)
{
    double a_abs = ABS(row_vals[0]);
    double b_abs = ABS(row_vals[1]);

    /* Choose the column with the larger absolute value. */
    if (a_abs != b_abs)
    {
        return (a_abs > b_abs) ? 0 : 1;
    }

    /* Equal magnitude: the sparser column. */
    return (col_size0 < col_size1) ? 0 : 1;
}

/* Reserves the round's per-record arrays for 'need' records. False if an
   allocation fails. */
static inline bool dton_reserve_records(DtonWorkspace *dton_work, size_t need)
{
    if (need <= (size_t) dton_work->substs.cap)
    {
        return true;
    }
    bool ok = true;
    ok &= ps_grow(&dton_work->substs.recs, need, sizeof(DtonSubst));
    ok &= ps_grow(&dton_work->substs.order, need, sizeof(int));
    ok &= ps_grow(&dton_work->substs.succ, need, sizeof(int));
    ok &= ps_grow(&dton_work->substs.drop_priority, need, sizeof(int));
    ok &= ps_grow(&dton_work->substs.depth, need, sizeof(int));
    ok &= ps_grow(&dton_work->substs.stamp, need, sizeof(int));
    ok &= ps_grow(&dton_work->targets.list, need, sizeof(DtonTarget));
    if (!ok)
    {
        return false;
    }
    dton_work->substs.cap = (int) need;
    return true;
}

/* For every doubleton equality row, choose one column to substitute. If the chosen
   column is already claimed by another row, the row is deferred to the next round */
bool dton_claim(Problem *prob, DtonWorkspace *dton_work, int *deferred,
                int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const double *rhs = constraints->rhs;
    const int *row_sizes = constraints->state->row_sizes;
    const int *col_sizes = constraints->state->col_sizes;
    const iVec *dton_rows = constraints->state->dton_rows;

    dton_work->substs.n_recs = 0;

    int i, ii, col0, col1, subst, stay, k;

    /* Every worklist row claims at most one record, and the targets are
       bounded by the records. */
    if (!dton_reserve_records(dton_work, (size_t) dton_rows->len))
    {
        return false;
    }

    for (ii = 0; ii < dton_rows->len; ++ii)
    {
        i = dton_rows->data[ii];

        /* a row that used to be a doubleton might have been modified */
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

        /* the col this row wants to substitute is already claimed by another row */
        if (dton_work->substs.col_subst[k] >= 0)
        {
            deferred[(*n_deferred)++] = i;
            continue;
        }

        /* store information about the substitution */
        assert(dton_work->substs.n_recs < dton_work->substs.cap);
        DtonSubst *rec = dton_work->substs.recs + dton_work->substs.n_recs;
        rec->k = k;
        rec->owner = i;
        rec->j = cols[stay];
        rec->aik = vals[subst];
        rec->aij = vals[stay];
        rec->dir_mult = -rec->aij / rec->aik;
        rec->dir_shift = rhs[i] / rec->aik;
        dton_work->substs.col_subst[k] = dton_work->substs.n_recs++;
    }
    return true;
}

/* Map substituted columns onto the round's survivors. Cycles are broken by
   ignoring the substitution done by the row with largest index. The row
   for the ignored substitution is appended to 'deferred'. */
void dton_compose(DtonWorkspace *dton_work, int *deferred, int *n_deferred)
{
    DtonSubst *recs = dton_work->substs.recs;
    int *col_subst = dton_work->substs.col_subst;
    int *succ = dton_work->substs.succ;
    int *depth = dton_work->substs.depth;
    int *order = dton_work->substs.order;
    int *drop_priority = dton_work->substs.drop_priority;
    int *stamp = dton_work->substs.stamp;
    int n_recs = dton_work->substs.n_recs;

    /* succ[i] is the index of the substitution that eliminates the stay column
       of substitution i, or -1 if no substitution eliminates the stay column */
    for (int i = 0; i < n_recs; ++i)
    {
        succ[i] = col_subst[recs[i].j];
        drop_priority[i] = recs[i].owner;
    }

    /* break cycles by ignoring the substitution done by the row with largest idx */
    int *dropped = deferred + *n_deferred;
    int n_dropped = compute_chain_depths(n_recs, succ, drop_priority, depth, order,
                                         dropped, stamp, dton_work->acc.touched);

    // -----------------------------------------------------------------------------
    //  map every eliminated col onto the surviving col at the end of its chain
    // -----------------------------------------------------------------------------
    int n_kept = n_recs - n_dropped;
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

    // -----------------------------------------------------------------------
    // reset ignored substitutions and remove the ignored records
    // -----------------------------------------------------------------------
    for (int ii = 0; ii < n_dropped; ++ii)
    {
        int idx = dropped[ii];
        col_subst[recs[idx].k] = -1;
        dropped[ii] = recs[idx].owner;
    }
    *n_deferred += n_dropped;

    int *new_index = stamp;
    int out = 0;
    for (int idx = 0; idx < n_recs; ++idx)
    {
        if (depth[idx] == CHAINS_DROPPED)
        {
            continue;
        }
        recs[out] = recs[idx];
        depth[out] = depth[idx];
        col_subst[recs[out].k] = out;
        new_index[idx] = out++;
    }
    dton_work->substs.n_recs = out;
    for (int ii = 0; ii < n_kept; ++ii)
    {
        order[ii] = new_index[order[ii]];
    }
}

/* Transfers the bounds of x_k onto x_j through the doubleton equality row i:
   aij * x_j + aik * x_k = rhs. Returns INFEASIBLE when a transferred bound
   contradicts the bounds of x_j, UNCHANGED otherwise. */
static PresolveStatus dton_transfer_bounds_link(Constraints *constraints, int i,
                                                int j, double aij, double aik,
                                                double rhs, double lb_k, double ub_k,
                                                ColTag col_tag_k)
{
    assert(aij != 0 && aik != 0);
    bool same_sign = (aik * aij > 0.0);
    PresolveStatus status = UNCHANGED;

    if (!HAS_TAG(col_tag_k, C_TAG_LB_INF))
    {
        double bound = (rhs - aik * lb_k) / aij;
        status = same_sign
                     ? update_ub(constraints, j, bound, i HUGE_BOUND_IS_NOT_OK)
                     : update_lb(constraints, j, bound, i HUGE_BOUND_IS_NOT_OK);
        RETURN_IF_INFEASIBLE(status);
    }

    if (!HAS_TAG(col_tag_k, C_TAG_UB_INF))
    {
        double bound = (rhs - aik * ub_k) / aij;
        status = same_sign
                     ? update_lb(constraints, j, bound, i HUGE_BOUND_IS_NOT_OK)
                     : update_ub(constraints, j, bound, i HUGE_BOUND_IS_NOT_OK);
    }
    return status;
}

/* Transfers the bounds of each eliminated column onto its direct stay
   column, one link at a time in descending chain depth, so each transfer sees
   the bounds tightened by the deeper links. A composed transfer would break
   dual complementarity on chains. Returns INFEASIBLE at the first
   contradicting transfer, UNCHANGED otherwise. */
PresolveStatus dton_transfer_bounds(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const double *rhs = constraints->rhs;
    const Bound *bounds = constraints->bounds;
    const ColTag *col_tags = constraints->col_tags;
    const DtonSubst *recs = dton_work->substs.recs;
    const int *order = dton_work->substs.order;

    for (int ii = 0; ii < dton_work->substs.n_recs; ++ii)
    {
        const DtonSubst *rec = recs + order[ii];
        int i = rec->owner;
        int k = rec->k;
        int j = rec->j;
        assert(ii == 0 || dton_work->substs.depth[order[ii - 1]] >=
                              dton_work->substs.depth[order[ii]]);

        PresolveStatus status =
            dton_transfer_bounds_link(constraints, i, j, rec->aij, rec->aik, rhs[i],
                                      bounds[k].lb, bounds[k].ub, col_tags[k]);
        RETURN_IF_INFEASIBLE(status);
    }
    return UNCHANGED;
}

/* Postsolve records, three per eliminated column k with owner row i
   and survivor s: ADDED_ROWS(i, column k, aik), SUB_COL_DTON(k, s, mult,
   shift) and DELETED_ROW(i, c_k / aik). Replayed in reverse they give x_k and
       y_i = (c_k - sum_{r != i} a_rk y_r) / a_ik,
   where the sum includes the owner rows whose stay column is k. Those records
   are deeper, so the columns are emitted by ascending depth and replayed
   deepest first. */
void dton_record(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT; /* pre-round transpose */
    const double *c = prob->obj->c;     /* pre-round objective */
    PostsolveInfo *info = constraints->state->postsolve_info;

    /* ascending depth: order lists the substitutions by descending depth */
    for (int ii = dton_work->substs.n_recs - 1; ii >= 0; --ii)
    {
        const DtonSubst *rec = dton_work->substs.recs + dton_work->substs.order[ii];
        int k = rec->k;
        int i = rec->owner;
        double aik = rec->aik;

        const int *rows = AT->i + AT->p[k].start;
        const double *col_vals = AT->x + AT->p[k].start;
        int len = AT->p[k].end - AT->p[k].start;

#ifndef NDEBUG
        /* Every row of the pre-round column is active. */
        for (int jj = 0; jj < len; ++jj)
        {
            assert(!HAS_TAG(constraints->row_tags[rows[jj]], R_TAG_INACTIVE));
        }
#endif

        save_retrieval_added_rows(info, i, rows, col_vals, (size_t) len, aik);
        save_retrieval_sub_col_dton(info, k, rec->target, rec->mult, rec->shift);
        save_retrieval_deleted_row(info, i, c[k] / aik);
    }
}
