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

#include "Activity.h"
#include "Binary_search.h"
#include "Bounds.h"
#include "Chains.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "DtonsEq_internal.h"
#include "Locks.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "Postsolver.h"
#include "Problem.h"
#include "RowColViews.h"
#include "SimpleReductions.h"
#include "State.h"
#include "Tags.h"
#include "Workspace.h"
#include "radix_sort.h"
#include <limits.h>

/* accumulator flag bits: the row originally had this column / the
   substitution contributed to it */
#define DTON_ACC_EXISTED ((uint8_t) 1)
#define DTON_ACC_MODIFIED ((uint8_t) 2)

DtonWorkspace *dton_ws_new(size_t n_rows, size_t n_cols)
{
    DtonWorkspace *ws = (DtonWorkspace *) ps_calloc(1, sizeof(DtonWorkspace));
    RETURN_PTR_IF_NULL(ws, NULL);
    ws->m = n_rows;
    ws->n = n_cols;

    // the record and row arrays are reserved per round (dton_claim,
    // dton_reserve_rows) and the borrowed arrays are set by dton_ws_attach
    ws->acc.value = (double *) ps_malloc(n_cols, sizeof(double));
    ws->acc.flags = (uint8_t *) ps_malloc(n_cols, sizeof(uint8_t));
    ws->at.cap = (int *) ps_malloc(n_cols, sizeof(int));
    ws->log.alloc = 1024;
    ws->log.col = (int *) ps_malloc((size_t) ws->log.alloc, sizeof(int));
    ws->log.row = (int *) ps_malloc((size_t) ws->log.alloc, sizeof(int));
    ws->log.val = (double *) ps_malloc((size_t) ws->log.alloc, sizeof(double));
    ws->log.row2 = (int *) ps_malloc((size_t) ws->log.alloc, sizeof(int));
    ws->log.val2 = (double *) ps_malloc((size_t) ws->log.alloc, sizeof(double));

    if (!ws->acc.value || !ws->acc.flags || !ws->at.cap || !ws->log.col ||
        !ws->log.row || !ws->log.val || !ws->log.row2 || !ws->log.val2)
    {
        dton_ws_free(ws);
        return NULL;
    }

    ws->tuning.rebuild_dirty_frac = 0.25;
    return ws;
}

void dton_ws_free(DtonWorkspace *ws)
{
    RETURN_IF_NULL(ws);
    PS_FREE(ws->substs.recs);
    PS_FREE(ws->substs.order);
    PS_FREE(ws->substs.succ);
    PS_FREE(ws->substs.drop_priority);
    PS_FREE(ws->substs.depth);
    PS_FREE(ws->substs.stamp);
    PS_FREE(ws->targets.list);
    PS_FREE(ws->targets.old_size);
    PS_FREE(ws->log.start);
    PS_FREE(ws->acc.value);
    PS_FREE(ws->acc.flags);
    PS_FREE(ws->rows.list);
    PS_FREE(ws->rows.idx);
    PS_FREE(ws->rows.aux);
    PS_FREE(ws->at.cap);
    PS_FREE(ws->log.col);
    PS_FREE(ws->log.row);
    PS_FREE(ws->log.val);
    PS_FREE(ws->log.row2);
    PS_FREE(ws->log.val2);
    PS_FREE(ws);
}

/* Points the borrowed per-column arrays at the presolver's scratch and
   initializes them. */
void dton_ws_attach(DtonWorkspace *ws, Work *work)
{
    ws->substs.col_subst = work->iwork1_max_nrows_ncols;
    ws->targets.col_to_target = work->iwork2_max_nrows_ncols;
    for (size_t k = 0; k < ws->n; ++k)
    {
        ws->substs.col_subst[k] = -1;
        ws->targets.col_to_target[k] = -1;
    }
    sparse_accumulator_init(&ws->acc, work->radix_aux, ws->acc.value, ws->acc.flags,
                            work->iwork_n_cols, ws->n);
}

/* Decides which of the two entries of a doubleton equality row to
   substitute: the larger coefficient, so that neither the substitution nor
   the postsolve divides by a tiny pivot. Equal magnitudes: singleton column
   first, then the sparser column. Returns the entry position (0 or 1). */
int dton_choose_subst(const double *row_vals, int col_size0, int col_size1)
{
    double a_abs = ABS(row_vals[0]);
    double b_abs = ABS(row_vals[1]);

    // the substitution and the postsolve divide by the substituted
    // coefficient, so it is the larger one; the fill-in multiplier is then
    // at most 1 in magnitude
    if (a_abs != b_abs)
    {
        return (a_abs > b_abs) ? 0 : 1;
    }
    // equal magnitude: a singleton column costs no fill-in, then the sparser
    if (col_size0 == 1 && col_size1 != 1)
    {
        return 0;
    }
    if (col_size0 != 1 && col_size1 == 1)
    {
        return 1;
    }
    return (col_size0 < col_size1) ? 0 : 1;
}

/* Phase 1, no mutation: walks state->dton_rows and claims the substituted
   column of every eliminable row into ws (the column-indexed map plus a
   packed DtonSubst record per accepted claim). A row whose chosen column is
   already claimed is appended to 'deferred' and retried next round; stale
   worklist entries are dropped. False if the round's records cannot be
   reserved (nothing is claimed then). */
bool dton_claim(Problem *prob, DtonWorkspace *ws, int *deferred, int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const double *rhs = constraints->rhs;
    const ColTag *col_tags = constraints->col_tags;
    const int *row_sizes = constraints->state->row_sizes;
    const int *col_sizes = constraints->state->col_sizes;
    const iVec *dton_rows = constraints->state->dton_rows;

    ws->substs.n = 0;

    // exact reservation for the round: every worklist row claims at most one
    // record, and the targets and log segments are bounded by the records
    size_t need = (size_t) dton_rows->len;
    if (need > (size_t) ws->substs.cap)
    {
        if (need > (size_t) INT_MAX)
        {
            return false;
        }
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
    }

    for (int ii = 0; ii < dton_rows->len; ++ii)
    {
        int i = dton_rows->data[ii];

        // stale entries (row shrunk or deactivated since it was listed)
        if (row_sizes[i] < 2)
        {
            continue;
        }

        assert(row_sizes[i] == 2);
        assert(HAS_TAG(constraints->row_tags[i], R_TAG_EQ) &&
               !HAS_TAG(constraints->row_tags[i], R_TAG_INACTIVE));
        assert(constraints->lhs[i] == rhs[i] && !IS_ABS_INF(rhs[i]));
        assert(A->p[i].end - A->p[i].start == 2);

        const int *cols = A->i + A->p[i].start;
        const double *vals = A->x + A->p[i].start;
        int c0 = cols[0];
        int c1 = cols[1];

        // a column may have been deactivated since the row was listed
        if (HAS_TAG((ColTag) (col_tags[c0] | col_tags[c1]), C_TAG_INACTIVE))
        {
            continue;
        }

        int slot = dton_choose_subst(vals, col_sizes[c0], col_sizes[c1]);

        int k = cols[slot];

        // claim conflict: the column is already owned by an earlier row
        // this round; the loser stays active and is retried next round
        // (no orientation flip)
        if (ws->substs.col_subst[k] >= 0)
        {
            deferred[(*n_deferred)++] = i;
            continue;
        }

        assert(ws->substs.n < ws->substs.cap);
        DtonSubst *rec = ws->substs.recs + ws->substs.n;
        rec->k = k;
        rec->owner = i;
        rec->j = cols[1 - slot];
        rec->dir_mult = -vals[1 - slot] / vals[slot];
        rec->dir_shift = rhs[i] / vals[slot];
        ws->substs.col_subst[k] = ws->substs.n++;
    }
    return true;
}

/* Phase 2, no mutation: composes the per-link maps into maps onto the
   round's final survivors and computes chain depths. The substitutions form a
   functional graph, record -> record of its stay column, that compute_chain_depths
   turns into depths; a cycle is broken by un-eliminating the record with the
   largest owner row on it, and that owner row is appended to 'deferred'. On
   return ws->substs.recs holds exactly the still-eliminated substitutions in
   claim order, ws->substs.depth their chain depths, ws->substs.maxdepth the
   round's maximum chain depth, and ws->substs.order the records by descending
   depth. */
void dton_compose(DtonWorkspace *ws, int *deferred, int *n_deferred)
{
    DtonSubst *recs = ws->substs.recs;
    int *col_subst = ws->substs.col_subst;
    int *succ = ws->substs.succ;
    int *depth = ws->substs.depth;
    int *order = ws->substs.order;
    int n = ws->substs.n;

    for (int idx = 0; idx < n; ++idx)
    {
        succ[idx] = col_subst[recs[idx].j]; // -1 when the stay column survives
        ws->substs.drop_priority[idx] = recs[idx].owner;
    }

    // the dropped records go straight to the deferred list, as owner rows;
    // the accumulator's touched list is free to serve as the walk path
    int *dropped = deferred + *n_deferred;
    int n_dropped =
        compute_chain_depths(n, succ, ws->substs.drop_priority, depth, order,
                             dropped, ws->substs.stamp, ws->acc.touched);
    for (int d = 0; d < n_dropped; ++d)
    {
        int idx = dropped[d];
        col_subst[recs[idx].k] = -1;
        dropped[d] = recs[idx].owner;
    }
    *n_deferred += n_dropped;

    // compose along the chains, shallow links first, so every record composes
    // with an already-final map: x_k = dir_mult * x_child + dir_shift with
    // x_child = mult_child * x_survivor + shift_child
    int n_kept = n - n_dropped;
    for (int o = n_kept - 1; o >= 0; --o)
    {
        int idx = order[o];
        DtonSubst *rec = recs + idx;
        if (depth[idx] == 0)
        {
            rec->target = rec->j; // a survivor, or an un-eliminated column
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

    // compact the records in claim order; col_subst and order follow them
    // (stamp is free now and maps old to new indices)
    int *map = ws->substs.stamp;
    int out = 0;
    ws->substs.maxdepth = 0;
    for (int idx = 0; idx < n; ++idx)
    {
        if (depth[idx] == CHAINS_DROPPED)
        {
            map[idx] = -1;
            continue;
        }
        recs[out] = recs[idx];
        depth[out] = depth[idx];
        col_subst[recs[out].k] = out;
        ws->substs.maxdepth = MAX(ws->substs.maxdepth, depth[out]);
        map[idx] = out++;
    }
    ws->substs.n = out;
    for (int o = 0; o < n_kept; ++o)
    {
        order[o] = map[order[o]];
    }
}

/* Transfers the bounds of the substituted variable onto the variable that
   stays, for the doubleton equality row i: aij * x_stay + aik * x_subst = rhs.
   lb_subst and ub_subst are the bounds of the substituted variable. Emits the
   bound-change postsolve records via update_lb/update_ub. Returns INFEASIBLE
   when a transferred bound contradicts the stay column's bounds (the row and
   the two columns' bounds are then inconsistent), UNCHANGED otherwise. */
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

            if (HAS_TAG(constraints->col_tags[stay], C_TAG_FIXED))
            {
                return UNCHANGED;
            }
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

            if (HAS_TAG(constraints->col_tags[stay], C_TAG_FIXED))
            {
                return UNCHANGED;
            }
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

/* Phase 3: transfers the bounds of each still-eliminated column onto its
   DIRECT partner via dton_transfer_bounds_link, one link at a time in
   DESCENDING chain depth, so a link's transfer sees the bounds already
   tightened by the deeper links (a composed shortcut breaks dual
   complementarity on chains). Runs entirely before any matrix mutation,
   so the owner rows and the substituted columns' bounds and tags are read
   live. Returns INFEASIBLE at the first contradicting transfer (the
   tightenings made so far are valid and stay), UNCHANGED otherwise. */
static PresolveStatus dton_transfer_bounds(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const double *rhs = constraints->rhs;
    const Bound *bounds = constraints->bounds;
    const ColTag *col_tags = constraints->col_tags;
    const DtonSubst *recs = ws->substs.recs;
    const int *order = ws->substs.order;

    // descending depth so a link's transfer sees the bounds and tags already
    // tightened by the deeper links; hence the substituted column's state is
    // read live
    for (int o = 0; o < ws->substs.n; ++o)
    {
        const DtonSubst *rec = recs + order[o];
        int i = rec->owner;
        int k = rec->k;
        int j = rec->j;
        assert(o == 0 ||
               ws->substs.depth[order[o - 1]] >= ws->substs.depth[order[o]]);
        assert(A->p[i].end - A->p[i].start == 2);

        // the owner row is untouched until the apply sweep
        const int *cols = A->i + A->p[i].start;
        const double *vals = A->x + A->p[i].start;
        int slot = (cols[0] == k) ? 0 : 1;
        assert(cols[slot] == k && cols[1 - slot] == j);
        double aik = vals[slot];
        double aij = vals[1 - slot];

        PresolveStatus status =
            dton_transfer_bounds_link(constraints, i, aij, aik, rhs[i], bounds[k].lb,
                                      bounds[k].ub, j, col_tags[k]);
        RETURN_IF_INFEASIBLE(status);

        // update_lb/update_ub cannot fix or deactivate columns, so the claim
        // set stays valid; if column fixing is ever added to them, the
        // eliminator needs an un-elimination cascade
        assert(!HAS_TAG(col_tags[j], C_TAG_INACTIVE));
        assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
    }
    return UNCHANGED;
}

/* Phase 3b: postsolve records, three per still-eliminated column k with
   owner row i and survivor s: ADDED_ROWS(i, pre-round column k, aik),
   SUB_COL_DTON(k, s, mult, shift), DELETED_ROW(i, ck / aik). Replayed in
   reverse they give x_k from the survivor and
       y_i = (c_k - sum_{r != i} a_rk y_r) / a_ik,
   stationarity for x_k written with the column and c_k of the round start.
   The rows r include the owner rows of this round whose stay column is k;
   those are strictly deeper, so the columns are emitted in ascending depth
   and replayed deepest first. Runs before any mutation (the columns are
   read from the pre-round AT, c_k before the objective update). */
void dton_record(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT; // pre-round transpose
    const double *c = prob->obj->c;     // pre-round objective
    PostsolveInfo *info = constraints->state->postsolve_info;

    // ascending depth (the order lists the records by descending depth)
    for (int o = ws->substs.n - 1; o >= 0; --o)
    {
        const DtonSubst *rec = ws->substs.recs + ws->substs.order[o];
        int k = rec->k;
        int i = rec->owner;

        // a_ik from the owner row, still intact
        const int *cols = A->i + A->p[i].start;
        const double *vals = A->x + A->p[i].start;
        assert(A->p[i].end - A->p[i].start == 2);
        assert(cols[0] == k || cols[1] == k);
        double aik = (cols[0] == k) ? vals[0] : vals[1];

        const int *rows = AT->i + AT->p[k].start;
        const double *col_vals = AT->x + AT->p[k].start;
        size_t len = (size_t) (AT->p[k].end - AT->p[k].start);

#ifndef NDEBUG
        // every row of the pre-round column is active: explorers flush their
        // deleted rows from AT before returning and this round has not
        // deactivated its owner rows yet
        for (size_t p = 0; p < len; ++p)
        {
            assert(!HAS_TAG(constraints->row_tags[rows[p]], R_TAG_INACTIVE));
        }
#endif

        save_retrieval_added_rows(info, i, rows, col_vals, len, aik);
        save_retrieval_sub_col_dton(info, k, rec->target, rec->mult, rec->shift);
        save_retrieval_deleted_row(info, i, c[k] / aik);
    }
}

/* Appends a tuple to the sweep's change log (an upsert, or a delete when
   'val' is zero), doubling the arrays as needed. On allocation failure the
   round is flagged so the refresh falls back to a full rebuild. */
static inline void dton_log_push(DtonWorkspace *ws, int col, int row, double val)
{
    if (ws->log.len == ws->log.alloc)
    {
        size_t cap = (size_t) ws->log.alloc * 2;
        bool ok = ws->log.alloc <= INT_MAX / 2;
        if (ok)
        {
            ok = ps_grow(&ws->log.col, cap, sizeof(int)) && ok;
            ok = ps_grow(&ws->log.row, cap, sizeof(int)) && ok;
            ok = ps_grow(&ws->log.val, cap, sizeof(double)) && ok;
            ok = ps_grow(&ws->log.row2, cap, sizeof(int)) && ok;
            ok = ps_grow(&ws->log.val2, cap, sizeof(double)) && ok;
        }
        if (!ok)
        {
            ws->log.overflow = true;
            return;
        }
        ws->log.alloc = (int) cap;
    }
    ws->log.col[ws->log.len] = col;
    ws->log.row[ws->log.len] = row;
    ws->log.val[ws->log.len] = val;
    ws->log.len++;
}

/* Reserves the row list for the round (one slot per entry of the eliminated
   columns in the pre-round AT). False if the allocation fails; nothing has
   been mutated at that point. */
static bool dton_reserve_rows(Problem *prob, DtonWorkspace *ws)
{
    const Matrix *AT = prob->constraints->AT;
    size_t need = 0;
    for (int idx = 0; idx < ws->substs.n; ++idx)
    {
        int k = ws->substs.recs[idx].k;
        need += (size_t) (AT->p[k].end - AT->p[k].start);
    }
    if (need <= (size_t) ws->rows.cap)
    {
        return true;
    }
    if (need > (size_t) INT_MAX)
    {
        return false;
    }
    bool ok = true;
    ok = ps_grow(&ws->rows.list, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->rows.idx, need, sizeof(int)) && ok;
    ok = ps_grow(&ws->rows.aux, need, sizeof(int)) && ok;
    if (!ok)
    {
        return false;
    }
    ws->rows.cap = (int) need;
    return true;
}

/* Finishes a swept row q whose entries now occupy [start, start + new_len): shifts
   the finite sides by the substituted constants, updates the sizes, and feeds the
   row worklists (the old_len guards keep deferred rows from being pushed twice). */
static void dton_sweep_row_finish(Problem *prob, int q, int old_len, int new_len,
                                  double chg, int *deferred, int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    RowTag *row_tags = constraints->row_tags;
    int start = A->p[q].start;

    assert(new_len <= old_len);
    A->p[q].end = start + new_len;
    constraints->state->row_sizes[q] = new_len;
    A->nnz -= (size_t) (old_len - new_len);
    DEBUG(ASSERT_INCREASING_I(A->i + start, (size_t) new_len););
    DEBUG(ASSERT_NO_ZEROS_D(A->x + start, (size_t) new_len););

    if (chg != 0.0)
    {
        if (!HAS_TAG(row_tags[q], R_TAG_LHS_INF))
        {
            constraints->lhs[q] -= chg;
        }
        if (!HAS_TAG(row_tags[q], R_TAG_RHS_INF))
        {
            constraints->rhs[q] -= chg;
        }
    }

    switch (new_len)
    {
        case 0:
            assert(!iVec_contains(constraints->state->empty_rows, q));
            iVec_append(constraints->state->empty_rows, q);
            assert(!HAS_TAG(row_tags[q], R_TAG_INACTIVE));
            break;
        case 1:
            if (old_len != 1)
            {
                assert(!iVec_contains(constraints->state->ston_rows, q));
                iVec_append(constraints->state->ston_rows, q);
            }
            break;
        case 2:
            if (HAS_TAG(row_tags[q], R_TAG_EQ) && old_len != 2)
            {
                deferred[(*n_deferred)++] = q;
            }
            break;
        default:
            break;
    }
}

/* Rewrites row q in place through the sparse accumulator. Every entry of an
   eliminated column lands on its composed target; untouched entries are kept
   verbatim; a touched value below ZERO_TOL is dropped as a cancellation. The
   activity is recomputed from scratch. */
static void dton_sweep_row(Problem *prob, DtonWorkspace *ws, int q, int *deferred,
                           int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    const ColTag *col_tags = constraints->col_tags;
    int start = A->p[q].start;
    int end = A->p[q].end;
    int old_len = end - start;
    SparseAccumulator *acc = &ws->acc;
    double chg = 0.0;

    sparse_accumulator_clear(acc);
    for (int p = start; p < end; ++p)
    {
        int c = A->i[p];
        double v = A->x[p];
        int ri = ws->substs.col_subst[c];
        if (ri >= 0)
        {
            const DtonSubst *r = ws->substs.recs + ri;
            sparse_accumulator_add(acc, r->target, v * r->mult, DTON_ACC_MODIFIED);
            chg += v * r->shift;
        }
        else
        {
            sparse_accumulator_add(acc, c, v, DTON_ACC_EXISTED);
        }
    }

    // rows must stay column-sorted
    sparse_accumulator_sort(acc);

    // rebuild in place (never grows: every eliminated entry maps onto exactly
    // one target, and merging only shrinks)
    int w = start;
    for (int t = 0; t < acc->n_touched; ++t)
    {
        int oc = acc->touched[t];
        double val = acc->value[oc];
        bool keep =
            (acc->flags[oc] == DTON_ACC_EXISTED) ? true : ABS(val) > ZERO_TOL;
        bool is_target = ws->targets.col_to_target[oc] >= 0;
        if (keep)
        {
            A->i[w] = oc;
            A->x[w] = val;
            ++w;
            if (is_target)
            {
                dton_log_push(ws, oc, q, val);
            }
        }
        else if (is_target && (acc->flags[oc] & DTON_ACC_EXISTED))
        {
            dton_log_push(ws, oc, q, 0.0); // delete
        }
    }
    int new_len = w - start;

    // full activity recompute for the rewritten row; Activity_init resets the
    // status, which is restored. A row listed in updated_activities must have
    // status ADDED and a computable side (verify_row_states), so a NOT_ADDED
    // row is promoted when it gains one.
    if (new_len > 0)
    {
        Activity *act = constraints->state->activities + q;
        uint8_t old_status = act->status;
        Activity_init(act, A->x + start, A->i + start, new_len, constraints->bounds,
                      col_tags);
        act->status = old_status;
        if (act->status == NOT_ADDED && (act->n_inf_min == 0 || act->n_inf_max == 0))
        {
            act->status = ADDED;
            iVec_append(constraints->state->updated_activities, q);
        }
    }

    dton_sweep_row_finish(prob, q, old_len, new_len, chg, deferred, n_deferred);
}

/* Phase 4: applies the round. Deactivates the owner rows, rewrites every
   affected row in place through the composed maps (the ZERO_TOL drop applies
   only to entries the substitution changed), shifts finite row sides,
   updates the affected rows' activities, feeds the row worklists (new
   doubleton-equality candidates go to 'deferred'), updates the objective and
   deactivates the eliminated columns. */
void dton_apply(Problem *prob, DtonWorkspace *ws, int *deferred, int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT; // pre-round transpose, read-only
    RowTag *row_tags = constraints->row_tags;
    ColTag *col_tags = constraints->col_tags;
    int *row_sizes = constraints->state->row_sizes;
    int *col_sizes = constraints->state->col_sizes;
    Objective *obj = prob->obj;

    // the composed targets: the only active columns whose AT content changes
    // this round. Their pre-round sizes feed the size transitions of phase 5.
    ws->targets.n = 0;
    ws->log.len = 0;
    ws->log.overflow = false;
    for (int idx = 0; idx < ws->substs.n; ++idx)
    {
        int T = ws->substs.recs[idx].target;
        if (ws->targets.col_to_target[T] < 0)
        {
            ws->targets.col_to_target[T] = ws->targets.n;
            ws->targets.list[ws->targets.n] = T;
            ws->targets.old_size[ws->targets.n] = col_sizes[T];
            ws->targets.n++;
        }
    }

    // deactivate the owner rows first: substitution turns them into 0 = 0,
    // and the sweep must skip them
    for (int idx = 0; idx < ws->substs.n; ++idx)
    {
        int i = ws->substs.recs[idx].owner;
        assert(row_sizes[i] == 2);
        assert(!HAS_TAG(row_tags[i], R_TAG_INACTIVE));
        RESET_TAG(row_tags[i], R_TAG_INACTIVE);
        row_sizes[i] = SIZE_INACTIVE_ROW;
        A->p[i].end = A->p[i].start;
        A->nnz -= 2;
    }

    // the affected rows: every active row of an eliminated column, from the
    // pre-round AT (space reserved by dton_reserve_rows), sorted ascending so
    // the change log is row-ascending per target
    int n_rows = 0;
    for (int idx = 0; idx < ws->substs.n; ++idx)
    {
        int k = ws->substs.recs[idx].k;
        for (int p = AT->p[k].start; p < AT->p[k].end; ++p)
        {
            int r = AT->i[p];
            if (!HAS_TAG(row_tags[r], R_TAG_INACTIVE))
            {
                assert(n_rows < ws->rows.cap);
                ws->rows.list[n_rows] = r;
                ws->rows.idx[n_rows] = n_rows;
                n_rows++;
            }
        }
    }
    ws->rows.n = n_rows;
    radix_sort_by_key(ws->rows.idx, (size_t) n_rows, ws->rows.list, ws->rows.aux);

    // sweep each affected row once (a row with several eliminated entries is
    // listed once per entry)
    int prev = -1;
    for (int g = 0; g < n_rows; ++g)
    {
        int q = ws->rows.list[ws->rows.idx[g]];
        if (q != prev)
        {
            dton_sweep_row(prob, ws, q, deferred, n_deferred);
            prev = q;
        }
    }

    // composed one-shot objective update in descending depth (a fixed order
    // keeps the summation into c[target] reproducible); targets are never
    // eliminated, so no c[k] read is clobbered by a c[target] write
    for (int o = 0; o < ws->substs.n; ++o)
    {
        const DtonSubst *r = ws->substs.recs + ws->substs.order[o];
        if (obj->c[r->k] != 0.0)
        {
            sub_var_in_obj_dton(obj, r->k, r->target, r->mult, r->shift);
        }
    }

    // deactivate the eliminated columns (after the activity updates, whose
    // asserts reject inactive columns)
    for (int idx = 0; idx < ws->substs.n; ++idx)
    {
        int k = ws->substs.recs[idx].k;
        assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
        col_tags[k] = C_TAG_INACTIVE;
        col_sizes[k] = SIZE_INACTIVE_COL;
    }
}

/* Full rebuild: the fallback of the refresh and the compaction of the tail
   (exported for unit tests). Transposes A into AT's existing allocation, which
   always has room since nnz never grows; the tail is whatever the allocation
   has left after the rows. */
void dton_rebuild_AT(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    Matrix *AT = constraints->AT;
    transpose_into(constraints->A, AT, constraints->state->work->iwork_n_cols);
    ws->at_valid = false;

    // only the targets' sizes can have changed (the killed columns are
    // already SIZE_INACTIVE_COL, every other active column is untouched)
    int *col_sizes = constraints->state->col_sizes;
    for (int t = 0; t < ws->targets.n; ++t)
    {
        int T = ws->targets.list[t];
        col_sizes[T] = AT->p[T].end - AT->p[T].start;
    }
}

/* Applies the round to AT column by column. Returns false if the tail ran
   out; the caller then rebuilds. A partial refresh leaves nothing the rebuild
   depends on: it only reads A. */
static bool dton_at_refresh(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    Matrix *AT = constraints->AT;
    const Matrix *A = constraints->A;
    const RowTag *row_tags = constraints->row_tags;
    int *col_sizes = constraints->state->col_sizes;

    if (!ws->at_valid)
    {
        row_slots_init(&ws->at, AT);
        ws->at_valid = true;
    }

    // stable counting sort of the log by target; the cursor array reuses the
    // accumulator's touched list, which is free now
    int *cursor = ws->acc.touched;
    for (int t = 0; t <= ws->targets.n; ++t)
    {
        ws->log.start[t] = 0;
    }
    for (int l = 0; l < ws->log.len; ++l)
    {
        ws->log.start[ws->targets.col_to_target[ws->log.col[l]] + 1]++;
    }
    for (int t = 0; t < ws->targets.n; ++t)
    {
        ws->log.start[t + 1] += ws->log.start[t];
        cursor[t] = ws->log.start[t];
    }
    for (int l = 0; l < ws->log.len; ++l)
    {
        int pos = cursor[ws->targets.col_to_target[ws->log.col[l]]]++;
        ws->log.row2[pos] = ws->log.row[l];
        ws->log.val2[pos] = ws->log.val[l];
    }

    // killed columns
    for (int idx = 0; idx < ws->substs.n; ++idx)
    {
        int k = ws->substs.recs[idx].k;
        AT->p[k].end = AT->p[k].start;
    }

    // each target merges its log segment (rows ascending) into its old
    // column; the owner rows of the round are inactive and dropped
    for (int t = 0; t < ws->targets.n; ++t)
    {
        int T = ws->targets.list[t];
        int ls = ws->log.start[t];
        int le = ws->log.start[t + 1];
        if (!matrix_update_row(AT, &ws->at, T, ws->log.row2 + ls, ws->log.val2 + ls,
                               le - ls, row_tags, R_TAG_INACTIVE))
        {
            return false; // tail full: the caller rebuilds
        }
        col_sizes[T] = AT->p[T].end - AT->p[T].start;
    }

    AT->nnz = A->nnz;
    return true;
}

/* Phase 5: brings A transpose up to date. The eliminated columns are
   emptied; each composed target (the only other columns whose content
   changed) becomes its old column minus the deactivated rows, merged with
   the sweep's log entries for it, in place when it fits its slot and
   otherwise moved to the tail. Falls back to a full rebuild into the existing
   allocation when the dirty content exceeds a quarter of nnz, the tail is
   full or its waste exceeds 2 nnz, or the change log could not grow. Then pushes the
   targets' size transitions to the empty/singleton column worklists and recounts
   their locks. */
static void dton_refresh_AT(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT;
    int *col_sizes = constraints->state->col_sizes;

    // dirty content of the round, known before anything touches AT
    size_t dirty = (size_t) ws->log.len;
    for (int t = 0; t < ws->targets.n; ++t)
    {
        int T = ws->targets.list[t];
        dirty += (size_t) (AT->p[T].end - AT->p[T].start);
    }
    bool waste =
        ws->at_valid && (size_t) (ws->at.tail_next - ws->at.tail_base) > 2 * A->nnz;
    bool rebuild =
        ws->log.overflow ||
        (double) dirty > ws->tuning.rebuild_dirty_frac * (double) A->nnz || waste;

    if (!rebuild && !dton_at_refresh(prob, ws))
    {
        rebuild = true;
    }
    if (rebuild)
    {
        dton_rebuild_AT(prob, ws);
    }
    ws->tuning.last_round_rebuilt = rebuild;
    AT = constraints->AT;
    assert(AT->nnz == A->nnz);

    for (int t = 0; t < ws->targets.n; ++t)
    {
        int T = ws->targets.list[t];
        ws->targets.col_to_target[T] = -1;
        assert(!HAS_TAG(constraints->col_tags[T], C_TAG_INACTIVE));
        assert(col_sizes[T] == AT->p[T].end - AT->p[T].start);

        // size-transition pushes (a target both gains fill-in and loses
        // its owner-row entries, so both directions occur)
        int new_size = col_sizes[T];
        int old_size = ws->targets.old_size[t];

        if (new_size == 0 && old_size != 0)
        {
            assert(!iVec_contains(constraints->state->empty_cols, T));
            iVec_append(constraints->state->empty_cols, T);
        }
        else if (new_size == 1 && old_size != 1)
        {
            assert(!iVec_contains(constraints->state->ston_cols, T));
            iVec_append(constraints->state->ston_cols, T);
        }

        // targets are the only active columns whose entries changed this
        // round, so only their locks need recounting
        RowView col_view =
            new_rowview(AT->x + AT->p[T].start, AT->i + AT->p[T].start,
                        col_sizes + T, AT->p + T, NULL, NULL, NULL, T);
        count_locks_one_column(&col_view, constraints->state->col_locks + T,
                               constraints->row_tags);
    }
    ws->targets.n = 0;

    ws->log.len = 0;
    ws->log.overflow = false;
}

PresolveStatus remove_dton_eq_rows(Problem *prob)
{
    Constraints *constraints = prob->constraints;
    State *state = constraints->state;
    Work *work = state->work;

    assert(state->ston_rows->len == 0);
    assert(state->empty_rows->len == 0);
    assert(state->empty_cols->len == 0);
    DEBUG(verify_problem_up_to_date(constraints));

    if (work->dton == NULL)
    {
        work->dton = dton_ws_new(constraints->m, constraints->n);
        // on allocation failure, skip the reduction (nothing was mutated)
        RETURN_PTR_IF_NULL(work->dton, UNCHANGED);
    }

    DtonWorkspace *ws = work->dton;
    dton_ws_attach(ws, work);

    // double ptr in case the appends realloc the vector
    iVec **dton_rows = &state->dton_rows;

    while ((*dton_rows)->len > 0)
    {
        int *deferred = work->iwork_n_rows;
        int n_deferred = 0;

        DEBUG(verify_no_duplicates_sort(*dton_rows));

        // a failed reservation (records or rows) gives the round up before
        // anything is mutated
        bool progress = dton_claim(prob, ws, deferred, &n_deferred);
        dton_compose(ws, deferred, &n_deferred);
        progress = progress && ws->substs.n > 0;
        if (progress && !dton_reserve_rows(prob, ws))
        {
            for (int idx = 0; idx < ws->substs.n; ++idx)
            {
                ws->substs.col_subst[ws->substs.recs[idx].k] = -1;
            }
            ws->substs.n = 0;
            progress = false;
        }

        if (progress)
        {
            // an infeasible transfer ends the presolve; nothing else has been
            // mutated, and the claim state is released for hygiene
            if (dton_transfer_bounds(prob, ws) == INFEASIBLE)
            {
                for (int idx = 0; idx < ws->substs.n; ++idx)
                {
                    ws->substs.col_subst[ws->substs.recs[idx].k] = -1;
                }
                ws->substs.n = 0;
                return INFEASIBLE;
            }
            dton_record(prob, ws);
            dton_apply(prob, ws, deferred, &n_deferred);
            dton_refresh_AT(prob, ws);

            // reset the claim state; the eliminated columns are inactive
            // and can never be claimed again
            for (int idx = 0; idx < ws->substs.n; ++idx)
            {
                ws->substs.col_subst[ws->substs.recs[idx].k] = -1;
            }
            ws->substs.n = 0;
        }

        // the worklist becomes the deferred rows plus the sweep's new
        // doubleton candidates
        iVec_clear_no_resize(*dton_rows);
        if (n_deferred > 0)
        {
            DEBUG(verify_no_duplicates_sort_ptr(deferred, (size_t) n_deferred));
            iVec_append_array(*dton_rows, deferred, (size_t) n_deferred);
        }

        if (!progress)
        {
            break; // only deferrals: retrying now would spin
        }
    }

    // bounds transfers cannot fix columns today; kept in case that changes
    if (state->fixed_cols_to_delete->len > 0)
    {
        delete_fixed_cols_from_problem(prob);
        delete_inactive_cols_from_A_and_AT(constraints);
    }

    DEBUG(verify_problem_up_to_date(constraints));

    return UNCHANGED;
}
