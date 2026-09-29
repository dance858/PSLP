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
#include "Constraints.h"
#include "Debugger.h"
#include "DtonsEq_internal.h"
#include "Locks.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "Problem.h"
#include "RowColViews.h"
#include "State.h"
#include "Tags.h"
#include "Workspace.h"
#include "radix_sort.h"
#include <limits.h>
#include <string.h>

/* accumulator flag bits: the row originally had this column / the
   substitution contributed to it */
#define DTON_ACC_EXISTED ((uint8_t) 1)
#define DTON_ACC_MODIFIED ((uint8_t) 2)

/* Appends a tuple to the sweep's change log (an upsert, or a delete when
   'val' is zero), doubling the arrays up to max_len as needed. On allocation
   failure the round is flagged so the refresh falls back to a full rebuild. */
static inline void dton_log_push(DtonWorkspace *ws, int col, int row, double val)
{
    if (ws->log.len == ws->log.alloc)
    {
        assert(ws->log.alloc < ws->log.max_len);
        size_t cap = MIN((size_t) ws->log.alloc * 2, (size_t) ws->log.max_len);
        bool ok = true;
        ok = ps_grow(&ws->log.col, cap, sizeof(int)) && ok;
        ok = ps_grow(&ws->log.row, cap, sizeof(int)) && ok;
        ok = ps_grow(&ws->log.val, cap, sizeof(double)) && ok;
        ok = ps_grow(&ws->log.row2, cap, sizeof(int)) && ok;
        ok = ps_grow(&ws->log.val2, cap, sizeof(double)) && ok;
        if (!ok)
        {
            ws->log.incomplete = true;
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
   columns in the pre-round AT). False if the allocation fails. */
bool dton_reserve_rows(Problem *prob, DtonWorkspace *ws)
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
    assert(need <= (size_t) INT_MAX);
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
            iVec_append(constraints->state->empty_rows, q);
            break;
        case 1:
            if (old_len != 1)
            {
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
   eliminated column lands on its composed target. Untouched entries are kept
   verbatim. A touched value below ZERO_TOL is dropped as a cancellation. The
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
    for (int ii = start; ii < end; ++ii)
    {
        int c = A->i[ii];
        double v = A->x[ii];
        int record_index = ws->substs.col_subst[c];
        if (record_index >= 0)
        {
            const DtonSubst *r = ws->substs.recs + record_index;
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
    int write_pos = start;
    for (int ii = 0; ii < acc->n_touched; ++ii)
    {
        int oc = acc->touched[ii];
        double val = acc->value[oc];
        bool keep = acc->flags[oc] == DTON_ACC_EXISTED || ABS(val) > ZERO_TOL;
        bool is_target = ws->targets.col_to_target[oc] >= 0;
        if (keep)
        {
            A->i[write_pos] = oc;
            A->x[write_pos] = val;
            ++write_pos;
            if (is_target)
            {
                dton_log_push(ws, oc, q, val);
            }
        }
        else if (acc->flags[oc] & DTON_ACC_EXISTED)
        {
            assert(is_target);
            dton_log_push(ws, oc, q, 0.0); // delete
        }
    }
    int new_len = write_pos - start;

    // full activity recompute for the rewritten row, keeping its status. A
    // NOT_ADDED row that gains a computable side is promoted.
    if (new_len > 0)
    {
        Activity *act = constraints->state->activities + q;
        uint8_t old_status = act->status;
        Activity_init(act, A->x + start, A->i + start, new_len, constraints->bounds,
                      col_tags);
        act->status = old_status;
        if (old_status == NOT_ADDED && (act->n_inf_min == 0 || act->n_inf_max == 0))
        {
            act->status = ADDED;
            iVec_append(constraints->state->updated_activities, q);
        }
    }

    dton_sweep_row_finish(prob, q, old_len, new_len, chg, deferred, n_deferred);
}

/* Phase 4: applies the round to A, the row sides, activities, worklists and
   objective, and deactivates the owner rows and eliminated columns. */
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

    // collect the composed targets, the only active columns whose content
    // changes this round, with their pre-round sizes for the size transitions
    // of phase 5. Deactivate the owner rows first: substitution turns them into
    // 0 = 0, and the sweep must skip them
    ws->targets.n = 0;
    ws->log.len = 0;
    ws->log.max_len = (int) A->nnz;
    ws->log.incomplete = false;
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

        int i = ws->substs.recs[idx].owner;
        assert(row_sizes[i] == 2);
        assert(!HAS_TAG(row_tags[i], R_TAG_INACTIVE));
        RESET_TAG(row_tags[i], R_TAG_INACTIVE);
        row_sizes[i] = SIZE_INACTIVE_ROW;
        A->p[i].end = A->p[i].start;
        A->nnz -= 2;
    }

    // the affected rows: every active row of an eliminated column, from the
    // pre-round AT, sorted ascending so the change log is row-ascending per
    // target
    int n_rows = 0;
    for (int idx = 0; idx < ws->substs.n; ++idx)
    {
        int k = ws->substs.recs[idx].k;
        for (int jj = AT->p[k].start; jj < AT->p[k].end; ++jj)
        {
            int r = AT->i[jj];
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
    for (int ii = 0; ii < n_rows; ++ii)
    {
        int q = ws->rows.list[ws->rows.idx[ii]];
        if (q != prev)
        {
            dton_sweep_row(prob, ws, q, deferred, n_deferred);
            prev = q;
        }
    }

    // composed one-shot objective update in descending depth (a fixed order
    // keeps the summation into c[target] reproducible). Targets are never
    // eliminated, so no c[k] read is clobbered by a c[target] write
    for (int ii = 0; ii < ws->substs.n; ++ii)
    {
        const DtonSubst *r = ws->substs.recs + ws->substs.order[ii];
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

/* Full rebuild: the fallback of the refresh and the compaction of the tail.
   Transposes A into AT's existing allocation, which always has room since nnz
   never grows. The tail is whatever the allocation has left after the rows. */
void dton_rebuild_AT(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    Matrix *AT = constraints->AT;
    transpose_into(constraints->A, AT, constraints->state->work->iwork_n_cols);
    ws->at_valid = false;

    // only the targets' sizes can have changed (the eliminated columns are
    // already SIZE_INACTIVE_COL, every other active column is untouched)
    int *col_sizes = constraints->state->col_sizes;
    for (int ii = 0; ii < ws->targets.n; ++ii)
    {
        int T = ws->targets.list[ii];
        col_sizes[T] = AT->p[T].end - AT->p[T].start;
    }
}

/* Applies the round to AT column by column. Returns false if the tail ran
   out. The caller then rebuilds. A partial refresh leaves nothing the rebuild
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

    // stable counting sort of the log by target
    int *cursor = ws->acc.touched;
    memset(ws->log.start, 0, (size_t) (ws->targets.n + 1) * sizeof(int));
    for (int ii = 0; ii < ws->log.len; ++ii)
    {
        ws->log.start[ws->targets.col_to_target[ws->log.col[ii]] + 1]++;
    }
    for (int ii = 0; ii < ws->targets.n; ++ii)
    {
        ws->log.start[ii + 1] += ws->log.start[ii];
        cursor[ii] = ws->log.start[ii];
    }
    for (int ii = 0; ii < ws->log.len; ++ii)
    {
        int pos = cursor[ws->targets.col_to_target[ws->log.col[ii]]]++;
        ws->log.row2[pos] = ws->log.row[ii];
        ws->log.val2[pos] = ws->log.val[ii];
    }

    // empty the eliminated columns
    for (int idx = 0; idx < ws->substs.n; ++idx)
    {
        int k = ws->substs.recs[idx].k;
        AT->p[k].end = AT->p[k].start;
    }

    // each target merges its log segment (rows ascending) into its old
    // column. The owner rows of the round are inactive and dropped
    for (int ii = 0; ii < ws->targets.n; ++ii)
    {
        int T = ws->targets.list[ii];
        int log_start = ws->log.start[ii];
        int log_end = ws->log.start[ii + 1];
        if (!matrix_update_row(AT, &ws->at, T, ws->log.row2 + log_start,
                               ws->log.val2 + log_start, log_end - log_start,
                               row_tags, R_TAG_INACTIVE))
        {
            return false; // tail full: the caller rebuilds
        }
        col_sizes[T] = AT->p[T].end - AT->p[T].start;
    }

    AT->nnz = A->nnz;
    return true;
}

/* Phase 5: brings A transpose up to date. The eliminated columns are emptied
   and each target is merged with its log segment, in place or moved to the
   tail. Falls back to a full rebuild when the refresh does not pay off or
   cannot finish. Then updates the targets' size worklists and locks. */
void dton_refresh_AT(Problem *prob, DtonWorkspace *ws)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT;
    int *col_sizes = constraints->state->col_sizes;

    // dirty content of the round, known before anything touches AT
    size_t dirty = (size_t) ws->log.len;
    for (int ii = 0; ii < ws->targets.n; ++ii)
    {
        int T = ws->targets.list[ii];
        dirty += (size_t) (AT->p[T].end - AT->p[T].start);
    }
    bool waste =
        ws->at_valid && (size_t) (ws->at.tail_next - ws->at.tail_base) > 2 * A->nnz;
    bool rebuild = ws->log.incomplete ||
                   (double) dirty > ws->rebuild_dirty_frac * (double) A->nnz ||
                   waste || !dton_at_refresh(prob, ws);
    if (rebuild)
    {
        dton_rebuild_AT(prob, ws);
    }
    ws->last_round_rebuilt = rebuild;
    assert(AT->nnz == A->nnz);

    for (int ii = 0; ii < ws->targets.n; ++ii)
    {
        int T = ws->targets.list[ii];
        ws->targets.col_to_target[T] = -1;
        assert(!HAS_TAG(constraints->col_tags[T], C_TAG_INACTIVE));
        assert(col_sizes[T] == AT->p[T].end - AT->p[T].start);

        // size-transition pushes
        int new_size = col_sizes[T];
        int old_size = ws->targets.old_size[ii];

        assert(old_size > 0);
        if (new_size == 0)
        {
            iVec_append(constraints->state->empty_cols, T);
        }
        else if (new_size == 1 && old_size != 1)
        {
            iVec_append(constraints->state->ston_cols, T);
        }

        // recount the targets' locks
        RowView col_view =
            new_rowview(AT->x + AT->p[T].start, AT->i + AT->p[T].start,
                        col_sizes + T, AT->p + T, NULL, NULL, NULL, T);
        count_locks_one_column(&col_view, constraints->state->col_locks + T,
                               constraints->row_tags);
    }
}
