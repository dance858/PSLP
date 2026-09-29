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

/* Accumulator flag bits: the row originally had this column, or a
   substitution contributed to it. */
#define DTON_ACC_EXISTED ((uint8_t) 1)
#define DTON_ACC_MODIFIED ((uint8_t) 2)

/* Appends a tuple to the sweep's change log (an upsert, or a delete when
   'val' is zero), doubling the arrays up to max_len as needed. On allocation
   failure the round is flagged so the update falls back to a full rebuild. */
static inline void dton_log_push(DtonLog *change_log, int target_index, int row,
                                 double val)
{
    if (change_log->len == change_log->cap)
    {
        assert(change_log->cap < change_log->max_len);
        size_t cap = MIN((size_t) change_log->cap * 2, (size_t) change_log->max_len);
        bool ok = true;
        ok &= ps_grow(&change_log->target_index, cap, sizeof(int));
        ok &= ps_grow(&change_log->row, cap, sizeof(int));
        ok &= ps_grow(&change_log->val, cap, sizeof(double));
        ok &= ps_grow(&change_log->sorted_row, cap, sizeof(int));
        ok &= ps_grow(&change_log->sorted_val, cap, sizeof(double));
        if (!ok)
        {
            change_log->incomplete = true;
            return;
        }
        change_log->cap = (int) cap;
    }
    change_log->target_index[change_log->len] = target_index;
    change_log->row[change_log->len] = row;
    change_log->val[change_log->len] = val;
    change_log->len++;
}

/* Reserves the row list for the round (one slot per entry of the eliminated
   columns in the pre-round AT). False if the allocation fails. */
bool dton_reserve_rows(Problem *prob, DtonWorkspace *dton_work)
{
    const Matrix *AT = prob->constraints->AT;
    int subst_cols_nnz = 0;

    for (int ii = 0; ii < dton_work->substs.n; ++ii)
    {
        int k = dton_work->substs.recs[ii].k;
        subst_cols_nnz += AT->p[k].end - AT->p[k].start;
    }
    if (subst_cols_nnz <= dton_work->rows.cap)
    {
        return true;
    }

    bool ok = true;
    ok &= ps_grow(&dton_work->rows.list, (size_t) subst_cols_nnz, sizeof(int));
    ok &= ps_grow(&dton_work->rows.perm, (size_t) subst_cols_nnz, sizeof(int));
    ok &=
        ps_grow(&dton_work->rows.sort_scratch, (size_t) subst_cols_nnz, sizeof(int));
    if (!ok)
    {
        return false;
    }
    dton_work->rows.cap = subst_cols_nnz;
    return true;
}

/* Finishes a swept row whose entries now occupy [start, start + new_len): shifts
   the finite sides by the substituted constants, updates the sizes, and feeds the
   row worklists (the old_len guards keep deferred rows from being pushed twice). */
static void dton_sweep_row_finish(Problem *prob, int row, int old_len, int new_len,
                                  double side_shift, int *deferred, int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    RowTag *row_tags = constraints->row_tags;
    int start = A->p[row].start;

    assert(new_len <= old_len);
    A->p[row].end = start + new_len;
    constraints->state->row_sizes[row] = new_len;
    A->nnz -= (size_t) (old_len - new_len);
    DEBUG(ASSERT_INCREASING_I(A->i + start, (size_t) new_len););
    DEBUG(ASSERT_NO_ZEROS_D(A->x + start, (size_t) new_len););

    if (side_shift != 0.0)
    {
        if (!HAS_TAG(row_tags[row], R_TAG_LHS_INF))
        {
            constraints->lhs[row] -= side_shift;
        }
        if (!HAS_TAG(row_tags[row], R_TAG_RHS_INF))
        {
            constraints->rhs[row] -= side_shift;
        }
    }

    switch (new_len)
    {
        case 0:
            iVec_append(constraints->state->empty_rows, row);
            break;
        case 1:
            if (old_len != 1)
            {
                iVec_append(constraints->state->ston_rows, row);
            }
            break;
        case 2:
            if (HAS_TAG(row_tags[row], R_TAG_EQ) && old_len != 2)
            {
                deferred[(*n_deferred)++] = row;
            }
            break;
        default:
            break;
    }
}

/* Rewrites a row in place through the sparse accumulator. Every entry of an
   eliminated column lands on its composed target. Untouched entries are kept
   verbatim. A touched value below ZERO_TOL is dropped as a cancellation. The
   activity is recomputed from scratch. */
static void dton_sweep_row(Problem *prob, DtonWorkspace *dton_work, int row,
                           int *deferred, int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    const ColTag *col_tags = constraints->col_tags;
    int start = A->p[row].start;
    int end = A->p[row].end;
    int old_len = end - start;
    SparseAccumulator *acc = &dton_work->acc;
    double side_shift = 0.0;

    sparse_accumulator_clear(acc);
    for (int ii = start; ii < end; ++ii)
    {
        int col = A->i[ii];
        double coeff = A->x[ii];
        int record_index = dton_work->substs.col_subst[col];
        if (record_index >= 0)
        {
            const DtonSubst *rec = dton_work->substs.recs + record_index;
            sparse_accumulator_add(acc, rec->target, coeff * rec->mult,
                                   DTON_ACC_MODIFIED);
            side_shift += coeff * rec->shift;
        }
        else
        {
            sparse_accumulator_add(acc, col, coeff, DTON_ACC_EXISTED);
        }
    }

    /* Rows must stay column-sorted. */
    sparse_accumulator_sort(acc);

    /* Rebuild the row in place. It never grows: every eliminated entry maps
       onto exactly one target, and merging only shrinks. */
    int write_pos = start;
    for (int ii = 0; ii < acc->n_touched; ++ii)
    {
        int col = acc->touched[ii];
        double val = acc->value[col];
        bool keep = acc->flags[col] == DTON_ACC_EXISTED || ABS(val) > ZERO_TOL;
        int target_index = dton_work->targets.col_to_target[col];
        bool is_target = target_index >= 0;
        if (keep)
        {
            A->i[write_pos] = col;
            A->x[write_pos] = val;
            ++write_pos;
            if (is_target)
            {
                dton_log_push(&dton_work->log, target_index, row, val);
            }
        }
        else if (acc->flags[col] & DTON_ACC_EXISTED)
        {
            assert(is_target);
            dton_log_push(&dton_work->log, target_index, row, 0.0); /* Delete. */
        }
    }
    int new_len = write_pos - start;

    /* Recompute the row's activity from scratch, keeping its status. A
       NOT_ADDED row that gains a computable side is promoted. */
    if (new_len > 0)
    {
        Activity *act = constraints->state->activities + row;
        uint8_t old_status = act->status;
        Activity_init(act, A->x + start, A->i + start, new_len, constraints->bounds,
                      col_tags);
        act->status = old_status;
        if (old_status == NOT_ADDED && (act->n_inf_min == 0 || act->n_inf_max == 0))
        {
            act->status = ADDED;
            iVec_append(constraints->state->updated_activities, row);
        }
    }

    dton_sweep_row_finish(prob, row, old_len, new_len, side_shift, deferred,
                          n_deferred);
}

/* Phase 4: applies the round to A, the row sides, activities, worklists and
   objective, and deactivates the owner rows and eliminated columns. */
void dton_apply(Problem *prob, DtonWorkspace *dton_work, int *deferred,
                int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT; /* pre-round transpose, read-only */
    RowTag *row_tags = constraints->row_tags;
    ColTag *col_tags = constraints->col_tags;
    int *row_sizes = constraints->state->row_sizes;
    int *col_sizes = constraints->state->col_sizes;
    Objective *obj = prob->obj;
    const DtonSubsts *substs = &dton_work->substs;
    DtonTargets *targets = &dton_work->targets;
    DtonRows *rows = &dton_work->rows;

    /* Start an empty change log. */
    dton_work->log.len = 0;
    dton_work->log.max_len = (int) A->nnz;
    dton_work->log.incomplete = false;

    /* Collect the composed targets, the only active columns whose content
       changes this round, with their pre-round sizes. */
    targets->n = 0;
    for (int idx = 0; idx < substs->n; ++idx)
    {
        int target = substs->recs[idx].target;
        if (targets->col_to_target[target] < 0)
        {
            targets->col_to_target[target] = targets->n;
            targets->list[targets->n] = target;
            targets->old_size[targets->n] = col_sizes[target];
            targets->n++;
        }
    }

    /* Deactivate the owner rows. Substitution turns them into 0 = 0, and the
       sweep must skip them. */
    for (int idx = 0; idx < substs->n; ++idx)
    {
        int i = substs->recs[idx].owner;
        assert(row_sizes[i] == 2);
        assert(!HAS_TAG(row_tags[i], R_TAG_INACTIVE));
        RESET_TAG(row_tags[i], R_TAG_INACTIVE);
        row_sizes[i] = SIZE_INACTIVE_ROW;
        A->p[i].end = A->p[i].start;
        A->nnz -= 2;
    }

    /* The affected rows: every active row of an eliminated column, from the
       pre-round AT, sorted ascending so the change log is row-ascending per
       target. */
    int n_rows = 0;
    for (int idx = 0; idx < substs->n; ++idx)
    {
        int k = substs->recs[idx].k;
        for (int jj = AT->p[k].start; jj < AT->p[k].end; ++jj)
        {
            int row = AT->i[jj];
            if (!HAS_TAG(row_tags[row], R_TAG_INACTIVE))
            {
                assert(n_rows < rows->cap);
                rows->list[n_rows] = row;
                rows->perm[n_rows] = n_rows;
                n_rows++;
            }
        }
    }
    radix_sort_by_key(rows->perm, (size_t) n_rows, rows->list, rows->sort_scratch);

    /* Sweep each affected row once. A row with several eliminated entries is
       listed once per entry. */
    int prev_row = -1;
    for (int ii = 0; ii < n_rows; ++ii)
    {
        int row = rows->list[rows->perm[ii]];
        if (row != prev_row)
        {
            dton_sweep_row(prob, dton_work, row, deferred, n_deferred);
            prev_row = row;
        }
    }

    /* Substitute the eliminated columns in the objective. Targets are never
       eliminated, so no c[k] read follows a c[target] write. */
    for (int ii = 0; ii < substs->n; ++ii)
    {
        const DtonSubst *rec = substs->recs + substs->order[ii];
        if (obj->c[rec->k] != 0.0)
        {
            sub_var_in_obj_dton(obj, rec->k, rec->target, rec->mult, rec->shift);
        }
    }

    /* Deactivate the eliminated columns, after the activity updates, whose
       asserts reject inactive columns. */
    for (int idx = 0; idx < substs->n; ++idx)
    {
        int k = substs->recs[idx].k;
        assert(!HAS_TAG(col_tags[k], C_TAG_INACTIVE));
        RESET_TAG(col_tags[k], C_TAG_INACTIVE);
        col_sizes[k] = SIZE_INACTIVE_COL;
    }
}

/* Full rebuild of AT from A. */
void dton_rebuild_AT(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    Matrix *AT = constraints->AT;
    transpose_into(constraints->A, AT, constraints->state->work->iwork_n_cols);
    dton_work->at_valid = false;

    /* Update the targets' column sizes. */
    int *col_sizes = constraints->state->col_sizes;
    for (int ii = 0; ii < dton_work->targets.n; ++ii)
    {
        int target = dton_work->targets.list[ii];
        col_sizes[target] = AT->p[target].end - AT->p[target].start;
    }
}

/* Sorts the change log by target into sorted_row/sorted_val, keeping rows ascending
   within a target. Target t's entries end up in [start[t], start[t + 1]).
   'cursor' is scratch of n_targets ints. */
static inline void dton_log_sort_by_target(DtonLog *change_log, int n_targets,
                                           int *cursor)
{
    memset(change_log->start, 0, (size_t) (n_targets + 1) * sizeof(int));
    for (int ii = 0; ii < change_log->len; ++ii)
    {
        change_log->start[change_log->target_index[ii] + 1]++;
    }
    for (int ii = 0; ii < n_targets; ++ii)
    {
        change_log->start[ii + 1] += change_log->start[ii];
        cursor[ii] = change_log->start[ii];
    }
    for (int ii = 0; ii < change_log->len; ++ii)
    {
        int pos = cursor[change_log->target_index[ii]]++;
        change_log->sorted_row[pos] = change_log->row[ii];
        change_log->sorted_val[pos] = change_log->val[ii];
    }
}

/* Applies the round to AT column by column. Returns false if the tail ran
   out. The caller then rebuilds. A partial merge leaves nothing the rebuild
   depends on: it only reads A. */
static bool dton_merge_into_AT(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    Matrix *AT = constraints->AT;
    const Matrix *A = constraints->A;
    const RowTag *row_tags = constraints->row_tags;
    int *col_sizes = constraints->state->col_sizes;
    DtonLog *change_log = &dton_work->log;
    const int *targets = dton_work->targets.list;
    int n_targets = dton_work->targets.n;

    if (!dton_work->at_valid)
    {
        row_slots_init(&dton_work->at, AT);
        dton_work->at_valid = true;
    }

    dton_log_sort_by_target(change_log, n_targets, dton_work->acc.touched);

    /* Empty the eliminated columns. */
    for (int idx = 0; idx < dton_work->substs.n; ++idx)
    {
        int k = dton_work->substs.recs[idx].k;
        AT->p[k].end = AT->p[k].start;
    }

    /* Each target merges its log segment (rows ascending) into its old column.
       The owner rows of the round are inactive and dropped. */
    for (int ii = 0; ii < n_targets; ++ii)
    {
        int target = targets[ii];
        int log_start = change_log->start[ii];
        int log_end = change_log->start[ii + 1];
        if (!matrix_update_row(AT, &dton_work->at, target,
                               change_log->sorted_row + log_start,
                               change_log->sorted_val + log_start,
                               log_end - log_start, row_tags, R_TAG_INACTIVE))
        {
            return false; /* Tail full: the caller rebuilds. */
        }
        col_sizes[target] = AT->p[target].end - AT->p[target].start;
    }

    AT->nnz = A->nnz;
    return true;
}

/* Phase 5: brings A transpose up to date. The eliminated columns are emptied
   and each target is merged with its log segment, in place or moved to the
   tail. Falls back to a full rebuild when merging does not pay off or
   cannot finish. Then updates the targets' size worklists and locks. */
void dton_update_AT(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT;
    int *col_sizes = constraints->state->col_sizes;

    /* Dirty content of the round, known before anything touches AT. */
    double dirty = dton_work->log.len;
    for (int ii = 0; ii < dton_work->targets.n; ++ii)
    {
        int target = dton_work->targets.list[ii];
        dirty += AT->p[target].end - AT->p[target].start;
    }
    bool tail_bloated =
        dton_work->at_valid &&
        (size_t) (dton_work->at.tail_next - dton_work->at.tail_base) > 2 * A->nnz;
    bool rebuild = dton_work->log.incomplete ||
                   dirty > dton_work->rebuild_dirty_frac * (double) A->nnz ||
                   tail_bloated;
    if (!rebuild)
    {
        rebuild = !dton_merge_into_AT(prob, dton_work);
    }
    if (rebuild)
    {
        dton_rebuild_AT(prob, dton_work);
    }
    dton_work->last_round_rebuilt = rebuild;
    assert(AT->nnz == A->nnz);

    for (int ii = 0; ii < dton_work->targets.n; ++ii)
    {
        int target = dton_work->targets.list[ii];
        dton_work->targets.col_to_target[target] = -1;
        assert(!HAS_TAG(constraints->col_tags[target], C_TAG_INACTIVE));
        assert(col_sizes[target] == AT->p[target].end - AT->p[target].start);

        /* Push the size transitions. */
        int new_size = col_sizes[target];
        int old_size = dton_work->targets.old_size[ii];

        assert(old_size > 0);
        if (new_size == 0)
        {
            iVec_append(constraints->state->empty_cols, target);
        }
        else if (new_size == 1 && old_size != 1)
        {
            iVec_append(constraints->state->ston_cols, target);
        }

        /* Recount the target's locks. */
        RowView col_view = new_rowview(
            AT->x + AT->p[target].start, AT->i + AT->p[target].start,
            col_sizes + target, AT->p + target, NULL, NULL, NULL, target);
        count_locks_one_column(&col_view, constraints->state->col_locks + target,
                               constraints->row_tags);
    }
}
