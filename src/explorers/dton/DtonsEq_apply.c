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

/* Accumulator flag bits: the row originally had this column, or a
   substitution contributed to it. */
#define DTON_ACC_EXISTED ((uint8_t) 1)
#define DTON_ACC_MODIFIED ((uint8_t) 2)

/* Appends a tuple to a target's segment of the change log (an upsert, or a
   delete when 'val' is zero). Does nothing in a round that rebuilds. */
static inline void dton_log_push(DtonWorkspace *dton_work, int target_index, int row,
                                 double val)
{
    if (!dton_work->merge_round)
    {
        return;
    }
    DtonTarget *target = dton_work->targets.list + target_index;
    int pos = target->log_end++;
    assert(pos < (target_index + 1 < dton_work->targets.n ? target[1].log_start
                                                          : dton_work->log.cap));
    dton_work->log.row[pos] = row;
    dton_work->log.val[pos] = val;
}

/* Decides whether the round merges its changes into AT or rebuilds it, and
   lays out the change log of a merge round. On entry each target's log_end
   holds the room its segment needs. */
static void dton_log_start(DtonWorkspace *dton_work, size_t nnz)
{
    DtonTargets *targets = &dton_work->targets;
    DtonLog *change_log = &dton_work->log;

    /* Dirty content of the round. */
    double dirty = 0.0;
    for (int ii = 0; ii < targets->n; ++ii)
    {
        dirty += targets->list[ii].log_end;
    }
    bool tail_bloated =
        dton_work->at_valid &&
        (size_t) (dton_work->at.tail_next - dton_work->at.tail_base) > 2 * nnz;
    dton_work->merge_round =
        !tail_bloated && dirty <= dton_work->rebuild_dirty_frac * (double) nnz;
    if (!dton_work->merge_round)
    {
        return;
    }

    /* A failed allocation turns the round into a rebuild. */
    if (dirty > change_log->cap)
    {
        if (!ps_grow(&change_log->row, (size_t) dirty, sizeof(int)) ||
            !ps_grow(&change_log->val, (size_t) dirty, sizeof(double)))
        {
            dton_work->merge_round = false;
            return;
        }
        change_log->cap = (int) dirty;
    }

    int log_start = 0;
    for (int ii = 0; ii < targets->n; ++ii)
    {
        DtonTarget *target = targets->list + ii;
        int room = target->log_end;
        target->log_start = log_start;
        target->log_end = log_start;
        log_start += room;
    }
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
                dton_log_push(dton_work, target_index, row, val);
            }
        }
        else if (acc->flags[col] & DTON_ACC_EXISTED)
        {
            assert(is_target);
            dton_log_push(dton_work, target_index, row, 0.0); /* Delete. */
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
    uint64_t *swept_rows = dton_work->swept_rows;

    /* Collect the composed targets, the only active columns whose content
       changes this round, with their pre-round sizes. A target's log segment
       needs room for one tuple per entry of its old column. */
    targets->n = 0;
    for (int idx = 0; idx < substs->n; ++idx)
    {
        int col = substs->recs[idx].target;
        if (targets->col_to_target[col] < 0)
        {
            targets->col_to_target[col] = targets->n;
            DtonTarget *target = targets->list + targets->n++;
            target->col = col;
            target->old_size = col_sizes[col];
            target->log_end = col_sizes[col];
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

    /* Mark the rows to sweep: every active row of an eliminated column, from
       the pre-round AT. Each marked entry adds one tuple of room to the log
       segment of the column's target. */
    for (int idx = 0; idx < substs->n; ++idx)
    {
        const DtonSubst *rec = substs->recs + idx;
        DtonTarget *target = targets->list + targets->col_to_target[rec->target];
        for (int jj = AT->p[rec->k].start; jj < AT->p[rec->k].end; ++jj)
        {
            int row = AT->i[jj];
            if (!HAS_TAG(row_tags[row], R_TAG_INACTIVE))
            {
                swept_rows[row / 64] |= (uint64_t) 1 << (row % 64);
                target->log_end++;
            }
        }
    }

    dton_log_start(dton_work, A->nnz);

    /* Sweep the marked rows in ascending order, so every log segment is
       row-ascending, and clear the marks. */
    int n_words = (dton_work->m + 63) / 64;
    for (int ii = 0; ii < n_words; ++ii)
    {
        uint64_t word = swept_rows[ii];
        swept_rows[ii] = 0;
        for (int row = 64 * ii; word != 0; ++row, word >>= 1)
        {
            if (word & 1)
            {
                dton_sweep_row(prob, dton_work, row, deferred, n_deferred);
            }
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
        int col = dton_work->targets.list[ii].col;
        col_sizes[col] = AT->p[col].end - AT->p[col].start;
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
    const DtonLog *change_log = &dton_work->log;

    if (!dton_work->at_valid)
    {
        row_slots_init(&dton_work->at, AT);
        dton_work->at_valid = true;
    }

    /* Empty the eliminated columns. */
    for (int idx = 0; idx < dton_work->substs.n; ++idx)
    {
        int k = dton_work->substs.recs[idx].k;
        AT->p[k].end = AT->p[k].start;
    }

    /* Each target merges its log segment (rows ascending) into its old column.
       The owner rows of the round are inactive and dropped. */
    for (int ii = 0; ii < dton_work->targets.n; ++ii)
    {
        const DtonTarget *target = dton_work->targets.list + ii;
        if (!matrix_update_row(
                AT, &dton_work->at, target->col, change_log->row + target->log_start,
                change_log->val + target->log_start,
                target->log_end - target->log_start, row_tags, R_TAG_INACTIVE))
        {
            return false; /* Tail full: the caller rebuilds. */
        }
        col_sizes[target->col] = AT->p[target->col].end - AT->p[target->col].start;
    }

    AT->nnz = A->nnz;
    return true;
}

/* Phase 5: brings A transpose up to date. In a merge round the eliminated
   columns are emptied and each target is merged with its log segment, in
   place or moved to the tail. Otherwise, or when the merge cannot finish, A
   transpose is rebuilt. Then updates the targets' size worklists and locks. */
void dton_update_AT(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT;
    int *col_sizes = constraints->state->col_sizes;

    if (dton_work->merge_round)
    {
        dton_work->merge_round = dton_merge_into_AT(prob, dton_work);
    }
    if (!dton_work->merge_round)
    {
        dton_rebuild_AT(prob, dton_work);
    }
    assert(AT->nnz == constraints->A->nnz);

    for (int ii = 0; ii < dton_work->targets.n; ++ii)
    {
        int col = dton_work->targets.list[ii].col;
        dton_work->targets.col_to_target[col] = -1;
        assert(!HAS_TAG(constraints->col_tags[col], C_TAG_INACTIVE));
        assert(col_sizes[col] == AT->p[col].end - AT->p[col].start);

        /* Push the size transitions. */
        int new_size = col_sizes[col];
        int old_size = dton_work->targets.list[ii].old_size;

        assert(old_size > 0);
        if (new_size == 0)
        {
            iVec_append(constraints->state->empty_cols, col);
        }
        else if (new_size == 1 && old_size != 1)
        {
            iVec_append(constraints->state->ston_cols, col);
        }

        /* Recount the target's locks. */
        RowView col_view =
            new_rowview(AT->x + AT->p[col].start, AT->i + AT->p[col].start,
                        col_sizes + col, AT->p + col, NULL, NULL, NULL, col);
        count_locks_one_column(&col_view, constraints->state->col_locks + col,
                               constraints->row_tags);
    }
}
