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

/* Phases 4 and 5 of a round. Phase 4 (dton_apply) substitutes the eliminated
   columns in A: every entry of an eliminated column moves onto the surviving
   column it is expressed in. Phase 5 (dton_update_AT) brings AT up to date, by
   merging the changes logged in phase 4 into those surviving columns or by
   rebuilding AT from A. */

/* Accumulator flag: a substitution contributed to the column. */
#define DTON_ACC_MODIFIED ((uint8_t) 1)

/* We maintain a log (dton_work->log) that lists the entries of surviving
   columns that the sweep changed, so that AT can be updated without a rebuild.
   This function records that a surviving column's entry in 'row' is now 'val',
   or was removed when 'val' is zero. */
static inline void dton_log_push(DtonWorkspace *dton_work, int target_index, int row,
                                 double val)
{
    if (!dton_work->merge_round)
    {
        return;
    }
    DtonTarget *target = dton_work->targets.list + target_index;
    int pos = target->log_end++;
    assert(pos < dton_work->log.cap);
    dton_work->log.row[pos] = row;
    dton_work->log.val[pos] = val;
}

/* True if merging the round's changes into AT beats rebuilding it.
   'affected_nnz' is the number of entries in the eliminated columns and in the
   surviving columns they are mapped onto. */
static bool dton_merge_pays_off(const DtonWorkspace *dton_work, double affected_nnz,
                                size_t nnz)
{
    /* rebuild when the round changes a large part of the matrix */
    bool rebuild = affected_nnz > dton_work->rebuild_frac * (double) nnz;

    /* rebuild when the tail of AT has grown large */
    if (dton_work->AT_slots_valid)
    {
        int tail_len = dton_work->AT_slots.tail_next - dton_work->AT_slots.tail_base;
        rebuild = rebuild || (size_t) tail_len > 2 * nnz;
    }

    return !rebuild;
}

/* Gives each changed surviving column its segment of the change log. On entry
   its log_end holds the number of tuples the segment must have room for. False
   if the log cannot be grown. */
static bool dton_log_layout(DtonWorkspace *dton_work, int log_len)
{
    DtonLog *change_log = &dton_work->log;
    if (log_len > change_log->cap)
    {
        if (!ps_grow(&change_log->row, (size_t) log_len, sizeof(int)) ||
            !ps_grow(&change_log->val, (size_t) log_len, sizeof(double)))
        {
            return false;
        }
        change_log->cap = log_len;
    }

    int log_start = 0;
    for (int ii = 0; ii < dton_work->targets.n_targets; ++ii)
    {
        DtonTarget *target = dton_work->targets.list + ii;
        int room = target->log_end;
        target->log_start = log_start;
        target->log_end = log_start;
        log_start += room;
    }
    return true;
}

/* Updates what depends on a swept row whose entries now occupy
   [start, start + new_len): its size, activity, sides and worklists. */
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

    /* Recompute the row's activity from scratch, keeping its status. A
       NOT_ADDED row that gains a computable side is promoted. */
    if (new_len > 0)
    {
        Activity *act = constraints->state->activities + row;
        uint8_t old_status = act->status;
        Activity_init(act, A->x + start, A->i + start, new_len, constraints->bounds,
                      constraints->col_tags);
        act->status = old_status;
        if (old_status == NOT_ADDED && (act->n_inf_min == 0 || act->n_inf_max == 0))
        {
            act->status = ADDED;
            iVec_append(constraints->state->updated_activities, row);
        }
    }

    /* Shift the finite sides by the constants of the substitutions. */
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

    /* Push the row onto the worklist of its new size. A singleton or doubleton
       row whose size did not change is already on its worklist. */
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

/* Rewrites a row of A in place. Every entry of an eliminated column moves
   onto the surviving column it is expressed in. A modified entry below
   ZERO_TOL is dropped. The other entries are kept unchanged. */
static void dton_sweep_row(Problem *prob, DtonWorkspace *dton_work, int row,
                           int *deferred, int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    int start = A->p[row].start;
    int end = A->p[row].end;
    int old_len = end - start;
    SparseAccumulator *acc = &dton_work->acc;
    double side_shift = 0.0;

    /* substitute any eliminated column with its survivor */
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
            sparse_accumulator_add(acc, col, coeff, 0);
        }
    }

    /* sort the row so the column indices within the row are increasing */
    sparse_accumulator_sort(acc);

    /* write the nonzeros of the modified row back to A and log modified entries */
    int write_pos = start;
    for (int ii = 0; ii < acc->n_touched; ++ii)
    {
        int col = acc->touched[ii];
        double val = acc->value[col];
        bool modified = acc->flags[col] == DTON_ACC_MODIFIED;
        bool keep = !modified || ABS(val) > ZERO_TOL;
        if (keep)
        {
            A->i[write_pos] = col;
            A->x[write_pos] = val;
            ++write_pos;
        }

        /* log modified entries for updating AT later */
        if (modified)
        {
            dton_log_push(dton_work, dton_work->targets.col_to_target[col], row,
                          keep ? val : 0.0);
        }
    }
    int new_len = write_pos - start;

    dton_sweep_row_finish(prob, row, old_len, new_len, side_shift, deferred,
                          n_deferred);
}

/* Collects the surviving columns that the round changes, decides whether the
   round merges its changes into AT or rebuilds it, and lays out the change log
   of a merge round. */
static void dton_collect_targets(Problem *prob, DtonWorkspace *dton_work)
{
    const int *col_sizes = prob->constraints->state->col_sizes;
    const DtonSubsts *substs = &dton_work->substs;
    DtonTargets *targets = &dton_work->targets;

    /* Collect the columns that eliminated columns are mapped onto (the only
       surviving columns that change this round). */
    double affected_nnz = 0.0;
    int log_len = 0;
    targets->n_targets = 0;
    for (int idx = 0; idx < substs->n_recs; ++idx)
    {
        const DtonSubst *rec = substs->recs + idx;
        int target_index = targets->col_to_target[rec->target];
        if (target_index < 0)
        {
            target_index = targets->n_targets++;
            targets->col_to_target[rec->target] = target_index;
            targets->list[target_index].col = rec->target;
            targets->list[target_index].log_end = 0;
            affected_nnz += col_sizes[rec->target];
        }
        targets->list[target_index].log_end += col_sizes[rec->k];
        affected_nnz += col_sizes[rec->k];
        log_len += col_sizes[rec->k];
    }

    /* A failed log allocation turns the round into a rebuild. */
    dton_work->merge_round =
        dton_merge_pays_off(dton_work, affected_nnz, prob->constraints->A->nnz) &&
        dton_log_layout(dton_work, log_len);
}

/* Applies the round to A, the lhs and rhs, activities, worklists and objective */
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
    uint64_t *swept_rows = dton_work->swept_rows;

    /* mark doubleton rows used for substitution as inactive */
    for (int idx = 0; idx < substs->n_recs; ++idx)
    {
        int i = substs->recs[idx].owner;
        assert(row_sizes[i] == 2 && !HAS_TAG(row_tags[i], R_TAG_INACTIVE));
        RESET_TAG(row_tags[i], R_TAG_INACTIVE);
        row_sizes[i] = SIZE_INACTIVE_ROW;
        A->p[i].end = A->p[i].start;
    }

    /* update the number of nonzeros in A */
    A->nnz -= 2 * (size_t) substs->n_recs;

    /* collect the surviving columns that the eliminated columns are mapped onto */
    dton_collect_targets(prob, dton_work);

    /* mark rows that need to be updated (ie., rows with eliminated columns) */
    for (int idx = 0; idx < substs->n_recs; ++idx)
    {
        int k = substs->recs[idx].k;
        for (int jj = AT->p[k].start; jj < AT->p[k].end; ++jj)
        {
            int row = AT->i[jj];
            if (!HAS_TAG(row_tags[row], R_TAG_INACTIVE))
            {
                bitmap_set(swept_rows, row);
            }
        }
    }

    /* update marked rows (this affects A but not AT) */
    int n_words = bitmap_words(dton_work->n_rows);
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

    /* substitute the eliminated columns in the objective and mark them inactive */
    for (int ii = 0; ii < substs->n_recs; ++ii)
    {
        const DtonSubst *rec = substs->recs + substs->order[ii];
        if (obj->c[rec->k] != 0.0)
        {
            sub_var_in_obj_dton(obj, rec->k, rec->target, rec->mult, rec->shift);
        }

        assert(!HAS_TAG(col_tags[rec->k], C_TAG_INACTIVE));
        RESET_TAG(col_tags[rec->k], C_TAG_INACTIVE);
        col_sizes[rec->k] = SIZE_INACTIVE_COL;
    }
}

/* Full rebuild of AT from A. */
void dton_rebuild_AT(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    Matrix *AT = constraints->AT;
    transpose_into(constraints->A, AT, constraints->state->work->iwork_n_cols);
    dton_work->AT_slots_valid = false;
}

/* Merges the change log into AT, one column at a time. Returns false if the
   tail of AT runs out of room, and the caller then rebuilds. A partial merge
   is harmless because the rebuild only reads A. */
static bool dton_merge_into_AT(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    Matrix *AT = constraints->AT;
    const Matrix *A = constraints->A;
    const RowTag *row_tags = constraints->row_tags;
    const DtonLog *change_log = &dton_work->log;

    /* recompute how much room each column of AT has, if a rebuild changed it. */
    if (!dton_work->AT_slots_valid)
    {
        row_slots_init(&dton_work->AT_slots, AT);
        dton_work->AT_slots_valid = true;
    }

    /* mark eliminated columns as empty in AT */
    for (int idx = 0; idx < dton_work->substs.n_recs; ++idx)
    {
        int k = dton_work->substs.recs[idx].k;
        AT->p[k].end = AT->p[k].start;
    }

    /* Every surviving column that the round changed has its own segment of the
       log. Merge each segment into its column. Entries in the inactive owner
       rows are dropped. */
    for (int ii = 0; ii < dton_work->targets.n_targets; ++ii)
    {
        const DtonTarget *target = dton_work->targets.list + ii;
        if (!matrix_update_row(AT, &dton_work->AT_slots, target->col,
                               change_log->row + target->log_start,
                               change_log->val + target->log_start,
                               target->log_end - target->log_start, row_tags,
                               R_TAG_INACTIVE))
        {
            /* if the merge fails for one column (ie., one row of AT), we must
               rebuild AT from scratch */
            return false;
        }
    }

    AT->nnz = A->nnz;
    return true;
}

/* Updates the size, worklists and locks of each changed surviving column.
   AT must be up to date. */
static void dton_update_changed_cols(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT;
    int *col_sizes = constraints->state->col_sizes;

    for (int ii = 0; ii < dton_work->targets.n_targets; ++ii)
    {
        int col = dton_work->targets.list[ii].col;
        dton_work->targets.col_to_target[col] = -1;
        assert(!HAS_TAG(constraints->col_tags[col], C_TAG_INACTIVE));

        int old_size = col_sizes[col];
        int new_size = AT->p[col].end - AT->p[col].start;
        col_sizes[col] = new_size;
        assert(old_size > 0);

        /* check if the surviving column is now empty or a singleton */
        if (new_size == 0)
        {
            iVec_append(constraints->state->empty_cols, col);
        }
        else if (new_size == 1 && old_size != 1)
        {
            iVec_append(constraints->state->ston_cols, col);
        }

        /* Recount the column's locks. */
        RowView col_view =
            new_rowview(AT->x + AT->p[col].start, AT->i + AT->p[col].start,
                        col_sizes + col, AT->p + col, NULL, NULL, NULL, col);
        count_locks_one_column(&col_view, constraints->state->col_locks + col,
                               constraints->row_tags);
    }
}

/* Brings AT up to date (so it is consistent with A). Then updates the size and
 * locks of each changed column. */
void dton_update_AT(Problem *prob, DtonWorkspace *dton_work)
{
    /* try to merge the changes into AT */
    if (dton_work->merge_round)
    {
        dton_work->merge_round = dton_merge_into_AT(prob, dton_work);
    }

    /* if the merge round failed, we need to rebuild AT from scratch */
    if (!dton_work->merge_round)
    {
        dton_rebuild_AT(prob, dton_work);
    }
    assert(prob->constraints->AT->nnz == prob->constraints->A->nnz);

    dton_update_changed_cols(prob, dton_work);
}
