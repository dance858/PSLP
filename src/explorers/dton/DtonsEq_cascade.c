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

#include "Constraints.h"
#include "DtonsEq_internal.h"
#include "Matrix.h"
#include "Numerics.h"
#include "Problem.h"
#include "State.h"
#include "Tags.h"
#include <string.h>

/* Phase 2b of a round. A row that the sweep of phase 4 would turn into a
   doubleton equality is claimed in this round instead of the next one. */

/* Column k is eliminated this round. Every row holding k loses one more
   column, so the row's count of claimed columns grows by one. A row that is
   left with exactly two unclaimed columns and is an equality joins the
   candidates. The count never shrinks, so a row joins at most once. */
static void dton_walk_column(DtonWorkspace *dton_work,
                             const Constraints *constraints, int k, int *n_cand)
{
    const Matrix *AT = constraints->AT; /* pre-round transpose */
    const int *row_sizes = constraints->state->row_sizes;
    const RowTag *row_tags = constraints->row_tags;
    int *claimed_in_row = dton_work->claimed_in_row;

    for (int jj = AT->p[k].start; jj < AT->p[k].end; ++jj)
    {
        int s = AT->i[jj];
        assert(!HAS_TAG(row_tags[s], R_TAG_INACTIVE));
        claimed_in_row[s]++;
        if (row_sizes[s] - claimed_in_row[s] == 2 && HAS_TAG(row_tags[s], R_TAG_EQ))
        {
            dton_work->cand_rows[(*n_cand)++] = s;
        }
    }
}

/* Claims row s if the round's substitutions turn it into a doubleton equality
   whose chosen column can be eliminated this round. Returns the index of the new
   record, or -1 when the row waits for the next round. */
static int dton_cascade_claim(Problem *prob, DtonWorkspace *dton_work, int s)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const int *col_sizes = constraints->state->col_sizes;
    DtonSubsts *substs = &dton_work->substs;
    const int *col_to_target = dton_work->targets.col_to_target;
    SparseAccumulator *acc = &dton_work->acc;

    assert(constraints->lhs[s] == constraints->rhs[s] &&
           !IS_ABS_INF(constraints->rhs[s]));

    /* the row as the sweep will leave it: a doubleton keeps two entries */
    double side_shift = dton_accumulate_row(dton_work, A, s);
    int cols[2];
    int n_kept = 0;
    for (int ii = 0; ii < acc->n_touched; ++ii)
    {
        int col = acc->touched[ii];
        bool modified = acc->flags[col] == DTON_ACC_MODIFIED;
        if (!modified || ABS(acc->value[col]) > ZERO_TOL)
        {
            if (n_kept < 2)
            {
                cols[n_kept] = col;
            }
            n_kept++;
        }
    }
    if (n_kept != 2)
    {
        return -1;
    }
    if (cols[0] > cols[1])
    {
        int tmp = cols[0];
        cols[0] = cols[1];
        cols[1] = tmp;
    }
    double vals[2] = {acc->value[cols[0]], acc->value[cols[1]]};
    double rhs = constraints->rhs[s] - side_shift;

    bool in_P[2], sub[2];
    double Pkk[2];
    dton_P_flags(prob, cols[0], cols[1], in_P, sub, Pkk);
    int subst = dton_choose_subst_qp(vals, col_sizes[cols[0]], col_sizes[cols[1]],
                                     in_P[0], in_P[1], sub[0], sub[1]);

    /* The chosen column must be unclaimed and no record's surviving column, since
       the postsolve reads it from the pre-round transpose. */
    if (subst < 0 || substs->col_subst[cols[subst]] >= 0 ||
        col_to_target[cols[subst]] >= 0)
    {
        return -1;
    }
    int k = cols[subst];
    assert(acc->flags[k] == 0); /* pristine coefficient */

    if (substs->n_recs == substs->cap &&
        !dton_reserve_records(dton_work, 2 * (size_t) substs->cap))
    {
        return -1;
    }
    /* the row's claimed columns were claimed earlier, so the record ranks above
       every record before it */
    int idx = dton_push_record(dton_work, s, k, cols[1 - subst], vals[subst],
                               vals[1 - subst], rhs, Pkk[subst]);
    dton_work->targets.col_to_target[cols[1 - subst]] = 0;
    assert(substs->n_first > 0);
    int below = idx == substs->n_first ? substs->depth[substs->order[0]]
                                       : substs->depth[idx - 1];
    substs->depth[idx] = below + 1;
    return idx;
}

/* The surviving columns are marked in col_to_target while the cascade runs.
   claimed_in_row keeps the counts of the round for dton_apply. */
void dton_cascade(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const int *row_sizes = constraints->state->row_sizes;
    DtonSubsts *substs = &dton_work->substs;
    int *col_to_target = dton_work->targets.col_to_target;
    int n_first = substs->n_recs;
    int n_cand = 0;

    for (int idx = 0; idx < n_first; ++idx)
    {
        col_to_target[substs->recs[idx].target] = 0;
        dton_walk_column(dton_work, constraints, substs->recs[idx].k, &n_cand);
    }

    while (n_cand > 0)
    {
        int s = dton_work->cand_rows[--n_cand];
        /* the count may have grown since the row was stacked */
        if (row_sizes[s] - dton_work->claimed_in_row[s] != 2)
        {
            continue;
        }
        int idx = dton_cascade_claim(prob, dton_work, s);
        if (idx >= 0)
        {
            dton_walk_column(dton_work, constraints, substs->recs[idx].k, &n_cand);
        }
    }

    /* reset the marks through the records. dton_apply zeroes the counts when
       it walks the same columns. */
    for (int idx = 0; idx < substs->n_recs; ++idx)
    {
        col_to_target[substs->recs[idx].target] = -1;
    }

    /* the cascade records lead the order, latest first */
    int n_cascade = substs->n_recs - n_first;
    memmove(substs->order + n_cascade, substs->order,
            (size_t) n_first * sizeof(int));
    for (int ii = 0; ii < n_cascade; ++ii)
    {
        substs->order[ii] = substs->n_recs - 1 - ii;
    }
}
