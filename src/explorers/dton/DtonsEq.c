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
#include "Debugger.h"
#include "DtonsEq_internal.h"
#include "Memory_wrapper.h"
#include "Problem.h"
#include "State.h"
#include "Workspace.h"

#define REBUILD_FRAC 0.25

/* allocation but no initialization */
DtonWorkspace *dton_workspace_new(size_t n_rows, size_t n_cols)
{
    DtonWorkspace *dton_work = (DtonWorkspace *) ps_calloc(1, sizeof(DtonWorkspace));
    RETURN_PTR_IF_NULL(dton_work, NULL);
    dton_work->n_rows = (int) n_rows;
    dton_work->n_cols = (int) n_cols;
    dton_work->rebuild_frac = REBUILD_FRAC;

    dton_work->swept_rows = (uint64_t *) ps_calloc(
        (size_t) bitmap_words(dton_work->n_rows), sizeof(uint64_t));
    dton_work->acc.value = (double *) ps_malloc(n_cols, sizeof(double));
    dton_work->acc.flags = (uint8_t *) ps_malloc(n_cols, sizeof(uint8_t));
    dton_work->AT_slots.cap = (int *) ps_malloc(n_cols, sizeof(int));

    if (!dton_work->swept_rows || !dton_work->acc.value || !dton_work->acc.flags ||
        !dton_work->AT_slots.cap)
    {
        dton_workspace_free(dton_work);
        return NULL;
    }

    return dton_work;
}

void dton_workspace_free(DtonWorkspace *dton_work)
{
    RETURN_IF_NULL(dton_work);
    PS_FREE(dton_work->substs.recs);
    PS_FREE(dton_work->substs.depth);
    PS_FREE(dton_work->substs.order);
    PS_FREE(dton_work->substs.succ);
    PS_FREE(dton_work->substs.drop_priority);
    PS_FREE(dton_work->substs.stamp);
    PS_FREE(dton_work->targets.list);
    PS_FREE(dton_work->swept_rows);
    PS_FREE(dton_work->acc.value);
    PS_FREE(dton_work->acc.flags);
    PS_FREE(dton_work->AT_slots.cap);
    PS_FREE(dton_work->log.row);
    PS_FREE(dton_work->log.val);
    PS_FREE(dton_work);
}

/* no allocation but initialization */
void dton_workspace_init(DtonWorkspace *dton_work, Work *work)
{
    int n_cols = dton_work->n_cols;
    int *col_subst = work->iwork1_max_nrows_ncols;
    int *col_to_target = work->iwork2_max_nrows_ncols;
    for (int k = 0; k < n_cols; ++k)
    {
        col_subst[k] = -1;
        col_to_target[k] = -1;
    }
    dton_work->substs.col_subst = col_subst;
    dton_work->targets.col_to_target = col_to_target;
    sparse_accumulator_init(&dton_work->acc, work->radix_aux, dton_work->acc.value,
                            dton_work->acc.flags, work->iwork_n_cols,
                            (size_t) dton_work->n_cols);
}

PresolveStatus remove_dton_eq_rows(Problem *prob)
{
    Constraints *constraints = prob->constraints;
    State *state = constraints->state;
    Work *work = state->work;
    DtonWorkspace *dton_work = work->dton;
    dton_workspace_init(dton_work, work);
    iVec *dton_rows = state->dton_rows;
    int *deferred = work->iwork_n_rows;

    assert(state->ston_rows->len == 0);
    assert(state->empty_rows->len == 0);
    assert(state->empty_cols->len == 0);
    DEBUG(verify_problem_up_to_date(constraints));

    while (dton_rows->len > 0)
    {
        int n_deferred = 0;

        DEBUG(verify_no_duplicates_sort(dton_rows));

        /* determine which columns to subtitute */
        bool progress = dton_claim(prob, dton_work, deferred, &n_deferred);
        dton_compose(dton_work, deferred, &n_deferred);
        progress = progress && dton_work->substs.n_recs > 0;

        if (progress)
        {
            /* transfer bounds for substituted columns */
            if (dton_transfer_bounds(prob, dton_work) == INFEASIBLE)
            {
                return INFEASIBLE;
            }

            /* save postsolve information */
            dton_record(prob, dton_work);

            /* XXX some comment here */
            dton_apply(prob, dton_work, deferred, &n_deferred);
            dton_update_AT(prob, dton_work);

            /* clear the substituted columns for the next round. */
            for (int idx = 0; idx < dton_work->substs.n_recs; ++idx)
            {
                dton_work->substs.col_subst[dton_work->substs.recs[idx].k] = -1;
            }
            dton_work->substs.n_recs = 0;
        }

        /* refresh the list of dton_rows for the next round */
        iVec_clear_no_resize(dton_rows);
        if (n_deferred > 0)
        {
            DEBUG(verify_no_duplicates_sort_ptr(deferred, (size_t) n_deferred));
            iVec_append_array(dton_rows, deferred, (size_t) n_deferred);
        }
    }

    DEBUG(verify_problem_up_to_date(constraints));

    return UNCHANGED;
}
