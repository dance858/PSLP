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

DtonWorkspace *dton_ws_new(size_t n_rows, size_t n_cols)
{
    DtonWorkspace *ws = (DtonWorkspace *) ps_calloc(1, sizeof(DtonWorkspace));
    RETURN_PTR_IF_NULL(ws, NULL);
    ws->m = n_rows;
    ws->n = n_cols;

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

    ws->rebuild_dirty_frac = 0.25;
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

PresolveStatus remove_dton_eq_rows(Problem *prob)
{
    Constraints *constraints = prob->constraints;
    State *state = constraints->state;
    Work *work = state->work;

    assert(state->ston_rows->len == 0);
    assert(state->empty_rows->len == 0);
    assert(state->empty_cols->len == 0);
    DEBUG(verify_problem_up_to_date(constraints));

    assert(work->dton != NULL);

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
            progress = false;
        }

        if (progress)
        {
            // an infeasible transfer ends the presolve
            if (dton_transfer_bounds(prob, ws) == INFEASIBLE)
            {
                return INFEASIBLE;
            }
            dton_record(prob, ws);
            dton_apply(prob, ws, deferred, &n_deferred);
            dton_refresh_AT(prob, ws);

            // reset the claim state. The eliminated columns are inactive
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

    DEBUG(verify_problem_up_to_date(constraints));

    return UNCHANGED;
}
