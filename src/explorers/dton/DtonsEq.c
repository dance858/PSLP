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
    DtonWorkspace *dton_work = (DtonWorkspace *) ps_calloc(1, sizeof(DtonWorkspace));
    RETURN_PTR_IF_NULL(dton_work, NULL);
    dton_work->m = (int) n_rows;
    dton_work->n = (int) n_cols;
    dton_work->rebuild_dirty_frac = 0.25;

    dton_work->swept_rows =
        (uint64_t *) ps_calloc((n_rows + 63) / 64, sizeof(uint64_t));
    dton_work->acc.value = (double *) ps_malloc(n_cols, sizeof(double));
    dton_work->acc.flags = (uint8_t *) ps_malloc(n_cols, sizeof(uint8_t));
    dton_work->at.cap = (int *) ps_malloc(n_cols, sizeof(int));

    if (!dton_work->swept_rows || !dton_work->acc.value || !dton_work->acc.flags ||
        !dton_work->at.cap)
    {
        dton_ws_free(dton_work);
        return NULL;
    }

    return dton_work;
}

void dton_ws_free(DtonWorkspace *dton_work)
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
    PS_FREE(dton_work->at.cap);
    PS_FREE(dton_work->log.row);
    PS_FREE(dton_work->log.val);
    PS_FREE(dton_work);
}

void dton_ws_attach(DtonWorkspace *dton_work, Work *work)
{
    int n = dton_work->n;
    int *col_subst = work->iwork1_max_nrows_ncols;
    int *col_to_target = work->iwork2_max_nrows_ncols;
    for (int k = 0; k < n; ++k)
    {
        col_subst[k] = -1;
        col_to_target[k] = -1;
    }
    dton_work->substs.col_subst = col_subst;
    dton_work->targets.col_to_target = col_to_target;
    sparse_accumulator_init(&dton_work->acc, work->radix_aux, dton_work->acc.value,
                            dton_work->acc.flags, work->iwork_n_cols,
                            (size_t) dton_work->n);
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

    DtonWorkspace *dton_work = work->dton;
    dton_ws_attach(dton_work, work);

    iVec *dton_rows = state->dton_rows;
    int *deferred = work->iwork_n_rows;

    while (dton_rows->len > 0)
    {
        int n_deferred = 0;

        DEBUG(verify_no_duplicates_sort(dton_rows));

        /* A failed reservation gives the round up before anything is
           mutated. */
        bool progress = dton_claim(prob, dton_work, deferred, &n_deferred);
        dton_compose(dton_work, deferred, &n_deferred);
        progress = progress && dton_work->substs.n > 0;

        if (progress)
        {
            /* An infeasible transfer ends the presolve. */
            if (dton_transfer_bounds(prob, dton_work) == INFEASIBLE)
            {
                return INFEASIBLE;
            }
            dton_record(prob, dton_work);
            dton_apply(prob, dton_work, deferred, &n_deferred);
            dton_update_AT(prob, dton_work);

            /* Clear the substitutions for the next round. */
            for (int idx = 0; idx < dton_work->substs.n; ++idx)
            {
                dton_work->substs.col_subst[dton_work->substs.recs[idx].k] = -1;
            }
            dton_work->substs.n = 0;
        }

        /* The worklist becomes the deferred rows plus the sweep's new doubleton
           candidates. */
        iVec_clear_no_resize(dton_rows);
        if (n_deferred > 0)
        {
            DEBUG(verify_no_duplicates_sort_ptr(deferred, (size_t) n_deferred));
            iVec_append_array(dton_rows, deferred, (size_t) n_deferred);
        }

        if (!progress)
        {
            break; /* Only deferrals: retrying now would spin. */
        }
    }

    DEBUG(verify_problem_up_to_date(constraints));

    return UNCHANGED;
}
