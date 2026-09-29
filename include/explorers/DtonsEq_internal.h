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

#ifndef DTONS_EQ_INTERNAL_H
#define DTONS_EQ_INTERNAL_H

/* Workspace and round kernels of the doubleton eliminator. Internal to
   DtonsEq.c; included by the unit tests. Each phase is documented at its
   definition. */

#include "DtonsEq.h"
#include "RowSlots.h"
#include "SparseAccumulator.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
struct Problem;
struct Work;

/* One substitution of the round: row 'owner' eliminates column k. Appended in
   claim order; ws->substs.col_subst[k] indexes it (-1 when k is not claimed). */
typedef struct DtonSubst
{
    int k;           /* substituted column */
    int owner;       /* owner row */
    int j;           /* direct stay column */
    double dir_mult; /* x_k = dir_shift + dir_mult * x_j */
    double dir_shift;
    int target;  /* composed survivor of the round */
    double mult; /* x_k = shift + mult * x_target */
    double shift;
} DtonSubst;

/* Phases 1-2: the round's substitutions in claim order, compacted after
   composition to the still-eliminated columns. Reserved exactly per round
   (one slot per worklist row) by dton_claim. */
typedef struct DtonSubsts
{
    DtonSubst *recs; /* [cap] */
    int *depth;      /* [cap] chain depth after composition */
    int *order;      /* [cap] record indices by descending chain depth */
    int *succ; /* [cap] compute_chain_depths input: record of the stay column */
    int *drop_priority; /* [cap] compute_chain_depths input: owner row */
    int *stamp;         /* [cap] compute_chain_depths scratch */
    int cap;            /* also sizes DtonTargets.list/old_size and DtonLog.start */
    int n;
    int maxdepth;   /* max chain depth this round */
    int *col_subst; /* [n] record index of a claimed column, -1 otherwise
                       (borrowed: iwork1_max_nrows_ncols) */
} DtonSubsts;

/* Phase 4: the rows of the eliminated columns, read from the pre-round AT
   (one entry per column entry, so with duplicates) and sorted (idx is the
   permutation, aux the sort scratch). Reserved per round with the exact
   count by dton_reserve_rows. */
typedef struct DtonRows
{
    int *list;
    int *idx;
    int *aux;
    int cap;
    int n;
} DtonRows;

/* Phases 4-5: the unique composed targets with their pre-round sizes, for
   the size transitions (empty/singleton column worklists) and the lock
   recount. */
typedef struct DtonTargets
{
    int *list;     /* [DtonSubsts.cap] */
    int *old_size; /* [DtonSubsts.cap] parallel to list */
    int n;
    int *col_to_target; /* [n] index into list, -1 otherwise; reset through list
                        (borrowed: iwork2_max_nrows_ncols) */
} DtonTargets;

/* Phases 4-5: change log of the sweep, one (target, row, value) tuple for
   every target entry a swept row ends up with, and value 0 for a present
   target entry that was dropped (A holds no exact zeros, so 0 marks a
   delete); rows ascending. The refresh sorts it by target into row2/val2
   (stable, so rows stay ascending) and applies it onto the old column. */
typedef struct DtonLog
{
    int *col;
    int *row;
    double *val;
    int *row2;
    double *val2;
    int len;
    int alloc;
    bool overflow; /* an append failed: the round falls back to a rebuild */
    int *start;    /* [DtonSubsts.cap + 1] segment starts of the sorted log per
                      target */
} DtonLog;

/* Tunables and diagnostics, visible to the tests; nothing depends on them
   for correctness. */
typedef struct DtonTuning
{
    double rebuild_dirty_frac; /* rebuild when dirty content > frac * nnz (0.25) */
    bool last_round_rebuilt;   /* the last round took the rebuild path */
} DtonTuning;

/* Scratch for the doubleton eliminator (work->dton). Only 'at' persists across
   calls; the rest is per round or per row and is reset through the round's own
   lists. The per-round arrays are reserved by the first round and the borrowed
   arrays are unset until dton_ws_attach. */
typedef struct DtonWorkspace
{
    size_t m; /* rows of A at allocation */
    size_t n; /* cols of A at allocation */
    DtonSubsts substs;
    DtonRows rows;
    /* Phase 4: the row being rewritten. stamp is borrowed from radix_aux and
       touched from iwork_n_cols. touched also serves as the compose path and
       the refresh cursors. */
    SparseAccumulator acc;
    DtonTargets targets;
    RowSlots at;   /* slots of A transpose's rows, persist across calls */
    bool at_valid; /* 'at' matches the layout of A transpose */
    DtonLog log;
    DtonTuning tuning;
} DtonWorkspace;

/* Points the borrowed per-column arrays at the presolver's shared scratch and
   initializes them (col_subst and col_to_target to -1, cstamp to 0). Once per
   call of the eliminator, before any kernel runs. */
void dton_ws_attach(DtonWorkspace *ws, struct Work *work);

/* Which entry (0 or 1) of a doubleton equality row to substitute. */
int dton_choose_subst(const double *row_vals, int col_size0, int col_size1);

/* Phase 1: claims the substituted column of every eliminable row. False if
   the round's records could not be reserved. */
bool dton_claim(struct Problem *prob, DtonWorkspace *ws, int *deferred,
                int *n_deferred);

/* Phase 2: composes the substitution chains and breaks cycles. */
void dton_compose(DtonWorkspace *ws, int *deferred, int *n_deferred);

/* Phase 3b: emits the postsolve records of the round. */
void dton_record(struct Problem *prob, DtonWorkspace *ws);

/* Phase 4: applies the round to A, the row sides, worklists and objective. */
void dton_apply(struct Problem *prob, DtonWorkspace *ws, int *deferred,
                int *n_deferred);

/* The full-transpose fallback of the A transpose refresh. */
void dton_rebuild_AT(struct Problem *prob, DtonWorkspace *ws);

#endif
