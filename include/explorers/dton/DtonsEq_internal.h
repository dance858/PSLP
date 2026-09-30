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

/* Workspace and round kernels of the doubleton eliminator. Internal to the
   DtonsEq*.c files and included by the unit tests. */

#include "DtonsEq.h"
#include "RowSlots.h"
#include "SparseAccumulator.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
struct Problem;
struct Work;

/* One substitution of the round: row 'owner' eliminates column k. Appended in
   claim order. substs.col_subst[k] indexes it (-1 when k is not claimed). */
typedef struct DtonSubst
{
    int k;      /* substituted column */
    int owner;  /* owner row */
    int j;      /* direct stay column */
    int target; /* composed survivor of the round */
    double aik; /* owner row: aij * x_j + aik * x_k = rhs */
    double aij;
    double dir_mult; /* x_k = dir_shift + dir_mult * x_j */
    double dir_shift;
    double mult; /* x_k = shift + mult * x_target */
    double shift;
} DtonSubst;

/* Phases 1-2: the round's substitutions in claim order, compacted after
   composition to the still-eliminated columns. */
typedef struct DtonSubsts
{
    DtonSubst *recs; /* [cap] */
    int *depth;      /* [cap] chain depth after composition */
    int *order;      /* [cap] record indices by descending chain depth */
    int *succ; /* [cap] compute_chain_depths input: record of the stay column */
    int *drop_priority; /* [cap] compute_chain_depths input: owner row */
    int *stamp;         /* [cap] compute_chain_depths scratch */
    int cap;            /* also sizes DtonTargets.list */
    int n_recs;
    int *col_subst; /* [n_cols] record index of a claimed column, -1 otherwise. Reset
                       through recs. */
} DtonSubsts;

/* Phases 4-5: a composed target of the round, with its pre-round size for the
   size transitions (empty/singleton column worklists). Its part of the change
   log is [log_start, log_end). */
typedef struct DtonTarget
{
    int col;
    int old_size;
    int log_start;
    int log_end;
} DtonTarget;

typedef struct DtonTargets
{
    DtonTarget *list; /* [DtonSubsts.cap] the unique targets */
    int n_targets;
    int *col_to_target; /* [n_cols] index into list, -1 otherwise. Reset through
                           list. */
} DtonTargets;

/* Phases 4-5: change log of the sweep, one (row, value) tuple for every
   target entry a swept row ends up with, and value 0 for a present target
   entry that was dropped (A holds no exact zeros, so 0 marks a delete). Every
   target has its own segment with rows ascending, which the merge applies
   onto the old column. */
typedef struct DtonLog
{
    int *row;
    double *val;
    int cap;
} DtonLog;

/* Scratch for the doubleton eliminator (work->dton). Only 'AT' persists across
   calls. The rest is per round or per row and is reset through the round's own
   lists. */
typedef struct DtonWorkspace
{
    int n_rows; /* rows of A at allocation */
    int n_cols; /* cols of A at allocation */
    DtonSubsts substs;
    DtonTargets targets;
    uint64_t *swept_rows; /* [ceil(n_rows / 64)] bitmap of the rows the round sweeps.
                             Zero between rounds. */
    SparseAccumulator acc;
    DtonLog log;
    RowSlots AT;      /* slots of A transpose's rows */
    bool AT_valid;    /* 'AT' matches the layout of A transpose */
    bool merge_round; /* the round logs its changes and merges them into A
                         transpose. Otherwise A transpose is rebuilt. */
    /* For the tests. Nothing depends on it for correctness. */
    double rebuild_dirty_frac; /* rebuild when dirty content > frac * nnz (0.25) */
} DtonWorkspace;

/* Number of 64-bit words in a bitmap of n_bits bits. */
static inline int bitmap_words(int n_bits)
{
    return (n_bits + 63) / 64;
}

/* Points the borrowed per-column arrays at the presolver's shared scratch and
   initializes them. Once per call of the eliminator, before any kernel runs. */
void dton_workspace_init(DtonWorkspace *dton_work, struct Work *work);

/* Which entry (0 or 1) of a doubleton equality row to substitute. */
int dton_choose_subst(const double *row_vals, int col_size0, int col_size1);

/* Phase 1: claims the substituted column of every eliminable row. False if
   the round's records could not be reserved. */
bool dton_claim(struct Problem *prob, DtonWorkspace *dton_work, int *deferred,
                int *n_deferred);

/* Phase 2: composes the substitution chains and breaks cycles. */
void dton_compose(DtonWorkspace *dton_work, int *deferred, int *n_deferred);

/* Phase 3: transfers the bounds of the eliminated columns onto their stay
   columns. INFEASIBLE if a transfer contradicts a bound. */
PresolveStatus dton_transfer_bounds(struct Problem *prob, DtonWorkspace *dton_work);

/* Phase 3b: emits the postsolve records of the round. */
void dton_record(struct Problem *prob, DtonWorkspace *dton_work);

/* Phase 4: applies the round to A, the row sides, worklists and objective. */
void dton_apply(struct Problem *prob, DtonWorkspace *dton_work, int *deferred,
                int *n_deferred);

/* The full-transpose fallback of dton_update_AT. */
void dton_rebuild_AT(struct Problem *prob, DtonWorkspace *dton_work);

/* Phase 5: brings A transpose up to date. */
void dton_update_AT(struct Problem *prob, DtonWorkspace *dton_work);

#endif
