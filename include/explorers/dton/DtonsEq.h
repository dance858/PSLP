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

#ifndef DTONS_EQ_H
#define DTONS_EQ_H

#include "PSLP_status.h"
#include <stddef.h>
struct Problem;
struct DtonWorkspace;

/* Spare entries reserved after the rows of A transpose, where the eliminator
   relocates columns that grow. */
static inline size_t dton_extra_memory(size_t nnz)
{
    return nnz / 16 > 1024 ? nnz / 16 : 1024;
}

/* Allocates the eliminator's workspace (work->dton). Returns NULL if any
   allocation fails. */
struct DtonWorkspace *dton_ws_new(size_t n_rows, size_t n_cols);

/* Frees the eliminator's workspace (work->dton). */
void dton_ws_free(struct DtonWorkspace *dton_work);

/* Eliminates the doubleton equality rows on state->dton_rows, round by round.
   Returns INFEASIBLE when a transferred bound contradicts the bounds of the
   column that stays, UNCHANGED otherwise. */
PresolveStatus remove_dton_eq_rows(struct Problem *prob);

#endif
