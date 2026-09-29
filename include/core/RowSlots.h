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

#ifndef CORE_ROWSLOTS_H
#define CORE_ROWSLOTS_H

#include "Matrix.h"
#include <stdbool.h>
#include <stdint.h>

/* Room reserved for each row of a matrix whose rows may grow. A row that
   outgrows its slot moves to a fresh slot in the spare tail
   [tail_base, n_alloc), filled from tail_next.

   Rows are not reordered when one moves. The matrix then has its rows out of
   index order and gaps between them, so a row must be read through its own
   range [p[r].start, p[r].end). Only rebuilding the matrix packs the rows in
   index order again. */
typedef struct RowSlots
{
    int *cap; /* [m] room reserved for each row, owned by the caller */
    int tail_base;
    int tail_next;
} RowSlots;

/* Sets the slots from the layout of M, whose rows must be contiguous and in
   index order with p[m].start marking their end (as built by transpose). */
void row_slots_init(RowSlots *slots, const Matrix *M);

/* Replaces row r of M by the merge of its entries and n updates (cols
   ascending). A zero value deletes the entry, any other value inserts or
   overwrites it. Deleting an absent entry does nothing. Entries whose column
   c has tags[c] & drop_bits are dropped (tags may be NULL). Writes in place
   when the result fits the row's slot, otherwise into a fresh slot in the
   tail, and updates M->nnz. The merge is built in the free tail, so it needs
   room there for the result even when the row stays in its slot. Returns
   false, leaving the row untouched, if the tail has no room. */
bool matrix_update_row(Matrix *M, RowSlots *slots, int r, const int *cols,
                       const double *vals, int n, const uint8_t *tags,
                       uint8_t drop_bits);

#endif // CORE_ROWSLOTS_H
