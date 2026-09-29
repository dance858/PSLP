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

#include "RowSlots.h"
#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include <string.h>

void row_slots_init(RowSlots *slots, const Matrix *M)
{
    for (int r = 0; r < (int) M->m; ++r)
    {
        slots->cap[r] = M->p[r + 1].start - M->p[r].start;
        assert(slots->cap[r] >= M->p[r].end - M->p[r].start);
    }
    slots->tail_base = M->p[M->m].start; // the end of the rows
    slots->tail_next = slots->tail_base;
}

static inline bool dropped(const uint8_t *tags, uint8_t drop_bits, int c)
{
    return tags != NULL && (tags[c] & drop_bits) != 0;
}

bool matrix_update_row(Matrix *M, RowSlots *slots, int row, const int *cols,
                       const double *vals, int n, const uint8_t *tags,
                       uint8_t drop_bits)
{
    int start = M->p[row].start;
    int old = start;
    int old_end = M->p[row].end;
    int k = 0;

    // merge the old entries and the updates into the free tail
    int merged_start = slots->tail_next;
    int merged_end = merged_start;
    while (old < old_end || k < n)
    {
        int c_old = (old < old_end) ? M->i[old] : INT_MAX;
        int c_upd = (k < n) ? cols[k] : INT_MAX;
        int c;
        double v;
        if (c_old < c_upd)
        {
            c = c_old;
            v = M->x[old++];
            if (dropped(tags, drop_bits, c))
            {
                continue;
            }
        }
        else
        {
            c = c_upd;
            v = vals[k++];
            if (c_old == c_upd)
            {
                old++; // overwritten
            }
            if (v == 0.0)
            {
                continue; // deleted
            }
        }
        if ((size_t) merged_end == M->n_alloc)
        {
            return false; // tail full: only free space was written
        }
        M->i[merged_end] = c;
        M->x[merged_end] = v;
        merged_end++;
    }

    // a result that fits the row's slot goes back into it; otherwise the
    // merged entries already are the row's new slot
    int len = merged_end - merged_start;
    if (len <= slots->cap[row])
    {
        memcpy(M->i + start, M->i + merged_start, (size_t) len * sizeof(int));
        memcpy(M->x + start, M->x + merged_start, (size_t) len * sizeof(double));
        merged_start = start;
    }
    else
    {
        slots->tail_next = merged_start + len;
        slots->cap[row] = len;
    }
    M->p[row].start = merged_start;
    M->p[row].end = merged_start + len;
    M->nnz = M->nnz + (size_t) len - (size_t) (old_end - start);
    return true;
}
