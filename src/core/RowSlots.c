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

bool matrix_update_row(Matrix *M, RowSlots *slots, int r, const int *cols,
                       const double *vals, int n, const uint8_t *tags,
                       uint8_t drop_bits)
{
    int start = M->p[r].start;
    int o = start;
    int oe = M->p[r].end;
    int u = 0;

    // merge the old entries and the updates into the free tail
    int d0 = slots->tail_next;
    int d = d0;
    while (o < oe || u < n)
    {
        int c_old = (o < oe) ? M->i[o] : INT_MAX;
        int c_upd = (u < n) ? cols[u] : INT_MAX;
        int c;
        double v;
        if (c_old < c_upd)
        {
            c = c_old;
            v = M->x[o++];
            if (dropped(tags, drop_bits, c))
            {
                continue;
            }
        }
        else
        {
            c = c_upd;
            v = vals[u++];
            if (c_old == c_upd)
            {
                o++; // overwritten
            }
            if (v == 0.0)
            {
                continue; // deleted
            }
        }
        if ((size_t) d == M->n_alloc)
        {
            return false; // tail full: only free space was written
        }
        M->i[d] = c;
        M->x[d] = v;
        d++;
    }

    // a result that fits the row's slot goes back into it; otherwise the
    // merged entries already are the row's new slot
    int len = d - d0;
    if (len <= slots->cap[r])
    {
        memcpy(M->i + start, M->i + d0, (size_t) len * sizeof(int));
        memcpy(M->x + start, M->x + d0, (size_t) len * sizeof(double));
        d0 = start;
    }
    else
    {
        slots->tail_next = d0 + len;
        slots->cap[r] = len;
    }
    M->p[r].start = d0;
    M->p[r].end = d0 + len;
    M->nnz = M->nnz + (size_t) len - (size_t) (oe - start);
    return true;
}
