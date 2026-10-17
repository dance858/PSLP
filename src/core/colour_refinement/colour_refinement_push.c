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

/* The push.

   After a half-round of side 'from', every line of a non-singleton class of
   the other side 'to' has the signature against the classes 'from' had
   before that half-round. Only the changed lines of 'from' moved, so the
   signatures of 'to' are brought up to date by replacing, for every entry of a
   changed line, the hash of its old class by the hash of its new class. The
   lines reached are the frontier of 'to': the only lines whose signature can
   differ from that of their class.

   A line reached for the first time gets bit 0 of its signature set, so it
   enters the frontier once. The candidate passes clear it. */

#include "colour_refinement_internal.h"
#include <string.h>

/* A push below this many changed lines runs on the calling thread. */
#define SMALL_PUSH (1 << 14)
#define PUSH_LINES_PER_CHUNK 256

typedef struct PushPass
{
    const RefineSide *from;
    RefineSide *to;
    RefineWorkspace *work;
    int n_chunks;
} PushPass;

/* Appends the 'count' buffered frontier lines of a thread to the frontier. */
static void flush_push_buffer(PushPass *pass, int thread, int count)
{
    if (count == 0)
    {
        return;
    }
    int start = ps_atomic_fetch_add_int(&pass->to->n_frontier, count);
    memcpy(pass->to->frontier + start, pass->work->push_buffers[thread],
           (size_t) count * sizeof(int));
}

static void push_chunk(void *context, int chunk, int thread)
{
    PushPass *pass = (PushPass *) context;
    const RefineSide *from = pass->from;
    uint64_t *signature = pass->to->signature;
    int *buffer = pass->work->push_buffers[thread];

    /* the lines the thread buffered in its earlier chunks */
    int count = pass->work->push_counts[thread];
    int start = (int) ((int64_t) from->n_changed * chunk / pass->n_chunks);
    int end = (int) ((int64_t) from->n_changed * (chunk + 1) / pass->n_chunks);

    for (int ii = start; ii < end; ++ii)
    {
        int line = from->changed[ii];
        int old_class = from->changed_old_class[ii];
        int new_class = from->class_id[line];
        int line_end = from->p[line].end;
        for (int idx = from->p[line].start; idx < line_end; ++idx)
        {
            int other = from->i[idx];
            uint64_t bits = value_bits(from->x[idx]);
            uint64_t change =
                entry_hash(new_class, bits) - entry_hash(old_class, bits);
            uint64_t old = ps_atomic_fetch_add_u64(&signature[other], change);

            /* the first thread to set the mark adds the line to the frontier */
            if (has_mark(old))
            {
                continue;
            }
            uint64_t marked =
                ps_atomic_fetch_or_u64(&signature[other], SIGNATURE_MARK);
            if (!has_mark(marked))
            {
                buffer[count++] = other;
                if (count == PUSH_BUFFER_SIZE)
                {
                    flush_push_buffer(pass, thread, count);
                    count = 0;
                }
            }
        }
    }
    pass->work->push_counts[thread] = count;
}

void refine_push(const RefineSide *from, RefineSide *to, RefineWorkspace *work,
                 ThreadPool *pool)
{
    to->n_frontier = 0;
    work->work_done += from->changed_entries;

    if (from->n_changed < SMALL_PUSH || work->n_threads == 1)
    {
        uint64_t *signature = to->signature;
        int *frontier = to->frontier;
        int n_changed = from->n_changed;
        int n_frontier = 0;
        for (int ii = 0; ii < n_changed; ++ii)
        {
            int line = from->changed[ii];
            int old_class = from->changed_old_class[ii];
            int new_class = from->class_id[line];
            int line_end = from->p[line].end;
            for (int idx = from->p[line].start; idx < line_end; ++idx)
            {
                int other = from->i[idx];
                uint64_t bits = value_bits(from->x[idx]);
                uint64_t change =
                    entry_hash(new_class, bits) - entry_hash(old_class, bits);
                uint64_t old = signature[other];
                signature[other] = (old + change) | SIGNATURE_MARK;
                if (!has_mark(old))
                {
                    frontier[n_frontier++] = other;
                }
            }
        }
        to->n_frontier = n_frontier;
        return;
    }

    PushPass pass = {from, to, work, from->n_changed / PUSH_LINES_PER_CHUNK + 1};
    memset(work->push_counts, 0, (size_t) work->n_threads * sizeof(int));
    thread_pool_run(pool, push_chunk, &pass, pass.n_chunks);
    for (int thread = 0; thread < work->n_threads; ++thread)
    {
        flush_push_buffer(&pass, thread, work->push_counts[thread]);
    }
}
