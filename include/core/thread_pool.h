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

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <stdint.h>

/* 1 if the pool can have worker threads (pthreads), 0 on Windows. */
#if defined(_WIN32) || defined(_WIN64)
#define THREAD_POOL_HAS_WORKERS 0
#else
#define THREAD_POOL_HAS_WORKERS 1
#endif

/* Thread pool. A pass is split into chunks and 'function' is called once per
   chunk. 'thread' is the index of the thread running the chunk (0 is the
   calling thread), which a pass can use to pick per-thread scratch. Chunks are
   handed out in ascending order, so a chunk may wait for earlier chunks to
   finish. Without workers every chunk runs on the calling thread. */
typedef void (*ThreadPoolChunkFunction)(void *context, int chunk, int thread);

typedef struct ThreadPool ThreadPool;

/* Pool with 'n_workers' worker threads (0 for none). NULL on failure. */
ThreadPool *thread_pool_new(int n_workers);
void thread_pool_free(ThreadPool *pool);

/* Workers plus the calling thread. */
int thread_pool_n_threads(const ThreadPool *pool);

/* Runs the chunks of one pass and returns when all of them are done. */
void thread_pool_run(ThreadPool *pool, ThreadPoolChunkFunction function,
                     void *context, int n_chunks);

/* Chunks of a pass: at most eight per thread, each with at least
   'items_per_chunk' items. */
int thread_pool_chunk_count(int n_threads, int n_items, int items_per_chunk);

/* bounds[0..n_chunks] split [0, n_items) into even parts. */
void thread_pool_chunk_bounds(int n_items, int n_chunks, int *bounds);

/* bounds[0..n_chunks] split [0, n_items) so that the non-decreasing prefix
   (n_items + 1 entries) grows evenly from chunk to chunk. */
void thread_pool_chunk_bounds_by_prefix(const int *prefix, int n_items, int n_chunks,
                                        int *bounds);

/* Turns per-chunk bucket counts, counts[chunk * n_buckets + bucket], into the
   write positions of a scatter and writes the bucket offsets (n_buckets + 1
   entries) when 'offsets' is given. Returns the total. */
int thread_pool_bucket_offsets(int *counts, int n_chunks, int n_buckets,
                               int *offsets);

/* Atomic operations for chunks that write to shared data. Loads are acquire
   and stores release; the fetch operations are relaxed. Without workers they
   are plain operations. */
#if !THREAD_POOL_HAS_WORKERS

static inline uint64_t ps_atomic_fetch_add_u64(uint64_t *p, uint64_t v)
{
    uint64_t old = *p;
    *p = old + v;
    return old;
}

static inline uint64_t ps_atomic_fetch_or_u64(uint64_t *p, uint64_t v)
{
    uint64_t old = *p;
    *p = old | v;
    return old;
}

static inline int ps_atomic_fetch_add_int(int *p, int v)
{
    int old = *p;
    *p = old + v;
    return old;
}

static inline int ps_atomic_load_int(const int *p)
{
    return *p;
}

static inline void ps_atomic_store_int(int *p, int v)
{
    *p = v;
}

#else

static inline uint64_t ps_atomic_fetch_add_u64(uint64_t *p, uint64_t v)
{
    return __atomic_fetch_add(p, v, __ATOMIC_RELAXED);
}

static inline uint64_t ps_atomic_fetch_or_u64(uint64_t *p, uint64_t v)
{
    return __atomic_fetch_or(p, v, __ATOMIC_RELAXED);
}

static inline int ps_atomic_fetch_add_int(int *p, int v)
{
    return __atomic_fetch_add(p, v, __ATOMIC_RELAXED);
}

static inline int ps_atomic_load_int(const int *p)
{
    return __atomic_load_n(p, __ATOMIC_ACQUIRE);
}

static inline void ps_atomic_store_int(int *p, int v)
{
    __atomic_store_n(p, v, __ATOMIC_RELEASE);
}

#endif

#endif // THREAD_POOL_H
