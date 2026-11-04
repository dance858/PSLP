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

/* Thread pool (pthreads). Workers sleep on a condition variable between passes.
   On Windows the pool has no workers. */

#include "thread_pool.h"
#include "Memory_wrapper.h"
#include <stdbool.h>
#include <stddef.h>

#if THREAD_POOL_HAS_WORKERS
#include <pthread.h>
#include <unistd.h>
#endif

struct ThreadPool
{
    int n_workers;
    ThreadPoolChunkFunction function;
    void *context;
    int n_chunks;
    int next_chunk;
#if THREAD_POOL_HAS_WORKERS
    int n_started;
    int n_finished;
    pthread_t *threads;
    pthread_mutex_t mutex;
    pthread_cond_t wake;
    pthread_cond_t done;
    int generation;
    int stop;
#endif
};

static void run_chunks(ThreadPool *pool, int thread)
{
    for (;;)
    {
        int chunk = ps_atomic_fetch_add_int(&pool->next_chunk, 1);
        if (chunk >= pool->n_chunks)
        {
            break;
        }
        pool->function(pool->context, chunk, thread);
    }
}

#if THREAD_POOL_HAS_WORKERS

static void *worker_main(void *arg)
{
    ThreadPool *pool = (ThreadPool *) arg;
    int seen = 0;

    pthread_mutex_lock(&pool->mutex);
    /* Thread indices go to the workers in the order they start. */
    int thread = ++pool->n_started;
    for (;;)
    {
        while (pool->generation == seen && !pool->stop)
        {
            pthread_cond_wait(&pool->wake, &pool->mutex);
        }
        if (pool->stop)
        {
            break;
        }
        seen = pool->generation;
        pthread_mutex_unlock(&pool->mutex);

        run_chunks(pool, thread);

        pthread_mutex_lock(&pool->mutex);
        if (++pool->n_finished == pool->n_workers)
        {
            pthread_cond_signal(&pool->done);
        }
    }
    pthread_mutex_unlock(&pool->mutex);

    return NULL;
}

static void stop_workers(ThreadPool *pool)
{
    pthread_mutex_lock(&pool->mutex);
    pool->stop = 1;
    pthread_cond_broadcast(&pool->wake);
    pthread_mutex_unlock(&pool->mutex);

    for (int ii = 0; ii < pool->n_workers; ++ii)
    {
        pthread_join(pool->threads[ii], NULL);
    }
}

/* Starts the workers. Threads that fail to start are left out of the pool. */
static bool start_workers(ThreadPool *pool, int n_workers)
{
    pool->threads = (pthread_t *) ps_malloc((size_t) n_workers, sizeof(pthread_t));
    if (!pool->threads)
    {
        return false;
    }
    pthread_mutex_init(&pool->mutex, NULL);
    pthread_cond_init(&pool->wake, NULL);
    pthread_cond_init(&pool->done, NULL);

    for (int ii = 0; ii < n_workers; ++ii)
    {
        if (pthread_create(&pool->threads[ii], NULL, worker_main, pool) != 0)
        {
            break;
        }
        pool->n_workers = ii + 1;
    }
    return true;
}

static void free_workers(ThreadPool *pool)
{
    if (pool->threads)
    {
        stop_workers(pool);
        pthread_mutex_destroy(&pool->mutex);
        pthread_cond_destroy(&pool->wake);
        pthread_cond_destroy(&pool->done);
    }
    PS_FREE(pool->threads);
}

/* Runs a pass on the workers and the calling thread. */
static void run_parallel(ThreadPool *pool)
{
    pthread_mutex_lock(&pool->mutex);
    pool->n_finished = 0;
    pool->generation++;
    pthread_cond_broadcast(&pool->wake);
    pthread_mutex_unlock(&pool->mutex);

    run_chunks(pool, 0);

    pthread_mutex_lock(&pool->mutex);
    while (pool->n_finished < pool->n_workers)
    {
        pthread_cond_wait(&pool->done, &pool->mutex);
    }
    pthread_mutex_unlock(&pool->mutex);
}

#endif

ThreadPool *thread_pool_new(int n_workers)
{
    ThreadPool *pool = (ThreadPool *) ps_calloc(1, sizeof(ThreadPool));
    if (!pool)
    {
        return NULL;
    }
#if THREAD_POOL_HAS_WORKERS
    if (n_workers > 0 && !start_workers(pool, n_workers))
    {
        thread_pool_free(pool);
        return NULL;
    }
#else
    (void) n_workers;
#endif
    return pool;
}

void thread_pool_free(ThreadPool *pool)
{
    if (!pool)
    {
        return;
    }
#if THREAD_POOL_HAS_WORKERS
    free_workers(pool);
#endif
    PS_FREE(pool);
}

int thread_pool_n_threads(const ThreadPool *pool)
{
    return pool->n_workers + 1;
}

int thread_pool_online_cpus(void)
{
#if THREAD_POOL_HAS_WORKERS
    long online = sysconf(_SC_NPROCESSORS_ONLN);
    return online < 1 ? 1 : (int) online;
#else
    return 1;
#endif
}

int thread_pool_resolve_n_threads(int requested)
{
    int n_threads = requested;
    if (n_threads <= 0)
    {
        int online = thread_pool_online_cpus();
        n_threads =
            online < THREAD_POOL_AUTO_THREADS ? online : THREAD_POOL_AUTO_THREADS;
    }
    if (n_threads > THREAD_POOL_MAX_THREADS)
    {
        n_threads = THREAD_POOL_MAX_THREADS;
    }
    return n_threads < 1 ? 1 : n_threads;
}

void thread_pool_run(ThreadPool *pool, ThreadPoolChunkFunction function,
                     void *context, int n_chunks)
{
    if (n_chunks <= 0)
    {
        return;
    }
    pool->function = function;
    pool->context = context;
    pool->n_chunks = n_chunks;
    pool->next_chunk = 0;
#if THREAD_POOL_HAS_WORKERS
    if (pool->n_workers > 0 && n_chunks > 1)
    {
        run_parallel(pool);
        return;
    }
#endif
    run_chunks(pool, 0);
}

int thread_pool_chunk_count(int n_threads, int n_items, int items_per_chunk)
{
    int by_threads = 8 * n_threads;
    int by_items = n_items / items_per_chunk + 1;
    return by_threads < by_items ? by_threads : by_items;
}

void thread_pool_chunk_bounds(int n_items, int n_chunks, int *bounds)
{
    for (int ii = 0; ii <= n_chunks; ++ii)
    {
        bounds[ii] = (int) ((int64_t) n_items * ii / n_chunks);
    }
}

void thread_pool_chunk_bounds_by_prefix(const int *prefix, int n_items, int n_chunks,
                                        int *bounds)
{
    bounds[0] = 0;
    for (int ii = 1; ii < n_chunks; ++ii)
    {
        int goal = (int) ((int64_t) prefix[n_items] * ii / n_chunks);
        int lo = 0, hi = n_items;
        while (lo < hi)
        {
            int mid = lo + (hi - lo) / 2;
            if (prefix[mid] < goal)
            {
                lo = mid + 1;
            }
            else
            {
                hi = mid;
            }
        }
        bounds[ii] = lo;
    }
    bounds[n_chunks] = n_items;
}

int thread_pool_bucket_offsets(int *counts, int n_chunks, int n_buckets,
                               int *offsets)
{
    int position = 0;
    for (int bucket = 0; bucket < n_buckets; ++bucket)
    {
        if (offsets)
        {
            offsets[bucket] = position;
        }
        for (int chunk = 0; chunk < n_chunks; ++chunk)
        {
            int *slot =
                counts + (size_t) chunk * (size_t) n_buckets + (size_t) bucket;
            int count = *slot;
            *slot = position;
            position += count;
        }
    }
    if (offsets)
    {
        offsets[n_buckets] = position;
    }
    return position;
}
