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

#include "radix_sort.h"
#include "thread_pool.h"
#include <stdint.h>
#include <string.h>

/* Chunks of the parallel sort: the threads of the pool, at most this many. */
#define RADIX_SORT_CHUNK_CAP 4

static void insertion_sort_rows(int *rows, size_t n, const int *sparsity_IDs,
                                const int *coeff_hashes)
{
    for (size_t i = 1; i < n; i++)
    {
        int key = rows[i];
        int key_sp = sparsity_IDs[key];
        int key_ch = coeff_hashes[key];
        size_t j = i;
        while (j > 0)
        {
            int prev_sp = sparsity_IDs[rows[j - 1]];
            int prev_ch = coeff_hashes[rows[j - 1]];
            if (prev_sp < key_sp || (prev_sp == key_sp && prev_ch < key_ch))
            {
                break;
            }
            rows[j] = rows[j - 1];
            j--;
        }
        rows[j] = key;
    }
}

void radix_sort_rows(int *rows, size_t n, const int *sparsity_IDs,
                     const int *coeff_hashes, int *aux)
{
    if (n < 256)
    {
        insertion_sort_rows(rows, n, sparsity_IDs, coeff_hashes);
        return;
    }

    size_t counts[256];
    int *src = rows, *dst = aux;

    // Phase 1: 4 passes on coeff_hashes (secondary key, sorted first in LSD)
    // Phase 2: 4 passes on sparsity_IDs (primary key, sorted last in LSD)
    for (int phase = 0; phase < 2; phase++)
    {
        const int *keys = (phase == 0) ? coeff_hashes : sparsity_IDs;

        for (int pass = 0; pass < 4; pass++)
        {
            int shift = pass * 8;

            // histogram
            memset(counts, 0, 256 * sizeof(size_t));
            for (size_t i = 0; i < n; i++)
            {
                unsigned byte = ((uint32_t) keys[src[i]] >> shift) & 0xFF;
                counts[byte]++;
            }

            // skip pass if all values fall in one bucket (never skip pass 0)
            if (!(phase == 0 && pass == 0))
            {
                int skip = 0;
                for (int b = 0; b < 256; b++)
                {
                    if (counts[b] == n)
                    {
                        skip = 1;
                        break;
                    }
                }
                if (skip)
                {
                    continue;
                }
            }

            // prefix sum
            size_t total = 0;
            for (int b = 0; b < 256; b++)
            {
                size_t c = counts[b];
                counts[b] = total;
                total += c;
            }

            // scatter: backward for pass 0 (reverse stability), forward otherwise
            if (phase == 0 && pass == 0)
            {
                // backward scatter reverses equal-byte elements
                for (size_t i = n; i > 0; i--)
                {
                    unsigned byte = ((uint32_t) keys[src[i - 1]] >> shift) & 0xFF;
                    dst[counts[byte]++] = src[i - 1];
                }
            }
            else
            {
                for (size_t i = 0; i < n; i++)
                {
                    unsigned byte = ((uint32_t) keys[src[i]] >> shift) & 0xFF;
                    dst[counts[byte]++] = src[i];
                }
            }

            // swap src and dst
            int *tmp = src;
            src = dst;
            dst = tmp;
        }
    }

    // if result ended up in aux, copy back to rows
    if (src != rows)
    {
        memcpy(rows, src, n * sizeof(int));
    }
}

// --- Single-key radix sort ---

static void insertion_sort_by_key(int *indices, size_t n, const int *keys)
{
    for (size_t i = 1; i < n; i++)
    {
        int idx = indices[i];
        int idx_key = keys[idx];
        size_t j = i;
        while (j > 0 && keys[indices[j - 1]] > idx_key)
        {
            indices[j] = indices[j - 1];
            j--;
        }
        indices[j] = idx;
    }
}

void radix_sort_by_key(int *indices, size_t n, const int *keys, int *aux)
{
    if (n < 256)
    {
        insertion_sort_by_key(indices, n, keys);
        return;
    }

    size_t counts[256];
    int *src = indices, *dst = aux;

    for (int pass = 0; pass < 4; pass++)
    {
        int shift = pass * 8;

        // histogram
        memset(counts, 0, 256 * sizeof(size_t));
        for (size_t i = 0; i < n; i++)
        {
            unsigned byte = ((uint32_t) keys[src[i]] >> shift) & 0xFF;
            counts[byte]++;
        }

        // skip pass if all values fall in one bucket
        if (pass > 0)
        {
            int skip = 0;
            for (int b = 0; b < 256; b++)
            {
                if (counts[b] == n)
                {
                    skip = 1;
                    break;
                }
            }
            if (skip)
            {
                continue;
            }
        }

        // prefix sum
        size_t total = 0;
        for (int b = 0; b < 256; b++)
        {
            size_t c = counts[b];
            counts[b] = total;
            total += c;
        }

        // scatter (forward — stable)
        for (size_t i = 0; i < n; i++)
        {
            unsigned byte = ((uint32_t) keys[src[i]] >> shift) & 0xFF;
            dst[counts[byte]++] = src[i];
        }

        // swap src and dst
        int *tmp = src;
        src = dst;
        dst = tmp;
    }

    // if result ended up in aux, copy back to indices
    if (src != indices)
    {
        memcpy(indices, src, n * sizeof(int));
    }
}

// --- Parallel radix sort infrastructure ---

typedef struct
{
    int *rows;
    size_t n;
    const int *sparsity_IDs;
    const int *coeff_hashes;
    int *aux;
} SortChunkArg;

// Pool callback: sorts one chunk. The thread index is not needed because each
// chunk has its own slice of rows and aux.
static void sort_chunk(void *context, int chunk, int thread)
{
    (void) thread;
    SortChunkArg *arg = (SortChunkArg *) context + chunk;
    radix_sort_rows(arg->rows, arg->n, arg->sparsity_IDs, arg->coeff_hashes,
                    arg->aux);
}

// Compare two row indices by (sparsity_ID, coeff_hash) as uint32_t.
// Returns <0 if a comes first, >0 if b comes first, 0 if equal.
static inline int merge_compare(int a, int b, const int *sparsity_IDs,
                                const int *coeff_hashes)
{
    uint32_t sa = (uint32_t) sparsity_IDs[a];
    uint32_t sb = (uint32_t) sparsity_IDs[b];
    if (sa != sb) return (sa < sb) ? -1 : 1;
    uint32_t ca = (uint32_t) coeff_hashes[a];
    uint32_t cb = (uint32_t) coeff_hashes[b];
    if (ca != cb) return (ca < cb) ? -1 : 1;
    return 0;
}

// Merges the sorted chunks into dst. The sequential sort leaves equal keys in
// reverse input order, so the last chunk wins ties: the result does not depend
// on the number of chunks.
static void merge_chunks(int *dst, int *const *chunk_ptrs, const size_t *chunk_sizes,
                         int n_chunks, const int *sparsity_IDs,
                         const int *coeff_hashes, size_t total_n)
{
    size_t pos[RADIX_SORT_MAX_CHUNKS] = {0};

    for (size_t out = 0; out < total_n; out++)
    {
        int best = -1;
        for (int k = 0; k < n_chunks; k++)
        {
            if (pos[k] >= chunk_sizes[k]) continue;
            if (best < 0 ||
                merge_compare(chunk_ptrs[k][pos[k]], chunk_ptrs[best][pos[best]],
                              sparsity_IDs, coeff_hashes) <= 0)
            {
                best = k;
            }
        }
        dst[out] = chunk_ptrs[best][pos[best]++];
    }
}

void radix_sort_rows_chunked(int *rows, size_t n, const int *sparsity_IDs,
                             const int *coeff_hashes, int *aux, ThreadPool *pool,
                             int n_chunks)
{
    if (n_chunks > RADIX_SORT_MAX_CHUNKS)
    {
        n_chunks = RADIX_SORT_MAX_CHUNKS;
    }
    if (n_chunks > 1 && (size_t) n_chunks > n)
    {
        n_chunks = (int) n;
    }
    if (!pool || n_chunks <= 1)
    {
        radix_sort_rows(rows, n, sparsity_IDs, coeff_hashes, aux);
        return;
    }

    // Divide into n_chunks chunks. The first (n % n_chunks) get one extra element.
    size_t base = n / (size_t) n_chunks;
    size_t extra = n % (size_t) n_chunks;

    SortChunkArg args[RADIX_SORT_MAX_CHUNKS];
    size_t chunk_sizes[RADIX_SORT_MAX_CHUNKS];
    int *chunk_ptrs[RADIX_SORT_MAX_CHUNKS];
    size_t offset = 0;

    for (int k = 0; k < n_chunks; k++)
    {
        chunk_sizes[k] = base + ((size_t) k < extra ? 1 : 0);
        chunk_ptrs[k] = rows + offset;
        args[k].rows = rows + offset;
        args[k].n = chunk_sizes[k];
        args[k].sparsity_IDs = sparsity_IDs;
        args[k].coeff_hashes = coeff_hashes;
        args[k].aux = aux + offset;
        offset += chunk_sizes[k];
    }

    // Sort the chunks on the pool; the calling thread takes chunks too.
    thread_pool_run(pool, sort_chunk, args, n_chunks);

    // Merge into aux, then copy back.
    merge_chunks(aux, chunk_ptrs, chunk_sizes, n_chunks, sparsity_IDs, coeff_hashes,
                 n);
    memcpy(rows, aux, n * sizeof(int));
}

void parallel_radix_sort_rows(int *rows, size_t n, const int *sparsity_IDs,
                              const int *coeff_hashes, int *aux, ThreadPool *pool)
{
    if (!pool || n < PARALLEL_SORT_THRESHOLD)
    {
        radix_sort_rows(rows, n, sparsity_IDs, coeff_hashes, aux);
        return;
    }
    int n_chunks = thread_pool_n_threads(pool);
    if (n_chunks > RADIX_SORT_CHUNK_CAP)
    {
        n_chunks = RADIX_SORT_CHUNK_CAP;
    }
    radix_sort_rows_chunked(rows, n, sparsity_IDs, coeff_hashes, aux, pool,
                            n_chunks);
}
