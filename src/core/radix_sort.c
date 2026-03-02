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
#include "pslp_thread.h"
#include <stdint.h>
#include <string.h>

#define NUM_SORT_THREADS 4
#define PARALLEL_SORT_THRESHOLD 100000

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

static void insertion_sort_by_key(int *indices, size_t n,
                                  const int *keys)
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

void radix_sort_by_key(int *indices, size_t n,
                       const int *keys, int *aux)
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
            unsigned byte =
                ((uint32_t) keys[src[i]] >> shift) & 0xFF;
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
            unsigned byte =
                ((uint32_t) keys[src[i]] >> shift) & 0xFF;
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

static void *sort_chunk_func(void *arg)
{
    SortChunkArg *a = (SortChunkArg *) arg;
    radix_sort_rows(a->rows, a->n, a->sparsity_IDs, a->coeff_hashes, a->aux);
    return NULL;
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

// 4-way merge of sorted chunks into dst.
static void merge_4way(int *dst, int *chunk_ptrs[NUM_SORT_THREADS],
                       size_t chunk_sizes[NUM_SORT_THREADS], const int *sparsity_IDs,
                       const int *coeff_hashes, size_t total_n)
{
    size_t pos[NUM_SORT_THREADS] = {0, 0, 0, 0};

    for (size_t out = 0; out < total_n; out++)
    {
        int best = -1;
        for (int k = 0; k < NUM_SORT_THREADS; k++)
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

void parallel_radix_sort_rows(int *rows, size_t n, const int *sparsity_IDs,
                              const int *coeff_hashes, int *aux)
{
    if (n < PARALLEL_SORT_THRESHOLD)
    {
        radix_sort_rows(rows, n, sparsity_IDs, coeff_hashes, aux);
        return;
    }

    // Divide into NUM_SORT_THREADS chunks.
    // First (n % NUM_SORT_THREADS) chunks get one extra element.
    size_t base = n / NUM_SORT_THREADS;
    size_t extra = n % NUM_SORT_THREADS;

    SortChunkArg args[NUM_SORT_THREADS];
    size_t chunk_sizes[NUM_SORT_THREADS];
    int *chunk_ptrs[NUM_SORT_THREADS];
    size_t offset = 0;

    for (int k = 0; k < NUM_SORT_THREADS; k++)
    {
        chunk_sizes[k] = base + (k < (int) extra ? 1 : 0);
        chunk_ptrs[k] = rows + offset;
        args[k].rows = rows + offset;
        args[k].n = chunk_sizes[k];
        args[k].sparsity_IDs = sparsity_IDs;
        args[k].coeff_hashes = coeff_hashes;
        args[k].aux = aux + offset;
        offset += chunk_sizes[k];
    }

    // Launch worker threads for chunks 1..3, main thread sorts chunk 0.
    ps_thread_t threads[NUM_SORT_THREADS - 1];
    for (int k = 1; k < NUM_SORT_THREADS; k++)
    {
        ps_thread_create(&threads[k - 1], NULL, sort_chunk_func, &args[k]);
    }
    sort_chunk_func(&args[0]);

    for (int k = 1; k < NUM_SORT_THREADS; k++)
    {
        ps_thread_join(&threads[k - 1], NULL);
    }

    // 4-way merge into aux, then copy back.
    merge_4way(aux, chunk_ptrs, chunk_sizes, sparsity_IDs, coeff_hashes, n);
    memcpy(rows, aux, n * sizeof(int));
}
