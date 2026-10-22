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

#ifndef RADIX_SORT_H
#define RADIX_SORT_H

#include "thread_pool.h"
#include <stddef.h>

// The parallel sort splits inputs of at least this many rows. Measured on
// mittelmann, miplib, gams and qap instances (2026-10-18): below a million
// lines the chunks gain nothing, above they save about a tenth of the sort.
#define PARALLEL_SORT_THRESHOLD 1000000
// The most chunks of a chunked sort: one per thread of the largest pool.
#define RADIX_SORT_MAX_CHUNKS THREAD_POOL_MAX_THREADS

// LSD radix sort on composite key (sparsity_ID, coeff_hash).
// Sorts row indices in-place. aux must have space for n ints.
void radix_sort_rows(int *rows, size_t n, const int *sparsity_IDs,
                     const int *coeff_hashes, int *aux);

// Sorts the chunks on the pool and merges them. One chunk, or a NULL pool, is
// the sequential sort. The result does not depend on the number of chunks.
void radix_sort_rows_chunked(int *rows, size_t n, const int *sparsity_IDs,
                             const int *coeff_hashes, int *aux,
                             struct ThreadPool *pool, int n_chunks);

// Parallel version: from PARALLEL_SORT_THRESHOLD rows on, one chunk per thread
// of the pool, up to a cap.
void parallel_radix_sort_rows(int *rows, size_t n, const int *sparsity_IDs,
                              const int *coeff_hashes, int *aux,
                              struct ThreadPool *pool);

// Single-key LSD radix sort. Sorts indices by keys[indices[i]]
// ascending (unsigned order). Stable. aux must have space for n ints.
void radix_sort_by_key(int *indices, size_t n, const int *keys, int *aux);

#endif /* RADIX_SORT_H */
