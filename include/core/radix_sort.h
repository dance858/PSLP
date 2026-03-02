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

#include <stddef.h>

// LSD radix sort on composite key (sparsity_ID, coeff_hash).
// Sorts row indices in-place. aux must have space for n ints.
void radix_sort_rows(int *rows, size_t n, const int *sparsity_IDs,
                     const int *coeff_hashes, int *aux);

// Parallel version: splits into 4 chunks, sorts in parallel, merges.
// Falls back to sequential radix_sort_rows for n < 100000.
void parallel_radix_sort_rows(int *rows, size_t n, const int *sparsity_IDs,
                              const int *coeff_hashes, int *aux);

// Single-key LSD radix sort. Sorts indices by keys[indices[i]]
// ascending (unsigned order). Stable. aux must have space for n ints.
void radix_sort_by_key(int *indices, size_t n, const int *keys, int *aux);

#endif /* RADIX_SORT_H */
