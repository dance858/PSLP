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

#ifndef CORE_SPARSEACCUMULATOR_H
#define CORE_SPARSEACCUMULATOR_H

#include <stddef.h>
#include <stdint.h>

/* Sums values into the columns of one sparse vector at a time. Each column
   holds a value and flag bits that are OR-ed together. A column counts as
   touched once something is added to it in the current vector, and only
   touched columns are read. Starting a new vector clears it without writing
   to the n columns. */
typedef struct SparseAccumulator
{
    int *stamp;     /* [n] version that last touched the column */
    double *value;  /* [n] */
    uint8_t *flags; /* [n] */
    int *touched;   /* [n] touched columns in first-touch order */
    int n_touched;
    int version;
} SparseAccumulator;

/* Sets the arrays, all owned by the caller, and zeroes stamp. Call
   sparse_accumulator_clear before the first add. */
void sparse_accumulator_init(SparseAccumulator *acc, int *stamp, double *value,
                             uint8_t *flags, int *touched, size_t n);

/* Starts a new vector with no touched columns. */
void sparse_accumulator_clear(SparseAccumulator *acc);

void sparse_accumulator_add(SparseAccumulator *acc, int col, double val,
                            uint8_t flag);

/* Sorts the touched columns ascending. Insertion sort, meant for short
   vectors. */
void sparse_accumulator_sort(SparseAccumulator *acc);

#endif // CORE_SPARSEACCUMULATOR_H
