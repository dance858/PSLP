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

#include "SparseAccumulator.h"
#include <string.h>

void sparse_accumulator_init(SparseAccumulator *acc, int *stamp, double *value,
                             uint8_t *flags, int *touched, size_t n)
{
    acc->stamp = stamp;
    acc->value = value;
    acc->flags = flags;
    acc->touched = touched;
    acc->n_touched = 0;
    acc->version = 0;
    memset(stamp, 0, n * sizeof(int));
}

void sparse_accumulator_clear(SparseAccumulator *acc)
{
    acc->version++;
    acc->n_touched = 0;
}

void sparse_accumulator_add(SparseAccumulator *acc, int col, double val,
                            uint8_t flag)
{
    if (acc->stamp[col] != acc->version)
    {
        acc->stamp[col] = acc->version;
        acc->value[col] = 0.0;
        acc->flags[col] = 0;
        acc->touched[acc->n_touched++] = col;
    }
    acc->value[col] += val;
    acc->flags[col] |= flag;
}

void sparse_accumulator_sort(SparseAccumulator *acc)
{
    int *touched = acc->touched;
    for (int a = 1; a < acc->n_touched; ++a)
    {
        int key = touched[a];
        int b = a - 1;
        while (b >= 0 && touched[b] > key)
        {
            touched[b + 1] = touched[b];
            --b;
        }
        touched[b + 1] = key;
    }
}
