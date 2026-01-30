/*
 * Copyright 2025 Daniel Cederberg
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

#include "utils.h"
#include "Numerics.h"

void dPtr_shrink(double *ptr, const int *map, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        if (map[i] != -1)
        {
            ptr[map[i]] = ptr[i];
        }
    }
}

void iPtr_shrink(int *ptr, const int *map, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        if (map[i] != -1)
        {
            ptr[map[i]] = ptr[i];
        }
    }
}

void rowTagPtr_shrink(RowTag *ptr, const int *map, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        if (map[i] != -1)
        {
            ptr[map[i]] = ptr[i];
        }
    }
}

void colTagPtr_shrink(ColTag *ptr, const int *map, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        if (map[i] != -1)
        {
            ptr[map[i]] = ptr[i];
        }
    }
}

void shrink_idx_vector(iVec *vec, const int *map)
{
    size_t len = vec->len;
    size_t curr = 0;

    for (size_t i = 0; i < len; ++i)
    {
        if (map[vec->data[i]] != -1)
        {
            vec->data[curr] = map[vec->data[i]];
            curr++;
        }
    }

    vec->len = curr;
}

double get_max_abs(const double *vals, size_t len)
{
    double max_abs = 0.0;
    for (size_t i = 0; i < len; ++i)
    {
        double abs_val = ABS(vals[i]);
        max_abs = MAX(max_abs, abs_val);
    }
    return max_abs;
}
