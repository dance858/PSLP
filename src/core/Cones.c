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

#include "Cones.h"

#include <limits.h>
#include <math.h>

static const ConeTypeInfo cone_types[] = {
    [PSLP_CONE_SOC] = {1, INT_MAX},
    [PSLP_CONE_ROTATED] = {2, INT_MAX},
};

const ConeTypeInfo *cone_type_info(ConeType type)
{
    if (type != PSLP_CONE_SOC && type != PSLP_CONE_ROTATED)
    {
        return NULL;
    }
    return cone_types + type;
}

bool cones_valid(const ConeBlock *cones, size_t n_cones, size_t m, const double *rhs)
{
    if (n_cones == 0)
    {
        return true;
    }
    if (cones == NULL || cones[0].start < 0)
    {
        return false;
    }

    size_t next_start = (size_t) cones[0].start;
    for (size_t ii = 0; ii < n_cones; ++ii)
    {
        const ConeBlock *block = cones + ii;
        const ConeTypeInfo *info = cone_type_info(block->type);
        if (info == NULL || block->dim < info->dim_min ||
            block->dim > info->dim_max || block->start < 0 ||
            (size_t) block->start != next_start)
        {
            return false;
        }
        next_start += (size_t) block->dim;
    }
    if (next_start != m)
    {
        return false;
    }

    /* The cone rows are the rows from the first block to the end. */
    for (size_t i = (size_t) cones[0].start; i < m; ++i)
    {
        if (!isfinite(rhs[i]))
        {
            return false;
        }
    }
    return true;
}
