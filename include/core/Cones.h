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

#ifndef CORE_CONES_H
#define CORE_CONES_H

#include "PSLP_API.h"
#include "PSLP_stats.h"

#include <stdbool.h>
#include <stddef.h>

/* cone-type-specific information */
typedef struct
{
    int dim_min;
    int dim_max; /* INT_MAX when the dimension is not limited */
} ConeTypeInfo;

const ConeTypeInfo *cone_type_info(ConeType type);

/* validate the cones of the problem (returns true if everything is consistent) */
bool cones_valid(const ConeBlock *cones, size_t n_cones, size_t m,
                 const double *rhs);

struct Problem;

/* Tries the cone-specific reductions on every block, then deletes the rows
   they marked. Returns INFEASIBLE or UNCHANGED. */
PresolveStatus process_cones(struct Problem *prob, const Settings *stgs,
                             const PresolveStats *stats);

#endif /* CORE_CONES_H */
