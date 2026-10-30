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
#include "Cones_soc.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "Problem.h"
#include "glbopts.h"

#include <limits.h>
#include <math.h>

static const ConeTypeInfo cone_types[] = {
    [PSLP_CONE_SOC] = {1, INT_MAX},
    [PSLP_CONE_ROTATED] = {2, INT_MAX},
    [PSLP_CONE_EXP] = {3, 3},
};

const ConeTypeInfo *cone_type_info(ConeType type)
{
    if (type != PSLP_CONE_SOC && type != PSLP_CONE_ROTATED && type != PSLP_CONE_EXP)
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

PresolveStatus process_cones(Problem *prob, const Settings *stgs,
                             const PresolveStats *stats)
{
    Constraints *constraints = prob->constraints;
    DEBUG(verify_problem_up_to_date(constraints));

    bool rotated_to_P = stgs->rotated_to_P && stats->Pnnz_original == 0;

    for (size_t ii = 0; ii < constraints->n_cones; ++ii)
    {
        const ConeBlock *block = constraints->cones + ii;

        /* A block keeps its start and dim until constraints_clean. Its heads
           only go with the whole block, so an inactive head row is a block
           deleted in an earlier round. */
        if (HAS_TAG(constraints->row_tags[block->start], R_TAG_INACTIVE))
        {
            continue;
        }

        /* EXP blocks are passed through unchanged */
        if (block->type == PSLP_CONE_SOC || block->type == PSLP_CONE_ROTATED)
        {
            RETURN_IF_INFEASIBLE(process_soc_block(prob, block, rotated_to_P));
        }
    }

    delete_inactive_rows(constraints);
    return UNCHANGED;
}
