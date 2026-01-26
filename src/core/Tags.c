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

#include "Tags.h"
#include "Memory_wrapper.h"
#include "glbopts.h"

RowTag *new_rowtags(double *lhs, double *rhs, int n_rows)
{
    RowTag *row_tags = (RowTag *) ps_calloc((size_t) n_rows, sizeof(RowTag));
    RETURN_PTR_IF_NULL(row_tags, NULL);

    for (int i = 0; i < n_rows; ++i)
    {
        if (IS_POS_INF(rhs[i]))
        {
            rhs[i] = INF;
            UPDATE_TAG(row_tags[i], R_TAG_RHS_INF);
        }

        if (IS_NEG_INF(lhs[i]))
        {
            lhs[i] = -INF;
            UPDATE_TAG(row_tags[i], R_TAG_LHS_INF);
        }
        else if (lhs[i] == rhs[i])
        {
            UPDATE_TAG(row_tags[i], R_TAG_EQ);
        }
    }

    return row_tags;
}
