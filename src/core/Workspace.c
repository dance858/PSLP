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

#include "Workspace.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include <stdio.h>

#define INT_VEC_INITIALIZATION 25

Work *new_work(size_t n_rows, size_t n_cols)
{
    Work *work = (Work *) ps_malloc(1, sizeof(Work));
    RETURN_PTR_IF_NULL(work, NULL);

    work->iwork_n_cols = (int *) ps_calloc(n_cols, sizeof(int));
    work->iwork_n_rows = (int *) ps_calloc(n_rows, sizeof(int));
    work->iwork1_max_nrows_ncols =
        (int *) ps_calloc(MAX(n_rows, n_cols), sizeof(int));
    work->iwork2_max_nrows_ncols =
        (int *) ps_calloc(MAX(n_rows, n_cols), sizeof(int));
    work->int_vec = iVec_new(INT_VEC_INITIALIZATION);
    work->mappings = (Mapping *) ps_malloc(1, sizeof(Mapping));
    work->mappings->cols = (int *) ps_malloc(n_cols, sizeof(int));
    work->mappings->rows = (int *) ps_malloc(n_rows, sizeof(int));

    if (!work->iwork_n_cols || !work->iwork_n_rows ||
        !work->iwork1_max_nrows_ncols || !work->iwork2_max_nrows_ncols ||
        !work->int_vec || !work->mappings || !work->mappings->cols ||
        !work->mappings->rows)
    {
        free_work(work);
        return NULL;
    }

    return work;
}

void free_work(Work *work)
{
    if (work == NULL)
    {
        return;
    }

    PS_FREE(work->iwork_n_cols);
    PS_FREE(work->iwork_n_rows);
    PS_FREE(work->iwork1_max_nrows_ncols);
    PS_FREE(work->iwork2_max_nrows_ncols);
    iVec_free(work->int_vec);
    PS_FREE(work->mappings->cols);
    PS_FREE(work->mappings->rows);
    PS_FREE(work->mappings);
    PS_FREE(work);
}
