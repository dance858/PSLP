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

#include "State.h"
#include "Activity.h"
#include "Locks.h"
#include "Numerics.h"
#include "PSLP_API.h"

#include "Workspace.h"
#include "utils.h"

#define MIN_INIT_STON_ROWS 100
#define INIT_STON_ROWS_FRACTION 0.01
#define MIN_INIT_STON_COLS 100
#define INIT_STON_COLS_FRACTION 0.01
#define MIN_INIT_EMPTY_ROWS 50
#define INIT_EMPTY_ROWS_FRACTION 0.001
#define MIN_INIT_EMPTY_COLS 50
#define INIT_EMPTY_COLS_FRACTION 0.001
#define MIN_INIT_FIXED_COLS_TO_DELETE 1000
#define INIT_FIXED_COLS_TO_DELETE_FRACTION 0.01
#define MIN_INIT_SUB_COLS_TO_DELETE 1000
#define INIT_SUB_COLS_TO_DELETE_FRACTION 0.01
#define MIN_INIT_ROWS_TO_DELETE 1000
#define INIT_ROWS_TO_DELETE_FRACTION 0.01
#define MIN_INIT_DTON_ROWS 100
#define INIT_DTON_ROWS_FRACTION 0.01
#define MIN_INIT_UPDATED_ACTIVITIES 1000
#define INIT_UPDATED_ACTIVITIES_FRACTION 0.01

State *new_state(int *row_sizes, int *col_sizes, Lock *col_locks, int n_rows,
                 int n_cols, Activity *activities, Work *work,
                 const RowTag *row_tags)
{
    State *data = (State *) ps_malloc(1, sizeof(State));
    RETURN_PTR_IF_NULL(data, NULL);
    data->col_locks = col_locks;
    data->row_sizes = row_sizes;
    data->col_sizes = col_sizes;
    data->activities = activities;
    data->work = work;

    // --------------------------------------------------------------------------
    //                    allocate a bunch of vectors
    // --------------------------------------------------------------------------
    data->ston_rows = iVec_new(
        MAX(MIN_INIT_STON_ROWS, (size_t) (n_rows * INIT_STON_ROWS_FRACTION)));
    data->ston_cols = iVec_new(
        MAX(MIN_INIT_STON_COLS, (size_t) (n_cols * INIT_STON_COLS_FRACTION)));
    data->dton_rows = iVec_new(
        MAX(MIN_INIT_DTON_ROWS, (size_t) (n_rows * INIT_DTON_ROWS_FRACTION)));
    data->empty_cols = iVec_new(
        MAX(MIN_INIT_EMPTY_COLS, (size_t) (n_cols * INIT_EMPTY_COLS_FRACTION)));
    data->empty_rows = iVec_new(
        MAX(MIN_INIT_EMPTY_ROWS, (size_t) (n_rows * INIT_EMPTY_ROWS_FRACTION)));
    data->updated_activities =
        iVec_new(MAX(MIN_INIT_UPDATED_ACTIVITIES,
                     (size_t) (n_rows * INIT_UPDATED_ACTIVITIES_FRACTION)));
    data->fixed_cols_to_delete =
        iVec_new(MAX(MIN_INIT_FIXED_COLS_TO_DELETE,
                     (size_t) (n_cols * INIT_FIXED_COLS_TO_DELETE_FRACTION)));
    data->sub_cols_to_delete =
        iVec_new(MAX(MIN_INIT_SUB_COLS_TO_DELETE,
                     (size_t) (n_cols * INIT_SUB_COLS_TO_DELETE_FRACTION)));
    data->rows_to_delete = iVec_new(MAX(
        MIN_INIT_ROWS_TO_DELETE, (size_t) (n_rows * INIT_ROWS_TO_DELETE_FRACTION)));
    data->postsolve_info = postsolve_info_new(n_rows, n_cols);

    if (!data->ston_rows || !data->ston_cols || !data->dton_rows ||
        !data->empty_cols || !data->empty_rows || !data->fixed_cols_to_delete ||
        !data->sub_cols_to_delete || !data->rows_to_delete || !data->postsolve_info)
    {
        free_state(data);
        return NULL;
    }

    // --------------------------------------------------------------------------
    //           record empty, singleton, and doubleton equality rows
    // --------------------------------------------------------------------------
    for (int i = 0; i < n_rows; ++i)
    {
        switch (row_sizes[i])
        {
            case 0:
                assert(!iVec_contains(data->empty_rows, i));
                iVec_append(data->empty_rows, i);
                break;
            case 1:
                assert(!iVec_contains(data->ston_rows, i));
                iVec_append(data->ston_rows, i);
                break;
            case 2:
                if (HAS_TAG(row_tags[i], R_TAG_EQ))
                {
                    assert(!iVec_contains(data->ston_rows, i));
                    iVec_append(data->dton_rows, i);
                }
                break;
        }
    }

    // --------------------------------------------------------------------------
    //                    record empty and singleton columns
    // --------------------------------------------------------------------------
    for (int i = 0; i < n_cols; ++i)
    {
        switch (col_sizes[i])
        {
            case 0:
                iVec_append(data->empty_cols, i);
                break;
            case 1:
                iVec_append(data->ston_cols, i);
                break;
        }
    }

    // --------------------------------------------------------------------------
    //                    record activities to update
    // --------------------------------------------------------------------------
    for (int i = 0; i < n_rows; ++i)
    {
        if (activities[i].n_inf_min == 0 || activities[i].n_inf_max == 0)
        {
            activities[i].status = ADDED;
            iVec_append(data->updated_activities, i);
        }
    }

    return data;
}

void free_state(State *data)
{
    if (!data)
    {
        return;
    }

    iVec_free(data->ston_rows);
    iVec_free(data->ston_cols);
    iVec_free(data->dton_rows);
    iVec_free(data->empty_cols);
    iVec_free(data->empty_rows);
    iVec_free(data->updated_activities);
    iVec_free(data->fixed_cols_to_delete);
    iVec_free(data->sub_cols_to_delete);
    iVec_free(data->rows_to_delete);
    postsolve_info_free(data->postsolve_info);
    PS_FREE(data);
}

static void shrink_locks(Lock *ptr, const int *map, int len)
{
    for (int i = 0; i < len; ++i)
    {
        if (map[i] != -1)
        {
            ptr[map[i]] = ptr[i];
        }
    }
}

static void shrink_activities(Activity *ptr, const int *map, int len)
{
    for (int i = 0; i < len; ++i)
    {
        if (map[i] != -1)
        {
            ptr[map[i]] = ptr[i];
        }
    }
}

void clean_state(State *data, const Mapping *maps, int n_rows, int n_cols)
{
    iPtr_shrink(data->row_sizes, maps->rows, n_rows);
    shrink_activities(data->activities, maps->rows, n_rows);
    shrink_locks(data->col_locks, maps->cols, n_cols);
    iPtr_shrink(data->col_sizes, maps->cols, n_cols);

    shrink_idx_vector(data->ston_rows, maps->rows);
    shrink_idx_vector(data->dton_rows, maps->rows);
    shrink_idx_vector(data->empty_rows, maps->rows);
    shrink_idx_vector(data->updated_activities, maps->rows);
    shrink_idx_vector(data->ston_cols, maps->cols);
    shrink_idx_vector(data->empty_cols, maps->cols);
}
