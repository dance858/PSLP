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

#ifndef CORE_INTERNALDATA_H
#define CORE_INTERNALDATA_H

#include "Postsolver.h"
#include "Tags.h"
#include "iVec.h"
#include <stdbool.h>

// forward declarations
struct Mapping;
struct Activity;
struct Work;
struct Lock;

typedef struct State
{
    iVec *ston_rows;
    iVec *ston_cols;
    iVec *dton_rows;
    iVec *empty_cols;
    iVec *empty_rows;
    iVec *updated_activities;
    struct Lock *col_locks;
    struct Activity *activities;
    struct Work *work;
    PostsolveInfo *postsolve_info;

    // number of non-zeros within each row and column
    int *row_sizes;
    int *col_sizes;

    /* row sizes of P, SIZE_INACTIVE_COL for removed variables, NULL for an LP */
    int *P_row_sizes;

    iVec *fixed_cols_to_delete;
    iVec *sub_cols_to_delete;
    iVec *rows_to_delete;
} State;

/* Constructor and destructor */
State *new_state(int *row_sizes, int *col_sizes, int *P_row_sizes,
                 struct Lock *col_locks, size_t n_rows, size_t n_cols,
                 struct Activity *activities, struct Work *work,
                 const RowTag *row_tags);
void free_state(State *data);

/* True if column k has a non-empty row in P. Singleton and parallel columns
   skip such a variable; the doubleton eliminator and the dual fix handle it
   through the P-aware paths (see Problem.h). For an inactive column the row
   size is SIZE_INACTIVE_COL and this returns true; callers test
   C_TAG_INACTIVE first. Always false for an LP. */
static inline bool col_in_P(const State *data, int k)
{
    return data->P_row_sizes != NULL && data->P_row_sizes[k] != 0;
}

/* True if empty_cols holds a column with an empty P row, which
   remove_empty_cols always removes. Columns with entries in P stay in the
   worklist; a diagonal-only one is also removable but not detected here. */
static inline bool has_removable_empty_cols(const State *data)
{
    for (size_t i = 0; i < data->empty_cols->len; ++i)
    {
        if (!col_in_P(data, data->empty_cols->data[i]))
        {
            return true;
        }
    }

    return false;
}

/* Maintains consistency of all internal data structures when columns and
   removes are deleted */
void clean_state(State *data, const struct Mapping *maps, size_t n_rows,
                 size_t n_cols);

#endif // CORE_INTERNALDATA_H
