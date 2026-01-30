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

#ifndef CORE_INTERNALDATA_H
#define CORE_INTERNALDATA_H

#include "Postsolver.h"
#include "Tags.h"
#include "iVec.h"

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

    iVec *fixed_cols_to_delete;
    iVec *sub_cols_to_delete;
    iVec *rows_to_delete;
} State;

/* Constructor and destructor */
State *new_state(int *row_sizes, int *col_sizes, struct Lock *col_locks,
                 PSLP_uint n_rows, PSLP_uint n_cols, struct Activity *activities,
                 struct Work *work, const RowTag *row_tags);
void free_state(State *data);

/* Maintains consistency of all internal data structures when columns and
   removes are deleted */
void clean_state(State *data, const struct Mapping *maps, PSLP_uint n_rows,
                 PSLP_uint n_cols);

#endif // CORE_INTERNALDATA_H
