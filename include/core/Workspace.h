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

#ifndef CORE_WORKSPACE_H
#define CORE_WORKSPACE_H

#include "Bounds.h"
#include "iVec.h"

typedef struct Mapping
{
    int *rows;
    int *cols;
} Mapping;

/* The int_vec in the following struct is used within
    1. parallel_rows to store bin_starts
*/
typedef struct Work
{
    int *iwork_n_cols;

    // iwork_n_rows are used within parallel_rows
    int *iwork_n_rows;

    // used within parallel_rows and parallel_cols
    // iwork1_max_nrows_ncols is also used within dual propagation
    int *iwork1_max_nrows_ncols;
    int *iwork2_max_nrows_ncols;

    // int_vec is used within parallel_rows and dual propagation
    iVec *int_vec;

    Mapping *mappings;
} Work;

Work *new_work(int n_rows, int n_col);
void free_work(Work *work);

#endif // CORE_WORKSPACE_H