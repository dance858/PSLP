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

#ifndef CORE_WORKSPACE_H
#define CORE_WORKSPACE_H

#include "Bounds.h"
#include "iVec.h"
#include <stdbool.h>

typedef struct Mapping
{
    int *rows;
    int *cols;
} Mapping;

struct DtonWorkspace;
struct FoldMaps;
struct Solution;

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

    // radix_aux is used as auxiliary buffer for radix sort in parallel_rows/cols
    // TODO: do we really need this? or can we reuse space?
    int *radix_aux;

    Mapping *mappings;

    // the fold, kept for postsolve; NULL when the problem was not folded
    struct FoldMaps *fold_maps;
    struct Solution *folded_sol; // scratch of the folded dimensions

    // allocated when the doubleton explorer is on
    struct DtonWorkspace *dton;
} Work;

Work *new_work(size_t n_rows, size_t n_cols, bool dton_eq);
void free_work(Work *work);

/* Frees the arrays that only presolve uses. The mappings stay. */
void free_work_scratch(Work *work);

/* A solution of the given dimensions (one spare entry so an empty problem
   still gets an allocation), or NULL when an allocation fails. */
struct Solution *new_solution(size_t n_rows, size_t n_cols);
void free_solution(struct Solution *sol);

#endif // CORE_WORKSPACE_H
