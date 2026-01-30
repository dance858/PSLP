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

#ifndef PSLP_STATS_H
#define PSLP_STATS_H

typedef struct PresolveStats
{
    size_t n_rows_original;
    size_t n_cols_original;
    size_t nnz_original;
    size_t n_rows_reduced;
    size_t n_cols_reduced;
    size_t nnz_reduced;

    /* reduction stats */
    size_t nnz_removed_trivial;
    size_t nnz_removed_fast;
    size_t nnz_removed_primal_propagation;
    size_t nnz_removed_parallel_rows;
    size_t nnz_removed_parallel_cols;

    /* time stats */
    double time_init;
    double time_fast_reductions;
    double time_medium_reductions;
    double time_primal_propagation;
    double time_parallel_rows;
    double time_parallel_cols;
    double time_presolve;
    double time_postsolve;

} PresolveStats;

#endif // PSLP_STATS_H
