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
    int n_rows_original;
    int n_cols_original;
    int nnz_original;
    int n_rows_reduced;
    int n_cols_reduced;
    int nnz_reduced;

    /* reduction stats */
    int nnz_removed_trivial;
    int nnz_removed_fast;
    int nnz_removed_primal_propagation;
    int nnz_removed_parallel_rows;
    int nnz_removed_parallel_cols;

    /* time stats */
    double ps_time_init;
    double ps_time_fast;
    double ps_time_medium;
    double ps_time_primal_propagation;
    double ps_time_parallel_rows;
    double ps_time_parallel_cols;
    double presolve_total_time;
    double ps_time_post_solve;

} PresolveStats;

#endif // PSLP_STATS_H
