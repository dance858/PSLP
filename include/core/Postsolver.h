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

#ifndef CORE_POSTSOLVER_H
#define CORE_POSTSOLVER_H

#include <stdbool.h>
#include <stdint.h>

#include "PSLP_sol.h"
#include "Tags.h"

struct u16Vec;
struct dVec;
struct iVec;
struct Constraints;

typedef uint16_t ReductionType;

enum ReductionTypes
{
    // required for primal postsolve
    FIXED_COL = 0,
    FIXED_COL_INF = 1 << 0,
    SUB_COL = 1 << 1,
    PARALLEL_COL = 1 << 2,

    // required for dual postsolve
    DELETED_ROW = 1 << 3,
    ADDED_ROW = 1 << 4,
    ADDED_ROWS = 1 << 5,
    LHS_CHANGE = 1 << 6,
    RHS_CHANGE = 1 << 7,
    EQ_TO_INEQ = 1 << 8,
    BOUND_CHANGE_NO_ROW = 1 << 9,
    BOUND_CHANGE_THE_ROW = 1 << 10,
};

typedef struct PostsolveInfo
{
    size_t n_cols_reduced;
    size_t n_rows_reduced;

    // contains the type of reduction
    struct u16Vec *type;

    // contains the start index of the information required to undo a
    // reduction
    struct iVec *starts;

    // indices and vals contain the information required to undo a reduction
    struct iVec *indices;
    struct dVec *vals;

    // maps from original problem to reduced problem. If col_map[4] = 2 it means
    // that the 2nd column in the reduced problem corresponds to the 4th column
    // in the original problem. If col_map[4] = -1 it means that the 4th column
    // in the original problem was removed.
    const int *col_map;
    const int *row_map;
} PostsolveInfo;

PostsolveInfo *postsolve_info_new(size_t n_rows, size_t n_cols);
void postsolve_info_free(PostsolveInfo *info);
void postsolver_update(PostsolveInfo *info, size_t n_cols_reduced,
                       size_t n_rows_reduced, const int *col_map,
                       const int *row_map);
void postsolver_run(const PostsolveInfo *info, Solution *sol, const double *x,
                    const double *y, const double *z);

/* Saves the information required to retrieve variable xk that was fixed
   to val. To recover the dual variable we need zk = ck - ak^T y
   * info->vals stores    [val, ck, ak[0], ak[1], .., dots, ak[len - 1]
   * info->indices stores [col, dummy, rows[0], rows[1], ... rows[len - 1]].
*/
void save_retrieval_fixed_col(PostsolveInfo *info, int col, double val, double ck,
                              const double *vals, const int *rows, size_t len);

/* Saves the information required to retrieve variable xk that was fixed
   to either +INF or -INF.
    * info->indices stores [sign(xk), k, row_len1, cols1, row_len2,
                            cols2, ..., row_len_nrows, cols_nrows]
    where sign(xk) = 1 if xk is fixed to +INF and -1 if xk is fixed to -INF.

    * info->vals stores [nrows, bound, rhs / lhs 1, coeffs1, rhs / lhs 2,
                         coeffs2, ..., rhs / lhs nrows, coeffs_nrows]
*/
void save_retrieval_fixed_col_inf(PostsolveInfo *info, int col, int pos_inf,
                                  const struct Constraints *constraints,
                                  double bound);

/* This function saves the information required to retrieve variable xk
   that was substituted from the problem using equality constraint i:
   aik xk + sum_{j != k} aij xj = rhs.
   info->vals stores    [rhs, vals[0], vals[1], ... vals[len - 1], ck].
   info->indices stores [k  , cols[0], cols[1], ... cols[len - 1], i].
 */
void save_retrieval_sub_col(PostsolveInfo *info, int col, int *cols, double *coeffs,
                            size_t len, double rhs, int i, double ck);

/* This function saves the information required to retrieve variable xj
   and xk that were replaced with a new variable x_new = xj + ratio * xk
   due to parallel column reduction.
    info->vals stores    [lb_j, ub_j, lb_k, ub_k, ratio].
    info->indices stores [j, k, cTag_j, cTag_k, dummy_value].
*/
void save_retrieval_parallel_col(PostsolveInfo *info, double ub_j, double lb_j,
                                 double lb_k, double ub_k, double ratio, int j,
                                 int k, ColTag cTag_j, ColTag cTag_k);

/* This function saves the value on yi when row i is deleted */
void save_retrieval_deleted_row(PostsolveInfo *info, int row, double val);

/* This function saves the information required to retrieve yi when
   row i is added to row j so the new row j becomes aj = aj + ratio * ai
*/
void save_retrieval_added_row(PostsolveInfo *info, int i, int j, double ratio);

/* This function saves the information required to retrieve yi and yj
   when the rhs or lhs of row i is changed because of row j. Here we assume that
   ai = ratio * aj.

   * info->vals stores [new_side, ratio, vals[0], vals[1], ... vals[len - 1]].
   * info->indices stores [i, j, cols[0], cols[1], ... cols[len - 1]].

   where vals and cols correspond to row i.
*/
void save_retrieval_rhs_or_lhs_change(PostsolveInfo *info, int i, const double *vals,
                                      const int *cols, size_t len, double new_side,
                                      int j, double ratio, bool is_lhs_change);

/* This function saves the information required to retrieve yi when
   row i has been added to many other rows. In this case we have
   yi = yi - \sum_{j \neq i} (ajk / aik) yj.

   * info->vals stores [aik, vals[0], vals[1], ... vals[len - 1]].
   * info->indices stores [i, rows[0], rows[1], ... rows[len - 1]].
*/
void save_retrieval_added_rows(PostsolveInfo *info, int i, const int *rows,
                               const double *vals, size_t len, double aik);

/* This function saves the information required to undo the effect of a bound
   change. Suppose we use row 'i' to update one bound on variable 'j'.
   If the implied bound on variable 'j' is active at the optimal solution
   of the reduced problem, we set yi = yi + zj / aij. For every variable k
   appearing in row 'i', we update zk = zk - (aik / aij) * zj.

   * info->vals stores [implied_bound, original_bound].
   * info->indices stores [j, is_original_lower_bound]
*/
void save_retrieval_bound_change_no_row(PostsolveInfo *info, int j,
                                        double implied_bound,
                                        double original_other_bound,
                                        int is_original_other_bound_lower_bound);

/* This function saves the actual row that was used to derive several bound changes.
   info->vals stores [(double) num_of_bound_changes), vals]
   info->indices stores [i, cols]
*/
void save_retrieval_bound_change_the_row(PostsolveInfo *info, int i, const int *cols,
                                         const double *vals, size_t len,
                                         int num_of_bound_changes);

/* This function saves the information required to retrieve yi when
   equality row i has been transformed into an inequality by eliminating
   column k. The postsolve is
   yi = (ck / aik) + yi.

   * info->vals stores [ck / aik].
   * info->indices stores [i].
*/
void save_retrieval_eq_to_ineq(PostsolveInfo *info, int row, double val);

#endif // CORE_POSTSOLVER_H
