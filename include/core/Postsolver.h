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

#ifndef CORE_POSTSOLVER_H
#define CORE_POSTSOLVER_H

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>

#include "PSLP_sol.h"
#include "Tags.h"
#include "glbopts.h"

#define COL_NOT_RETRIEVED INF
#define ROW_NOT_RETRIEVED INF
#define DUMMY_VALUE -382749

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

    // only required for mapping a solution to the reduced problem
    PARALLEL_ROW = 1 << 11,
    SIDE_RELAXED = 1 << 12,
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
void postsolver_run_primal_infeas_ray(const PostsolveInfo *info, Solution *sol,
                                      const double *y, const double *z);
void postsolver_run_dual_infeas_ray(const PostsolveInfo *info, Solution *sol,
                                    const double *x);

/* Maps a primal-dual point (x, y) of the original problem to the reduced
   problem by replaying the recorded reductions in forward order (the inverse of
   'postsolver_run', which replays them backwards). 'col_map' and 'row_map' are
   the original-to-reduced index maps, and 'x_work' / 'y_work' are scratch
   buffers of the original dimensions that must not alias x / y. Passing
   x_red = NULL skips the primal part and y_red = NULL skips the dual part. The
   reduced dual slack z is not produced here; it should be computed from the
   reduced problem data as z_red = c_red - A_red^T y_red. */
void postsolver_map_to_reduced(const PostsolveInfo *info, const int *col_map,
                               const int *row_map, size_t n_cols_orig,
                               size_t n_rows_orig, const double *x, const double *y,
                               double *x_work, double *y_work, double *x_red,
                               double *y_red);

void retrieve_deleted_row(Solution *sol, int row, double val);
void retrieve_added_row(Solution *sol, int i, int j, double ratio);
void retrieve_added_rows(Solution *sol, int i, const int *rows, const double *vals,
                         int len, double aik);

/* -------------------------------------------------------------------------
   Records
   -------------------------------------------------------------------------
   Reduction i is stored as the slices indices[starts[i] .. starts[i + 1]) and
   vals[starts[i] .. starts[i + 1]) of the two parallel arrays, so both slices
   have the same length, called the record length below. Each record type has
   a save_retrieval_* function that writes it and a decode_* function that
   reads it back into a *Record struct. The decoders take the two slices (i.e.
   indices + starts[i], vals + starts[i]) and the record length, and are the
   only place that knows the layouts; consumers must not index into a record
   themselves. Placeholder entries hold DUMMY_VALUE and are checked. The
   pointers in a decoded record point into the record and stay valid as long
   as the PostsolveInfo is not modified.
   ------------------------------------------------------------------------- */

/* Saves the information required to retrieve variable xk that was fixed
   to val. To recover the dual variable we need zk = ck - ak^T y
   * info->vals stores    [val, ck, ak[0], ak[1], .., dots, ak[len - 1]
   * info->indices stores [col, dummy, rows[0], rows[1], ... rows[len - 1]].
   len is 0 for an empty column.
*/
void save_retrieval_fixed_col(PostsolveInfo *info, int col, double val, double ck,
                              const double *vals, const int *rows, size_t len);

typedef struct FixedColRecord
{
    int col;
    double val;
    double ck;
    int len;
    const int *rows;
    const double *vals;
} FixedColRecord;

static inline FixedColRecord decode_fixed_col(const int *indices, const double *vals,
                                              int record_len)
{
    FixedColRecord r;
    assert(record_len >= 2);
    assert(indices[1] == DUMMY_VALUE);
    r.col = indices[0];
    r.val = vals[0];
    r.ck = vals[1];
    r.len = record_len - 2;
    r.rows = indices + 2;
    r.vals = vals + 2;
    return r;
}

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

/* The rows are stored one after the other, starting at row_indices /
   row_vals: row r occupies row_indices[c] = row_len, row_indices[c + 1 ..
   c + row_len] = cols and row_vals[c] = side, row_vals[c + 1 .. c + row_len]
   = coeffs, after which the next row starts at c + row_len + 1. */
typedef struct FixedColInfRecord
{
    bool pos_inf;
    int col;
    int n_rows;
    double bound;
    const int *row_indices;
    const double *row_vals;
} FixedColInfRecord;

static inline FixedColInfRecord
decode_fixed_col_inf(const int *indices, const double *vals, int record_len)
{
    FixedColInfRecord r;
    assert(record_len >= 2);
    assert(indices[0] == 1 || indices[0] == -1);
    r.pos_inf = (indices[0] > 0);
    r.col = indices[1];
    r.n_rows = (int) vals[0];
    r.bound = vals[1];
    r.row_indices = indices + 2;
    r.row_vals = vals + 2;

#ifndef NDEBUG
    // the rows must fill the record exactly
    int counter = 2;
    for (int i = 0; i < r.n_rows; ++i)
    {
        counter += indices[counter] + 1;
    }
    assert(counter == record_len);
#else
    (void) record_len;
#endif

    return r;
}

/* This function saves the information required to retrieve variable xk
   that was substituted from the problem using equality constraint i:
   aik xk + sum_{j != k} aij xj = rhs.
   info->vals stores    [rhs, vals[0], vals[1], ... vals[len - 1], ck].
   info->indices stores [k  , cols[0], cols[1], ... cols[len - 1], i].
   (vals, cols) is the whole row i including column k, so len >= 2.
 */
void save_retrieval_sub_col(PostsolveInfo *info, int col, const int *cols,
                            const double *coeffs, size_t len, double rhs, int i,
                            double ck);

typedef struct SubColRecord
{
    int k;
    double rhs;
    int len;
    const int *cols;
    const double *vals;
    int row;
    double ck;
} SubColRecord;

static inline SubColRecord decode_sub_col(const int *indices, const double *vals,
                                          int record_len)
{
    SubColRecord r;
    assert(record_len >= 4);
    r.k = indices[0];
    r.rhs = vals[0];
    r.len = record_len - 2;
    r.cols = indices + 1;
    r.vals = vals + 1;
    r.row = indices[1 + r.len];
    r.ck = vals[1 + r.len];
    return r;
}

/* This function saves the information required to retrieve variable xj
   and xk that were replaced with a new variable x_new = xj + ratio * xk
   due to parallel column reduction. (Note that the parameter order is
   ub_j, lb_j while the record stores lb_j first.)
    info->vals stores    [lb_j, ub_j, lb_k, ub_k, ratio].
    info->indices stores [j, k, cTag_j, cTag_k, dummy_value].
*/
void save_retrieval_parallel_col(PostsolveInfo *info, double ub_j, double lb_j,
                                 double lb_k, double ub_k, double ratio, int j,
                                 int k, ColTag cTag_j, ColTag cTag_k);

typedef struct ParallelColRecord
{
    int j;
    int k;
    ColTag cTag_j;
    ColTag cTag_k;
    double lb_j;
    double ub_j;
    double lb_k;
    double ub_k;
    double ratio;
} ParallelColRecord;

static inline ParallelColRecord
decode_parallel_col(const int *indices, const double *vals, int record_len)
{
    ParallelColRecord r;
    assert(record_len == 5);
    (void) record_len; // only used by the assert
    assert(indices[4] == DUMMY_VALUE);
    r.j = indices[0];
    r.k = indices[1];
    r.cTag_j = (ColTag) indices[2];
    r.cTag_k = (ColTag) indices[3];
    r.lb_j = vals[0];
    r.ub_j = vals[1];
    r.lb_k = vals[2];
    r.ub_k = vals[3];
    r.ratio = vals[4];
    return r;
}

/* This function saves the value on yi when row i is deleted.
 * info->vals stores    [val].
 * info->indices stores [row].
 */
void save_retrieval_deleted_row(PostsolveInfo *info, int row, double val);

typedef struct DeletedRowRecord
{
    int row;
    double val;
} DeletedRowRecord;

static inline DeletedRowRecord decode_deleted_row(const int *indices,
                                                  const double *vals, int record_len)
{
    DeletedRowRecord r;
    assert(record_len == 1);
    (void) record_len; // only used by the assert
    r.row = indices[0];
    r.val = vals[0];
    return r;
}

/* This function saves the information required to retrieve yi when
   row i is added to row j so the new row j becomes aj = aj + ratio * ai.
   * info->vals stores    [ratio, dummy].
   * info->indices stores [i, j].
*/
void save_retrieval_added_row(PostsolveInfo *info, int i, int j, double ratio);

typedef struct AddedRowRecord
{
    int i;
    int j;
    double ratio;
} AddedRowRecord;

static inline AddedRowRecord decode_added_row(const int *indices, const double *vals,
                                              int record_len)
{
    AddedRowRecord r;
    assert(record_len == 2);
    (void) record_len; // only used by the assert
    assert(vals[1] == DUMMY_VALUE);
    r.i = indices[0];
    r.j = indices[1];
    r.ratio = vals[0];
    return r;
}

/* This function saves the information required to retrieve yi and yj
   when the rhs or lhs of row i is changed because of row j. Here we assume that
   ai = ratio * aj.

   * info->vals stores [new_side, ratio, vals[0], vals[1], ... vals[len - 1]].
   * info->indices stores [i, j, cols[0], cols[1], ... cols[len - 1]].

   where vals and cols correspond to row i (len >= 2).
*/
void save_retrieval_rhs_or_lhs_change(PostsolveInfo *info, int i, const double *vals,
                                      const int *cols, size_t len, double new_side,
                                      int j, double ratio, bool is_lhs_change);

typedef struct SideChangeRecord
{
    int i;
    int j;
    double new_side;
    double ratio;
    int len;
    const int *cols;
    const double *vals;
} SideChangeRecord;

static inline SideChangeRecord decode_side_change(const int *indices,
                                                  const double *vals, int record_len)
{
    SideChangeRecord r;
    assert(record_len >= 4);
    r.i = indices[0];
    r.j = indices[1];
    r.new_side = vals[0];
    r.ratio = vals[1];
    r.len = record_len - 2;
    r.cols = indices + 2;
    r.vals = vals + 2;
    return r;
}

/* This function saves the information required to retrieve yi when
   row i has been added to many other rows. In this case we have
   yi = yi - \sum_{j \neq i} (ajk / aik) yj.

   * info->vals stores [aik, vals[0], vals[1], ... vals[len - 1]].
   * info->indices stores [i, rows[0], rows[1], ... rows[len - 1]].
   (len >= 1)
*/
void save_retrieval_added_rows(PostsolveInfo *info, int i, const int *rows,
                               const double *vals, size_t len, double aik);

typedef struct AddedRowsRecord
{
    int i;
    double aik;
    int len;
    const int *rows;
    const double *vals;
} AddedRowsRecord;

static inline AddedRowsRecord decode_added_rows(const int *indices,
                                                const double *vals, int record_len)
{
    AddedRowsRecord r;
    assert(record_len >= 2);
    r.i = indices[0];
    r.aik = vals[0];
    r.len = record_len - 1;
    r.rows = indices + 1;
    r.vals = vals + 1;
    return r;
}

/* This function saves the information required to undo the effect of a bound
   change. Suppose we use row 'i' to update one bound on variable 'j'.
   If the implied bound on variable 'j' is active at the optimal solution
   of the reduced problem, we set yi = yi + zj / aij. For every variable k
   appearing in row 'i', we update zk = zk - (aik / aij) * zj.

   * info->vals stores [implied_bound, original_other_bound].
   * info->indices stores [j, is_original_other_bound_lower_bound]
*/
void save_retrieval_bound_change_no_row(PostsolveInfo *info, int j,
                                        double implied_bound,
                                        double original_other_bound,
                                        int is_original_other_bound_lower_bound);

typedef struct BoundChangeNoRowRecord
{
    int j;
    int is_original_other_bound_lower_bound;
    double implied_bound;
    double original_other_bound;
} BoundChangeNoRowRecord;

static inline BoundChangeNoRowRecord
decode_bound_change_no_row(const int *indices, const double *vals, int record_len)
{
    BoundChangeNoRowRecord r;
    assert(record_len == 2);
    (void) record_len; // only used by the assert
    r.j = indices[0];
    r.is_original_other_bound_lower_bound = indices[1];
    r.implied_bound = vals[0];
    r.original_other_bound = vals[1];
    return r;
}

/* This function saves the actual row that was used to derive several bound
   changes (the BOUND_CHANGE_NO_ROW records that precede this one).
   info->vals stores    [(double) num_of_bound_changes, vals[0], .., vals[len - 1]]
   info->indices stores [i, cols[0], .., cols[len - 1]]
   (len >= 1)
*/
void save_retrieval_bound_change_the_row(PostsolveInfo *info, int i, const int *cols,
                                         const double *vals, size_t len,
                                         int num_of_bound_changes);

typedef struct BoundChangeTheRowRecord
{
    int i;
    int num_of_bound_changes;
    int len;
    const int *cols;
    const double *vals;
} BoundChangeTheRowRecord;

static inline BoundChangeTheRowRecord
decode_bound_change_the_row(const int *indices, const double *vals, int record_len)
{
    BoundChangeTheRowRecord r;
    assert(record_len >= 2);
    r.i = indices[0];
    r.num_of_bound_changes = (int) vals[0];
    r.len = record_len - 1;
    r.cols = indices + 1;
    r.vals = vals + 1;
    return r;
}

/* This function saves that parallel row j (with aj = ai / ratio) was removed
   in favour of row i. Postsolve does not need this (yj is recovered through
   LHS_CHANGE / RHS_CHANGE or is zero), but the forward map to the reduced
   problem uses it to transfer the multiplier of row j to row i.

   * info->vals stores [ratio, dummy].
   * info->indices stores [i, j].
*/
void save_retrieval_parallel_row(PostsolveInfo *info, int i, int j, double ratio);

typedef struct ParallelRowRecord
{
    int i;
    int j;
    double ratio;
} ParallelRowRecord;

static inline ParallelRowRecord
decode_parallel_row(const int *indices, const double *vals, int record_len)
{
    ParallelRowRecord r;
    assert(record_len == 2);
    (void) record_len; // only used by the assert
    assert(vals[1] == DUMMY_VALUE);
    r.i = indices[0];
    r.j = indices[1];
    r.ratio = vals[0];
    return r;
}

/* This function saves the information required to retrieve yi when
   equality row i has been transformed into an inequality by eliminating
   column k. The postsolve is
   yi = (ck / aik) + yi.

   'sign' is the sign that yi may have after the reduction: -1 if the row kept
   its rhs (yi <= 0), +1 if it kept its lhs (yi >= 0). Postsolve does not need
   it; the forward map to the reduced problem projects yi onto it.

   * info->vals stores [ck / aik, dummy].
   * info->indices stores [i, sign].
*/
void save_retrieval_eq_to_ineq(PostsolveInfo *info, int row, double val, int sign);

typedef struct EqToIneqRecord
{
    int row;
    int sign;
    double val;
} EqToIneqRecord;

static inline EqToIneqRecord decode_eq_to_ineq(const int *indices,
                                               const double *vals, int record_len)
{
    EqToIneqRecord r;
    assert(record_len == 2);
    (void) record_len; // only used by the assert
    assert(vals[1] == DUMMY_VALUE);
    assert(indices[1] == 1 || indices[1] == -1);
    r.row = indices[0];
    r.sign = indices[1];
    r.val = vals[0];
    return r;
}

/* This function saves that a finite side of row i was dropped as redundant
   (check_activities in SimpleReductions.c). Postsolve does not need it; the
   forward map to the reduced problem projects yi onto the sign the surviving
   side admits: -1 if the row kept its rhs (yi <= 0), +1 if it kept its lhs
   (yi >= 0).

   * info->vals stores [dummy, dummy].
   * info->indices stores [i, sign].
*/
void save_retrieval_side_relaxed(PostsolveInfo *info, int row, int sign);

typedef struct SideRelaxedRecord
{
    int row;
    int sign;
} SideRelaxedRecord;

static inline SideRelaxedRecord
decode_side_relaxed(const int *indices, const double *vals, int record_len)
{
    SideRelaxedRecord r;
    assert(record_len == 2);
    (void) record_len; // only used by the assert
    assert(vals[0] == DUMMY_VALUE && vals[1] == DUMMY_VALUE);
    (void) vals;
    assert(indices[1] == 1 || indices[1] == -1);
    r.row = indices[0];
    r.sign = indices[1];
    return r;
}

#endif // CORE_POSTSOLVER_H
