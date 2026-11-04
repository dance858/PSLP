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

#include "Bounds.h"
#include "Chains.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "DtonsEq_internal.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "Postsolver.h"
#include "Problem.h"
#include "State.h"
#include "Tags.h"

/* Decides which of the two entries of a doubleton equality row to
   substitute. Returns the index (0 or 1). We prefer the sparser column,
   unless the coefficient of the sparser column is much smaller than the
   other one */
int dton_choose_subst(const double *row_vals, int col_size0, int col_size1)
{
    double a_abs = ABS(row_vals[0]);
    double b_abs = ABS(row_vals[1]);

    /* equal sizes: the column with the larger absolute value */
    if (col_size0 == col_size1)
    {
        return (a_abs > b_abs) ? 0 : 1;
    }

    int sparse = (col_size0 < col_size1) ? 0 : 1;
    double sparse_abs = (sparse == 0) ? a_abs : b_abs;
    double dense_abs = (sparse == 0) ? b_abs : a_abs;

    /* numerical stability */
    if (sparse_abs < DTON_PIVOT_RATIO * dense_abs)
    {
        return 1 - sparse;
    }
    return sparse;
}

/* Which entry (0 or 1) of a doubleton equality row to substitute, or -1 when
   the row cannot be eliminated. A column absent from P is always preferred. A
   column in P is substituted only when the other column is in P too and its
   row of P is a positive diagonal alone that folds into the other column's
   diagonal without fill-in (sub0/sub1). */
int dton_choose_subst_qp(const double *row_vals, int col_size0, int col_size1,
                         bool in_P0, bool in_P1, bool sub0, bool sub1)
{
    if (in_P0 != in_P1)
    {
        return in_P0 ? 1 : 0;
    }
    if (!in_P0 || (sub0 && sub1))
    {
        return dton_choose_subst(row_vals, col_size0, col_size1);
    }
    if (!sub0 && !sub1)
    {
        return -1;
    }
    return sub0 ? 0 : 1;
}

bool dton_reserve_records(DtonWorkspace *dton_work, size_t n_records)
{
    if (n_records <= (size_t) dton_work->substs.cap)
    {
        return true;
    }
    bool ok = true;
    ok &= ps_grow(&dton_work->substs.recs, n_records, sizeof(DtonSubst));
    ok &= ps_grow(&dton_work->substs.order, n_records, sizeof(int));
    ok &= ps_grow(&dton_work->substs.succ, n_records, sizeof(int));
    ok &= ps_grow(&dton_work->substs.drop_priority, n_records, sizeof(int));
    ok &= ps_grow(&dton_work->substs.depth, n_records, sizeof(int));
    ok &= ps_grow(&dton_work->substs.stamp, n_records, sizeof(int));
    ok &= ps_grow(&dton_work->targets.list, n_records, sizeof(DtonTarget));
    if (!ok)
    {
        return false;
    }
    dton_work->substs.cap = (int) n_records;
    return true;
}

void dton_P_flags(const Problem *prob, int col0, int col1, bool *in_P, bool *sub,
                  double *Pkk)
{
    const State *state = prob->constraints->state;
    const int *P_row_sizes = state->P_row_sizes;
    const Matrix *P = prob->obj->P;

    in_P[0] = col_in_P(state, col0);
    in_P[1] = col_in_P(state, col1);
    sub[0] = false;
    sub[1] = false;
    Pkk[0] = 0.0;
    Pkk[1] = 0.0;
    if (in_P[0] && in_P[1])
    {
        /* only symmetry is checked on input, so a diagonal is not given */
        sub[0] = diagonal_only_P_entry(P, P_row_sizes, col0, Pkk) &&
                 P_row_has_diagonal(P, P_row_sizes, col1);
        sub[1] = diagonal_only_P_entry(P, P_row_sizes, col1, Pkk + 1) &&
                 P_row_has_diagonal(P, P_row_sizes, col0);
    }
}

int dton_push_record(DtonWorkspace *dton_work, int i, int k, int j, double aik,
                     double aij, double rhs, double Pkk)
{
    assert(dton_work->substs.n_recs < dton_work->substs.cap);
    int idx = dton_work->substs.n_recs++;
    DtonSubst *rec = dton_work->substs.recs + idx;
    rec->k = k;
    rec->owner = i;
    rec->j = j;
    rec->aik = aik;
    rec->aij = aij;
    rec->dir_mult = -aij / aik;
    rec->dir_shift = rhs / aik;
    rec->target = j;
    rec->mult = rec->dir_mult;
    rec->shift = rec->dir_shift;
    rec->rhs = rhs;
    rec->Pkk = Pkk;
    dton_work->substs.col_subst[k] = idx;
    return idx;
}

/* For every doubleton equality row, choose one column to substitute. If the chosen
   column is already claimed by another row, the row is deferred to the next round */
bool dton_claim(Problem *prob, DtonWorkspace *dton_work, int *deferred,
                int *n_deferred)
{
    Constraints *constraints = prob->constraints;
    const Matrix *A = constraints->A;
    const double *rhs = constraints->rhs;
    const State *state = constraints->state;
    const int *row_sizes = state->row_sizes;
    const int *col_sizes = state->col_sizes;
    const iVec *dton_rows = state->dton_rows;

    dton_work->substs.n_recs = 0;

    int i, ii, col0, col1, subst, stay, k;

    /* Every worklist row claims at most one record, and the changed surviving
       columns are bounded by the records. */
    if (!dton_reserve_records(dton_work, (size_t) dton_rows->len))
    {
        return false;
    }

    for (ii = 0; ii < dton_rows->len; ++ii)
    {
        i = dton_rows->data[ii];

        /* a row that used to be a doubleton might have been modified */
        if (row_sizes[i] < 2)
        {
            continue;
        }

        assert(row_sizes[i] == 2 && A->p[i].end - A->p[i].start == 2 &&
               constraints->lhs[i] == rhs[i] && !IS_ABS_INF(rhs[i]) &&
               HAS_TAG(constraints->row_tags[i], R_TAG_EQ) &&
               !HAS_TAG(constraints->row_tags[i], R_TAG_INACTIVE));

        const int *cols = A->i + A->p[i].start;
        const double *vals = A->x + A->p[i].start;
        col0 = cols[0];
        col1 = cols[1];

        assert(!HAS_TAG(constraints->col_tags[col0] | constraints->col_tags[col1],
                        C_TAG_INACTIVE));

        bool in_P[2], sub[2];
        double Pkk[2];
        dton_P_flags(prob, col0, col1, in_P, sub, Pkk);
        subst = dton_choose_subst_qp(vals, col_sizes[col0], col_sizes[col1], in_P[0],
                                     in_P[1], sub[0], sub[1]);

        /* the row cannot be eliminated, or the col this row wants to
           substitute is already claimed by another row */
        if (subst < 0 || dton_work->substs.col_subst[cols[subst]] >= 0)
        {
            deferred[(*n_deferred)++] = i;
            continue;
        }
        stay = 1 - subst;
        k = cols[subst];

        /* store information about the substitution */
        dton_push_record(dton_work, i, k, cols[stay], vals[subst], vals[stay],
                         rhs[i], Pkk[subst]);
    }
    return true;
}

/* Map substituted columns onto the round's survivors. Cycles are broken by
   ignoring the substitution done by the row with largest index. The row
   for the ignored substitution is appended to 'deferred'. */
void dton_compose(DtonWorkspace *dton_work, int *deferred, int *n_deferred)
{
    DtonSubst *recs = dton_work->substs.recs;
    int *col_subst = dton_work->substs.col_subst;
    int *succ = dton_work->substs.succ;
    int *depth = dton_work->substs.depth;
    int *order = dton_work->substs.order;
    int *drop_priority = dton_work->substs.drop_priority;
    int *stamp = dton_work->substs.stamp;
    int n_recs = dton_work->substs.n_recs;

    /* succ[i] is the index of the substitution that eliminates the stay column
       of substitution i, or -1 if no substitution eliminates the stay column */
    for (int i = 0; i < n_recs; ++i)
    {
        succ[i] = col_subst[recs[i].j];
        drop_priority[i] = recs[i].owner;
    }

    /* break cycles by ignoring the substitution done by the row with largest idx */
    int *dropped = deferred + *n_deferred;
    int n_dropped = compute_chain_depths(n_recs, succ, drop_priority, depth, order,
                                         dropped, stamp, dton_work->acc.touched);

    // -----------------------------------------------------------------------------
    //  map every eliminated col onto the surviving col at the end of its chain
    // -----------------------------------------------------------------------------
    int n_kept = n_recs - n_dropped;
    for (int ii = n_kept - 1; ii >= 0; --ii)
    {
        int idx = order[ii];
        DtonSubst *rec = recs + idx;
        if (depth[idx] == 0)
        {
            rec->target = rec->j;
            rec->mult = rec->dir_mult;
            rec->shift = rec->dir_shift;
        }
        else
        {
            const DtonSubst *child = recs + succ[idx];
            rec->target = child->target;
            rec->mult = rec->dir_mult * child->mult;
            rec->shift = rec->dir_mult * child->shift + rec->dir_shift;
        }
    }

    // -----------------------------------------------------------------------
    // reset ignored substitutions and remove the ignored records
    // -----------------------------------------------------------------------
    for (int ii = 0; ii < n_dropped; ++ii)
    {
        int idx = dropped[ii];
        col_subst[recs[idx].k] = -1;
        dropped[ii] = recs[idx].owner;
    }
    *n_deferred += n_dropped;

    int *new_index = stamp;
    int out = 0;
    for (int idx = 0; idx < n_recs; ++idx)
    {
        if (depth[idx] == CHAINS_DROPPED)
        {
            continue;
        }
        recs[out] = recs[idx];
        depth[out] = depth[idx];
        col_subst[recs[out].k] = out;
        new_index[idx] = out++;
    }
    dton_work->substs.n_recs = out;
    dton_work->substs.n_first = out;
    for (int ii = 0; ii < n_kept; ++ii)
    {
        order[ii] = new_index[order[ii]];
    }
}

/* Transfers the bounds of a record's eliminated column x_k onto its stay
   column x_j through the row as the record sees it, aij * x_j + aik * x_k =
   rhs (row owner of A for a first-pass record), given in column order for the
   postsolve record. Returns INFEASIBLE when a transferred bound contradicts the
   bounds of x_j, UNCHANGED otherwise. */
static PresolveStatus dton_transfer_bounds_rec(Constraints *constraints,
                                               const DtonSubst *rec)
{
    int i = rec->owner;
    int j = rec->j;
    int k = rec->k;
    double aij = rec->aij;
    double aik = rec->aik;
    const Bound *bound_k = constraints->bounds + k;
    ColTag col_tag_k = constraints->col_tags[k];
    int cols[2] = {k < j ? k : j, k < j ? j : k};
    double vals[2] = {k < j ? aik : aij, k < j ? aij : aik};

    assert(aij != 0 && aik != 0);
    bool same_sign = (aik * aij > 0.0);
    PresolveStatus status = UNCHANGED;

    if (!HAS_TAG(col_tag_k, C_TAG_LB_INF))
    {
        double bound = (rec->rhs - aik * bound_k->lb) / aij;
        status = same_sign ? update_ub_with_row(constraints, j, bound, i, cols, vals,
                                                2 HUGE_BOUND_IS_NOT_OK)
                           : update_lb_with_row(constraints, j, bound, i, cols, vals,
                                                2 HUGE_BOUND_IS_NOT_OK);
        RETURN_IF_INFEASIBLE(status);
    }

    if (!HAS_TAG(col_tag_k, C_TAG_UB_INF))
    {
        double bound = (rec->rhs - aik * bound_k->ub) / aij;
        status = same_sign ? update_lb_with_row(constraints, j, bound, i, cols, vals,
                                                2 HUGE_BOUND_IS_NOT_OK)
                           : update_ub_with_row(constraints, j, bound, i, cols, vals,
                                                2 HUGE_BOUND_IS_NOT_OK);
    }
    return status;
}

/* Transfers the bounds of each first-pass eliminated column onto its direct
   stay column, one link at a time in descending chain depth, so each transfer
   sees the bounds tightened by the deeper links. A composed transfer would
   break dual complementarity on chains. A record of dton_cascade transfers in
   dton_record instead, right before its own records. Returns INFEASIBLE at the
   first contradicting transfer, UNCHANGED otherwise. */
PresolveStatus dton_transfer_bounds(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const DtonSubst *recs = dton_work->substs.recs;
    const int *order = dton_work->substs.order;

    for (int ii = 0; ii < dton_work->substs.n_recs; ++ii)
    {
        assert(ii == 0 || dton_work->substs.depth[order[ii - 1]] >=
                              dton_work->substs.depth[order[ii]]);
        if (order[ii] >= dton_work->substs.n_first)
        {
            continue;
        }
        PresolveStatus status =
            dton_transfer_bounds_rec(constraints, recs + order[ii]);
        RETURN_IF_INFEASIBLE(status);
    }
    return UNCHANGED;
}

/* Postsolve records, three per eliminated column k with owner row i
   and survivor s: ADDED_ROWS(i, column k, aik), SUB_COL_DTON(k, s, mult,
   shift) and DELETED_ROW(i, c_k / aik). Replayed in reverse they give x_k and
       y_i = (c_k - sum_{r != i} a_rk y_r) / a_ik,
   where the sum includes the owner rows whose stay column is k and the owner
   rows claimed in dton_cascade that held k. Those records are deeper, so the
   columns are emitted by ascending depth and replayed deepest first. A column in P
   (row k of P is the diagonal P_kk alone) adds P_kk x_k / a_ik to y_i through a
   DELETED_ROW_ZK record, saved first and hence replayed last. */
PresolveStatus dton_record(Problem *prob, DtonWorkspace *dton_work)
{
    Constraints *constraints = prob->constraints;
    const Matrix *AT = constraints->AT; /* pre-round transpose */
    const double *c = prob->obj->c;     /* pre-round objective */
    PostsolveInfo *info = constraints->state->postsolve_info;

    /* ascending depth: order lists the substitutions by descending depth */
    for (int ii = dton_work->substs.n_recs - 1; ii >= 0; --ii)
    {
        int idx = dton_work->substs.order[ii];
        const DtonSubst *rec = dton_work->substs.recs + idx;
        int k = rec->k;
        int i = rec->owner;
        double aik = rec->aik;

        /* a cascade record transfers its bounds here, see dton_transfer_bounds */
        if (idx >= dton_work->substs.n_first)
        {
            PresolveStatus status = dton_transfer_bounds_rec(constraints, rec);
            RETURN_IF_INFEASIBLE(status);
        }

        const int *rows = AT->i + AT->p[k].start;
        const double *col_vals = AT->x + AT->p[k].start;
        int len = AT->p[k].end - AT->p[k].start;

#ifndef NDEBUG
        /* Every row of the pre-round column is active. */
        for (int jj = 0; jj < len; ++jj)
        {
            assert(!HAS_TAG(constraints->row_tags[rows[jj]], R_TAG_INACTIVE));
        }
#endif

        if (rec->Pkk != 0.0)
        {
            save_retrieval_deleted_row_zk(info, i, k, aik, rec->Pkk);
        }
        save_retrieval_added_rows(info, i, rows, col_vals, (size_t) len, aik);
        save_retrieval_sub_col_dton(info, k, rec->target, rec->mult, rec->shift);
        save_retrieval_deleted_row(info, i, c[k] / aik);
    }
    return UNCHANGED;
}
