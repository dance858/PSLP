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

#include "Problem.h"
#include "Binary_search.h"
#include "Constraints.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "RowColViews.h"
#include "State.h"
#include "Workspace.h"
#include "glbopts.h"
#include "utils.h"

Problem *new_problem(Constraints *constraints, Objective *obj)
{
    Problem *problem = (Problem *) ps_malloc(1, sizeof(Problem));
    RETURN_PTR_IF_NULL(problem, NULL);

    problem->constraints = constraints;
    problem->obj = obj;

    return problem;
}

void free_problem(Problem *problem)
{
    RETURN_IF_NULL(problem);
    RETURN_IF_NULL(problem->constraints);
    State *data = problem->constraints->state;

    if (data)
    {
        free_work(data->work);
        PS_FREE(data->col_locks);
        PS_FREE(data->activities);
        PS_FREE(data->row_sizes);
        PS_FREE(data->col_sizes);
        PS_FREE(data->P_row_sizes);
        free_state(data);
    }

    free_constraints(problem->constraints);
    free_objective(problem->obj);
    PS_FREE(problem);
}

Objective *objective_new(double *c, Matrix *P)
{
    Objective *obj = (Objective *) ps_malloc(1, sizeof(Objective));
    RETURN_PTR_IF_NULL(obj, NULL);
    obj->c = c;
    obj->offset = 0.0;
    obj->P = P;

    return obj;
}

void free_objective(Objective *obj)
{
    RETURN_IF_NULL(obj);
    PS_FREE(obj->c);
    free_matrix(obj->P);
    PS_FREE(obj);
}

static void objective_clean(Objective *obj, const int *P_row_sizes,
                            const int *col_map, size_t n_cols_old, size_t n_cols_new)
{
    dPtr_shrink(obj->c, col_map, n_cols_old);

    if (obj->P == NULL)
    {
        return;
    }

    remove_extra_space(obj->P, P_row_sizes, col_map, n_cols_new);
}

void mark_P_row_inactive(Problem *prob, int k)
{
    if (prob->obj->P == NULL)
    {
        return;
    }

    int *P_row_sizes = prob->constraints->state->P_row_sizes;
    assert(P_row_sizes[k] == 0);
    assert(prob->obj->P->p[k].start == prob->obj->P->p[k].end);
    P_row_sizes[k] = SIZE_INACTIVE_COL;
}

bool diagonal_only_P_entry(const Matrix *P, const int *P_row_sizes, int k,
                           double *Pkk)
{
    if (P_row_sizes[k] != 1)
    {
        return false;
    }

    int start = P->p[k].start;
    if (P->i[start] != k || P->x[start] <= 0.0)
    {
        return false;
    }

    *Pkk = P->x[start];
    return true;
}

bool P_row_has_diagonal(const Matrix *P, const int *P_row_sizes, int k)
{
    return sorted_find(P->i + P->p[k].start, P_row_sizes[k], k) >= 0;
}

/* c + upd, snapped to zero when the sum is rounding noise relative to its
   terms. An exact cancellation leaves a residue whose sign is noise, and
   downstream sign tests (e.g. the unboundedness check in process_colston_ineq)
   must never see it. */
static double add_cancelling(double c, double upd)
{
    double sum = c + upd;
    return (ABS(sum) <= CANCEL_TOL_REL * MAX(ABS(c), ABS(upd))) ? 0.0 : sum;
}

void fix_var_in_obj(Problem *prob, int col, double value)
{
    Objective *obj = prob->obj;
    Matrix *P = obj->P;
    obj->offset += obj->c[col] * value;
    if (P == NULL)
    {
        return;
    }

    int *P_row_sizes = prob->constraints->state->P_row_sizes;
    int start = P->p[col].start;
    int end = P->p[col].end;
    assert(P_row_sizes[col] == end - start);

    for (int j = start; j < end; ++j)
    {
        int l = P->i[j];
        double Pkl = P->x[j];

        if (l == col)
        {
            obj->offset += 0.5 * Pkl * value * value;
            continue;
        }

        obj->c[l] = add_cancelling(obj->c[l], Pkl * value);

        // P is symmetric, so P_lk exists in row l
        RowView row_l = new_rowview(P->x + P->p[l].start, P->i + P->p[l].start,
                                    P_row_sizes + l, P->p + l, NULL, NULL, NULL, l);
        remove_coeff(&row_l, col);
        P->nnz -= 1;
    }

    P->nnz -= (size_t) (end - start);
    P->p[col].end = start;
    P_row_sizes[col] = SIZE_INACTIVE_COL;
}

void problem_clean(Problem *prob)
{
    Mapping *maps = prob->constraints->state->work->mappings;
    size_t n_cols_old = prob->constraints->n;
    size_t n_rows_old = prob->constraints->m;
    constraints_clean(prob->constraints, maps);
    objective_clean(prob->obj, prob->constraints->state->P_row_sizes, maps->cols,
                    n_cols_old, prob->constraints->n);
    clean_state(prob->constraints->state, maps, n_rows_old, n_cols_old);
}

void sub_var_in_obj(Problem *prob, const double *vals, const int *cols, int len,
                    int k, double aik, double rhs)
{
    Objective *obj = prob->obj;
    mark_P_row_inactive(prob, k);

    if (obj->c[k] == 0.0)
    {
        return;
    }

    double ratio = obj->c[k] / aik;
    for (int i = 0; i < len; ++i)
    {
        obj->c[cols[i]] = add_cancelling(obj->c[cols[i]], -ratio * vals[i]);
    }

    obj->offset += rhs * ratio;
}

void sub_var_in_obj_dton(Objective *obj, int subst, int stay, double mult,
                         double shift)
{
    obj->c[stay] += mult * obj->c[subst];
    obj->offset += shift * obj->c[subst];
    obj->c[subst] = 0.0;
}

void sub_diag_var_in_obj_dton(Problem *prob, int k, int t, double mult, double shift)
{
    Objective *obj = prob->obj;
    Matrix *P = obj->P;
    int *P_row_sizes = prob->constraints->state->P_row_sizes;
    assert(P_row_sizes[k] == 1 && P->i[P->p[k].start] == k);
    double Pkk = P->x[P->p[k].start];

    /* (1/2) P_kk x_k^2 + c_k x_k with x_k = shift + mult x_t */
    obj->c[t] += mult * (obj->c[k] + Pkk * shift);
    obj->offset += shift * (obj->c[k] + 0.5 * Pkk * shift);
    obj->c[k] = 0.0;

    int start_t = P->p[t].start;
    int rel = sorted_find(P->i + start_t, P_row_sizes[t], t);
    assert(rel >= 0);
    P->x[start_t + rel] += mult * mult * Pkk;

    /* column k of P is the diagonal alone, so no other row changes */
    P->p[k].end = P->p[k].start;
    P->nnz -= 1;
    P_row_sizes[k] = 0;
}
