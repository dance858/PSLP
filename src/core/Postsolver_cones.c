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

#include "Postsolver_cones.h"
#include "Constraints.h"
#include "Matrix.h"
#include "dVec.h"
#include "iVec.h"
#include "u16Vec.h"
#include <assert.h>

void save_retrieval_rotated_to_p(PostsolveInfo *info, const Constraints *constraints,
                                 double c_j, int j, int w, int h,
                                 const ConeBlock *block, int n_tails)
{
    const Matrix *A = constraints->A;
    const double *rhs = constraints->rhs;
    double a = A->x[A->p[h].start];
    double r_h = rhs[h];
    double omega = -rhs[w];

    u16Vec_append(info->type, ROTATED_TO_P);
    iVec_append(info->indices, j);
    iVec_append(info->indices, w);
    iVec_append(info->indices, h);
    iVec_append(info->indices, n_tails);
    dVec_append(info->vals, a);
    dVec_append(info->vals, r_h);
    dVec_append(info->vals, c_j);
    dVec_append(info->vals, omega);
    int n_saved = 0;
    for (int t = block->start + 2; t < block->start + block->dim; ++t)
    {
        if (HAS_TAG(constraints->row_tags[t], R_TAG_INACTIVE))
        {
            continue;
        }
        int p = A->p[t].start;
        iVec_append(info->indices, t);
        iVec_append(info->indices, A->i[p]);
        dVec_append(info->vals, A->x[p]);
        dVec_append(info->vals, rhs[t]);
        n_saved++;
    }
    assert(n_saved == n_tails);
    (void) n_saved; // only used by the assert
    iVec_append(info->starts, (int) info->indices->len);
    assert(info->starts->len == info->type->len + 1);
    assert(info->indices->len == info->vals->len);
}

void retrieve_rotated_to_p(Solution *sol, const RotatedToPRecord *r)
{
    double kappa = r->c_j / (r->a * r->omega);
    double s_h = 0.0;
    for (int t = 0; t < r->n_tails; ++t)
    {
        int k = r->tails[2 * t + 1];
        assert(sol->x[k] != COL_NOT_RETRIEVED);
        double s_t = r->tail_vals[2 * t] * sol->x[k] - r->tail_vals[2 * t + 1];
        s_h += s_t * s_t;
        sol->y[r->tails[2 * t]] = -kappa * s_t;
    }
    s_h /= 2.0 * r->omega;
    sol->x[r->j] = (r->r_h + s_h) / r->a;
    sol->z[r->j] = 0.0;
    sol->y[r->h] = r->c_j / r->a;
    sol->y[r->w] = kappa * s_h;
}
