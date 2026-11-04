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

#ifndef CORE_POSTSOLVER_CONES_H
#define CORE_POSTSOLVER_CONES_H

/* Postsolve records of the conic reductions. Their types are in the
   ReductionTypes enum of Postsolver.h and they are replayed by postsolver_run
   with the other records. */

#include "PSLP_API.h"
#include "Postsolver.h"

struct Constraints;

/* This function saves the removal of the singleton column x_j (cost c_j) by
   the rotated-to-P reduction of Cones_soc.c. Its entry a x_j is the head row h
   of a rotated block whose other head row w is the constant omega. The tail is
   the n_tails active rows of the block after its two heads (rows already
   inactive are skipped), each with one entry f_t x_k and side r_t. Postsolve
   recomputes x_j and the multipliers of the recorded rows. Everything but c_j
   is read from the constraints, which must still hold the block.

   * info->vals stores    [a, r_h, c_j, omega, f_0, r_0, f_1, r_1, ...].
   * info->indices stores [j, w, h, n_tails, t_0, k_0, t_1, k_1, ...] (t = tail
     row, k = its column).
*/
void save_retrieval_rotated_to_p(PostsolveInfo *info,
                                 const struct Constraints *constraints, double c_j,
                                 int j, int w, int h, const ConeBlock *block,
                                 int n_tails);

typedef struct RotatedToPRecord
{
    int j;
    int w;
    int h;
    int n_tails;
    double a;
    double r_h;
    double c_j;
    double omega;
    const int *tails;        /* (row, col) pairs */
    const double *tail_vals; /* (f, r) pairs */
} RotatedToPRecord;

static inline RotatedToPRecord
decode_rotated_to_p(const int *indices, const double *vals, int record_len)
{
    RotatedToPRecord r;
    r.j = indices[0];
    r.w = indices[1];
    r.h = indices[2];
    r.n_tails = indices[3];
    assert(record_len == 4 + 2 * r.n_tails);
    (void) record_len; // only used by the assert
    r.a = vals[0];
    r.r_h = vals[1];
    r.c_j = vals[2];
    r.omega = vals[3];
    r.tails = indices + 4;
    r.tail_vals = vals + 4;
    return r;
}

/* retrieves x_j, z_j and the multipliers of the block */
void retrieve_rotated_to_p(Solution *sol, const RotatedToPRecord *r);

#endif // CORE_POSTSOLVER_CONES_H
