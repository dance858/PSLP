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

#ifndef CORE_PROBLEM_H
#define CORE_PROBLEM_H

#include <stdbool.h>
#include <stddef.h>
struct Constraints;
struct Matrix;

/* Objective (1/2) x'Px + c'x + offset. P is the full symmetric matrix in CSR
   without extra row space. For an LP, P is NULL. */
typedef struct
{
    double *c;
    double offset;
    struct Matrix *P;
} Objective;

typedef struct Problem
{
    struct Constraints *constraints;
    Objective *obj;
} Problem;

/* Constructor and destructor */
Problem *new_problem(struct Constraints *constraints, Objective *obj);
void free_problem(Problem *problem);

/* Cleans the problem in the sense that all data associated with
   inactive rows and inactive columns is removed. Rows and columns
   are also re-indexed to take the removal of rows/columns into
   account. */
void problem_clean(Problem *problem);

/* Constructor and destructor. The objective takes ownership of c and P. */
Objective *objective_new(double *c, struct Matrix *P);
void free_objective(Objective *obj);

/* Fixes variable k = 'col' to 'value' in the objective:
     offset += c_k value + (1/2) P_kk value^2,   c_l += P_kl value (l != k),
   then removes row and column k of P and marks the row inactive. The caller
   saves the postsolve record first since it needs c_k and row k of P. */
void fix_var_in_obj(Problem *prob, int col, double value);

/* Substitutes variable 'k' from the objective using row 'i'
   (assumed to be an equality). Row 'i' is specified by
   (vals, cols, len, rhs). Row k of P must be empty and is marked inactive. */
void sub_var_in_obj(Problem *prob, const double *vals, const int *cols, int len,
                    int k, double aik, double rhs);

/* Marks row k of P (which must be empty) as removed. */
void mark_P_row_inactive(Problem *prob, int k);

/* True if row k of P is a single positive diagonal entry, returned in *Pkk.
   A size-one row is not necessarily a diagonal: check_P_input tolerates tiny
   off-diagonals in a row with no diagonal. */
bool diagonal_only_P_entry(const struct Matrix *P, const int *P_row_sizes, int k,
                           double *Pkk);

/* True if row k of P holds a diagonal entry. */
bool P_row_has_diagonal(const struct Matrix *P, const int *P_row_sizes, int k);

/* Substitutes x_subst = shift + mult * x_stay in the objective and zeroes
   c[subst]. */
void sub_var_in_obj_dton(Objective *obj, int subst, int stay, double mult,
                         double shift);

/* Substitutes x_k = shift + mult * x_t in the objective when row k of P is the
   positive diagonal Pkk alone (diagonal_only_P_entry) and t has a diagonal:
     c_t += mult (c_k + P_kk shift),  offset += shift (c_k + (1/2) P_kk shift),
     P_tt += mult^2 P_kk,
   then empties row k of P. The caller marks it inactive. */
void sub_diag_var_in_obj_dton(Problem *prob, int k, int t, double mult,
                              double shift);

#endif // CORE_PROBLEM_H
