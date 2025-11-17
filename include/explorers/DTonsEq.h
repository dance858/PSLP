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

#ifndef DTONS_EQ_H
#define DTONS_EQ_H

#include "PresolveStatus.h"
struct Matrix;
struct Problem;
struct PostsolveInfo;

/* This function removes doubleton equality rows from the problem. All internal
   data structures (such as row sizes, columns sizes, activities, bounds, lists
   of ston rows etc.) are updated accordingly. Doubleton rows that appear
   when another doubleton row is eliminated are also removed. After this
   function, the list of doubleton rows should be empty.

   When this function finishes, the problem is up to date, and no extra
   synchronization is needed.

   When we process a doubleton row, we first choose the variable that should be
   eliminated according to the following priorities.
   1. Integral ratios. (For example, if the first variable gives an integral
      ratio but not the second, we will substitute the first variable regardless
      of its column length etc.)
   2. Column lengths. (We substitute the variable with the smallest length.)
   3. We then have a simple check to reject large or small pivots.
   4. We then theck if there is sufficient space in AT for the substitution.
      If yes, we substitute the variable. If no, we reject the reduction.
      A possible modification is to check if it is possible to substitute the
      other variable if there wasn't enough space for substituting the other
      variable.

    A doubleton row is not eliminated if:
      * there is not sufficient space in A transpose to store the fill-in
        caused by the substitution.
      * the pivot is too large or too small.

   The current implementation is not that sophisticated. We might want to
   reject the reduction if the fill-in is too large, or add a more sophisticated
   rejection criteria for numerical stability. Perhaps we should only
   allow the substitution if the new bound is not too large?
 */
PresolveStatus remove_dton_eq_rows(struct Problem *prob, int max_shift);

typedef struct
{
    double old_coeff_stay;
    double old_coeff_subst;
    double new_coeff_stay;
} Old_and_new_coeff;

/* This function updates row 'q' in A when the variable 'k' is
   substituted using a doubleton row. The argument 'ratio' is aij / aik
   where the doubleton row that is substituted is aij * xj + aik * xk = rhs,
   and xj stays and xk is substituted.

   The function returns the new coefficient of the variable that stays.
   Note that we delete the substituted variable from the row before we insert
   the new variable, so we can substitute xk for xj even in rows with no extra
   space.

   This function is static inside DTonsEq.c, but we expose it for testing
   purposes.
*/
#ifdef TESTING
Old_and_new_coeff update_row_A_dton(struct Matrix *A, int i, int q, int j, int k,
                                    double aij, double aik, int *row_size,
                                    struct PostsolveInfo *postsolve_info);
#endif

#endif
