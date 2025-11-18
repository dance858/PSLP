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

#ifndef STONCOLS_HPP
#define STONCOLS_HPP

#include "PSLP_status.h"

// forward declaration
struct Problem;

/*  The following reductions are possible for a column singleton xk appearing
    in an equality i:
        1. Equality i is redundant if xk is (implied) free.
        2. If xk is implied free from either below or above it can be removed
           if the equality is changed to an inequality:

    Original constraint: aik*xk + \sum_{j \neq k} aij * xj = bi

       1. If aik > 0 and xk is implied free from above the new constraint is
          \sum_{j \neq k} aij * xj <= bi - ell_k * aik.
       2. If aik < 0 and xk is implied free from above the new constraint is
          \sum_{j \neq k} aij * xj >= bi - ell_k * aik.
       3. If aik > 0 and xk is implied free from below the new constraint is
          \sum_{j \neq k} aij * xj >= bi - u_k * aik.
       4. If aik < 0 and xk is implied free from below the new constraint is
          \sum_{j \neq k} aij * xj <= bi - u_k * aik.

       Why do we want to transform an equality constraint to an inequality,
       you may ask. This may remove locks and lead to further reductions etc.

 The following reductions are possible for a column singleton in an inequality
 constraint:
    1. Dual fix to lower bound.
    2. Dual fix to upper bound.
    3. Convert lower/upper inequality to equality constraint.
    4. The constraint is redundant if the column singleton is free.
*/

/* This function removes column singletons from the problem. When this function
   finishes, the problem is up to date, and no extra synchronization is needed.
   This also means that the internal data structures have been updated.

   After having deleted the column singletons, new column singletons may arise.
   This function deletes those new singletons as well, and it continues until
   no new column singletons are found. In other words, it eliminates chains
   of column singletons.

   The function returns UNBNDORINFEAS if the problem is unbounded or infeasible,
   or UNCHANGED (even if the problem was reduced).
*/
PresolveStatus remove_ston_cols(struct Problem *prob);

#endif // STONCOLS_HPP
