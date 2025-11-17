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

#ifndef LOCKS_H
#define LOCKS_H

#include <assert.h>
#include <stdbool.h>

#include "Tags.h"
struct Matrix;
struct RowView;

typedef struct Lock
{
    int up;
    int down;
} Lock;

/* Contructor that counts locks by traversing the rows of A */
Lock *new_locks(const struct Matrix *A, const RowTag *row_tags);
void free_locks(Lock *locks);

/* Updates locks when an equality constraint is transformed to an
   inequality constraint*/
void update_locks_eq_to_ineq(Lock *locks, const double *vals, const int *cols,
                             int len, bool to_upper);

/* Updates locks when an inequality constraint is transformed to an
   equality constraint*/
void update_locks_ineq_to_eq(Lock *locks, const double *vals, const int *cols,
                             int len, bool is_new_constraint_lhs);

// col should be a row view corresponding to the column
void count_locks_one_column(const struct RowView *col, Lock *lock,
                            const RowTag *row_tags);

static inline void update_lock_coeff_deletion(struct Lock *lock, double val,
                                              int isRhsInf, int isLhsInf)
{
    if (val > 0)
    {
        if (!isRhsInf)
        {
            lock->up--;
        }
        if (!isLhsInf)
        {
            lock->down--;
        }
    }
    else
    {
        assert(val != 0);
        if (!isRhsInf)
        {
            lock->down--;
        }
        if (!isLhsInf)
        {
            lock->up--;
        }
    }
}

#endif // LOCKS_H