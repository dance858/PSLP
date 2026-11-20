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

#include "Locks.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "RowColViews.h"
#include <assert.h>

void locks_init_cols(const Matrix *AT, const RowTag *row_tags, Lock *locks,
                     int start, int end)
{
    for (int col = start; col < end; ++col)
    {
        int up = 0, down = 0;
        int col_end = AT->p[col].end;
        for (int idx = AT->p[col].start; idx < col_end; ++idx)
        {
            RowTag tag = row_tags[AT->i[idx]];
            int rhs_finite = !HAS_TAG(tag, R_TAG_RHS_INF);
            int lhs_finite = !HAS_TAG(tag, R_TAG_LHS_INF);
            if (AT->x[idx] > 0)
            {
                up += rhs_finite;
                down += lhs_finite;
            }
            else
            {
                down += rhs_finite;
                up += lhs_finite;
            }
        }
        locks[col].up = up;
        locks[col].down = down;
    }
}

void free_locks(Lock *locks)
{
    PS_FREE(locks);
}

void update_locks_eq_to_ineq(Lock *locks, const double *vals, const int *cols,
                             int len, bool to_upper)
{
    // --------------------------------------------------------
    // if the new constraint is an upper inequality constraint
    // --------------------------------------------------------
    if (to_upper)
    {
        for (int i = 0; i < len; i++)
        {
            assert(vals[i] != 0);
            if (vals[i] > 0)
            {
                locks[cols[i]].down--;
            }
            else
            {
                locks[cols[i]].up--;
            }
        }
    }
    // --------------------------------------------------------
    // if the new constraint is a lower inequality constraint
    // --------------------------------------------------------
    else
    {
        for (int i = 0; i < len; i++)
        {
            assert(vals[i] != 0);
            if (vals[i] > 0)
            {
                locks[cols[i]].up--;
            }
            else
            {
                locks[cols[i]].down--;
            }
        }
    }
}

void update_locks_ineq_to_eq(Lock *locks, const double *vals, const int *cols,
                             int len, bool is_new_constraint_lhs)
{
    // --------------------------------
    // if we get a new LHS constraint
    // --------------------------------
    if (is_new_constraint_lhs)
    {
        for (int i = 0; i < len; i++)
        {
            assert(vals[i] != 0);
            if (vals[i] > 0)
            {
                locks[cols[i]].down++;
            }
            else
            {
                locks[cols[i]].up++;
            }
        }
    }
    // -------------------------------
    // if we get a new RHS constraint
    // -------------------------------
    else
    {
        for (int i = 0; i < len; i++)
        {
            assert(vals[i] != 0);
            if (vals[i] > 0)
            {
                locks[cols[i]].up++;
            }
            else
            {
                locks[cols[i]].down++;
            }
        }
    }
}

void count_locks_one_column(const RowView *col, Lock *lock, const RowTag *row_tags)
{
    lock->up = 0;
    lock->down = 0;

    const int *rows = col->cols;

    for (int i = 0; i < *col->len; ++i)
    {
        if (col->vals[i] > 0)
        {
            if (!HAS_TAG(row_tags[rows[i]], R_TAG_RHS_INF))
            {
                lock->up++;
            }

            if (!HAS_TAG(row_tags[rows[i]], R_TAG_LHS_INF))
            {
                lock->down++;
            }
        }
        else
        {
            assert(col->vals[i] < 0);
            if (!HAS_TAG(row_tags[rows[i]], R_TAG_RHS_INF))
            {
                lock->down++;
            }

            if (!HAS_TAG(row_tags[rows[i]], R_TAG_LHS_INF))
            {
                lock->up++;
            }
        }
    }
}
