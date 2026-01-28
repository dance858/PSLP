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

#include "Constraints.h"
#include "Bounds.h"
#include "Debugger.h"
#include "Locks.h"
#include "Matrix.h"
#include "State.h"
#include "Workspace.h"
#include "glbopts.h"
#include "utils.h"

Constraints *constraints_new(Matrix *A, Matrix *AT, double *lhs, double *rhs,
                             Bound *bounds, State *state, RowTag *row_tags,
                             ColTag *col_tags)
{
    Constraints *constraints = (Constraints *) ps_malloc(1, sizeof(Constraints));
    RETURN_PTR_IF_NULL(constraints, NULL);

    constraints->lhs = lhs;
    constraints->rhs = rhs;
    constraints->m = A->m;
    constraints->n = A->n;
    constraints->A = A;
    constraints->AT = AT;
    constraints->bounds = bounds;
    constraints->state = state;
    constraints->row_tags = row_tags;
    constraints->col_tags = col_tags;

    return constraints;
}

void free_constraints(Constraints *constraints)
{
    if (!constraints)
    {
        return;
    }

    PS_FREE(constraints->lhs);
    PS_FREE(constraints->rhs);
    free_matrix(constraints->A);
    free_matrix(constraints->AT);
    PS_FREE(constraints->bounds);
    PS_FREE(constraints->row_tags);
    PS_FREE(constraints->col_tags);
    PS_FREE(constraints);
}

void delete_inactive_rows(Constraints *constraints)
{
    iVec *rows_to_delete = constraints->state->rows_to_delete;

    if (rows_to_delete->len == 0)
    {
        return;
    }

    int i, j, row, col, shift;
    bool is_rhs_inf, is_lhs_inf;
    int *row_sizes = constraints->state->row_sizes;
    int *col_sizes = constraints->state->col_sizes;
    Lock *locks = constraints->state->col_locks;
    Matrix *A = constraints->A;
    Matrix *AT = constraints->AT;
    RowTag *row_tags = constraints->row_tags;
    RowRange *row_r = A->p;
    RowRange *col_r = AT->p;

    // ------------------------------------------------------------------------
    //                      Delete rows of A.
    //            (We also update column sizes and locks).
    // ------------------------------------------------------------------------
    for (i = 0; i < rows_to_delete->len; ++i)
    {
        row = rows_to_delete->data[i];
        assert(row_sizes[row] != SIZE_INACTIVE_ROW);

        // skip empty rows
        if (row_sizes[row] == 0)
        {
            row_sizes[row] = SIZE_INACTIVE_ROW;
            continue;
        }

        A->nnz -= (size_t) row_sizes[row];
        row_sizes[row] = SIZE_INACTIVE_ROW;
        is_rhs_inf = HAS_TAG(row_tags[row], R_TAG_RHS_INF);
        is_lhs_inf = HAS_TAG(row_tags[row], R_TAG_LHS_INF);

        for (j = row_r[row].start; j < row_r[row].end; ++j)
        {
            col = A->i[j];
            if (col_sizes[col] == SIZE_INACTIVE_COL)
            {
                continue;
            }

            col_sizes[col]--;
            update_lock_coeff_deletion(locks + col, A->x[j], is_rhs_inf, is_lhs_inf);
        }
        row_r[row].end = row_r[row].start;
    }

    // ------------------------------------------------------------------------
    // Process each row of AT and place coefficients contiguously in memory.
    // (Looping over rows of AT is much more cache friendly than looping over
    //  the columns of AT that change)
    // ------------------------------------------------------------------------
    iVec *empty_cols = constraints->state->empty_cols;
    iVec *ston_cols = constraints->state->ston_cols;
    const ColTag *col_tags = constraints->col_tags;
    for (col = 0; col < AT->m; ++col)
    {
        // skip rows of AT that haven't changed
        if (HAS_TAG(col_tags[col], C_TAG_INACTIVE) ||
            col_sizes[col] == col_r[col].end - col_r[col].start)
        {
            continue;
        }

        assert(col_sizes[col] != SIZE_INACTIVE_COL);

        // check for new empty and singleton columns,
        switch (col_sizes[col])
        {
            case 0:
                iVec_append(empty_cols, col);
                col_r[col].end = col_r[col].start;
                break;
            case 1:
                iVec_append(ston_cols, col);
                break;
        }

        // place coefficients contiguously in memory
        shift = 0;
        for (j = col_r[col].start; j < col_r[col].end; ++j)
        {
            if (row_sizes[AT->i[j]] == SIZE_INACTIVE_ROW)
            {
                shift++;
            }
            else
            {
                AT->x[j - shift] = AT->x[j];
                AT->i[j - shift] = AT->i[j];
            }
        }

        // update end of column
        col_r[col].end -= shift;
        assert(col_r[col].start + col_sizes[col] == col_r[col].end);
    }

    AT->nnz = A->nnz;
    iVec_clear_no_resize(rows_to_delete);
}

void delete_inactive_cols_from_A_and_AT(Constraints *constraints)
{
    int i, j, row, col, shift;
    iVec *fixed_cols_to_delete = constraints->state->fixed_cols_to_delete;
    iVec *sub_cols_to_delete = constraints->state->sub_cols_to_delete;
    int *row_sizes = constraints->state->row_sizes;
    int *col_sizes = constraints->state->col_sizes;
    Matrix *A = constraints->A;
    const Matrix *AT = constraints->AT;
    const RowTag *row_tags = constraints->row_tags;
    RowRange *row_r = A->p;
    RowRange *col_r = AT->p;

    // ------------------------------------------------------------------------
    //              Delete rows of A and update column sizes.
    // ------------------------------------------------------------------------
    for (i = 0; i < fixed_cols_to_delete->len; ++i)
    {
        col = fixed_cols_to_delete->data[i];
        col_sizes[col] = SIZE_INACTIVE_COL;

        for (j = col_r[col].start; j < col_r[col].end; ++j)
        {
            if (row_sizes[AT->i[j]] == SIZE_INACTIVE_ROW)
            {
                continue;
            }

            row_sizes[AT->i[j]]--;
        }
        col_r[col].end = col_r[col].start;
    }

    for (i = 0; i < sub_cols_to_delete->len; ++i)
    {
        col = sub_cols_to_delete->data[i];
        col_sizes[col] = SIZE_INACTIVE_COL;

        for (j = col_r[col].start; j < col_r[col].end; ++j)
        {
            if (row_sizes[AT->i[j]] == SIZE_INACTIVE_ROW)
            {
                continue;
            }

            row_sizes[AT->i[j]]--;
        }
        col_r[col].end = col_r[col].start;
    }

    // ------------------------------------------------------------------------
    // Process each row of A and place coefficients contiguously in memory.
    // ------------------------------------------------------------------------
    iVec *empty_rows = constraints->state->empty_rows;
    iVec *ston_rows = constraints->state->ston_rows;
    iVec *dton_rows = constraints->state->dton_rows;
    for (row = 0; row < A->m; ++row)
    {
        // skip rows of A that haven't changed
        if (HAS_TAG(row_tags[row], R_TAG_INACTIVE) ||
            row_sizes[row] == row_r[row].end - row_r[row].start)
        {
            continue;
        }

        assert(row_sizes[row] != SIZE_INACTIVE_ROW &&
               !HAS_TAG(row_tags[row], R_TAG_INACTIVE));

        // check for new empty, singleton and doubleton rows
        switch (row_sizes[row])
        {
            case 0:
                assert(!iVec_contains(empty_rows, row));
                iVec_append(empty_rows, row);
                A->nnz -= (size_t) (row_r[row].end - row_r[row].start);
                row_r[row].end = row_r[row].start;
                break;
            case 1:
                assert(!iVec_contains(ston_rows, row));
                iVec_append(ston_rows, row);
                break;
            case 2:
                if (HAS_TAG(row_tags[row], R_TAG_EQ))
                {
                    assert(!iVec_contains(dton_rows, row));
                    iVec_append(dton_rows, row);
                }
                break;
        }

        // place coefficients contiguously in memory
        shift = 0;
        for (j = row_r[row].start; j < row_r[row].end; ++j)
        {
            col = A->i[j];
            if (col_sizes[col] == SIZE_INACTIVE_COL)
            {
                shift++;
            }
            else
            {
                A->x[j - shift] = A->x[j];
                A->i[j - shift] = A->i[j];
            }
        }

        // update end of row
        row_r[row].end -= shift;
        A->nnz -= (size_t) shift;
        assert(row_r[row].start + row_sizes[row] == row_r[row].end);
    }

    constraints->AT->nnz = A->nnz;
    iVec_clear_no_resize(fixed_cols_to_delete);
    iVec_clear_no_resize(sub_cols_to_delete);
}

static void bounds_shrink(Bound *ptr, int *map, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        if (map[i] != -1)
        {
            ptr[map[i]] = ptr[i];
        }
    }
}

void constraints_clean(Constraints *constraints, Mapping *map, bool remove_all)
{
    const int *row_sizes = constraints->state->row_sizes;
    const int *col_sizes = constraints->state->col_sizes;

    remove_extra_space(constraints->A, row_sizes, col_sizes, remove_all, map->cols);
    remove_extra_space(constraints->AT, col_sizes, row_sizes, remove_all, map->rows);

    dPtr_shrink(constraints->rhs, map->rows, constraints->m);
    dPtr_shrink(constraints->lhs, map->rows, constraints->m);
    rowTagPtr_shrink(constraints->row_tags, map->rows, constraints->m);
    colTagPtr_shrink(constraints->col_tags, map->cols, constraints->n);
    bounds_shrink(constraints->bounds, map->cols, constraints->n);
    constraints->m = constraints->A->m;
    constraints->n = constraints->A->n;
}
