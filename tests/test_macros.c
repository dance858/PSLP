#include "test_macros.h"
#include "Bounds.h"
#include "Matrix.h"
#include "Numerics.h"
#include "minunit.h"
#include <stdio.h>

void init_bounds(Bound *bounds, double low_b, double upp_b, int len)
{
    for (int i = 0; i < len; i++)
    {
        bounds[i].lb = low_b;
        bounds[i].ub = upp_b;
    }
}

bool check_row_starts(const Matrix *A, const int *row_starts_correct)
{
    for (int i = 0; i < A->m + 1; ++i)
    {
        if (A->p[i].start != row_starts_correct[i])
        {
            printf("error, row_starts not equal at %d: %d != %d\n", i, A->p[i].start,
                   row_starts_correct[i]);
            return false;
        }
    }
    return true;
}

bool check_row_ends(const Matrix *A, const int *row_ends_correct)
{
    for (int i = 0; i < A->m + 1; ++i)
    {
        if (A->p[i].end != row_ends_correct[i])
        {
            printf("error, row_ends not equal at %d: %d != %d\n", i, A->p[i].end,
                   row_ends_correct[i]);
            return false;
        }
    }
    return true;
}

bool check_locks(const Lock *locks, const int *correct_up, const int *correct_down,
                 int size)
{
    for (int i = 0; i < size; ++i)
    {
        if (correct_up[i] != locks[i].up)
        {
            printf("Error up lock at %d: %d != %d\n", i, locks[i].up, correct_up[i]);
            return false;
        }
        if (correct_down[i] != locks[i].down)
        {
            printf("Error down lock at %d: %d != %d\n", i, locks[i].down,
                   correct_down[i]);
            return false;
        }
    }
    return true;
}

bool check_row_sizes(const Matrix *A, const int *row_sizes)
{
    for (int i = 0; i < A->m; i++)
    {
        int row_size = A->p[i].end - A->p[i].start;
        if (row_size != row_sizes[i])
        {
            return false;
        }
    }

    return true;
}

bool check_col_sizes(const Matrix *AT, const int *col_sizes)
{
    for (int i = 0; i < AT->m; i++)
    {
        int col_size = AT->p[i].end - AT->p[i].start;
        if (col_size != col_sizes[i])
        {
            return false;
        }
    }
    return true;
}

bool check_bounds(const Bound *bounds, const double *lbs, const double *ubs,
                  int size)
{
    for (int i = 0; i < size; i++)
    {
        if (bounds[i].lb != lbs[i] || bounds[i].ub != ubs[i])
        {
            return false;
        }
    }
    return true;
}

bool contains_int(const int *arr, int len, int val)
{
    for (int i = 0; i < len; ++i)
    {
        if (arr[i] == val)
        {
            return true;
        }
    }
    return false;
}

bool is_solution_correct(const double *x, const double *correct_x, const double *y,
                         const double *correct_y, const double *z,
                         const double *correct_z, int n_rows, int n_cols, double tol)
{
    for (int i = 0; i < n_cols; ++i)
    {
        if (ABS(x[i] - correct_x[i]) > tol)
        {
            printf("i = %d, x = %f, correct_x = %f\n", i, x[i], correct_x[i]);
            return false;
        }

        if (ABS(z[i] - correct_z[i]) > tol)
        {
            printf("i = %d, z = %f, correct_z = %f\n", i, z[i], correct_z[i]);
            return false;
        }
    }

    for (int i = 0; i < n_rows; ++i)
    {
        if (ABS(y[i] - correct_y[i]) > tol)
        {
            printf("i = %d, y = %f, correct_y = %f\n", i, y[i], correct_y[i]);
            return false;
        }
    }
    return true;
}