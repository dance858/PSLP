#ifndef TEST_MACROS_H
#define TEST_MACROS_H

#include "Bounds.h"
#include "Matrix.h"
#include "Numerics.h"

#define INIT_BOUNDS(bounds, low_b, upp_b, len)                                      \
    do                                                                              \
    {                                                                               \
        for (size_t i = 0; i < (len); i++)                                          \
        {                                                                           \
            (bounds)[i].lb = (low_b);                                               \
            (bounds)[i].ub = (upp_b);                                               \
        }                                                                           \
    } while (0)

#define CHECK_ROW_STARTS(A, row_starts_correct)                                     \
    do                                                                              \
    {                                                                               \
        int row_starts[A->m + 1];                                                   \
        for (int i = 0; i < A->m + 1; ++i)                                          \
        {                                                                           \
            row_starts[i] = A->p[i].start;                                          \
        }                                                                           \
        mu_assert("error, row_starts not equal",                                    \
                  ARRAYS_EQUAL(row_starts, row_starts_correct, A->m + 1));          \
    } while (0)

#define CHECK_ROW_ENDS(A, row_ends_correct)                                         \
    do                                                                              \
    {                                                                               \
        int row_ends[A->m + 1];                                                     \
        for (int i = 0; i < A->m + 1; ++i)                                          \
        {                                                                           \
            row_ends[i] = A->p[i].end;                                              \
        }                                                                           \
        mu_assert("error, row_starts not equal",                                    \
                  ARRAYS_EQUAL(row_ends, row_ends_correct, A->m + 1));              \
    } while (0)

#define CHECK_LOCKS(locks, correct_up, correct_down, size)                          \
    do                                                                              \
    {                                                                               \
        for (int i = 0; i < size; ++i)                                              \
        {                                                                           \
            mu_assert("Error up lock", correct_up[i] == locks[i].up);               \
            mu_assert("Error down lock", correct_down[i] == locks[i].down);         \
        }                                                                           \
    } while (0)

#define PRINT_ROW_STARTS(A)                                                         \
    do                                                                              \
    {                                                                               \
        for (size_t i = 0; i < A->m + 1; ++i)                                       \
        {                                                                           \
            printf("%d ", A->p[i].start);                                           \
        }                                                                           \
        printf("\n");                                                               \
    } while (0)

#define PRINT_ROW_ENDS(A)                                                           \
    do                                                                              \
    {                                                                               \
        for (size_t i = 0; i < A->m + 1; ++i)                                       \
        {                                                                           \
            printf("%d ", A->p[i].end);                                             \
        }                                                                           \
        printf("\n");                                                               \
    } while (0)

bool CHECK_ROW_SIZES(const Matrix *A, const int *row_sizes)
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

bool CHECK_COL_SIZES(const Matrix *AT, const int *col_sizes)
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

bool CHECK_BOUNDS(const Bound *bounds, const double *lbs, const double *ubs,
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

#endif // TEST_MACROS_H