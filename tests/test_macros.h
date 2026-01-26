#ifndef TEST_MACROS_H
#define TEST_MACROS_H

#include "Bounds.h"
#include "Locks.h"
#include "Matrix.h"
#include "Numerics.h"

void init_bounds(Bound *bounds, double low_b, double upp_b, int len);
bool check_row_starts(const Matrix *A, const int *row_starts_correct);
bool check_row_ends(const Matrix *A, const int *row_ends_correct);
bool check_locks(const Lock *locks, const int *correct_up, const int *correct_down,
                 int size);
bool CHECK_ROW_SIZES(const Matrix *A, const int *row_sizes);
bool CHECK_COL_SIZES(const Matrix *AT, const int *col_sizes);
bool CHECK_BOUNDS(const Bound *bounds, const double *lbs, const double *ubs,
                  int size);
bool contains_int(const int *arr, int len, int val);
bool is_solution_correct(const double *x, const double *correct_x, const double *y,
                         const double *correct_y, const double *z,
                         const double *correct_z, int n_rows, int n_cols,
                         double tol);

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

#endif // TEST_MACROS_H
