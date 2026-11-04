#ifndef TEST_ROWSLOTS_H
#define TEST_ROWSLOTS_H

#include "Matrix.h"
#include "RowSlots.h"
#include "minunit.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

int counter_rowslots = 0;

/* A 3 x 4 matrix with 5 entries and a spare tail of 6, rows
   {0: 1, 2: 2}, {1: 3}, {0: 4, 3: 5}, with its slots initialized. */
static Matrix *rowslots_matrix(RowSlots *slots, int *cap)
{
    Matrix *M = matrix_alloc(3, 4, 11);
    M->n_alloc = 11;
    int cols[] = {0, 2, 1, 0, 3};
    double vals[] = {1, 2, 3, 4, 5};
    int starts[] = {0, 2, 3, 5};
    for (int q = 0; q < 5; ++q)
    {
        M->i[q] = cols[q];
        M->x[q] = vals[q];
    }
    for (int r = 0; r < 3; ++r)
    {
        M->p[r].start = starts[r];
        M->p[r].end = starts[r + 1];
    }
    M->p[3].start = 5;
    M->p[3].end = 5;
    M->nnz = 5;
    slots->cap = cap;
    row_slots_init(slots, M);
    return M;
}

/* Row r holds exactly the entries (cols, vals). */
static bool rowslots_row_is(const Matrix *M, int r, const int *cols,
                            const double *vals, int len)
{
    if (M->p[r].end - M->p[r].start != len)
    {
        return false;
    }
    for (int q = 0; q < len; ++q)
    {
        if (M->i[M->p[r].start + q] != cols[q] || M->x[M->p[r].start + q] != vals[q])
        {
            return false;
        }
    }
    return true;
}

static char *test_rowslots_init()
{
    RowSlots s;
    int cap[3];
    Matrix *M = rowslots_matrix(&s, cap);
    mu_assert("caps", cap[0] == 2 && cap[1] == 1 && cap[2] == 2);
    mu_assert("tail", s.tail_base == 5 && s.tail_next == 5);
    free_matrix(M);
    return 0;
}

/* Row 0: overwrite column 0, insert column 1, column 2 dropped by its tag. The
   result has two entries and fits the slot. */
static char *test_rowslots_in_place()
{
    RowSlots s;
    int cap[3];
    Matrix *M = rowslots_matrix(&s, cap);
    uint8_t tags[] = {0, 0, 1, 0};
    int ucols[] = {0, 1};
    double uvals[] = {7, 8};
    mu_assert("update", matrix_update_row(M, &s, 0, ucols, uvals, 2, tags, 1));
    int cols[] = {0, 1};
    double vals[] = {7, 8};
    mu_assert("content", rowslots_row_is(M, 0, cols, vals, 2));
    mu_assert("in place", M->p[0].start == 0 && cap[0] == 2 && s.tail_next == 5);
    mu_assert("nnz", M->nnz == 5);
    free_matrix(M);
    return 0;
}

/* Row 2: a zero update deletes column 3; deleting the absent column 1 then
   changes nothing. */
static char *test_rowslots_delete()
{
    RowSlots s;
    int cap[3];
    Matrix *M = rowslots_matrix(&s, cap);
    int ucols[] = {3};
    double uvals[] = {0};
    mu_assert("delete", matrix_update_row(M, &s, 2, ucols, uvals, 1, NULL, 0));
    int cols[] = {0};
    double vals[] = {4};
    mu_assert("content", rowslots_row_is(M, 2, cols, vals, 1));
    mu_assert("nnz", M->nnz == 4);

    int absent[] = {1};
    mu_assert("absent", matrix_update_row(M, &s, 2, absent, uvals, 1, NULL, 0));
    mu_assert("unchanged", rowslots_row_is(M, 2, cols, vals, 1) && M->nnz == 4);
    free_matrix(M);
    return 0;
}

/* Row 1 gains columns 0 and 2: three entries do not fit its slot of one, so
   it moves to the tail. */
static char *test_rowslots_move_to_tail()
{
    RowSlots s;
    int cap[3];
    Matrix *M = rowslots_matrix(&s, cap);
    int ucols[] = {0, 2};
    double uvals[] = {1, 2};
    mu_assert("update", matrix_update_row(M, &s, 1, ucols, uvals, 2, NULL, 0));
    int cols[] = {0, 1, 2};
    double vals[] = {1, 3, 2};
    mu_assert("content", rowslots_row_is(M, 1, cols, vals, 3));
    mu_assert("moved", M->p[1].start == 5 && cap[1] == 3 && s.tail_next == 8);
    mu_assert("nnz", M->nnz == 7);
    int cols0[] = {0, 2}, cols2[] = {0, 3};
    double vals0[] = {1, 2}, vals2[] = {4, 5};
    mu_assert("other rows untouched", rowslots_row_is(M, 0, cols0, vals0, 2) &&
                                          rowslots_row_is(M, 2, cols2, vals2, 2));

    // the moved row is updated within its new slot and stays there
    int dcols[] = {2};
    double dvals[] = {0};
    mu_assert("update again", matrix_update_row(M, &s, 1, dcols, dvals, 1, NULL, 0));
    int cols1[] = {0, 1};
    double vals1[] = {1, 3};
    mu_assert("content again", rowslots_row_is(M, 1, cols1, vals1, 2));
    mu_assert("stays", M->p[1].start == 5 && cap[1] == 3 && s.tail_next == 8);
    free_matrix(M);
    return 0;
}

/* After row 1 moved, three tail entries are left; row 2 would grow to four
   entries, so the update fails and leaves row 2 untouched. */
static char *test_rowslots_tail_full()
{
    RowSlots s;
    int cap[3];
    Matrix *M = rowslots_matrix(&s, cap);
    int ucols[] = {0, 2};
    double uvals[] = {1, 2};
    mu_assert("first move", matrix_update_row(M, &s, 1, ucols, uvals, 2, NULL, 0));

    int gcols[] = {1, 2};
    double gvals[] = {6, 7};
    mu_assert("tail full", !matrix_update_row(M, &s, 2, gcols, gvals, 2, NULL, 0));
    int cols[] = {0, 3};
    double vals[] = {4, 5};
    mu_assert("row untouched", rowslots_row_is(M, 2, cols, vals, 2) &&
                                   M->p[2].start == 3 && cap[2] == 2);
    mu_assert("tail untouched", s.tail_next == 8 && M->nnz == 7);
    free_matrix(M);
    return 0;
}

/* The merge is built in the tail, so even an update that would stay in its
   slot fails once the tail is full, and leaves the row untouched. */
static char *test_rowslots_in_place_needs_room()
{
    RowSlots s;
    int cap[3];
    Matrix *M = rowslots_matrix(&s, cap);
    int cols1[] = {0, 2}, cols0[] = {1};
    double vals1[] = {1, 2}, vals0[] = {6};
    mu_assert("move row 1", matrix_update_row(M, &s, 1, cols1, vals1, 2, NULL, 0));
    mu_assert("move row 0", matrix_update_row(M, &s, 0, cols0, vals0, 1, NULL, 0));
    mu_assert("tail full", s.tail_next == 11 && M->nnz == 8);

    int dcols[] = {3};
    double dvals[] = {0};
    mu_assert("no room", !matrix_update_row(M, &s, 2, dcols, dvals, 1, NULL, 0));
    int cols[] = {0, 3};
    double vals[] = {4, 5};
    mu_assert("row untouched", rowslots_row_is(M, 2, cols, vals, 2) &&
                                   M->p[2].start == 3 && M->nnz == 8);
    free_matrix(M);
    return 0;
}

/* Without tags nothing is dropped: overwriting column 2 of row 0 keeps both
   entries. */
static char *test_rowslots_no_tags()
{
    RowSlots s;
    int cap[3];
    Matrix *M = rowslots_matrix(&s, cap);
    int ucols[] = {2};
    double uvals[] = {9};
    mu_assert("update", matrix_update_row(M, &s, 0, ucols, uvals, 1, NULL, 0));
    int cols[] = {0, 2};
    double vals[] = {1, 9};
    mu_assert("content", rowslots_row_is(M, 0, cols, vals, 2));
    free_matrix(M);
    return 0;
}

static const char *all_tests_rowslots()
{
    mu_run_test(test_rowslots_init, counter_rowslots);
    mu_run_test(test_rowslots_in_place, counter_rowslots);
    mu_run_test(test_rowslots_delete, counter_rowslots);
    mu_run_test(test_rowslots_move_to_tail, counter_rowslots);
    mu_run_test(test_rowslots_tail_full, counter_rowslots);
    mu_run_test(test_rowslots_in_place_needs_room, counter_rowslots);
    mu_run_test(test_rowslots_no_tags, counter_rowslots);
    return 0;
}

int test_rowslots()
{
    const char *result = all_tests_rowslots();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("row slots: TEST FAILED!\n");
    }
    else
    {
        printf("row slots: ALL TESTS PASSED\n");
    }
    printf("row slots: Tests run: %d\n", counter_rowslots);
    return result == 0;
}

#endif // TEST_ROWSLOTS_H
