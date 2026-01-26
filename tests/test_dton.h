#ifndef TEST_DTON_H
#define TEST_DTON_H

#include "DTonsEq.h"
#include "Debugger.h"
#include "PSLP_API.h"

#include "Problem.h"
#include "Workspace.h"
#include "debug_macros.h"
#include "minunit.h"
#include <stdio.h>

static int counter_dton = 0;

/*  test update_row_A_dton
    Substitution: x1 = 1 - 2x2
    x1  +     + x3 + x4 = 1  -> -2x2 + x3 + x4 = 0   (fill-in)
    x1  + x2  + x3 + x4 = 1  -> -x2  + x3 + x4 = 0   (inplace)
    2x1 + 4x2 + x3 + x4 = 1  ->        x3 + x4 = 1   (unexpected cancellation)
*/
static char *test_00_dton()
{
    double Ax[] = {1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1};
    int Ai[] = {0, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    int Ap[] = {0, 3, 7, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 4;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);

    int stay = 1;
    int subst = 0;
    double ratio = 2;
    int row_sizes[] = {3, 4, 4};
    PostsolveInfo *postsolve_info = postsolve_info_new(n_rows, n_cols);
    update_row_A_dton(A, 0, 0, stay, subst, 2, 1, row_sizes + 0, postsolve_info);
    update_row_A_dton(A, 0, 1, stay, subst, 2, 1, row_sizes + 1, postsolve_info);
    update_row_A_dton(A, 0, 2, stay, subst, 2, 1, row_sizes + 2, postsolve_info);

    // check that new row sizes are correct
    int row_sizes_correct[] = {3, 3, 2};
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes_correct, row_sizes, 3));

    int col_sizes[] = {SIZE_INACTIVE_COL, 2, 3, 3};
    int map[4] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    // check that new A is correct
    double Ax_correct[] = {-2, 1, 1, -1, 1, 1, 1, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 1, 2};
    int Ap_correct[] = {0, 3, 6, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));
    mu_assert("wrong nnz", A->nnz == 8);

    postsolve_info_free(postsolve_info);
    free_matrix(A);
    return 0;
}

/*  test update_row_A_dton
    Substitution: x3 = 1 - 2x4
    x1 + 2x2 + x3  +      + x5 + x6 = 1  ->  x1 + 2x2 - 2x4 + x5 + x6 =  0
   (fill-in) x1 +  x2 + x3  + x4   + x5 + x6 = 1  ->  x1 +  x2 -  x4 + x5 + x6 =
   0   (inplace) x1 +  x2 + 2x3 + 4x4  + x5 + x6 = 1  ->  x1 +  x2 +     + x5 +
   x6 = -1   (unexpected cancellation)
*/
static char *test_01_dton()
{
    double Ax[] = {1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1};
    int Ai[] = {0, 1, 2, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
    int Ap[] = {0, 5, 11, 17};
    int nnz = 17;
    int n_rows = 3;
    int n_cols = 6;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);

    int stay = 3;
    int subst = 2;
    double ratio = 2;
    int row_sizes[] = {5, 6, 6};
    PostsolveInfo *postsolve_info = postsolve_info_new(n_rows, n_cols);
    update_row_A_dton(A, 0, 0, stay, subst, 2, 1, row_sizes + 0, postsolve_info);
    update_row_A_dton(A, 0, 1, stay, subst, 2, 1, row_sizes + 1, postsolve_info);
    update_row_A_dton(A, 0, 2, stay, subst, 2, 1, row_sizes + 2, postsolve_info);

    // check that new row sizes are correct
    int row_sizes_correct[] = {5, 5, 4};
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes_correct, row_sizes, 3));

    int col_sizes[] = {3, 3, SIZE_INACTIVE_COL, 2, 3, 3};
    int map[6] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    // check that new A is correct
    double Ax_correct[] = {1, 2, -2, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1};
    int Ai_correct[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 3, 4};
    int Ap_correct[] = {0, 5, 10, 14};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 14));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 14));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));
    mu_assert("wrong nnz", A->nnz == 14);

    postsolve_info_free(postsolve_info);
    free_matrix(A);
    return 0;
}

/*  test update_row_A_dton
    Substitution: x5 = 1 - 2x6
    x1 + 2x2 + x3  + x4   + x5 +    = 1  ->  x1 + 2x2 + x3 + x4 +   - 2x6 =  0
   (fill-in) x1 +  x2 + x3  + x4   + x5 + x6 = 1  ->  x1 +  x2 + x3 + x4 +   -
   x6 =  0   (inplace) x1 +  x2 + 2x3 + 4x4  + x5 + x6 = 1  ->  x1 +  x2 + x3 +
   x4           = -1   (unexpected cancellation)
*/
static char *test_02_dton()
{
    double Ax[] = {1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
    int Ap[] = {0, 5, 11, 17};
    int nnz = 17;
    int n_rows = 3;
    int n_cols = 6;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);

    int stay = 5;
    int subst = 4;
    double ratio = 2;
    int row_sizes[] = {5, 6, 6};
    PostsolveInfo *postsolve_info = postsolve_info_new(n_rows, n_cols);
    update_row_A_dton(A, 0, 0, stay, subst, 2, 1, row_sizes + 0, postsolve_info);
    update_row_A_dton(A, 0, 1, stay, subst, 2, 1, row_sizes + 1, postsolve_info);
    update_row_A_dton(A, 0, 2, stay, subst, 2, 1, row_sizes + 2, postsolve_info);

    // check that new row sizes are correct
    int row_sizes_correct[] = {5, 5, 4};
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes_correct, row_sizes, 3));

    int col_sizes[] = {3, 3, 3, 3, SIZE_INACTIVE_COL, 2};
    int map[6] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    // check that new A is correct
    double Ax_correct[] = {1, 2, 1, 1, -2, 1, 1, 1, 1, -1, 1, 1, 1, 1};
    int Ai_correct[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3};
    int Ap_correct[] = {0, 5, 10, 14};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 14));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 14));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));
    mu_assert("wrong nnz", A->nnz == 14);

    postsolve_info_free(postsolve_info);
    free_matrix(A);
    return 0;
}

/*  test update_row_A_dton
    Substitution: x2 = 0.5 - 0.5x1
    x1  +     + x3 + x4 = 1  ->     x1  + x3 + x4 = 1        (unchanged)
    x1  + x2  + x3 + x4 = 1  ->  0.5x1  + x3 + x4 = 0.5      (inplace)
    2x1 + 4x2 + x3 + x4 = 1  ->           x3 + x4 = 1        (unexpected
   cancellation)
*/
static char *test_03_dton()
{
    double Ax[] = {1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1};
    int Ai[] = {0, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    int Ap[] = {0, 3, 7, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 4;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);

    int stay = 0;
    int subst = 1;
    double ratio = 0.5;
    int row_sizes[] = {3, 4, 4};
    PostsolveInfo *postsolve_info = postsolve_info_new(n_rows, n_cols);
    update_row_A_dton(A, 0, 1, stay, subst, 1, 2, row_sizes + 1, postsolve_info);
    update_row_A_dton(A, 0, 2, stay, subst, 1, 2, row_sizes + 2, postsolve_info);

    // check that new row sizes are correct
    int row_sizes_correct[] = {3, 3, 2};
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes_correct, row_sizes, 3));

    int col_sizes[] = {2, SIZE_INACTIVE_COL, 3, 3};
    int map[4] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    // check that new A is correct
    double Ax_correct[] = {1, 1, 1, 0.5, 1, 1, 1, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 1, 2};
    int Ap_correct[] = {0, 3, 6, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));
    mu_assert("wrong nnz", A->nnz == 8);

    postsolve_info_free(postsolve_info);
    free_matrix(A);
    return 0;
}

/*  test update_row_A_dton
    Substitution: x4 = 0.5 - 0.5x3
    x1 + 2x2 + x3  +      + x5 + x6 = 1  ->  x1 + 2x2 +    x3  +      + x5 + x6
   = 1     (unchanged) x1 +  x2 + x3  + x4   + x5 + x6 = 1  ->  x1 +  x2 + 0.5x3
   +      + x5 + x6 = 0.5   (inplace) x1 +  x2 + 2x3 + 4x4  + x5 + x6 = 1  -> x1
   +  x2 +               + x5 + x6 = -1     unexpected cancellation)
*/
static char *test_04_dton()
{
    double Ax[] = {1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1};
    int Ai[] = {0, 1, 2, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
    int Ap[] = {0, 5, 11, 17};
    int nnz = 17;
    int n_rows = 3;
    int n_cols = 6;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);

    int stay = 2;
    int subst = 3;
    double ratio = 0.5;
    int row_sizes[] = {5, 6, 6};
    PostsolveInfo *postsolve_info = postsolve_info_new(n_rows, n_cols);
    update_row_A_dton(A, 0, 1, stay, subst, 1, 2, row_sizes + 1, postsolve_info);
    update_row_A_dton(A, 0, 2, stay, subst, 1, 2, row_sizes + 2, postsolve_info);

    // check that new row sizes are correct
    int row_sizes_correct[] = {5, 5, 4};
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes_correct, row_sizes, 3));

    int col_sizes[] = {3, 3, 2, SIZE_INACTIVE_COL, 3, 3};
    int map[6] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    // check that new A is correct
    double Ax_correct[] = {1, 2, 1, 1, 1, 1, 1, 0.5, 1, 1, 1, 1, 1, 1};
    int Ai_correct[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 3, 4};
    int Ap_correct[] = {0, 5, 10, 14};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 14));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 14));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));
    mu_assert("wrong nnz", A->nnz == 14);

    postsolve_info_free(postsolve_info);
    free_matrix(A);
    return 0;
}

/*  test update_row_A_dton
    Substitution: x6 = 0.5 - 0.5x5
    x1 + 2x2 + x3  + x4   + x5 +    = 1  ->  x1 + 2x2 + x3  + x4   +  x5 +  = 1
   (unchanged-in) x1 +  x2 + x3  + x4   + x5 + x6 = 1  ->  x1 +  x2 + x3 +  x4 +
   0.5 x5   = 0.5   (inplace) x1 +  x2 + 2x3 + 4x4  + x5 + x6 = 1  ->  x1 +  x2
   + x3 +  x4            = -1    (unexpected cancellation)
*/
static char *test_05_dton()
{
    double Ax[] = {1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
    int Ap[] = {0, 5, 11, 17};
    int nnz = 17;
    int n_rows = 3;
    int n_cols = 6;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);

    int stay = 4;
    int subst = 5;
    double ratio = 0.5;
    int row_sizes[] = {5, 6, 6};
    PostsolveInfo *postsolve_info = postsolve_info_new(n_rows, n_cols);
    update_row_A_dton(A, 0, 1, stay, subst, 1, 2, row_sizes + 1, postsolve_info);
    update_row_A_dton(A, 0, 2, stay, subst, 1, 2, row_sizes + 2, postsolve_info);

    // check that new row sizes are correct
    int row_sizes_correct[] = {5, 5, 4};
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes_correct, row_sizes, 3));

    int col_sizes[] = {3, 3, 3, 3, 2, SIZE_INACTIVE_COL};
    int map[6] = {0};
    remove_extra_space(A, row_sizes, col_sizes, true, map);

    // check that new A is correct
    double Ax_correct[] = {1, 2, 1, 1, 1, 1, 1, 1, 1, 0.5, 1, 1, 1, 1};
    int Ai_correct[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3};
    int Ap_correct[] = {0, 5, 10, 14};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 14));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 14));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));
    mu_assert("wrong nnz", A->nnz == 14);

    postsolve_info_free(postsolve_info);
    free_matrix(A);
    return 0;
}

/*  test update_row_A_dton, WITH NO EXTRA SPACE
    Substitution: x1 = 1 - 2x2
    x1  +     + x3 + x4 = 1  -> -2x2 + x3 + x4 = 0   (fill-in)
    x1  + x2  + x3 + x4 = 1  -> -x2  + x3 + x4 = 0   (inplace)
    2x1 + 4x2 + x3 + x4 = 1  ->        x3 + x4 = 1   (unexpected cancellation)
*/
static char *test_06_dton()
{
    double Ax[] = {1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1};
    int Ai[] = {0, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    int Ap[] = {0, 3, 7, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 4;

    Matrix *A = matrix_new(Ax, Ai, Ap, n_rows, n_cols, nnz);
    int row_sizes[] = {3, 4, 4};
    int col_sizes[] = {3, 2, 3, 3};
    int map[4] = {0};

    remove_extra_space(A, row_sizes, col_sizes, true, map);

    int stay = 1;
    int subst = 0;
    double ratio = 2;
    PostsolveInfo *postsolve_info = postsolve_info_new(n_rows, n_cols);
    update_row_A_dton(A, 0, 0, stay, subst, 2, 1, row_sizes + 0, postsolve_info);
    update_row_A_dton(A, 0, 1, stay, subst, 2, 1, row_sizes + 1, postsolve_info);
    update_row_A_dton(A, 0, 2, stay, subst, 2, 1, row_sizes + 2, postsolve_info);

    // check that new row sizes are correct
    int row_sizes_correct[] = {3, 3, 2};
    mu_assert("error row_sizes", ARRAYS_EQUAL_INT(row_sizes_correct, row_sizes, 3));

    int col_sizes_new[] = {SIZE_INACTIVE_COL, 2, 3, 3};
    remove_extra_space(A, row_sizes, col_sizes_new, true, map);

    // check that new A is correct
    double Ax_correct[] = {-2, 1, 1, -1, 1, 1, 1, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 1, 2};
    int Ap_correct[] = {0, 3, 6, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));
    mu_assert("wrong nnz", A->nnz == 8);

    postsolve_info_free(postsolve_info);
    free_matrix(A);
    return 0;
}

/* dton row eliminated, bounds on the other variable are tightened.
min. [2 -1 -1 3 2]x
s.t. [2  1  -1  3  2]     [2]
     [1 -2   0  0  0] x = [-1]
     [1  4   0  3  1]     [4]
     0.4 <= x1 <= 0.6
     x2, x3, x4, x5 >= 0
*/
static char *test_1_dton()
{
    double Ax[] = {2, 1, -1, 3, 2, 1, -2, 1, 4, 3, 1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 3, 4};
    int Ap[] = {0, 5, 7, 11};
    int nnz = 11;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -1, 4};
    double rhs[] = {2, -1, 4};
    double lbs[] = {0.4, 0, 0, 0, 0};
    double ubs[] = {0.6, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error row size", CHECK_ROW_SIZES(A, constraints->state->row_sizes));
    mu_assert("error col size",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {5, -1, 3, 2, 6, 3, 1};
    int Ai_correct[] = {0, 1, 2, 3, 0, 2, 3};
    int Ap_correct[] = {0, 4, 7};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 7));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 7));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new coltags are correct
    ColTag col_tags_correct[] = {C_TAG_NONE, C_TAG_UB_INF, C_TAG_UB_INF,
                                 C_TAG_UB_INF};
    mu_assert("error col_tags",
              ARRAYS_EQUAL_COLTAG(col_tags_correct, constraints->col_tags, 4));

    // check that new variable bounds are correct
    double lbs_correct[] = {0.7, 0, 0, 0};
    double ubs_correct[] = {0.8, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {3, -1, 3, 2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == -2.0);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4, 5};
    double rhs_correct[] = {4, 5};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/*  dton row eliminated, unexpected cancellation.
    min. [2 -1 -1 3 2]x
    s.t. [2  1  -1  3  2]     [2]
         [1 -2   0  0  0] x = [-1]
         [1 -2   1  3  1]     [4]
         0.4 <= x1 <= 0.6
         x2, x3, x4, x5 >= 0
*/
static char *test_2_dton()
{
    double Ax[] = {2, 1, -1, 3, 2, 1, -2, 1, -2, 1, 3, 1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 2, 3, 4};
    int Ap[] = {0, 5, 7, 12};
    int nnz = 12;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -1, 4};
    double rhs[] = {2, -1, 4};
    double lbs[] = {0.4, 0, 0, 0, 0};
    double ubs[] = {0.6, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;

    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PS_FREE(stgs);
    free_presolver(presolver);

    /*
    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error", CHECK_ROW_SIZES(constraints->A,
                                       constraints->state->row_sizes));
    mu_assert("error", CHECK_COL_SIZES(constraints->AT,
                                       constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {5, -1, 3, 2, 1, 3, 1};
    int Ai_correct[] = {0, 1, 2, 3, 1, 2, 3};
    int Ap_correct[] = {0, 4, 7};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 7));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 7));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {0.7, 0, 0, 0};
    double ubs_correct[] = {0.8, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {3, -1, 3, 2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == -2.0);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4, 5};
    double rhs_correct[] = {4, 5};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    */
    // PS_FREE(stgs);
    // DEBUG(run_debugger(constraints, false));
    // free_presolver(presolver);
    return 0;
}

/*  dton row eliminated, unexpected cancellation resulting in zero column
    min. [2  -1 -1 3 2]x
    s.t. [2 -4  -1  3  2]     [2]
         [1 -2   0  0  0] x = [-1]
         [1 -2   1  3  1]     [4]
         1 <= x1 <= 2
         x3, x4, x5 >= 0
         x2 >= 1.1
*/
static char *test_3_dton()
{
    double Ax[] = {2, -4, -1, 3, 2, 1, -2, 1, -2, 1, 3, 1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 2, 3, 4};
    int Ap[] = {0, 5, 7, 12};
    int nnz = 12;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -1, 4};
    double rhs[] = {2, -1, 4};
    double lbs[] = {1, 1.1, 0, 0, 0};
    double ubs[] = {2, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {-1, 3, 2, 1, 3, 1};
    int Ai_correct[] = {1, 2, 3, 1, 2, 3};
    int Ap_correct[] = {0, 3, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {1.1, 0, 0, 0};
    double ubs_correct[] = {1.5, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {3, -1, 3, 2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error offset", prob->obj->offset == -2);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4, 5};
    double rhs_correct[] = {4, 5};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Two identical doubleton rows, feasible problem,
    min. [2  -1 -1 3 2]x
    s.t. [2 -5  -1  3  2]     [2]
         [1 -2   0  0  0] x = [-1]
         [1 -2   1  3  1]     [4]
         [1 -2   0  0  0]     [-1]
         1 <= x1 <= 2
         x3, x4, x5 >= 0
         x2 >= 1.1
*/
static char *test_004_dton()
{
    double Ax[] = {2, -5, -1, 3, 2, 1, -2, 1, -2, 1, 3, 1, 1, -2};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 2, 3, 4, 0, 1};
    int Ap[] = {0, 5, 7, 12, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, -1, 4, -1};
    double rhs[] = {2, -1, 4, -1};
    double lbs[] = {1, 1.1, 0, 0, 0};
    double ubs[] = {2, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {-1, -1, 3, 2, 1, 3, 1};
    int Ai_correct[] = {0, 1, 2, 3, 1, 2, 3};
    int Ap_correct[] = {0, 4, 7, 7};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 7));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 7));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {1.1, 0, 0, 0};
    double ubs_correct[] = {1.5, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {3, -1, 3, 2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error offset", prob->obj->offset == -2);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4, 5, 0};
    double rhs_correct[] = {4, 5, 0};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 3));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 3));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Two identical doubleton rows, feasible problem, empty column, empty row
    min. [2  -1 -1 3 2]x
    s.t. [2 -4  -1  3  2]     [2]
         [1 -2   0  0  0] x = [-1]
         [1 -2   1  3  1]     [4]
         [1 -2   0  0  0]     [-1]
         1 <= x1 <= 2
         x3, x4, x5 >= 0
         x2 >= 1.1
*/
static char *test_4_dton()
{
    double Ax[] = {2, -4, -1, 3, 2, 1, -2, 1, -2, 1, 3, 1, 1, -2};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 2, 3, 4, 0, 1};
    int Ap[] = {0, 5, 7, 12, 14};
    int nnz = 14;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, -1, 4, -1};
    double rhs[] = {2, -1, 4, -1};
    double lbs[] = {1, 1.1, 0, 0, 0};
    double ubs[] = {2, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {-1, 3, 2, 1, 3, 1};
    int Ai_correct[] = {1, 2, 3, 1, 2, 3};
    int Ap_correct[] = {0, 3, 6, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {1.1, 0, 0, 0};
    double ubs_correct[] = {1.5, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {3, -1, 3, 2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error obj", prob->obj->offset == -2);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4, 5, 0, 0};
    double rhs_correct[] = {4, 5, 0, 0};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 4));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 4));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*  Two identical doubleton rows, infeasible problem.
    min. [2  -1 -1 3 2]x
    s.t. [2 -4  -1  3  2]     [2]
         [1 -2   0  0  0] x = [-1]
         [1 -2   1  3  1]     [4]
         [1 -2   0  0  0]     [-1 + 1e-15]
         1 <= x1 <= 2
         x3, x4, x5 >= 0
         x2 >= 1.1
*/

static char *test_5_dton()
{
    double Ax[] = {2, 1, -1, 3, 2, 1, -2, 1, -2, 1, 3, 1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 2, 3, 4};
    int Ap[] = {0, 5, 7, 12};
    int nnz = 12;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -1, 4};
    double rhs[] = {2, -1, 4};
    double lbs[] = {0.4, 0, 0, 0, 0};
    double ubs[] = {0.6, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PS_FREE(stgs);
    // DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    mu_assert("error", 1 == 0);
    return 0;
}

/*  Two doubleton rows, presolver wants to substitute the same variable
         from both
    min. [2  -1 -1 3]x
    s.t. [1  0   2   0]     [2]
         [1  0   0   2] x = [2]
         [1  0   2  -4]     <= 4
         [0 -2   3   5]     <= 5
         [0  1  -6   7]     <= 7
         x >= 0
*/

static char *test_6_dton()
{
    double Ax[] = {1, 2, 1, 2, 1, 2, -4, -2, 3, 5, 1, -6, 7};
    int Ai[] = {0, 2, 0, 3, 0, 2, 3, 1, 2, 3, 1, 2, 3};
    int Ap[] = {0, 2, 4, 7, 10, 13};
    int nnz = 13;
    int n_rows = 5;
    int n_cols = 4;

    double lhs[] = {2, 2, -INF, -INF, -INF};
    double rhs[] = {2, 2, 4, 5, 7};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {-4, -2, 8, 1, 1};
    int Ai_correct[] = {1, 0, 1, 0, 1};
    int Ap_correct[] = {0, 1, 3, 5};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 5));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 5));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0};
    double ubs_correct[] = {INF, 1};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 2));

    // check that the objective function is correct
    double obj_correct[] = {-1, -2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 2));
    mu_assert("error offset", prob->obj->offset == 4);

    // check that lhs and rhs are correct
    double lhs_correct[] = {-INF, -INF, -INF};
    double rhs_correct[] = {2, 5, 7};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Two doubleton rows, eliminating one of them causes the other dton row
         to become a singleton row
    min. [2  -1 -1 3 2]x
    s.t. [2 -4  -1  3  2]     [2]
         [1 -2   0  0  0] x = [-1]
         [0 -2   1  3  1]     [4]
         [1 -1   0  0  0]     [1]
         0 <= x1 <= 5
         x2, x3, x4, x5 >= 0
*/
static char *test_7_dton()
{
    double Ax[] = {2, -4, -1, 3, 2, 1, -2, -2, 1, 3, 1, 1, -1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 1, 2, 3, 4, 0, 1};
    int Ap[] = {0, 5, 7, 11, 13};
    int nnz = 13;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, -1, 4, 1};
    double rhs[] = {2, -1, 4, 1};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {5, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {-1, 3, 2, -2, 1, 3, 1, 1};
    int Ai_correct[] = {1, 2, 3, 0, 1, 2, 3, 0};
    int Ap_correct[] = {0, 3, 7, 8};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 8));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 8));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {0.5, 0, 0, 0};
    double ubs_correct[] = {3, INF, INF, INF};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 4));

    // check that the objective function is correct
    double obj_correct[] = {3, -1, 3, 2};
    mu_assert("error obj", ARRAYS_EQUAL_DOUBLE(obj_correct, prob->obj->c, 4));
    mu_assert("error offset", prob->obj->offset == -2);

    // check that lhs and rhs are correct
    double lhs_correct[] = {4, 4, 2};
    double rhs_correct[] = {4, 4, 2};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 2));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 2));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* doubleton row with two fill-in in a column accepted
    min. [-2 -1 -1 -3]x
    s.t. [1  0   2   1]   =   [2]
         [0  1   0   1]   =   [3]
         [2  1   0   0]   <=  [5]
         [4  1   0   0]   <=  [6]
         [2  0   2   0] x <=  [7]
         [1  0   0   1]   <=  [9]
         [2  0   0   1]   <=  [10]
         [3  0   0   1]   <=  [13]
         [4  0   0   1]   <=  [14]
         x >= 0

*/
static char *test_8_dton()
{
    double Ax[] = {1, 2, 1, 1, 1, 2, 1, 4, 1, 2, 2, 1, 1, 2, 1, 3, 1, 4, 1};
    int Ai[] = {0, 2, 3, 1, 3, 0, 1, 0, 1, 0, 2, 0, 3, 0, 3, 0, 3, 0, 3};
    int Ap[] = {0, 3, 5, 7, 9, 11, 13, 15, 17, 19};
    int nnz = 19;
    int n_rows = 9;
    int n_cols = 4;

    double lhs[] = {2, 3, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {2, 3, 5, 6, 7, 9, 10, 13, 14};
    double lbs[4] = {0};
    double ubs[4] = {INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {1, 2, 1, 2, -1, 4, -1, 2, 2, 1, 1, 2, 1, 3, 1, 4, 1};
    int Ai_correct[] = {0, 1, 2, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2};
    int Ap_correct[] = {0, 3, 5, 7, 9, 11, 13, 15, 17};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 17));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 17));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0};
    double ubs_correct[] = {INF, INF, 3};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs_correct[] = {2, 2, 3, 7, 9, 10, 13, 14};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 8));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 8));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* doubleton row with three fill-in in a column accepted
    min. [-2 -1 -1 -3]x
    s.t. [1  0   2   1]   =   [2]
         [0  1   0   1]   =   [3]
         [2  1   0   0]   <=  [5]
         [4  1   0   0]   <=  [6]
         [2  1   2   0] x <=  [7]
         [1  0   0   1]   <=  [9]
         [2  0   0   1]   <=  [10]
         [3  0   0   1]   <=  [13]
         [4  0   0   1]   <=  [14]
         x >= 0
*/
static char *test_9_dton()
{
    double Ax[] = {1, 2, 1, 1, 1, 2, 1, 4, 1, 2, 1, 2, 1, 1, 2, 1, 3, 1, 4, 1};
    int Ai[] = {0, 2, 3, 1, 3, 0, 1, 0, 1, 0, 1, 2, 0, 3, 0, 3, 0, 3, 0, 3};
    int Ap[] = {0, 3, 5, 7, 9, 12, 14, 16, 18, 20};
    int nnz = 20;
    int n_rows = 9;
    int n_cols = 4;

    double lhs[] = {2, 3, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {2, 3, 5, 6, 7, 9, 10, 13, 14};
    double lbs[4] = {0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {-2, -1, -1, -3};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {1, 2, 1, 2, -1, 4, -1, 2, 2, -1, 1, 1, 2, 1, 3, 1, 4, 1};
    int Ai_correct[] = {0, 1, 2, 0, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 0, 2, 0, 2};
    int Ap_correct[] = {0, 3, 5, 7, 10, 12, 14, 16, 18};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 18));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 18));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, 0, 0};
    double ubs_correct[] = {INF, INF, 3};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs_correct[] = {2, 2, 3, 4, 9, 10, 13, 14};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 8));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 8));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* doubleton row with three fill-in in a column accepted
    min. [-2 -1 -1 -3]x
    s.t. [1  0   2   1]   =   [2]
         [0  1   0   1]   =   [3]
         [2  1   0   0]   <=  [5]
         [4  1   0   0]   <=  [6]
         [2  1   2   0] x <=  [7]
         [1  0   0   1]   <=  [9]
         [2  0   0   1]   <=  [10]
         [3  0   0   1]   <=  [13]
         [4  0   0   1]   <=  [14]
         x1, x2, x4 >= 0
         x3 free

*/
static char *test_10_dton()
{
    double Ax[] = {1, 2, 1, 1, 1, 2, 1, 4, 1, 2, 1, 2, 1, 1, 2, 1, 3, 1, 4, 1};
    int Ai[] = {0, 2, 3, 1, 3, 0, 1, 0, 1, 0, 1, 2, 0, 3, 0, 3, 0, 3, 0, 3};
    int Ap[] = {0, 3, 5, 7, 9, 12, 14, 16, 18, 20};
    int nnz = 20;
    int n_rows = 9;
    int n_cols = 4;

    double lhs[] = {2, 3, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {2, 3, 5, 6, 7, 9, 10, 13, 14};
    double lbs[] = {0, 0, -INF, 0};
    double ubs[] = {INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 10);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {1, 2, 1, 2, -1, 4, -1, 2, 2, -1, 1, 1, 2, 1, 3, 1, 4, 1};
    int Ai_correct[] = {0, 1, 2, 0, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 0, 2, 0, 2};
    int Ap_correct[] = {0, 3, 5, 7, 10, 12, 14, 16, 18};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 18));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 18));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    // check that new variable bounds are correct
    double lbs_correct[] = {0, -INF, 0};
    double ubs_correct[] = {INF, INF, 3};
    mu_assert("error bounds",
              CHECK_BOUNDS(constraints->bounds, lbs_correct, ubs_correct, 3));

    // check that lhs and rhs are correct
    double lhs_correct[] = {2, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs_correct[] = {2, 2, 3, 4, 9, 10, 13, 14};
    mu_assert("error lhs", ARRAYS_EQUAL_DOUBLE(lhs_correct, constraints->lhs, 8));
    mu_assert("error rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, constraints->rhs, 8));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*
    doubleton row with four fill-in rejected
    min. [-2 -1 -1 -3]x
    s.t. [1  0   2   1]   =   [2]
         [0  1   0   2]   =   [3]
         [2  1   0   0]   <=  [5]
         [4  1   0   0]   <=  [6]
         [2  1   0   0] x <=  [7]
         [1  1   0   0]   <=  [9]
         [2  0   0   1]   <=  [10]
         [3  0   0   1]   <=  [13]
         [4  0   0   1]   <=  [14]
         x1, x2, x4 >= 0
         x3 free
*/
static char *test_11_dton()
{
    double Ax[] = {1, 2, 1, 1, 2, 2, 1, 4, 1, 2, 1, 1, 1, 2, 1, 3, 1, 4, 1};
    int Ai[] = {0, 2, 3, 1, 3, 0, 1, 0, 1, 0, 1, 0, 1, 0, 3, 0, 3, 0, 3};
    int Ap[] = {0, 3, 5, 7, 9, 11, 13, 15, 17, 19};
    int nnz = 19;
    int n_rows = 9;
    int n_cols = 4;

    double lhs[] = {2, 3, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {2, 3, 5, 6, 7, 9, 10, 13, 14};
    double lbs[] = {0, 0, -INF, 0};
    double ubs[] = {INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    PresolveStatus status = remove_dton_eq_rows(prob, 0);
    // mu_assert("error status", status == UNCHANGED);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // print_matrix(A);

    // check that new A is correct
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax, A->x, nnz));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai, A->i, nnz));
    check_row_starts(A, Ap);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/*
    doubleton row with four fill-in accepted (one extra space is released
        by another reduction)
    min. [-2 -1 -1 -3]x
    s.t. [1  0   2   1]   =   [2]
         [0  1   0   2]   =   [3]
         [2  1   0   0]   <=  [5]
         [4  1   0   0]   <=  [6]
         [2  1   0   0] x <=  [7]
         [1  1   0   0]   <=  [9]
         [2  0   0   1]   <=  [10]
         [3  0   0   1]   <=  [13]
         [4  0   0   1]   <=  [14]
         x1, x2, x4 >= 0
         x3 free
*/
static char *test_12_dton()
{
    double Ax[] = {2, 1, -1, 3, 2, 1, -2, 1, -2, 1, 3, 1};
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 0, 1, 2, 3, 4};
    int Ap[] = {0, 5, 7, 12};
    int nnz = 12;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, -1, 4};
    double rhs[] = {2, -1, 4};
    double lbs[] = {0.4, 0, 0, 0, 0};
    double ubs[] = {0.6, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    PS_FREE(stgs);
    // DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    mu_assert("error", 1 == 1);
    return 0;
}

/*
     corner case that caused a nasty bug on a NETLIB problem
     (excessive fill-in in one column due to two doubleton rows
     with one common variable). We set max shift to 0, but both
     reductions should be allowed since we can eat up the space
     in column 1 without actually shifting any elements (since
     column 1 is inactive after the first dton reduction.)
     min. [-2 -1 -1 -3, -1]x
     s.t. [2,  -1,    0,   0,  0]     = 0
          [2,   0,   -1,   0,  0]     = 0
          [0,   1,    0,   1,  2]     <= 3
          [0,   1,    0,   3,  4]     <= 4
          [0,   0,    2,   5,  6] x   <= 5
          [0,   0,    2,   7,  8]     <= 6
          [0,   0,    1,   7,  8]     <= 7
          [2,   0,    0,   7,  8]     <= 8
          [3,   0,    0,   7,  8]     <= 9
          [4,   0,    0,   7,  8]     <= 10
          [5,   0,    0,   7,  8]     <= 11
          [6,   0,    0,   7,  8]     <= 12
          x1, x2, x3, x4 >= 0
*/
static char *test_13_dton()
{
    double Ax[] = {2, -1, 2, -1, 1, 1, 2, 1, 3, 4, 2, 5, 6, 2, 7, 8, 1,
                   7, 8,  2, 7,  8, 3, 7, 8, 4, 7, 8, 5, 7, 8, 6, 7, 8};
    int Ai[] = {0, 1, 0, 2, 1, 3, 4, 1, 3, 4, 2, 3, 4, 2, 3, 4, 2,
                3, 4, 0, 3, 4, 0, 3, 4, 0, 3, 4, 0, 3, 4, 0, 3, 4};
    int Ap[] = {0, 2, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34};
    int nnz = 34;
    int n_rows = 12;
    int n_cols = 5;

    double lhs[] = {0,    0,    -INF, -INF, -INF, -INF,
                    -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {0, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {-2, -1, -1, 3, -1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 0);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, 1, 2, 2, 3, 4, 4, 5, 6, 4, 7, 8, 2, 7, 8,
                           2, 7, 8, 3, 7, 8, 4, 7, 8, 5, 7, 8, 6, 7, 8};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
                        0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
    int Ap_correct[] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 30));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 30));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    PS_FREE(stgs);
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* Same as the previous test, but we change columns 0 and 1. We set
   max shift to 0. The second reduction is then rejected, because although there
   is free space on the left we can't eat it up. Very interesting.
   Is this how we want it to be? */
static char *test_14_dton()
{
    double Ax[] = {-1, 2, 2, -1, 1, 1, 2, 1, 3, 4, 2, 5, 6, 2, 7, 8, 1,
                   7,  8, 2, 7,  8, 3, 7, 8, 4, 7, 8, 5, 7, 8, 6, 7, 8};
    int Ai[] = {0, 1, 1, 2, 0, 3, 4, 0, 3, 4, 2, 3, 4, 2, 3, 4, 2,
                3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4};
    int Ap[] = {0, 2, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34};
    int nnz = 34;
    int n_rows = 12;
    int n_cols = 5;

    double lhs[] = {0,    0,    -INF, -INF, -INF, -INF,
                    -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {0, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {-2, -1, -1, 3, -1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 0);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, -1, 2, 1, 2, 2, 3, 4, 2, 5, 6, 2, 7, 8, 1, 7,
                           8, 2,  7, 8, 3, 7, 8, 4, 7, 8, 5, 7, 8, 6, 7, 8};
    int Ai_correct[] = {0, 1, 0, 2, 3, 0, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
                        3, 0, 2, 3, 0, 2, 3, 0, 2, 3, 0, 2, 3, 0, 2, 3};
    int Ap_correct[] = {0, 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 32));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 32));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);

    return 0;
}

/* Same as previous test, but with max shift set to 4.
   The reduction should be accepted because shifting the
   right row requires shifting four elements.  */
static char *test_15_dton()
{
    double Ax[] = {-1, 2, 2, -1, 1, 1, 2, 1, 3, 4, 2, 5, 6, 2, 7, 8, 1,
                   7,  8, 2, 7,  8, 3, 7, 8, 4, 7, 8, 5, 7, 8, 6, 7, 8};
    int Ai[] = {0, 1, 1, 2, 0, 3, 4, 0, 3, 4, 2, 3, 4, 2, 3, 4, 2,
                3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4};
    int Ap[] = {0, 2, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34};
    int nnz = 34;
    int n_rows = 12;
    int n_cols = 5;

    double lhs[] = {0,    0,    -INF, -INF, -INF, -INF,
                    -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {0, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {-2, -1, -1, 3, -1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 4);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, 1, 2, 2, 3, 4, 4, 5, 6, 4, 7, 8, 2, 7, 8,
                           2, 7, 8, 3, 7, 8, 4, 7, 8, 5, 7, 8, 6, 7, 8};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
                        0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
    int Ap_correct[] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 30));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 30));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    // check new AT
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* Same as previous test, but with max shift set to 3.
   The reduction should be rejected because shifting the
   right row requires shifting four elements.  */
static char *test_16_dton()
{
    double Ax[] = {-1, 2, 2, -1, 1, 1, 2, 1, 3, 4, 2, 5, 6, 2, 7, 8, 1,
                   7,  8, 2, 7,  8, 3, 7, 8, 4, 7, 8, 5, 7, 8, 6, 7, 8};
    int Ai[] = {0, 1, 1, 2, 0, 3, 4, 0, 3, 4, 2, 3, 4, 2, 3, 4, 2,
                3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4};
    int Ap[] = {0, 2, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34};
    int nnz = 34;
    int n_rows = 12;
    int n_cols = 5;

    double lhs[] = {0,    0,    -INF, -INF, -INF, -INF,
                    -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {0, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {INF, INF, INF, INF, INF};
    double c[] = {-2, -1, -1, 3, -1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    remove_dton_eq_rows(prob, 3);
    problem_clean(prob, true);

    mu_assert("error row_sizes",
              CHECK_ROW_SIZES(constraints->A, constraints->state->row_sizes));
    mu_assert("error col_sizes",
              CHECK_COL_SIZES(constraints->AT, constraints->state->col_sizes));

    // check that new A is correct
    double Ax_correct[] = {2, -1, 2, 1, 2, 2, 3, 4, 2, 5, 6, 2, 7, 8, 1, 7,
                           8, 2,  7, 8, 3, 7, 8, 4, 7, 8, 5, 7, 8, 6, 7, 8};
    int Ai_correct[] = {0, 1, 0, 2, 3, 0, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
                        3, 0, 2, 3, 0, 2, 3, 0, 2, 3, 0, 2, 3, 0, 2, 3};
    int Ap_correct[] = {0, 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 32));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 32));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    PS_FREE(stgs); // check new AT
    DEBUG(run_debugger(constraints, false));
    DEBUG(mu_assert("error AT", verify_A_and_AT_consistency(A, constraints->AT)));
    free_presolver(presolver);

    return 0;
}

/* We eliminate a doubleton row and substitute a variable in a row with
   no extra space. We get in-place fill-in in A, but that should be
   fine.

    min. [2  -1 -1 3 2]x
    s.t. [2  1   0  0  0]      [2]
         [0  1   1  3  0] x >= [5]
         [0 -2   1  3  1]      [4]
         [1 -1   0  1  1]      [1]
         0 <= x1 <= 5
         x2, x3, x4, x5 >= 0
*/
static char *test_17_dton()
{
    double Ax[] = {2, 1, 1, 1, 3, -2, 1, 3, 1, 1, -1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 3, 1, 2, 3, 4, 0, 1, 3, 4};
    int Ap[] = {0, 2, 5, 9, 13};
    int nnz = 13;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, 5, 4, 1};
    double rhs[] = {2, INF, 4, 1};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {5, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    State *data = constraints->state;
    remove_extra_space(A, data->row_sizes, data->col_sizes, true,
                       data->work->mappings->cols);

    remove_dton_eq_rows(prob, 0);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-2, 1, 3, 4, 1, 3, 1, 3, 1, 1};
    int Ai_correct[] = {0, 1, 2, 0, 1, 2, 3, 0, 2, 3};
    int Ap_correct[] = {0, 3, 7, 10};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 10));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 10));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* We eliminate a doubleton row without having any extra space in
   AT (and A). We get in-place fill-in in AT and A, but that
   should be fine.

    min. [2  -1 -1 3 2]x
    s.t. [2  1   0  0  0]      [2]
         [0  1   1  3  0] x >= [5]
         [1 -1   0  1  1]      [1]
         0 <= x1 <= 5
         x2, x3, x4, x5 >= 0
*/
static char *test_18_dton()
{
    double Ax[] = {2, 1, 1, 1, 3, 1, -1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 3, 0, 1, 3, 4};
    int Ap[] = {0, 2, 5, 9};
    int nnz = 9;
    int n_rows = 3;
    int n_cols = 5;

    double lhs[] = {2, 5, 1};
    double rhs[] = {2, INF, 1};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {5, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    Matrix *AT = constraints->AT;
    State *data = constraints->state;
    remove_extra_space(A, data->row_sizes, data->col_sizes, true,
                       data->work->mappings->cols);
    remove_extra_space(AT, data->col_sizes, data->row_sizes, true,
                       data->work->mappings->rows);

    remove_dton_eq_rows(prob, 0);
    problem_clean(prob, true);

    // check that new A is correct
    double Ax_correct[] = {-2, 1, 3, 3, 1, 1};
    int Ai_correct[] = {0, 1, 2, 0, 2, 3};
    int Ap_correct[] = {0, 3, 6};
    mu_assert("error Ax", ARRAYS_EQUAL_DOUBLE(Ax_correct, A->x, 6));
    mu_assert("error Ai", ARRAYS_EQUAL_INT(Ai_correct, A->i, 6));
    mu_assert("error row starts", check_row_starts(A, Ap_correct));

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

/* We eliminate a doubleton row and substitute a variable in a row with
   no extra space. We get in-place fill-in, but that should be fine.

   Same as above, but infeasible constraint. But this should affect
   anything, but we still get a weird bug! What happened was that
   we got a chain of doubleton rows, so the entire A matrix was
   eliminated. We don't run this test for now since it might not be
   a bug.

    min. [2  -1 -1 3 2]x
    s.t. [2  1   0  0  0]       [2]
         [0  1   1  3  0] x   = [-4]
         [0 -2   1  3  1]       [4]
         [1 -1   0  1  0]       [1]
         0 <= x1 <= 5
         x2, x3, x4, x5 >= 0
*/
static char *test_19_dton()
{
    double Ax[] = {2, 1, 1, 1, 3, -2, 1, 3, 1, 1, -1, 1};
    int Ai[] = {0, 1, 1, 2, 3, 1, 2, 3, 4, 0, 1, 3};
    int Ap[] = {0, 2, 5, 9, 12};
    int nnz = 12;
    int n_rows = 4;
    int n_cols = 5;

    double lhs[] = {2, -4, 4, 1};
    double rhs[] = {2, -4, 4, 1};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {5, INF, INF, INF, INF};
    double c[] = {2, -1, -1, 3, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;
    State *data = constraints->state;
    remove_extra_space(A, data->row_sizes, data->col_sizes, true,
                       data->work->mappings->cols);
    remove_dton_eq_rows(prob, 0);
    problem_clean(prob, true);

    PS_FREE(stgs);
    DEBUG(run_debugger(constraints, false));
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_dton()
{
    mu_run_test(test_00_dton, counter_dton);  // implemented
    mu_run_test(test_01_dton, counter_dton);  // implemented
    mu_run_test(test_02_dton, counter_dton);  // implemented
    mu_run_test(test_03_dton, counter_dton);  // implemented
    mu_run_test(test_04_dton, counter_dton);  // implemented
    mu_run_test(test_05_dton, counter_dton);  // implemented
    mu_run_test(test_06_dton, counter_dton);  // implemented
    mu_run_test(test_1_dton, counter_dton);   // implemented
    mu_run_test(test_2_dton, counter_dton);   // implemented
    mu_run_test(test_3_dton, counter_dton);   // implemented
    mu_run_test(test_004_dton, counter_dton); // implemented
    mu_run_test(test_4_dton, counter_dton);   // implemented
    // mu_run_test(test_5_dton, counter_dton);
    mu_run_test(test_6_dton, counter_dton);  // implemented
    mu_run_test(test_7_dton, counter_dton);  // implemented
    mu_run_test(test_8_dton, counter_dton);  // implemented
    mu_run_test(test_9_dton, counter_dton);  // implemented
    mu_run_test(test_10_dton, counter_dton); // implemented
    mu_run_test(test_11_dton, counter_dton); // implemented
    //   mu_run_test(test_12_dton, counter_dton);
    mu_run_test(test_13_dton, counter_dton); // implemented
    mu_run_test(test_14_dton, counter_dton); // implemented
    printf("before test 15\n");
    mu_run_test(test_15_dton, counter_dton); // implemented
    printf("after test 15\n");
    mu_run_test(test_16_dton, counter_dton); // implemented
    mu_run_test(test_17_dton, counter_dton); // implemented
    mu_run_test(test_18_dton, counter_dton); // implemented
    //    mu_run_test(test_19_dton, counter_dton); // implemented but we don't
    //    run it
    return 0;
}

int test_dton()
{
    const char *result = all_tests_dton();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("dton: TEST FAILED!\n");
    }
    else
    {
        printf("dton: ALL TESTS PASSED\n");
    }
    printf("dton: Tests run: %d\n", counter_dton);
    return result == 0;
}

#endif // TEST_DTON_H
