#ifndef TEST_QP_H
#define TEST_QP_H

#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "DtonsEq_internal.h"
#include "Matrix.h"
#include "Numerics.h"
#include "PSLP_API.h"
#include "PSLP_sol.h"
#include "PSLP_stats.h"
#include "Postsolver.h"
#include "Problem.h"
#include "SimpleReductions.h"
#include "State.h"
#include "Workspace.h"
#include "dVec.h"
#include "iVec.h"
#include "kkt.h"
#include "minunit.h"
#include "test_macros.h"
#include "u16Vec.h"
#include <stdio.h>

/* Tests for the quadratic objective term (1/2) x'Px. A variable that
   appears in P can be fixed (its row and column of P are folded away). It is
   substituted through a doubleton equality row only into another variable in
   P and only when its row of P holds nothing but a positive diagonal and the
   other variable of the row (no fill-in in P); otherwise the other column is
   substituted or the row stays.
   Test 1: structural validation of P in new_presolver_from_problem
   Test 2: fix_col folds the coupling terms of P into c and removes the
           row/column of P
   Test 3: singleton equality row fixing a variable that appears in P
           (postsolve recovers y from the gradient)
   Test 4: empty columns: fixed by the sign of c when absent from P, at
           clip(-c/Pkk) with only a diagonal, kept with off-diagonal entries
   Test 5: simple dual fix skips a variable that appears in P
   Test 6: doubleton equality row substitutes the variable that is absent
           from P
   Test 7: parallel columns are not merged when one of them appears in P
   Test 8: dual infeasibility ray after fixing a variable that appears in P
   Test 9: primal propagation fixing coupled variables (bound-change replay)
   Test 10: doubleton equality row where the variable in P is a column
            singleton (the other one is substituted and the singleton grows)
   Test 11: diagonal dual fix fixes variables whose minimiser is beyond the
            bound on the unlocked side
   Test 15: two coupled variables fixed by two singleton equality rows, in
            both orders (the first multiplier needs the later-fixed neighbour)
   Test 16: primal infeasibility ray after fixing a variable that appears in P
   Test 17: diagonal-only empty columns with infinite bounds: clipped to a
            finite bound, inside with both bounds infinite, and c = 0
   Test 18: doubleton equality row with both variables diagonal in P
            substitutes the sparser one (P_tt grows, postsolve recovers y)
   Test 19: dton_choose_subst_qp: the column absent from P, else the one
            that can be substituted without fill-in, else none
   Test 20: doubleton equality row with both variables coupled to a third
            one in P stays
   Test 21: the diagonal-only variable is substituted with the smaller
            coefficient and other rows (ADDED_ROWS enters its multiplier)
   Test 22: chain of two substitutions of diagonal-only variables in one
            round (P_tt accumulates, multipliers replay deepest first)
   Test 23: coupled dual fix of a column that is empty in A (CVXQP-like): the
            derivative is nonnegative over the box, fixed at lb; a coupled
            column with an undetermined sign stays
   Test 24: coupled dual fix at ub of a column with down-locks only
   Test 25: doubleton equality row with both variables coupled to each other
            only: the sparser one is substituted, P_jj, c_j and the offset
            change, the multiplier needs P_kj x_j (DELETED_ROW_ZK)
   Test 26: the same when the stay column keeps a coupling to a third variable
   Test 29: lp_fold leaves a QP unfolded; P = 0 through new_presolver_from_problem
            folds as an LP
*/

#define QP_TOL 1e-8
static int counter_qp = 0;

/* new_presolver_from_problem with the problem data as arguments */
static Presolver *qp_presolver(const double *Ax, const int *Ai, const int *Ap,
                               size_t m, size_t n, size_t nnz, const double *lhs,
                               const double *rhs, const double *lbs,
                               const double *ubs, const double *c, const double *Px,
                               const int *Pi, const int *Pp, size_t Pnnz,
                               const Settings *stgs)
{
    ProblemData problem = {.Ax = Ax,
                           .Ai = Ai,
                           .Ap = Ap,
                           .m = m,
                           .n = n,
                           .nnz = nnz,
                           .lhs = lhs,
                           .rhs = rhs,
                           .c = c,
                           .lbs = lbs,
                           .ubs = ubs,
                           .Px = Px,
                           .Pi = Pi,
                           .Pp = Pp,
                           .Pnnz = Pnnz};
    return new_presolver_from_problem(&problem, stgs);
}

/*  Input validation. The constraint part is x0 + x1 <= 1, -1 <= x <= 1. */
static char *test_1_qp()
{
    double Ax[] = {1, 1};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    double lhs[] = {-INF};
    double rhs[] = {1};
    double lbs[] = {-1, -1};
    double ubs[] = {1, 1};
    double c[] = {1, 1};
    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver;

    // asymmetric values
    {
        double Px[] = {1, 0.5, 0.7, 1};
        int Pi[] = {0, 1, 0, 1};
        int Pp[] = {0, 2, 4};
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp, 4, stgs);
        mu_assert("asymmetric P should be rejected", presolver == NULL);
    }

    // only one triangle given
    {
        double Px[] = {1, 0.5, 1};
        int Pi[] = {0, 1, 1};
        int Pp[] = {0, 2, 3};
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp, 3, stgs);
        mu_assert("upper triangular P should be rejected", presolver == NULL);
    }

    // unsorted row
    {
        double Px[] = {1, 2, 1, 2};
        int Pi[] = {1, 0, 0, 1};
        int Pp[] = {0, 2, 4};
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp, 4, stgs);
        mu_assert("unsorted P row should be rejected", presolver == NULL);
    }

    // duplicate column index in a row
    {
        double Px[] = {1, 1, 1};
        int Pi[] = {0, 0, 1};
        int Pp[] = {0, 2, 3};
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp, 3, stgs);
        mu_assert("duplicate P entry should be rejected", presolver == NULL);
    }

    // column index out of range
    {
        double Px[] = {1, 1};
        int Pi[] = {0, 2};
        int Pp[] = {0, 1, 2};
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp, 2, stgs);
        mu_assert("out of range P column should be rejected", presolver == NULL);
    }

    // row pointers: Pp[0] != 0, non-decreasing, Pp[n] == Pnnz
    {
        double Px[] = {1, 1, 1};
        int Pi[] = {0, 1, 1};
        int Pp_first[] = {1, 2, 3};
        int Pp_decreasing[] = {0, 2, 1};
        int Pp_last[] = {0, 1, 2};
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp_first, 3, stgs);
        mu_assert("Pp[0] != 0 should be rejected", presolver == NULL);
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp_decreasing, 1, stgs);
        mu_assert("decreasing Pp should be rejected", presolver == NULL);
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp_last, 3, stgs);
        mu_assert("Pp[n] != Pnnz should be rejected", presolver == NULL);
    }

    // a non-empty P needs all three arrays
    {
        double Px[] = {1, 1};
        int Pi[] = {0, 1};
        int Pp[] = {0, 1, 2};
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px,
                                 NULL, Pp, 2, stgs);
        mu_assert("missing Pi should be rejected", presolver == NULL);
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 NULL, 2, stgs);
        mu_assert("missing Pp should be rejected", presolver == NULL);
    }

    // valid P, and explicit zeros are dropped
    {
        double Px[] = {2, 0, 0, 2};
        int Pi[] = {0, 1, 0, 1};
        int Pp[] = {0, 2, 4};
        presolver = qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, Px, Pi,
                                 Pp, 4, stgs);
        mu_assert("valid P should be accepted", presolver != NULL);
        mu_assert("explicit zeros in P should be dropped",
                  presolver->prob->obj->P->nnz == 2);
        mu_assert("P row sizes",
                  presolver->prob->constraints->state->P_row_sizes[0] == 1 &&
                      presolver->prob->constraints->state->P_row_sizes[1] == 1);
        mu_assert("Pnnz_original counts the entries after dropping zeros",
                  presolver->stats->Pnnz_original == 2);
        free_presolver(presolver);
    }

    // no P at all is an LP
    {
        presolver = new_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, stgs);
        mu_assert("empty P should be accepted", presolver != NULL);
        mu_assert("an LP has no P",
                  presolver->prob->obj->P == NULL &&
                      presolver->prob->constraints->state->P_row_sizes == NULL);
        mu_assert("empty P stats", presolver->stats->Pnnz_original == 0);
        free_presolver(presolver);
    }

    PS_FREE(stgs);
    return 0;
}

/*  Direct fix_col with a coupling term.
    min. (1/2) x'Px + [1 -1 0.5] x,  P = [2 1 0; 1 2 0; 0 0 4]
    s.t. x0 + x1 + x2 <= 10, 0 <= x <= 5
    Fixing x0 = 2 gives offset = 1*2 + (1/2)*2*4 = 6, c1 = -1 + 1*2 = 1,
    and P = [2 0; 0 4] on (x1, x2). */
static char *test_2_qp()
{
    double Ax[] = {1, 1, 1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    double lhs[] = {-INF};
    double rhs[] = {10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {5, 5, 5};
    double c[] = {1, -1, 0.5};
    double Px[] = {2, 1, 1, 2, 4};
    int Pi[] = {0, 1, 0, 1, 2};
    int Pp[] = {0, 2, 4, 5};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, 1, 3, 3, lhs, rhs, lbs, ubs, c,
                                        Px, Pi, Pp, 5, stgs);
    mu_assert("init", presolver != NULL);
    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    State *state = constraints->state;
    Matrix *P = prob->obj->P;
    DEBUG(run_debugger_problem(prob, false));

    mu_assert("status", fix_col(prob, 0, 2.0) == REDUCED);

    // objective is updated eagerly
    mu_assert("offset", prob->obj->offset == 6.0);
    mu_assert("c after fold", prob->obj->c[1] == 1.0 && prob->obj->c[2] == 0.5);

    // row and column 0 of P are gone
    mu_assert("P nnz", P->nnz == 2);
    mu_assert("P row sizes", state->P_row_sizes[0] == SIZE_INACTIVE_COL &&
                                 state->P_row_sizes[1] == 1 &&
                                 state->P_row_sizes[2] == 1);
    mu_assert("P row 0 empty", P->p[0].start == P->p[0].end);
    mu_assert("P row 1", P->p[1].end - P->p[1].start == 1 &&
                             P->i[P->p[1].start] == 1 && P->x[P->p[1].start] == 2.0);

    // the postsolve record carries the folded c0 and the P row (2 entries)
    PostsolveInfo *info = state->postsolve_info;
    mu_assert("one record",
              info->type->len == 1 && info->type->data[0] == FIXED_COL);
    FixedColRecord rec = decode_fixed_col(info->indices->data, info->vals->data,
                                          info->starts->data[1]);
    mu_assert("record col/val/ck", rec.col == 0 && rec.val == 2.0 && rec.ck == 1.0);
    mu_assert("record a", rec.len == 1 && rec.rows[0] == 0 && rec.vals[0] == 1);
    mu_assert("record P", rec.P_len == 2 && rec.P_cols[0] == 0 &&
                              rec.P_vals[0] == 2 && rec.P_cols[1] == 1 &&
                              rec.P_vals[1] == 1);

    delete_fixed_cols_from_problem(prob);
    delete_inactive_cols_from_A_and_AT(constraints);
    mu_assert("rhs", constraints->rhs[0] == 8.0);
    DEBUG(run_debugger_problem(prob, false));

    problem_clean(prob);
    DEBUG(run_debugger_problem(prob, true));
    mu_assert("P dims after clean", P->m == 2 && P->n == 2 && P->nnz == 2);
    double Px_correct[] = {2, 4};
    int Pi_correct[] = {0, 1};
    int Pp_correct[] = {0, 1, 2};
    mu_assert("Px after clean", ARRAYS_EQUAL_DOUBLE(Px_correct, P->x, 2));
    mu_assert("Pi after clean", ARRAYS_EQUAL_INT(Pi_correct, P->i, 2));
    mu_assert("Pp after clean", check_row_starts(P, Pp_correct));
    mu_assert("c after clean", prob->obj->c[0] == 1.0 && prob->obj->c[1] == 0.5);
    mu_assert("P row sizes after clean",
              state->P_row_sizes[0] == 1 && state->P_row_sizes[1] == 1);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Singleton equality row fixing a variable that appears in P.
    min. (1/2) x'Px + [1 -3 1] x,  P = [2 1 0; 1 2 0; 0 0 0]
    s.t. x0 = 2,  x1 + x2 <= 3,  0 <= x <= 10
    Reduced (x1, x2): min x1^2 - x1 + x2 s.t. x1 + x2 <= 3, offset 6.
    Reduced solution x = (0.5, 0), y = 0, z = (0, 1). */
static char *test_3_qp()
{
    double Ax[] = {1, 1, 1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 1, 3};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {2, -INF};
    double rhs[] = {2, 3};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {1, -3, 1};
    double Px[] = {2, 1, 1, 2};
    int Pi[] = {0, 1, 0, 1};
    int Pp[] = {0, 2, 4, 4};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 3, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 4, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 1 && red->n == 2 && red->nnz == 2);
    mu_assert("reduced offset", red->obj_offset == 6.0);
    mu_assert("reduced c", red->c[0] == -1.0 && red->c[1] == 1.0);
    mu_assert("reduced P nnz", red->Pnnz == 1);
    mu_assert("reduced P", red->Pp[0] == 0 && red->Pp[1] == 1 && red->Pp[2] == 1 &&
                               red->Pi[0] == 0 && red->Px[0] == 2.0);
    mu_assert("stats", presolver->stats->Pnnz_original == 4 &&
                           presolver->stats->Pnnz_reduced == 1);

    double x[] = {0.5, 0};
    double y[] = {0};
    double z[] = {0, 1};
    postsolve(presolver, x, y, z);

    double correct_x[] = {2, 0.5, 0};
    double correct_y[] = {5.5, 0};
    double correct_z[] = {0, 0, 1};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Empty columns with and without a quadratic term.
    min. (1/2) x'Px + [-4 1 0 0 0] x,
         P = 2 I on (x0, x2, x3), plus P_22 = P_44 = 2, P_24 = P_42 = 1
    s.t. 2 x2 + x3 = 4, -1 <= x0 <= 1, 0 <= x1 <= 5, -10 <= x2, x3, x4 <= 10
    x0, x1 and x4 are empty in A. x1 is not in P and goes to its lower bound
    (LP logic). x0 has only a diagonal, so x0 = clip(-c0 / P00, -1, 1) = 1.
    x4 is coupled to x2 and stays. Reduced (x2, x3, x4):
    min x2^2 + x2 x4 + x4^2 + x3^2 - 3 s.t. 2 x2 + x3 = 4,
    solution (32, 12, -16) / 19 with y = 24 / 19. */
static char *test_4_qp()
{
    double Ax[] = {2, 1};
    int Ai[] = {2, 3};
    int Ap[] = {0, 2};
    int n_rows = 1, n_cols = 5;
    double lhs[] = {4};
    double rhs[] = {4};
    double lbs[] = {-1, 0, -10, -10, -10};
    double ubs[] = {1, 5, 10, 10, 10};
    double c[] = {-4, 1, 0, 0, 0};
    double Px[] = {2, 2, 1, 2, 1, 2};
    int Pi[] = {0, 2, 4, 3, 2, 4};
    int Pp[] = {0, 1, 1, 3, 4, 6};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 2, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 6, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 1 && red->n == 3 && red->nnz == 2);
    mu_assert("reduced offset", red->obj_offset == -3.0);
    mu_assert("reduced P nnz", red->Pnnz == 5);
    double Px_correct[] = {2, 1, 2, 1, 2};
    int Pi_correct[] = {0, 2, 1, 0, 2};
    int Pp_correct[] = {0, 2, 3, 5};
    mu_assert("reduced Px", ARRAYS_EQUAL_DOUBLE(Px_correct, red->Px, 5));
    mu_assert("reduced Pi", ARRAYS_EQUAL_INT(Pi_correct, red->Pi, 5));
    mu_assert("reduced Pp", ARRAYS_EQUAL_INT(Pp_correct, red->Pp, 4));

    // the kept column (x4, now column 2) is still listed as empty
    State *state = presolver->prob->constraints->state;
    mu_assert("kept empty col", state->empty_cols->len == 1 &&
                                    state->empty_cols->data[0] == 2 &&
                                    !has_removable_empty_cols(state));

    double x[] = {32.0 / 19, 12.0 / 19, -16.0 / 19};
    double y[] = {24.0 / 19};
    double z[] = {0, 0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {1, 0, 32.0 / 19, 12.0 / 19, -16.0 / 19};
    double correct_y[] = {24.0 / 19};
    double correct_z[] = {-2, 1, 0, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Simple dual fix must skip a variable that appears in P.
    min. x0^2 + x0 x2 + x2^2 + x0 + x1
    s.t. x0 + x1 <= 5, x1 + x2 <= 5, -10 <= x <= 10
    For an LP all three variables would be dual-fixed. Here only x1 (not in
    P) is fixed to its lower bound; the rows become redundant singletons and
    x0, x2 empty columns that appear in P, so they stay. Reduced:
    min x0^2 + x0 x2 + x2^2 + x0 - 10, solution (-2/3, 1/3). */
static char *test_5_qp()
{
    double Ax[] = {1, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2};
    int Ap[] = {0, 2, 4};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {-INF, -INF};
    double rhs[] = {5, 5};
    double lbs[] = {-10, -10, -10};
    double ubs[] = {10, 10, 10};
    double c[] = {1, 1, 0};
    double Px[] = {2, 1, 1, 2};
    int Pi[] = {0, 2, 0, 2};
    int Pp[] = {0, 2, 2, 4};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 4, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 4, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 0 && red->n == 2 && red->nnz == 0);
    mu_assert("reduced offset", red->obj_offset == -10.0);
    mu_assert("reduced c", red->c[0] == 1.0 && red->c[1] == 0.0);
    double Px_correct[] = {2, 1, 1, 2};
    int Pi_correct[] = {0, 1, 0, 1};
    mu_assert("reduced P", red->Pnnz == 4 &&
                               ARRAYS_EQUAL_DOUBLE(Px_correct, red->Px, 4) &&
                               ARRAYS_EQUAL_INT(Pi_correct, red->Pi, 4));
    mu_assert("P unchanged",
              presolver->stats->Pnnz_original == presolver->stats->Pnnz_reduced);

    double x[] = {-2.0 / 3, 1.0 / 3};
    double z[] = {0, 0};
    postsolve(presolver, x, NULL, z);

    double correct_x[] = {-2.0 / 3, -10, 1.0 / 3};
    double correct_y[] = {0, 0};
    double correct_z[] = {0, 1, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Doubleton equality row where one variable appears in P.
    min. x0^2 - x2 s.t. x0 + x1 = 1, x0 + x2 <= 4, x1 + x2 <= 4, 0 <= x <= 10
    Only dton_eq is enabled. x1 = 1 - x0 must be substituted (x0 is in P),
    which transfers the bound x0 <= 1 and turns the third row into
    -x0 + x2 <= 3. Reduced (x0, x2): min x0^2 - x2 s.t. x0 + x2 <= 4,
    -x0 + x2 <= 3, solution (0.5, 3.5) with y = (0, -1). */
static char *test_6_qp()
{
    double Ax[] = {1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 2, 1, 2};
    int Ap[] = {0, 2, 4, 6};
    int n_rows = 3, n_cols = 3;
    double lhs[] = {1, -INF, -INF};
    double rhs[] = {1, 4, 4};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {0, 0, -1};
    double Px[] = {2};
    int Pi[] = {0};
    int Pp[] = {0, 1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 6, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 1, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x1 substituted, x0 kept",
              col_map[0] == 0 && col_map[1] == -1 && col_map[2] == 1);
    mu_assert("reduced dims", red->m == 2 && red->n == 2 && red->nnz == 4);
    mu_assert("reduced P", red->Pnnz == 1 && red->Px[0] == 2.0 && red->Pi[0] == 0);
    mu_assert("reduced c", red->c[0] == 0.0 && red->c[1] == -1.0);
    mu_assert("reduced bounds", red->ubs[0] == 1.0);
    double rhs_correct[] = {4, 3};
    mu_assert("reduced rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, red->rhs, 2));

    double x[] = {0.5, 3.5};
    double y[] = {0, -1};
    double z[] = {0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {0.5, 0.5, 3.5};
    double correct_y[] = {1, 0, -1};
    double correct_z[] = {0, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Parallel columns are not merged when one of them appears in P.
    min. x0^2 - x0 - 2 x1 s.t. x0 + 2 x1 + x2 <= 1, x0 + 2 x1 - x2 <= 1,
    -1 <= x <= 1. Columns 0 and 1 are parallel with c1 = 2 c0, so an LP
    would merge them. Nothing can be reduced here; the solution is
    x = (0, 0.5, 0), y = (-0.5, -0.5). */
static char *test_7_qp()
{
    double Ax[] = {1, 2, 1, 1, 2, -1};
    int Ai[] = {0, 1, 2, 0, 1, 2};
    int Ap[] = {0, 3, 6};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {-INF, -INF};
    double rhs[] = {1, 1};
    double lbs[] = {-1, -1, -1};
    double ubs[] = {1, 1, 1};
    double c[] = {-1, -2, 0};
    double Px[] = {2};
    int Pi[] = {0};
    int Pp[] = {0, 1, 1, 1};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 6, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 1, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == UNCHANGED);
    mu_assert("reduced dims",
              presolver->reduced_prob->n == 3 && presolver->reduced_prob->m == 2);
    mu_assert("P unchanged",
              presolver->reduced_prob->Pnnz == 1 &&
                  presolver->stats->Pnnz_original == presolver->stats->Pnnz_reduced);

    double x[] = {0, 0.5, 0};
    double y[] = {-0.5, -0.5};
    double z[] = {0, 0, 0};
    postsolve(presolver, x, y, z);
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, x, presolver->sol->y, y,
                                  presolver->sol->z, z, n_rows, n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Dual infeasibility ray after a variable in P is fixed by a singleton
    equality row.
    min. (1/2) x'Px + [1 0 -1] x,  P = [2 1 0; 1 2 0; 0 0 0]
    s.t. x0 = 1, x1 - x2 <= 5, 0 <= x0, x1 <= 10, 0 <= x2 < INF
    Reduced (x1, x2): min x1^2 + x1 - x2 s.t. x1 - x2 <= 5, unbounded along
    d = (0, 1). The recovered ray must be (0, 0, 1). */
static char *test_8_qp()
{
    double Ax[] = {1, 1, -1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 1, 3};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {1, -INF};
    double rhs[] = {1, 5};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, INF};
    double c[] = {1, 0, -1};
    double Px[] = {2, 1, 1, 2};
    int Pi[] = {0, 1, 0, 1};
    int Pp[] = {0, 2, 4, 4};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 3, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 4, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 1 && red->n == 2);
    mu_assert("reduced c", red->c[0] == 1.0 && red->c[1] == -1.0);
    mu_assert("reduced P", red->Pnnz == 1 && red->Px[0] == 2.0 && red->Pi[0] == 0);

    // the reduced ray is a valid certificate for the reduced problem
    double d_red[] = {0, 1};
    mu_assert("reduced ray", is_dual_ray_certificate(red, d_red, QP_TOL));

    double d[3];
    postsolve_dual_infeas_ray(presolver, d_red, d);
    double correct_d[] = {0, 0, 1};
    mu_assert("ray values", ABS(d[0] - correct_d[0]) <= QP_TOL &&
                                ABS(d[1] - correct_d[1]) <= QP_TOL &&
                                ABS(d[2] - correct_d[2]) <= QP_TOL);
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("ray certificate", is_dual_ray_certificate(&orig, d, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Primal propagation fixing two coupled variables (the FIXED_COL records
    are replayed inside the bound-change walk of the postsolver).
    min. (1/2) x'Px + [1 1] x,  P = [2 1; 1 2]
    s.t. x0 + x1 <= 0, 0 <= x <= 10
    Both variables are forced to 0, so the reduced problem is empty and the
    solution is x = 0, y = 0, z = (1, 1). */
static char *test_9_qp()
{
    double Ax[] = {1, 1};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    int n_rows = 1, n_cols = 2;
    double lhs[] = {-INF};
    double rhs[] = {0};
    double lbs[] = {0, 0};
    double ubs[] = {10, 10};
    double c[] = {1, 1};
    double Px[] = {2, 1, 1, 2};
    int Pi[] = {0, 1, 0, 1};
    int Pp[] = {0, 2, 4};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->primal_propagation = true;
    stgs->finite_bound_tightening = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 2, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 4, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 0 && red->n == 0 && red->Pnnz == 0);
    mu_assert("reduced offset", red->obj_offset == 0.0);

    postsolve(presolver, NULL, NULL, NULL);

    double correct_x[] = {0, 0};
    double correct_y[] = {0};
    double correct_z[] = {1, 1};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Doubleton equality row whose variable in P is a column singleton, so the
    other variable (which appears in more rows) must be substituted and the
    singleton column gains entries (found on QPLIB_3708).
    min. x0^2 - x2 s.t. x0 + x1 = 1, x1 + x2 <= 4, x1 - x2 <= 2, 0 <= x <= 10
    Substituting x1 = 1 - x0 gives -x0 + x2 <= 3, -x0 - x2 <= 1, x0 <= 1;
    solution (0.5, 3.5) with y = (-1, 0). */
static char *test_10_qp()
{
    double Ax[] = {1, 1, 1, 1, 1, -1};
    int Ai[] = {0, 1, 1, 2, 1, 2};
    int Ap[] = {0, 2, 4, 6};
    int n_rows = 3, n_cols = 3;
    double lhs[] = {1, -INF, -INF};
    double rhs[] = {1, 4, 2};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {0, 0, -1};
    double Px[] = {2};
    int Pi[] = {0};
    int Pp[] = {0, 1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 6, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 1, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x1 substituted, x0 kept",
              col_map[0] == 0 && col_map[1] == -1 && col_map[2] == 1);
    mu_assert("reduced dims", red->m == 2 && red->n == 2 && red->nnz == 4);
    mu_assert("reduced P", red->Pnnz == 1 && red->Px[0] == 2.0 && red->Pi[0] == 0);
    double rhs_correct[] = {3, 1};
    mu_assert("reduced rhs", ARRAYS_EQUAL_DOUBLE(rhs_correct, red->rhs, 2));

    double x[] = {0.5, 3.5};
    double y[] = {-1, 0};
    double z[] = {0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {0.5, 0.5, 3.5};
    double correct_y[] = {1, -1, 0};
    double correct_z[] = {0, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Diagonal dual fix (PRIMALC-like structure).
    min. x0^2 + x1^2 + x0 + 4 x1 s.t. x0 + x1 <= 5, x0 + 2 x1 <= 8, 0 <= x <= 10
    No variable has down-locks and both minimisers -c0 / P00 = -0.5 and
    -c1 / P11 = -2 are below lb = 0, so both are fixed at 0 and the rows become
    redundant. Everything is removed. */
static char *test_11_qp()
{
    double Ax[] = {1, 1, 1, 2};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    int n_rows = 2, n_cols = 2;
    double lhs[] = {-INF, -INF};
    double rhs[] = {5, 8};
    double lbs[] = {0, 0};
    double ubs[] = {10, 10};
    double c[] = {1, 4};
    double Px[] = {2, 2};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 4, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 2, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 0 && red->n == 0 && red->Pnnz == 0);
    mu_assert("reduced offset", red->obj_offset == 0.0);

    postsolve(presolver, NULL, NULL, NULL);

    double correct_x[] = {0, 0};
    double correct_y[] = {0, 0};
    double correct_z[] = {1, 4};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Two coupled variables fixed by singleton equality rows, in both row
    orders. The multiplier of the row fixed first needs the neighbour that is
    recovered by the later FIXED_COL record.
    min. (1/2) x'Px + [1 -1 -4] x,  P = [2 1 0; 1 2 1; 0 1 2]
    s.t. x1 = 3, x0 = 2, x1 + x2 <= 10, 0 <= x <= 10
    After the fixes x2 is an empty column (the last row becomes a bound) with
    only a diagonal in P, so it is fixed at -c2 / P22 = 0.5 and the reduced
    problem is empty, offset 17.75.
    Original: x = (2, 3, 0.5), z = 0, y = (7.5, 8, 0) for the row order
    above and (8, 7.5, 0) with the two equality rows swapped. */
static char *test_15_qp()
{
    for (int swap = 0; swap < 2; ++swap)
    {
        double Ax[] = {1, 1, 1, 1};
        int Ai[] = {1, 0, 1, 2};
        int Ap[] = {0, 1, 2, 4};
        int n_rows = 3, n_cols = 3;
        double lhs[] = {3, 2, -INF};
        double rhs[] = {3, 2, 10};
        double correct_y[] = {7.5, 8, 0};
        if (swap)
        {
            Ai[0] = 0;
            Ai[1] = 1;
            lhs[0] = rhs[0] = 2;
            lhs[1] = rhs[1] = 3;
            correct_y[0] = 8;
            correct_y[1] = 7.5;
        }
        double lbs[] = {0, 0, 0};
        double ubs[] = {10, 10, 10};
        double c[] = {1, -1, -4};
        double Px[] = {2, 1, 1, 2, 1, 1, 2};
        int Pi[] = {0, 1, 0, 1, 2, 1, 2};
        int Pp[] = {0, 2, 5, 7};

        Settings *stgs = default_settings();
        set_settings_false(stgs);
        Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 4, lhs, rhs,
                                            lbs, ubs, c, Px, Pi, Pp, 7, stgs);
        mu_assert("init", presolver != NULL);
        mu_assert("status", run_presolver(presolver) == REDUCED);

        ProblemData *red = presolver->reduced_prob;
        mu_assert("reduced dims", red->m == 0 && red->n == 0 && red->Pnnz == 0);
        mu_assert("reduced offset", red->obj_offset == 17.75);

        double x[] = {0};
        double y[] = {0};
        double z[] = {0};
        postsolve(presolver, x, y, z);

        double correct_x[] = {2, 3, 0.5};
        double correct_z[] = {0, 0, 0};
        mu_assert("postsolve values",
                  is_solution_correct(
                      presolver->sol->x, correct_x, presolver->sol->y, correct_y,
                      presolver->sol->z, correct_z, n_rows, n_cols, QP_TOL));
        ProblemData orig =
            problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs,
                                lbs, ubs, c, Px, Pi, Pp, Pp[n_cols]);
        mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                      presolver->sol->z, QP_TOL));

        PS_FREE(stgs);
        free_presolver(presolver);
    }
    return 0;
}

/*  Primal infeasibility ray after a variable in P is fixed by a singleton
    equality row.
    min. (1/2) x'Px + [1 -1 -4] x,  P = [2 1 0; 1 2 0; 0 0 0]
    s.t. x0 = 2, x1 + x2 >= 5, x1 + x2 <= 4, 0 <= x <= 10
    The reduced rows contradict each other; the reduced certificate (-1, 1)
    must come back as (0, -1, 1). */
static char *test_16_qp()
{
    double Ax[] = {1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 2, 1, 2};
    int Ap[] = {0, 1, 3, 5};
    int n_rows = 3, n_cols = 3;
    double lhs[] = {2, 5, -INF};
    double rhs[] = {2, INF, 4};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {1, -1, -4};
    double Px[] = {2, 1, 1, 2};
    int Pi[] = {0, 1, 0, 1};
    int Pp[] = {0, 2, 4, 4};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 5, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 4, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 2 && red->n == 2);

    double y_red[] = {-1, 1};
    mu_assert("reduced certificate", is_primal_ray_certificate(red, y_red, QP_TOL));

    double y[3];
    postsolve_primal_infeas_ray(presolver, y_red, y);
    double correct_y[] = {0, -1, 1};
    mu_assert("ray values", ABS(y[0] - correct_y[0]) <= QP_TOL &&
                                ABS(y[1] - correct_y[1]) <= QP_TOL &&
                                ABS(y[2] - correct_y[2]) <= QP_TOL);
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("ray certificate", is_primal_ray_certificate(&orig, y, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Diagonal-only empty columns with infinite bounds.
    min. (1/2) x'Px + [-6 2 0 -1 -1] x,  P = diag(2, 4, 1, 0, 0)
    s.t. x3 + x4 <= 5, x0 <= 1, x1 free, x2 >= 1, 0 <= x3, x4 <= 10
    x0: minimiser 3, clipped to 1 (z0 = -4). x1: minimiser -0.5 inside, both
    bounds infinite. x2: c = 0 gives 0, clipped to its lower bound 1 (z2 = 1).
    Reduced (x3, x4): min -x3 - x4 s.t. x3 + x4 <= 5, offset -5, one optimal
    point is x = (5, 0), y = -1, z = 0. */
static char *test_17_qp()
{
    double Ax[] = {1, 1};
    int Ai[] = {3, 4};
    int Ap[] = {0, 2};
    int n_rows = 1, n_cols = 5;
    double lhs[] = {-INF};
    double rhs[] = {5};
    double lbs[] = {-INF, -INF, 1, 0, 0};
    double ubs[] = {1, INF, INF, 10, 10};
    double c[] = {-6, 2, 0, -1, -1};
    double Px[] = {2, 4, 1};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3, 3, 3};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 2, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 3, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 1 && red->n == 2 && red->Pnnz == 0);
    mu_assert("reduced c", red->c[0] == -1.0 && red->c[1] == -1.0);
    mu_assert("reduced offset", red->obj_offset == -5.0);
    mu_assert("no empty cols left",
              presolver->prob->constraints->state->empty_cols->len == 0);

    double x[] = {5, 0};
    double y[] = {-1};
    double z[] = {0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {1, -0.5, 1, 5, 0};
    double correct_y[] = {-1};
    double correct_z[] = {-4, 0, 1, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Doubleton equality row with both variables diagonal in P.
    min. x0^2 + x1^2 - x2 s.t. x0 + x1 = 1, x1 + x2 <= 4, 0 <= x <= 10
    Only dton_eq is enabled. Equal coefficients, so the sparser column x0 is
    substituted: x0 = 1 - x1 gives x1 <= 1 and
      min. 2 x1^2 - 2 x1 - x2 + 1 s.t. x1 + x2 <= 4,
    solution (0.25, 3.75) with y = -1. Postsolve: x0 = 0.75 and
    y0 = (c0 + P00 x0) / a00 = 1.5. */
static char *test_18_qp()
{
    double Ax[] = {1, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2};
    int Ap[] = {0, 2, 4};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {1, -INF};
    double rhs[] = {1, 4};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {0, 0, -1};
    double Px[] = {2, 2};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2, 2};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 4, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 2, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x0 substituted",
              col_map[0] == -1 && col_map[1] == 0 && col_map[2] == 1);
    mu_assert("reduced dims", red->m == 1 && red->n == 2 && red->nnz == 2);
    mu_assert("reduced P", red->Pnnz == 1 && red->Px[0] == 4.0 && red->Pi[0] == 0);
    mu_assert("reduced c", red->c[0] == -2.0 && red->c[1] == -1.0);
    mu_assert("reduced offset", red->obj_offset == 1.0);
    mu_assert("reduced bounds", red->lbs[0] == 0.0 && red->ubs[0] == 1.0);
    mu_assert("reduced rhs", red->rhs[0] == 4.0);

    double x[] = {0.25, 3.75};
    double y[] = {-1};
    double z[] = {0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {0.75, 0.25, 3.75};
    double correct_y[] = {1.5, -1};
    double correct_z[] = {0, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  dton_choose_subst_qp: without P it is dton_choose_subst; with one column
    in P the other one is substituted, whatever its coefficient; with both in
    P the column that can be substituted without fill-in (sub0 / sub1: a
    positive diagonal alone, or coupled to the other column only) is
    substituted (dton_choose_subst decides when both can be), and a row
    without one is left alone. */
static char *test_19_qp()
{
    double v12[] = {1.0, 2.0};
    mu_assert("no P: larger c1",
              dton_choose_subst_qp(v12, 3, 3, false, false, false, false) == 1);
    mu_assert("c1 in P: c0",
              dton_choose_subst_qp(v12, 3, 3, false, true, false, true) == 0);
    mu_assert("c0 in P: c1",
              dton_choose_subst_qp(v12, 1, 5, true, false, true, false) == 1);
    mu_assert("both in P, coupled",
              dton_choose_subst_qp(v12, 3, 3, true, true, false, false) == -1);
    mu_assert("both in P, c0 fill-free",
              dton_choose_subst_qp(v12, 3, 3, true, true, true, false) == 0);
    mu_assert("both in P, c1 fill-free",
              dton_choose_subst_qp(v12, 3, 3, true, true, false, true) == 1);
    mu_assert("both fill-free: larger c1",
              dton_choose_subst_qp(v12, 1, 5, true, true, true, true) == 1);

    double vtiny[] = {1e-9, 1.0};
    mu_assert("forced tiny pivot",
              dton_choose_subst_qp(vtiny, 3, 3, false, true, false, false) == 0);
    mu_assert("unforced tiny pivot avoided",
              dton_choose_subst_qp(vtiny, 3, 3, false, false, false, false) == 1);
    return 0;
}

/*  Doubleton equality row with both variables coupled to a third variable
    in P. Substituting either would create fill-in, so the row stays on the
    worklist and the eliminator must still terminate.
    min. x0^2 + x1^2 + x2^2 + x0 x2 + x1 x2 - x2
    s.t. x0 + x1 = 1, x1 + x2 <= 4, -10 <= x <= 10
    Only dton_eq is enabled; nothing is reduced. Solution (0.5, 0.5, 0) with
    y = (1, 0). */
static char *test_20_qp()
{
    double Ax[] = {1, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2};
    int Ap[] = {0, 2, 4};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {1, -INF};
    double rhs[] = {1, 4};
    double lbs[] = {-10, -10, -10};
    double ubs[] = {10, 10, 10};
    double c[] = {0, 0, -1};
    double Px[] = {2, 1, 2, 1, 1, 1, 2};
    int Pi[] = {0, 2, 1, 2, 0, 1, 2};
    int Pp[] = {0, 2, 4, 7};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 4, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 7, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == UNCHANGED);

    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims",
              red->m == 2 && red->n == 3 && red->nnz == 4 && red->Pnnz == 7);
    mu_assert("row kept", red->lhs[0] == 1.0 && red->rhs[0] == 1.0);

    double x[] = {0.5, 0.5, 0};
    double y[] = {1, 0};
    double z[] = {0, 0, 0};
    postsolve(presolver, x, y, z);
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, x, presolver->sol->y, y,
                                  presolver->sol->z, z, n_rows, n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Doubleton equality row where the diagonal-only variable x0 has the smaller
    coefficient and sits in a second row.
    min. x0^2 + x1^2 + x1 x2 + x2^2 + 0.5 x0 + 1.5 x1 - x3
    s.t. x0 + 2 x1 = 2, x0 + x3 <= 5, -10 <= x <= 10
    Only dton_eq is enabled. x0 = 2 - 2 x1 (x1 is coupled) gives
    -4 <= x1 <= 6, P11 = 2 + 4 * 2 = 10, c1 = 1.5 - 2 (0.5 + 2 * 2) = -7.5,
    offset 2 (0.5 + 2) = 5 and the row -2 x1 + x3 <= 3. Reduced solution
    (1, -0.5, 5) with y = -1. Postsolve: x0 = 0 and
    y0 = (c0 + P00 x0 - a10 y1) / a00 = 1.5. */
static char *test_21_qp()
{
    double Ax[] = {1, 2, 1, 1};
    int Ai[] = {0, 1, 0, 3};
    int Ap[] = {0, 2, 4};
    int n_rows = 2, n_cols = 4;
    double lhs[] = {2, -INF};
    double rhs[] = {2, 5};
    double lbs[] = {-10, -10, -10, -10};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {0.5, 1.5, 0, -1};
    double Px[] = {2, 2, 1, 1, 2};
    int Pi[] = {0, 1, 2, 1, 2};
    int Pp[] = {0, 1, 3, 5, 5};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 4, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 5, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x0 substituted", col_map[0] == -1 && col_map[1] == 0 &&
                                    col_map[2] == 1 && col_map[3] == 2);
    mu_assert("reduced dims", red->m == 1 && red->n == 3 && red->nnz == 2);
    double Px_correct[] = {10, 1, 1, 2};
    int Pi_correct[] = {0, 1, 0, 1};
    int Pp_correct[] = {0, 2, 4, 4};
    mu_assert("reduced P", red->Pnnz == 4 &&
                               ARRAYS_EQUAL_DOUBLE(Px_correct, red->Px, 4) &&
                               ARRAYS_EQUAL_INT(Pi_correct, red->Pi, 4) &&
                               ARRAYS_EQUAL_INT(Pp_correct, red->Pp, 4));
    double c_correct[] = {-7.5, 0, -1};
    mu_assert("reduced c", ARRAYS_EQUAL_DOUBLE(c_correct, red->c, 3));
    mu_assert("reduced offset", red->obj_offset == 5.0);
    mu_assert("reduced bounds", red->lbs[0] == -4.0 && red->ubs[0] == 6.0);
    double Ax_correct[] = {-2, 1};
    mu_assert("reduced row", ARRAYS_EQUAL_DOUBLE(Ax_correct, red->Ax, 2) &&
                                 red->Ai[0] == 0 && red->Ai[1] == 2 &&
                                 red->rhs[0] == 3.0);

    double x[] = {1, -0.5, 5};
    double y[] = {-1};
    double z[] = {0, 0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {0, 1, -0.5, 5};
    double correct_y[] = {1.5, -1};
    double correct_z[] = {0, 0, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Two doubleton equality rows on diagonal-only variables form a chain in one
    round: 2 x0 + x1 = 1 substitutes x0 = 1/2 - x1/2 and 2 x1 + x2 = 2
    substitutes x1 = 1 - x2/2, so x0 = x2/4.
    min. x0^2 + x1^2 + x2^2 - 10.5 x2 - x3 s.t. the two rows, x2 + x3 <= 4,
    -10 <= x <= 10. Only dton_eq is enabled. Reduced (x2, x3):
    P22 = 2 + 2/16 + 2/4 = 2.625, c2 = -10.5 - 1 = -11.5, offset 1,
    solution (4, 0) with y = -1. Postsolve: x = (1, -1, 4, 0),
    y0 = (c0 + P00 x0) / 2 = 1 first, then y1 = (c1 + P11 x1 - y0) / 2 = -1.5. */
static char *test_22_qp()
{
    double Ax[] = {2, 1, 2, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 2, 3};
    int Ap[] = {0, 2, 4, 6};
    int n_rows = 3, n_cols = 4;
    double lhs[] = {1, 2, -INF};
    double rhs[] = {1, 2, 4};
    double lbs[] = {-10, -10, -10, -10};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {0, 0, -10.5, -1};
    double Px[] = {2, 2, 2};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3, 3};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 6, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 3, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x0, x1 substituted", col_map[0] == -1 && col_map[1] == -1 &&
                                        col_map[2] == 0 && col_map[3] == 1);
    mu_assert("reduced dims", red->m == 1 && red->n == 2 && red->nnz == 2);
    mu_assert("reduced P", red->Pnnz == 1 && red->Px[0] == 2.625 && red->Pi[0] == 0);
    mu_assert("reduced c", red->c[0] == -11.5 && red->c[1] == -1.0);
    mu_assert("reduced offset", red->obj_offset == 1.0);
    mu_assert("reduced row",
              red->rhs[0] == 4.0 && red->Ax[0] == 1.0 && red->Ax[1] == 1.0);

    double x[] = {4, 0};
    double y[] = {-1};
    double z[] = {0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {1, -1, 4, 0};
    double correct_y[] = {1, -1.5, -1};
    double correct_z[] = {0, 0, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Coupled dual fix of columns that are empty in A (CVXQP-like).
    min. x0^2 + x0 x1 + x1^2 + x1 x3 + x3^2 - 5.25 x1 - x2 - 5 x3
    s.t. x1 + x2 <= 4, 0.25 <= x0, x1, x3 <= 10, 0 <= x2 <= 10
    Only dual fix is enabled. x0 is empty in A and its derivative
    2 x0 + x1 >= 0.75 over the box, so x0 = 0.25: c1 = -5.25 + 0.25 = -5,
    offset 0.0625. x3 is empty in A too but -5 + x1 + 2 x3 changes sign over the
    box, so it stays; x1 has an up-lock and -5 + x3 + 2 x1 can be negative.
    Reduced (x1, x2, x3): solution (1, 3, 2) with y = -1. Postsolve:
    z0 = 2 x0 + x1 = 1.5. */
static char *test_23_qp()
{
    double Ax[] = {1, 1};
    int Ai[] = {1, 2};
    int Ap[] = {0, 2};
    int n_rows = 1, n_cols = 4;
    double lhs[] = {-INF};
    double rhs[] = {4};
    double lbs[] = {0.25, 0.25, 0, 0.25};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {0, -5.25, -1, -5};
    double Px[] = {2, 1, 1, 2, 1, 1, 2};
    int Pi[] = {0, 1, 0, 1, 3, 1, 3};
    int Pp[] = {0, 2, 5, 5, 7};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dual_fix = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 2, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 7, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x0 fixed", col_map[0] == -1 && col_map[1] == 0 && col_map[2] == 1 &&
                              col_map[3] == 2);
    mu_assert("reduced dims", red->m == 1 && red->n == 3 && red->nnz == 2);
    double Px_correct[] = {2, 1, 1, 2};
    int Pi_correct[] = {0, 2, 0, 2};
    int Pp_correct[] = {0, 2, 2, 4};
    mu_assert("reduced P", red->Pnnz == 4 &&
                               ARRAYS_EQUAL_DOUBLE(Px_correct, red->Px, 4) &&
                               ARRAYS_EQUAL_INT(Pi_correct, red->Pi, 4) &&
                               ARRAYS_EQUAL_INT(Pp_correct, red->Pp, 4));
    double c_correct[] = {-5, -1, -5};
    mu_assert("reduced c", ARRAYS_EQUAL_DOUBLE(c_correct, red->c, 3));
    mu_assert("reduced offset", red->obj_offset == 0.0625);

    double x[] = {1, 3, 2};
    double y[] = {-1};
    double z[] = {0, 0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {0.25, 1, 3, 2};
    double correct_y[] = {-1};
    double correct_z[] = {1.5, 0, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Coupled dual fix at the upper bound.
    min. x0^2 + x0 x1 + x1^2 - 4 x0 - 3 x1 - x2
    s.t. -x0 + x1 + x2 <= 2, x1 + 2 x2 <= 6, 0 <= x0, x1 <= 1, 0 <= x2 <= 10
    Only dual fix is enabled. x0 has a down-lock only and its derivative
    -4 + x1 + 2 x0 <= -1 over the box, so x0 = 1: c1 = -3 + 1 = -2,
    offset -4 + 1 = -3, first row x1 + x2 <= 3. x1 (up-locks, -2 + 2 x1 can
    be negative) and x2 stay. Reduced solution (0.5, 2.5) with y = (-1, 0).
    Postsolve: z0 = -4 + 2 x0 + x1 + y0 = -2.5. */
static char *test_24_qp()
{
    double Ax[] = {-1, 1, 1, 1, 2};
    int Ai[] = {0, 1, 2, 1, 2};
    int Ap[] = {0, 3, 5};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {-INF, -INF};
    double rhs[] = {2, 6};
    double lbs[] = {0, 0, 0};
    double ubs[] = {1, 1, 10};
    double c[] = {-4, -3, -1};
    double Px[] = {2, 1, 1, 2};
    int Pi[] = {0, 1, 0, 1};
    int Pp[] = {0, 2, 4, 4};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dual_fix = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 5, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 4, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x0 fixed", col_map[0] == -1 && col_map[1] == 0 && col_map[2] == 1);
    mu_assert("reduced dims", red->m == 2 && red->n == 2 && red->nnz == 4);
    mu_assert("reduced P", red->Pnnz == 1 && red->Px[0] == 2.0 && red->Pi[0] == 0);
    mu_assert("reduced c", red->c[0] == -2.0 && red->c[1] == -1.0);
    mu_assert("reduced offset", red->obj_offset == -3.0);
    mu_assert("reduced rhs", red->rhs[0] == 3.0 && red->rhs[1] == 6.0);

    double x[] = {0.5, 2.5};
    double y[] = {-1, 0};
    double z[] = {0, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {1, 0.5, 2.5};
    double correct_y[] = {-1, 0};
    double correct_z[] = {-2.5, 0, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  Doubleton equality row with both variables coupled to each other only.
    min. (1/2) x'Px + x2,  P = [2 1 0; 1 2 0; 0 0 0]
    s.t. x0 + x1 = 2, x0 + x2 <= 5, -10 <= x <= 10
    Only dton_eq is enabled. Row 1 of P holds the diagonal and column 0, so
    x1 = 2 - x0 can be substituted without fill-in (x1 is the sparser column):
    P00 = 2 + 2 - 2 = 2, c0 = -(0 + 2 * 2) + 1 * 2 = -2, offset 2 * 2 = 4,
    P01 removed, and the bounds of x1 transfer to -8 <= x0 <= 10.
    Reduced solution (x0, x2) = (1, -10), y = 0, z = (0, 1). Postsolve: x1 = 1
    and y0 = (P10 x0 + P11 x1) / a01 = 3. */
static char *test_25_qp()
{
    double Ax[] = {1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 2};
    int Ap[] = {0, 2, 4};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {2, -INF};
    double rhs[] = {2, 5};
    double lbs[] = {-10, -10, -10};
    double ubs[] = {10, 10, 10};
    double c[] = {0, 0, 1};
    double Px[] = {2, 1, 1, 2};
    int Pi[] = {0, 1, 0, 1};
    int Pp[] = {0, 2, 4, 4};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 4, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 4, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x1 substituted",
              col_map[0] == 0 && col_map[1] == -1 && col_map[2] == 1);
    mu_assert("reduced dims", red->m == 1 && red->n == 2 && red->nnz == 2);
    mu_assert("reduced P", red->Pnnz == 1 && red->Px[0] == 2.0 && red->Pi[0] == 0 &&
                               red->Pp[1] == 1 && red->Pp[2] == 1);
    mu_assert("reduced c", red->c[0] == -2.0 && red->c[1] == 1.0);
    mu_assert("reduced offset", red->obj_offset == 4.0);
    mu_assert("reduced bounds", red->lbs[0] == -8.0 && red->ubs[0] == 10.0);

    double x[] = {1, -10};
    double y[] = {0};
    double z[] = {0, 1};
    postsolve(presolver, x, y, z);

    double correct_x[] = {1, 1, -10};
    double correct_y[] = {3, 0};
    double correct_z[] = {0, 0, 1};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  As test 25, but the stay column x0 is also coupled to x3, so x0 cannot be
    substituted (fill-in) and x1 must be, whatever the column sizes.
    min. (1/2) x'Px + x2,  P = [2 1 0 1; 1 2 0 0; 0 0 0 0; 1 0 0 2]
    s.t. x0 + x1 = 2, x0 + x2 <= 5, x2 + x3 <= 4, -10 <= x <= 10
    Only dton_eq is enabled. x1 = 2 - x0: P00 = 2, P03 = P30 = 1, P33 = 2,
    c0 = -2, offset 4. Reduced solution (x0, x2, x3) = (4/3, -10, -2/3),
    y = 0, z = (0, 1, 0). Postsolve: x1 = 2/3, y0 = (4/3 + 2 * 2/3) / 1 = 8/3. */
static char *test_26_qp()
{
    double Ax[] = {1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 2, 2, 3};
    int Ap[] = {0, 2, 4, 6};
    int n_rows = 3, n_cols = 4;
    double lhs[] = {2, -INF, -INF};
    double rhs[] = {2, 5, 4};
    double lbs[] = {-10, -10, -10, -10};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {0, 0, 1, 0};
    double Px[] = {2, 1, 1, 1, 2, 1, 2};
    int Pi[] = {0, 1, 3, 0, 1, 0, 3};
    int Pp[] = {0, 3, 5, 5, 7};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 6, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 7, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    ProblemData *red = presolver->reduced_prob;
    const int *col_map = presolver->prob->constraints->state->work->mappings->cols;
    mu_assert("x1 substituted", col_map[0] == 0 && col_map[1] == -1 &&
                                    col_map[2] == 1 && col_map[3] == 2);
    mu_assert("reduced dims", red->m == 2 && red->n == 3 && red->nnz == 4);
    double Px_correct[] = {2, 1, 1, 2};
    int Pi_correct[] = {0, 2, 0, 2};
    int Pp_correct[] = {0, 2, 2, 4};
    mu_assert("reduced P", red->Pnnz == 4 &&
                               ARRAYS_EQUAL_DOUBLE(Px_correct, red->Px, 4) &&
                               ARRAYS_EQUAL_INT(Pi_correct, red->Pi, 4) &&
                               ARRAYS_EQUAL_INT(Pp_correct, red->Pp, 4));
    double c_correct[] = {-2, 1, 0};
    mu_assert("reduced c", ARRAYS_EQUAL_DOUBLE(c_correct, red->c, 3));
    mu_assert("reduced offset", red->obj_offset == 4.0);

    double x[] = {4.0 / 3.0, -10, -2.0 / 3.0};
    double y[] = {0, 0};
    double z[] = {0, 1, 0};
    postsolve(presolver, x, y, z);

    double correct_x[] = {4.0 / 3.0, 2.0 / 3.0, -10, -2.0 / 3.0};
    double correct_y[] = {8.0 / 3.0, 0, 0};
    double correct_z[] = {0, 0, 1, 0};
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, correct_x, presolver->sol->y,
                                  correct_y, presolver->sol->z, correct_z, n_rows,
                                  n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  No fold of a QP: three copies of  min -x0 - 2 x1  s.t.  x0 + x1 <= 1,
    0 <= x <= 1  fold to one copy as an LP (test_fold_presolver_lift), but with
    P00 = 2 on the first copy nothing folds, though lp_fold is enabled. KKT
    point: x = (0, 1) and z = (1, 0) per copy, y = -2. With P = 0 through
    new_presolver_from_problem the fold runs as for an LP. */
static char *test_29_qp()
{
    double Ax[] = {1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 2, 3, 4, 5};
    int Ap[] = {0, 2, 4, 6};
    int n_rows = 3, n_cols = 6;
    double lhs[] = {-INF, -INF, -INF};
    double rhs[] = {1, 1, 1};
    double lbs[] = {0, 0, 0, 0, 0, 0};
    double ubs[] = {1, 1, 1, 1, 1, 1};
    double c[] = {-1, -2, -1, -2, -1, -2};
    double Px[] = {2};
    int Pi[] = {0};
    int Pp[] = {0, 1, 1, 1, 1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->lp_fold = true;
    Presolver *presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 6, lhs, rhs, lbs,
                                        ubs, c, Px, Pi, Pp, 1, stgs);
    mu_assert("init", presolver != NULL);
    const PresolveStats *stats = presolver->stats;
    mu_assert("not folded", stats->n_rows_folded == 3 && stats->n_cols_folded == 6 &&
                                stats->nnz_removed_lp_fold == 0 &&
                                stats->Pnnz_original == 1);
    const Matrix *P = presolver->prob->obj->P;
    mu_assert("P unchanged", P->m == 6 && P->n == 6 && P->nnz == 1);
    mu_assert("presolve", run_presolver(presolver) != INFEASIBLE);
    ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims",
              red->m == 3 && red->n == 6 && red->nnz == 6 && red->Pnnz == 1);

    double x[] = {0, 1, 0, 1, 0, 1};
    double y[] = {-2, -2, -2};
    double z[] = {1, 0, 1, 0, 1, 0};
    mu_assert("reduced KKT", is_kkt_point(red, x, y, z, QP_TOL));
    postsolve(presolver, x, y, z);
    mu_assert("postsolve values",
              is_solution_correct(presolver->sol->x, x, presolver->sol->y, y,
                                  presolver->sol->z, z, n_rows, n_cols, QP_TOL));
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));
    free_presolver(presolver);

    // the same problem as an LP through the QP entry point folds to one copy
    presolver = qp_presolver(Ax, Ai, Ap, n_rows, n_cols, 6, lhs, rhs, lbs, ubs, c,
                             NULL, NULL, NULL, 0, stgs);
    mu_assert("init lp", presolver != NULL);
    mu_assert("folded", presolver->stats->n_rows_folded == 1 &&
                            presolver->stats->n_cols_folded == 2 &&
                            presolver->stats->nnz_removed_lp_fold == 4);
    mu_assert("no P for an LP", presolver->prob->obj->P == NULL);
    mu_assert("presolve lp", run_presolver(presolver) != INFEASIBLE);
    mu_assert("reduced lp", presolver->reduced_prob->m == 1 &&
                                presolver->reduced_prob->n == 2 &&
                                presolver->reduced_prob->Pnnz == 0);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/*  The objective offset of the input is carried over to the reduced problem.
    Same problem as test_3_qp, given with obj_offset = 2.5: the presolver adds
    6 when it fixes x0, so the reduced offset is 8.5. */
static char *test_30_qp()
{
    double Ax[] = {1, 1, 1};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 1, 3};
    int n_rows = 2, n_cols = 3;
    double lhs[] = {2, -INF};
    double rhs[] = {2, 3};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {1, -3, 1};
    double Px[] = {2, 1, 1, 2};
    int Pi[] = {0, 1, 0, 1};
    int Pp[] = {0, 2, 4, 4};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    ProblemData problem = {.Ax = Ax,
                           .Ai = Ai,
                           .Ap = Ap,
                           .m = (size_t) n_rows,
                           .n = (size_t) n_cols,
                           .nnz = 3,
                           .lhs = lhs,
                           .rhs = rhs,
                           .c = c,
                           .lbs = lbs,
                           .ubs = ubs,
                           .Px = Px,
                           .Pi = Pi,
                           .Pp = Pp,
                           .Pnnz = 4,
                           .obj_offset = 2.5};
    Presolver *presolver = new_presolver_from_problem(&problem, stgs);
    mu_assert("init", presolver != NULL);
    mu_assert("status", run_presolver(presolver) == REDUCED);

    const ProblemData *red = presolver->reduced_prob;
    mu_assert("reduced dims", red->m == 1 && red->n == 2 && red->nnz == 2);
    mu_assert("reduced offset", red->obj_offset == 8.5);

    double x[] = {0.5, 0};
    double y[] = {0};
    double z[] = {0, 1};
    postsolve(presolver, x, y, z);
    ProblemData orig =
        problem_from_csr_qp(Ax, Ai, Ap, n_rows, n_cols, Ap[n_rows], lhs, rhs, lbs,
                            ubs, c, Px, Pi, Pp, Pp[n_cols]);
    orig.obj_offset = 2.5;
    mu_assert("KKT", is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                                  presolver->sol->z, QP_TOL));

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

static const char *all_tests_qp()
{
    mu_run_test(test_1_qp, counter_qp);
    mu_run_test(test_2_qp, counter_qp);
    mu_run_test(test_3_qp, counter_qp);
    mu_run_test(test_4_qp, counter_qp);
    mu_run_test(test_5_qp, counter_qp);
    mu_run_test(test_6_qp, counter_qp);
    mu_run_test(test_7_qp, counter_qp);
    mu_run_test(test_8_qp, counter_qp);
    mu_run_test(test_9_qp, counter_qp);
    mu_run_test(test_10_qp, counter_qp);
    mu_run_test(test_11_qp, counter_qp);
    mu_run_test(test_15_qp, counter_qp);
    mu_run_test(test_16_qp, counter_qp);
    mu_run_test(test_17_qp, counter_qp);
    mu_run_test(test_18_qp, counter_qp);
    mu_run_test(test_19_qp, counter_qp);
    mu_run_test(test_20_qp, counter_qp);
    mu_run_test(test_21_qp, counter_qp);
    mu_run_test(test_22_qp, counter_qp);
    mu_run_test(test_23_qp, counter_qp);
    mu_run_test(test_24_qp, counter_qp);
    mu_run_test(test_25_qp, counter_qp);
    mu_run_test(test_26_qp, counter_qp);
    mu_run_test(test_29_qp, counter_qp);
    mu_run_test(test_30_qp, counter_qp);
    return 0;
}

int test_qp()
{
    const char *result = all_tests_qp();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("QP: TEST FAILED!\n");
    }
    else
    {
        printf("QP: ALL TESTS PASSED\n");
    }
    printf("QP: Tests run: %d\n", counter_qp);
    return result == 0;
}

#endif // TEST_QP_H
