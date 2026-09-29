#ifndef TEST_DTON_H
#define TEST_DTON_H

#include "Debugger.h"
#include "DtonsEq_internal.h"
#include "PSLP_API.h"

#include "Activity.h"
#include "Binary_search.h"
#include "Constraints.h"
#include "Matrix.h"
#include "PSLP_sol.h"
#include "Postsolver.h"
#include "Problem.h"
#include "State.h"
#include "Workspace.h"
#include "debug_macros.h"
#include "kkt.h"
#include "minunit.h"
#include "u16Vec.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int counter_dton = 0;

/* The claim record of column k (must be claimed in the current round). */
static const DtonSubst *dton_rec(const DtonWorkspace *dton_work, int k)
{
    assert(dton_work->substs.col_subst[k] >= 0);
    return dton_work->substs.recs + dton_work->substs.col_subst[k];
}

/* The chain depth of column k's record after composition. */
static int dton_depth(const DtonWorkspace *dton_work, int k)
{
    assert(dton_work->substs.col_subst[k] >= 0);
    return dton_work->substs.depth[dton_work->substs.col_subst[k]];
}

/* Workspace lifecycle: allocated lazily by the first elimination, sized by
   the problem, freed with the presolver (run under ASAN in CI). */
static char *test_dton_workspace()
{
    // x0 + x1 = 1 (doubleton equality), x0 + x1 + x2 <= 4
    double Ax[] = {1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 1, 2};
    int Ap[] = {0, 2, 5};
    int nnz = 5;
    int n_rows = 2;
    int n_cols = 3;

    double lhs[] = {1, -INF};
    double rhs[] = {1, 4};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    Work *work = presolver->prob->constraints->state->work;
    mu_assert("dton workspace must be allocated", work->dton != NULL);
    mu_assert("workspace m mismatch", work->dton->m == n_rows);
    mu_assert("workspace n mismatch", work->dton->n == n_cols);

    run_presolver(presolver);
    mu_assert("dton workspace must be freed after presolve", work->dton == NULL);
    free_presolver(presolver);

    stgs->dton_eq = false;
    presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);
    work = presolver->prob->constraints->state->work;
    mu_assert("dton workspace must not be allocated", work->dton == NULL);

    free_presolver(presolver);
    PS_FREE(stgs);
    dton_ws_free(NULL); // must tolerate NULL
    return 0;
}

static char *test_dton_choose_subst()
{
    // the larger coefficient is substituted, whatever the column sizes
    double v12[] = {1.0, 2.0};
    mu_assert("larger c1", dton_choose_subst(v12, 3, 3) == 1);
    mu_assert("larger c1 beats singleton c0", dton_choose_subst(v12, 1, 5) == 1);
    mu_assert("larger c1 beats sparser c0", dton_choose_subst(v12, 2, 3) == 1);
    double v21[] = {2.0, 1.0};
    mu_assert("larger c0", dton_choose_subst(v21, 3, 3) == 0);
    double vhuge[] = {1e9, 1.0};
    mu_assert("huge pivot substituted", dton_choose_subst(vhuge, 3, 3) == 0);
    double vtiny[] = {1e-9, 1.0};
    mu_assert("tiny pivot avoided", dton_choose_subst(vtiny, 3, 3) == 1);

    // equal magnitude: singleton column first, then the sparser column
    double v11[] = {1.0, -1.0};
    mu_assert("singleton c0", dton_choose_subst(v11, 1, 5) == 0);
    mu_assert("singleton c1", dton_choose_subst(v11, 5, 1) == 1);
    mu_assert("sparser c0", dton_choose_subst(v11, 2, 3) == 0);
    mu_assert("sparser c1", dton_choose_subst(v11, 3, 2) == 1);

    return 0;
}

/* Two doubleton rows choosing the same substituted column: the first
   claims it, the second is deferred (no orientation flip). */
static char *test_dton_claim_conflict()
{
    // r0: 3 x0 + x1 = 1, r1: 5 x0 + x2 = 2, r2: x1 + x2 <= 10 (pad)
    double Ax[] = {3, 1, 5, 1, 1, 1};
    int Ai[] = {0, 1, 0, 2, 1, 2};
    int Ap[] = {0, 2, 4, 6};
    double lhs[] = {1, 2, -INF};
    double rhs[] = {1, 2, 10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 3, 6, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    DtonWorkspace *dton_work = dton_ws_new(3, 3);
    dton_ws_attach(dton_work, presolver->prob->constraints->state->work);
    int deferred[8];
    int n_deferred = 0;
    mu_assert("claim",
              dton_claim(presolver->prob, dton_work, deferred, &n_deferred));

    mu_assert("one column claimed", dton_work->substs.n == 1);
    mu_assert("x0 claimed", dton_work->substs.recs[0].k == 0);
    mu_assert("x0 owned by r0", dton_rec(dton_work, 0)->owner == 0);
    mu_assert("record owner", dton_work->substs.recs[0].owner == 0);
    mu_assert("x0 stays into x1", dton_rec(dton_work, 0)->j == 1);
    mu_assert("record stay col", dton_work->substs.recs[0].j == 1);
    mu_assert("dir_mult", dton_rec(dton_work, 0)->dir_mult == -1.0 / 3.0);
    mu_assert("dir_shift", dton_rec(dton_work, 0)->dir_shift == 1.0 / 3.0);
    mu_assert("loser deferred", n_deferred == 1 && deferred[0] == 1);

    dton_compose(dton_work, deferred, &n_deferred);
    mu_assert("still one eliminated", dton_work->substs.n == 1);
    mu_assert("composed target", dton_rec(dton_work, 0)->target == 1);
    mu_assert("composed mult", dton_rec(dton_work, 0)->mult == -1.0 / 3.0);
    mu_assert("composed shift", dton_rec(dton_work, 0)->shift == 1.0 / 3.0);
    mu_assert("depth 0", dton_depth(dton_work, 0) == 0);
    mu_assert("compose defers nothing here", n_deferred == 1);

    dton_ws_free(dton_work);
    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* Chain x0 -> x1 -> x2: r0 eliminates x0 into x1 (singleton rule), r1
   eliminates x1 into x2 (larger coefficient); composition maps both onto
   the final survivor x2 with hand-computed maps and depths. */
static char *test_dton_chain_depth2()
{
    // r0: x0 + x1 = 1, r1: 2 x1 + x2 = 4, r2: x2 + x3 <= 10 (pad)
    double Ax[] = {1, 1, 2, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 2, 3};
    int Ap[] = {0, 2, 4, 6};
    double lhs[] = {1, 4, -INF};
    double rhs[] = {1, 4, 10};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {1, 1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 4, 6, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    DtonWorkspace *dton_work = dton_ws_new(3, 4);
    dton_ws_attach(dton_work, presolver->prob->constraints->state->work);
    int deferred[8];
    int n_deferred = 0;
    mu_assert("claim",
              dton_claim(presolver->prob, dton_work, deferred, &n_deferred));

    mu_assert("two columns claimed", dton_work->substs.n == 2);
    mu_assert("no deferrals", n_deferred == 0);
    mu_assert("x0 -> x1", dton_rec(dton_work, 0)->j == 1);
    mu_assert("x1 -> x2", dton_rec(dton_work, 1)->j == 2);

    dton_compose(dton_work, deferred, &n_deferred);

    // x1 = -0.5 x2 + 2, depth 0
    mu_assert("x1 target", dton_rec(dton_work, 1)->target == 2);
    mu_assert("x1 mult", dton_rec(dton_work, 1)->mult == -0.5);
    mu_assert("x1 shift", dton_rec(dton_work, 1)->shift == 2.0);
    mu_assert("x1 depth", dton_depth(dton_work, 1) == 0);

    // x0 = -x1 + 1 = 0.5 x2 - 1, depth 1
    mu_assert("x0 target", dton_rec(dton_work, 0)->target == 2);
    mu_assert("x0 mult", dton_rec(dton_work, 0)->mult == 0.5);
    mu_assert("x0 shift", dton_rec(dton_work, 0)->shift == -1.0);
    mu_assert("x0 depth", dton_depth(dton_work, 0) == 1);

    mu_assert("both survive composition", dton_work->substs.n == 2);
    mu_assert("compose defers nothing", n_deferred == 0);

    dton_ws_free(dton_work);
    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* A 2-cycle x0 -> x1 -> x0: the node with the larger owner row (x1,
   owned by r1) is un-eliminated and its row deferred; the remaining link
   composes onto the now-surviving x1. */
static char *test_dton_cycle_break()
{
    // r0: 2 x0 + x1 = 0 (claims x0), r1: x0 + 2 x1 = 0 (claims x1)
    double Ax[] = {2, 1, 1, 2};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    double lhs[] = {0, 0};
    double rhs[] = {0, 0};
    double lbs[] = {0, 0};
    double ubs[] = {10, 10};
    double c[] = {1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 2, 2, 4, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    DtonWorkspace *dton_work = dton_ws_new(2, 2);
    dton_ws_attach(dton_work, presolver->prob->constraints->state->work);
    int deferred[8];
    int n_deferred = 0;
    mu_assert("claim",
              dton_claim(presolver->prob, dton_work, deferred, &n_deferred));

    mu_assert("both columns claimed", dton_work->substs.n == 2);
    mu_assert("cycle x0 -> x1", dton_rec(dton_work, 0)->j == 1);
    mu_assert("cycle x1 -> x0", dton_rec(dton_work, 1)->j == 0);
    mu_assert("no claim deferrals", n_deferred == 0);

    dton_compose(dton_work, deferred, &n_deferred);

    mu_assert("one column survives the break", dton_work->substs.n == 1);
    mu_assert("x0 stays eliminated", dton_work->substs.recs[0].k == 0);
    mu_assert("record depth filled", dton_work->substs.depth[0] == 0);
    mu_assert("x1 un-eliminated", dton_work->substs.col_subst[1] < 0);
    mu_assert("broken owner deferred", n_deferred == 1 && deferred[0] == 1);
    mu_assert("x0 composes onto x1", dton_rec(dton_work, 0)->target == 1);
    mu_assert("x0 mult", dton_rec(dton_work, 0)->mult == -0.5);
    mu_assert("x0 shift", dton_rec(dton_work, 0)->shift == 0.0);
    mu_assert("x0 depth", dton_depth(dton_work, 0) == 0);

    dton_ws_free(dton_work);
    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* An ill-conditioned doubleton is eliminated through its large coefficient:
   x0 = (1 - x1) / 1e9, never x1 = 1 - 1e9 x0. */
static char *test_dton_pivot_large()
{
    // 1e9 x0 + x1 = 1
    double Ax[] = {1e9, 1};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    double lhs[] = {1};
    double rhs[] = {1};
    double lbs[] = {0, 0};
    double ubs[] = {10, 10};
    double c[] = {1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    DtonWorkspace *dton_work = dton_ws_new(1, 2);
    dton_ws_attach(dton_work, presolver->prob->constraints->state->work);
    int deferred[8];
    int n_deferred = 0;
    mu_assert("claim",
              dton_claim(presolver->prob, dton_work, deferred, &n_deferred));

    mu_assert("one claim", dton_work->substs.n == 1 && n_deferred == 0);
    mu_assert("x0 substituted",
              dton_rec(dton_work, 0)->k == 0 && dton_rec(dton_work, 0)->j == 1);
    mu_assert("multiplier -1e-9", dton_rec(dton_work, 0)->dir_mult == -1e-9);

    dton_ws_free(dton_work);
    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* Simplest full round: one doubleton whose substituted column appears in
   no other row. Bounds transfer onto the stay column, the owner row and
   the substituted column deactivate, the objective folds, and the stay
   column's size transition lands in ston_cols. */
static char *test_dton_eliminate_isolated()
{
    // r0: 2 x0 + x1 = 2 with x0 in [0.4, 0.6]  =>  x1 in [0.8, 1.2]
    // r1: x1 + x2 <= 10 (pad so x1 is not a singleton)
    double Ax[] = {2, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2};
    int Ap[] = {0, 2, 4};
    double lhs[] = {2, -INF};
    double rhs[] = {2, 10};
    double lbs[] = {0.4, 0, 0};
    double ubs[] = {0.6, 10, 10};
    double c[] = {1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 2, 3, 4, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;

    PresolveStatus status = remove_dton_eq_rows(prob);
    mu_assert("returns UNCHANGED", status == UNCHANGED);

    // bounds moved onto the stay column, its own bounds untouched
    mu_assert("stay lb tightened", constraints->bounds[1].lb == 0.8);
    mu_assert("stay ub tightened", constraints->bounds[1].ub == 1.2);
    mu_assert("subst bounds untouched",
              constraints->bounds[0].lb == 0.4 && constraints->bounds[0].ub == 0.6);

    // owner row and substituted column are gone
    mu_assert("owner row inactive",
              HAS_TAG(constraints->row_tags[0], R_TAG_INACTIVE));
    mu_assert("subst col inactive",
              HAS_TAG(constraints->col_tags[0], C_TAG_INACTIVE));
    mu_assert("nnz drops by the owner row", constraints->A->nnz == 2);
    mu_assert("row size inactive",
              constraints->state->row_sizes[0] == SIZE_INACTIVE_ROW);
    mu_assert("col size inactive",
              constraints->state->col_sizes[0] == SIZE_INACTIVE_COL);

    // objective: x0 = (2 - x1)/2 => c.x = 0.5 x1 + x2 + ... + 1
    mu_assert("obj stay coeff", prob->obj->c[1] == 0.5);
    mu_assert("obj subst zeroed", prob->obj->c[0] == 0.0);
    mu_assert("obj offset", prob->obj->offset == 1.0);

    // the stay column shrank 2 -> 1: singleton-column worklist
    mu_assert("stay col size", constraints->state->col_sizes[1] == 1);
    mu_assert("stay col pushed to ston_cols",
              iVec_contains(constraints->state->ston_cols, 1));

    mu_assert("worklist empty", constraints->state->dton_rows->len == 0);

    DEBUG(run_debugger(constraints, false));

    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* Substitution into other rows: fill-in in one row, an in-place merge
   that turns another row into a new doubleton candidate, side shifts
   respecting infinite sides, and a second round that consumes the new
   candidate. Hand-computed end state. */
static char *test_dton_apply_substitution()
{
    // r0: 2 x0 + x1 = 4          (round 1 eliminates x0 = 2 - 0.5 x1)
    // r1: 3 x0 + x2 + x3 <= 6    (fill-in of x1; rhs shifts by 6)
    // r2: 2 x0 + 3 x1 + 4 x4 = 9 (merge: becomes 2 x1 + 4 x4 = 5, a new
    //                             doubleton; round 2 eliminates x4)
    // r3: x1 + x3 <= 10          (untouched)
    double Ax[] = {2, 1, 3, 1, 1, 2, 3, 4, 1, 1};
    int Ai[] = {0, 1, 0, 2, 3, 0, 1, 4, 1, 3};
    int Ap[] = {0, 2, 5, 8, 10};
    double lhs[] = {4, -INF, 9, -INF};
    double rhs[] = {4, 6, 9, 10};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10, 10};
    double c[] = {1, 1, 1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 4, 5, 10, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    remove_dton_eq_rows(prob);

    // rows r0 and r2 eliminated, columns x0 and x4 eliminated
    mu_assert("r0 inactive", HAS_TAG(constraints->row_tags[0], R_TAG_INACTIVE));
    mu_assert("r2 inactive", HAS_TAG(constraints->row_tags[2], R_TAG_INACTIVE));
    mu_assert("x0 inactive", HAS_TAG(constraints->col_tags[0], C_TAG_INACTIVE));
    mu_assert("x4 inactive", HAS_TAG(constraints->col_tags[4], C_TAG_INACTIVE));
    mu_assert("final nnz", A->nnz == 5);

    // r1 = -1.5 x1 + x2 + x3 <= 0 (fill-in at the front, sorted)
    mu_assert("r1 length", constraints->state->row_sizes[1] == 3);
    mu_assert("r1 cols", A->i[A->p[1].start] == 1 && A->i[A->p[1].start + 1] == 2 &&
                             A->i[A->p[1].start + 2] == 3);
    mu_assert("r1 vals", A->x[A->p[1].start] == -1.5 &&
                             A->x[A->p[1].start + 1] == 1.0 &&
                             A->x[A->p[1].start + 2] == 1.0);
    mu_assert("r1 rhs shifted", constraints->rhs[1] == 0.0);
    mu_assert("r1 lhs stays -inf", constraints->lhs[1] == -INF);

    // r3 untouched
    mu_assert("r3 length", constraints->state->row_sizes[3] == 2);
    mu_assert("r3 rhs", constraints->rhs[3] == 10.0);

    // bounds chained across the two rounds: round 1 gives x1 <= 4, round
    // 2 tightens x1 <= 2.5 from x4's bounds through 2 x1 + 4 x4 = 5
    mu_assert("x1 lb", constraints->bounds[1].lb == 0.0);
    mu_assert("x1 ub", constraints->bounds[1].ub == 2.5);

    // objective: c[x1] = 1 - 0.5 - 0.5 = 0, offset = 2 + 1.25 = 3.25
    mu_assert("obj x1", prob->obj->c[1] == 0.0);
    mu_assert("obj x0 zeroed", prob->obj->c[0] == 0.0);
    mu_assert("obj x4 zeroed", prob->obj->c[4] == 0.0);
    mu_assert("obj offset", prob->obj->offset == 3.25);

    // col sizes against the rebuilt AT
    mu_assert("x1 size", constraints->state->col_sizes[1] == 2);
    mu_assert("x2 size", constraints->state->col_sizes[2] == 1);
    mu_assert("x3 size", constraints->state->col_sizes[3] == 2);

    mu_assert("worklist empty", constraints->state->dton_rows->len == 0);

    DEBUG(run_debugger(constraints, false));

    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* Cancellation policy: a substituted contribution that cancels an
   existing entry to (near) zero drops the entry, while a tiny coefficient
   the substitution never touched is kept verbatim. */
static char *test_dton_cancellation()
{
    // r0: x0 + 2 x1 = 0          (eliminates x1 = -0.5 x0)
    // r1: 0.5 x0 + x1 + x2 <= 5  (x0 entry cancels: 0.5 - 0.5 = 0 -> dropped)
    // r2: x1 + 1e-12 x3 + x4 <= 7 (tiny untouched x3 entry must survive)
    double Ax[] = {1, 2, 0.5, 1, 1, 1, 1e-12, 1};
    int Ai[] = {0, 1, 0, 1, 2, 1, 3, 4};
    int Ap[] = {0, 2, 5, 8};
    double lhs[] = {0, -INF, -INF};
    double rhs[] = {0, 5, 7};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10, 10};
    double c[] = {0, 0, 0, 0, 0};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 5, 8, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;
    Matrix *A = constraints->A;

    remove_dton_eq_rows(prob);

    mu_assert("x1 inactive", HAS_TAG(constraints->col_tags[1], C_TAG_INACTIVE));
    mu_assert("r0 inactive", HAS_TAG(constraints->row_tags[0], R_TAG_INACTIVE));

    // r1 collapsed to the single x2 entry and joined the singleton rows
    mu_assert("r1 shrank to 1", constraints->state->row_sizes[1] == 1);
    mu_assert("r1 keeps x2 only", A->i[A->p[1].start] == 2);
    mu_assert("r1 rhs unshifted", constraints->rhs[1] == 5.0);
    mu_assert("r1 in ston_rows", iVec_contains(constraints->state->ston_rows, 1));

    // r2: x1 mapped onto x0 (fill -0.5), the tiny untouched x3 entry stays
    mu_assert("r2 length", constraints->state->row_sizes[2] == 3);
    mu_assert("r2 cols", A->i[A->p[2].start] == 0 && A->i[A->p[2].start + 1] == 3 &&
                             A->i[A->p[2].start + 2] == 4);
    mu_assert("r2 fill value", A->x[A->p[2].start] == -0.5);
    mu_assert("r2 tiny entry kept", A->x[A->p[2].start + 1] == 1e-12);

    mu_assert("final nnz", A->nnz == 4);

    DEBUG(run_debugger(constraints, false));

    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* A chain that eliminates the entire matrix: three doubletons claim
   x0 -> x1, x2 -> x1,
   x3 -> x2, every row is an owner row, and the survivor column x1 ends
   empty. */
static char *test_dton_chain_empties_matrix()
{
    // r0: x0 + x1 = 1, r1: x1 + x2 = 1, r2: x2 + x3 = 1
    double Ax[] = {1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 2, 3};
    int Ap[] = {0, 2, 4, 6};
    double lhs[] = {1, 1, 1};
    double rhs[] = {1, 1, 1};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {1, 1, 1, 1};
    double c[] = {1, 1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 4, 6, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;

    remove_dton_eq_rows(prob);

    mu_assert("matrix fully eliminated", constraints->A->nnz == 0);
    for (int r = 0; r < 3; ++r)
    {
        mu_assert("all rows inactive",
                  HAS_TAG(constraints->row_tags[r], R_TAG_INACTIVE));
    }
    mu_assert("x0 inactive", HAS_TAG(constraints->col_tags[0], C_TAG_INACTIVE));
    mu_assert("x2 inactive", HAS_TAG(constraints->col_tags[2], C_TAG_INACTIVE));
    mu_assert("x3 inactive", HAS_TAG(constraints->col_tags[3], C_TAG_INACTIVE));
    mu_assert("survivor active", !HAS_TAG(constraints->col_tags[1], C_TAG_INACTIVE));
    mu_assert("survivor empty", constraints->state->col_sizes[1] == 0);
    mu_assert("survivor in empty_cols",
              iVec_contains(constraints->state->empty_cols, 1));

    // x0 = 1 - x1, x2 = 1 - x1, x3 = x1: objective folds to 2 + 0 * x1
    mu_assert("obj survivor", prob->obj->c[1] == 0.0);
    mu_assert("obj offset", prob->obj->offset == 2.0);

    mu_assert("worklist empty", constraints->state->dton_rows->len == 0);

    DEBUG(run_debugger(constraints, false));

    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* Within-round chain: bound transfer must run deepest link first and
   read the intermediate column's bounds LIVE. r0: x0 + x1 = 1 claims
   x0 -> x1 (singleton rule); r1: 2 x1 + x2 = 4 claims x1 -> x2 (larger
   coefficient; the r2 pad keeps x2 non-singleton). Depth-descending
   transfer: x0 in [0,1] gives x1 in [0,1], and the UPDATED x1 bounds give
   x2 = 4 - 2 x1 in [2, 4]. A stale snapshot of x1's claim-time bounds
   [-5,5] would give [-6, 14] instead. */
static char *test_dton_chain_bounds_live()
{
    double Ax[] = {1, 1, 2, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 2, 3};
    int Ap[] = {0, 2, 4, 6};
    double lhs[] = {1, 4, -INF};
    double rhs[] = {1, 4, 10};
    double lbs[] = {0, -5, -5, 0};
    double ubs[] = {1, 5, 5, 10};
    double c[] = {0, 0, 0, 0};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 4, 6, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;

    remove_dton_eq_rows(prob);

    // intermediate column x1 tightened by the deeper link (deactivation
    // does not reset bounds)
    mu_assert("x1 lb", constraints->bounds[1].lb == 0.0);
    mu_assert("x1 ub", constraints->bounds[1].ub == 1.0);

    // survivor x2 tightened from x1's LIVE bounds
    mu_assert("x2 lb", constraints->bounds[2].lb == 2.0);
    mu_assert("x2 ub", constraints->bounds[2].ub == 4.0);

    mu_assert("x0 inactive", HAS_TAG(constraints->col_tags[0], C_TAG_INACTIVE));
    mu_assert("x1 inactive", HAS_TAG(constraints->col_tags[1], C_TAG_INACTIVE));

    DEBUG(run_debugger(constraints, false));

    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* Tag-liveness variant: x1 starts with an INFINITE lower bound. The
   deeper link's update_lb(x1, 0) clears C_TAG_LB_INF; the shallower
   link's x2 upper-bound update happens only if it reads x1's tag LIVE
   (a stale claim-time snapshot tag would skip the branch entirely). */
static char *test_dton_chain_bounds_live_tags()
{
    double Ax[] = {1, 1, 2, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 2, 3};
    int Ap[] = {0, 2, 4, 6};
    double lhs[] = {1, 4, -INF};
    double rhs[] = {1, 4, 10};
    double lbs[] = {0, -INF, -5, 0};
    double ubs[] = {1, 5, 5, 10};
    double c[] = {0, 0, 0, 0};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 4, 6, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    Problem *prob = presolver->prob;
    Constraints *constraints = prob->constraints;

    remove_dton_eq_rows(prob);

    // deeper link: x1 gets lb 0 (clearing its inf tag) and ub 1
    mu_assert("x1 lb set", constraints->bounds[1].lb == 0.0);
    mu_assert("x1 ub", constraints->bounds[1].ub == 1.0);

    // shallower link must see the now-finite x1 lb: x2 ub = 4 - 2 * 0
    mu_assert("x2 ub from live tag", constraints->bounds[2].ub == 4.0);
    mu_assert("x2 lb", constraints->bounds[2].lb == 2.0);

    DEBUG(run_debugger(constraints, false));

    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* ------------------------------------------------------------------------
   Postsolve checks. The records are linear maps, so they are tested on a
   reduced point that satisfies reduced stationarity (z_red = c_red - A_red' y_red
   for arbitrary x_red, y_red): the postsolved point must satisfy stationarity
   on the original exactly and the removed equality rows must hold exactly.
   Optimality on real instances is covered by the map_checker sweep.
   ------------------------------------------------------------------------ */
#define DTON_POSTSOLVE_TOL 1e-9

static double dton_dual_residual(const double *Ax, const int *Ai, const int *Ap,
                                 int m, int n, const double *c, const double *y,
                                 const double *z)
{
    double *r = (double *) calloc((size_t) n, sizeof(double));
    for (int j = 0; j < n; ++j)
    {
        r[j] = c[j] - z[j];
    }
    for (int i = 0; i < m; ++i)
    {
        for (int p = Ap[i]; p < Ap[i + 1]; ++p)
        {
            r[Ai[p]] -= Ax[p] * y[i];
        }
    }
    double mx = 0.0;
    for (int j = 0; j < n; ++j)
    {
        mx = fmax(mx, fabs(r[j]));
    }
    free(r);
    return mx;
}

/* max |a_i x - rhs_i| over the equality rows that the presolve removed */
static double dton_removed_eq_residual(const double *Ax, const int *Ai,
                                       const int *Ap, int m, const double *lhs,
                                       const double *rhs, const int *rows_map,
                                       const double *x)
{
    double mx = 0.0;
    for (int i = 0; i < m; ++i)
    {
        if (rows_map[i] != -1 || lhs[i] != rhs[i])
        {
            continue;
        }
        double act = 0.0;
        for (int p = Ap[i]; p < Ap[i + 1]; ++p)
        {
            act += Ax[p] * x[Ai[p]];
        }
        mx = fmax(mx, fabs(act - rhs[i]));
    }
    return mx;
}

static void dton_stationary_reduced_point(const PresolvedProblem *red, double *x,
                                          double *y, double *z)
{
    for (size_t j = 0; j < red->n; ++j)
    {
        x[j] = 0.3 * (double) (j + 1);
        z[j] = red->c[j];
    }
    for (size_t i = 0; i < red->m; ++i)
    {
        y[i] = ((i % 2 == 0) ? 0.7 : -0.4) * (double) (i + 1);
        for (int p = red->Ap[i]; p < red->Ap[i + 1]; ++p)
        {
            z[red->Ai[p]] -= red->Ax[p] * y[i];
        }
    }
}

/* Presolves (all other explorers as set by 'all_explorers'), postsolves a
   stationary reduced point and checks the original point: the dual must stay
   stationary and every removed equality row must hold. Returns 0 or the
   failing check. */
static char *dton_check_postsolve(double *Ax, int *Ai, int *Ap, int m, int n,
                                  int nnz, double *lhs, double *rhs, double *lbs,
                                  double *ubs, double *c, bool all_explorers,
                                  double dirty_frac)
{
    Settings *stgs = default_settings(); // the presolver keeps a pointer to it
    if (all_explorers)
    {
        set_settings_true(stgs);
        stgs->parallel_cols = false;
    }
    else
    {
        set_settings_false(stgs);
        stgs->dton_eq = true;
    }
    stgs->verbose = false;
    Presolver *ps =
        new_presolver(Ax, Ai, Ap, m, n, nnz, lhs, rhs, lbs, ubs, c, stgs);
    if (ps == NULL)
    {
        return "presolver allocation failed";
    }
    // steer the transpose update: 1e9 never rebuilds (merge path), 0.0
    // always rebuilds (fallback path)
    Work *work = ps->prob->constraints->state->work;
    work->dton->rebuild_dirty_frac = dirty_frac;
    if (run_presolver(ps) != REDUCED)
    {
        return "presolve must reduce and stay feasible";
    }

    // the test must have exercised the elimination
    const PresolvedProblem *red = ps->reduced_prob;
    const u16Vec *types = ps->prob->constraints->state->postsolve_info->type;
    int n_dton = 0;
    for (size_t t = 0; t < types->len; ++t)
    {
        n_dton += (types->data[t] == SUB_COL_DTON);
    }
    if (n_dton == 0)
    {
        return "no SUB_COL_DTON record emitted";
    }

    double *x = (double *) calloc(red->n + 1, sizeof(double));
    double *y = (double *) calloc(red->m + 1, sizeof(double));
    double *z = (double *) calloc(red->n + 1, sizeof(double));
    dton_stationary_reduced_point(red, x, y, z);

    char *msg = 0;
    postsolve(ps, x, y, z);
    const Solution *sol = ps->sol;
    const int *rows_map = ps->prob->constraints->state->work->mappings->rows;
    if (dton_dual_residual(Ax, Ai, Ap, m, n, c, sol->y, sol->z) > DTON_POSTSOLVE_TOL)
    {
        msg = "postsolve breaks stationarity on the original";
    }
    else if (dton_removed_eq_residual(Ax, Ai, Ap, m, lhs, rhs, rows_map, sol->x) >
             DTON_POSTSOLVE_TOL)
    {
        msg = "postsolve violates a removed equality row";
    }

    free(x);
    free(y);
    free(z);
    free_presolver(ps);
    PS_FREE(stgs);
    return msg;
}

/* Single link, no chain: x0 = 2 - 0.5 x1 substituted into two rows. */
static char *test_dton_postsolve_single_link()
{
    // r0: 2 x0 + x1 = 4, r1: x0 + x1 + x2 <= 10, r2: 3 x0 + x2 + x3 >= 1
    double Ax[] = {2, 1, 1, 1, 1, 3, 1, 1};
    int Ai[] = {0, 1, 0, 1, 2, 0, 2, 3};
    int Ap[] = {0, 2, 5, 8};
    double lhs[] = {4, -INF, 1};
    double rhs[] = {4, 10, INF};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {1, 2, 3, 4};
    return dton_check_postsolve(Ax, Ai, Ap, 3, 4, 8, lhs, rhs, lbs, ubs, c, false,
                                1e9);
}

/* Depth-2 chain x0 -> x1 -> x2 with both eliminated columns in a kept row,
   so the owner-row multipliers need the deeper owner row and the kept row. */
static char *test_dton_postsolve_chain()
{
    // r0: x0 + x1 = 1, r1: 2 x1 + x2 = 4, r2: 2 x0 + x1 + x2 + x3 <= 10
    double Ax[] = {1, 1, 2, 1, 2, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 0, 1, 2, 3};
    int Ap[] = {0, 2, 4, 8};
    double lhs[] = {1, 4, -INF};
    double rhs[] = {1, 4, 10};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {1, 2, 3, 4};
    return dton_check_postsolve(Ax, Ai, Ap, 3, 4, 8, lhs, rhs, lbs, ubs, c, false,
                                1e9);
}

/* Tree: x0 -> x2 and x1 -> x2 (two heads), x2 -> x3; the column of x2 holds
   two owner rows of the same round. */
static char *test_dton_postsolve_tree()
{
    // r0: x0 + x2 = 1, r1: x1 + x2 = 1, r2: 2 x2 + x3 = 4,
    // r3: x0 + 2 x1 + 3 x2 + x3 <= 10
    double Ax[] = {1, 1, 1, 1, 2, 1, 1, 2, 3, 1};
    int Ai[] = {0, 2, 1, 2, 2, 3, 0, 1, 2, 3};
    int Ap[] = {0, 2, 4, 6, 10};
    double lhs[] = {1, 1, 4, -INF};
    double rhs[] = {1, 1, 4, 10};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {1, 2, 3, 4};
    return dton_check_postsolve(Ax, Ai, Ap, 4, 4, 10, lhs, rhs, lbs, ubs, c, false,
                                1e9);
}

/* Two rounds: r1 loses the claim on x0 in round one and is eliminated in
   round two after r0's substitution rewrote it. */
static char *test_dton_postsolve_two_rounds()
{
    // r0: 3 x0 + x1 = 1, r1: 5 x0 + x2 = 2, r2: x1 + x2 <= 10
    double Ax[] = {3, 1, 5, 1, 1, 1};
    int Ai[] = {0, 1, 0, 2, 1, 2};
    int Ap[] = {0, 2, 4, 6};
    double lhs[] = {1, 2, -INF};
    double rhs[] = {1, 2, 10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {1, 2, 3};
    return dton_check_postsolve(Ax, Ai, Ap, 3, 3, 6, lhs, rhs, lbs, ubs, c, false,
                                1e9);
}

/* Cycle: the broken owner row stays active, becomes a singleton row after
   the other substitution and is handled by the trivial explorers. */
static char *test_dton_postsolve_cycle()
{
    // r0: 2 x0 + x1 = 0, r1: x0 + 2 x1 = 0, r2: x0 + x1 + x2 <= 10
    double Ax[] = {2, 1, 1, 2, 1, 1, 1};
    int Ai[] = {0, 1, 0, 1, 0, 1, 2};
    int Ap[] = {0, 2, 4, 7};
    double lhs[] = {0, 0, -INF};
    double rhs[] = {0, 0, 10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {10, 10, 10};
    double c[] = {1, 2, 3};
    return dton_check_postsolve(Ax, Ai, Ap, 3, 3, 7, lhs, rhs, lbs, ubs, c, false,
                                1e9);
}

/* Exact cancellation of a substituted entry (dropped below ZERO_TOL). */
static char *test_dton_postsolve_cancellation()
{
    // r0: x0 + 2 x1 = 0, r1: 0.5 x0 + x1 + x2 <= 5, r2: x1 + x3 + x4 <= 7
    double Ax[] = {1, 2, 0.5, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 1, 2, 1, 3, 4};
    int Ap[] = {0, 2, 5, 8};
    double lhs[] = {0, -INF, -INF};
    double rhs[] = {0, 5, 7};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10, 10};
    double c[] = {1, 2, 3, 4, 5};
    return dton_check_postsolve(Ax, Ai, Ap, 3, 5, 8, lhs, rhs, lbs, ubs, c, false,
                                1e9);
}

/* The substituted column is free: nothing is transferred onto the stay
   column, whose bounds survive unchanged; the postsolve round trip holds. */
static char *test_dton_free_subst()
{
    // r0: x0 + x1 = 2 (x1 free and sparser -> substituted), r1: x0 + x2 <= 5,
    // r2: x0 + 2 x1 + x3 <= 7 (no cancellation: -x0 + x3 <= 3 after)
    double Ax[] = {1, 1, 1, 1, 1, 2, 1};
    int Ai[] = {0, 1, 0, 2, 0, 1, 3};
    int Ap[] = {0, 2, 4, 7};
    double lhs[] = {2, -INF, -INF};
    double rhs[] = {2, 5, 7};
    double lbs[] = {0, -INF, 0, 0};
    double ubs[] = {10, INF, 10, 10};
    double c[] = {1, 2, 3, 4};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 4, 7, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);
    mu_assert("must reduce", run_presolver(presolver) == REDUCED);
    const PresolvedProblem *red = presolver->reduced_prob;
    mu_assert("x1 eliminated", red->n == 3);
    mu_assert("stay lb unchanged", red->lbs[0] == 0.0);
    mu_assert("stay ub unchanged", red->ubs[0] == 10.0);
    free_presolver(presolver);
    PS_FREE(stgs);

    return dton_check_postsolve(Ax, Ai, Ap, 3, 4, 7, lhs, rhs, lbs, ubs, c, false,
                                1e9);
}

/* Primal infeasibility ray through a substitution. r0 eliminates
   x0 = 1 - x1 (x0 is a column singleton in r0? no: it also sits in r1, so the
   sparser rule picks it); the reduced rows r1: x2 + x3 >= 4 and
   r2: x1 + x2 + x3 <= 2 are infeasible together. */
static char *test_dton_primal_ray()
{
    // r0: x0 + x1 = 1, r1: x0 + x1 + x2 + x3 >= 5, r2: x1 + x2 + x3 <= 2
    double Ax[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 1, 2, 3, 1, 2, 3};
    int Ap[] = {0, 2, 6, 9};
    int m = 3, n = 4, nnz = 9;
    double lhs[] = {1, 5, -INF};
    double rhs[] = {1, INF, 2};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {1, 1, 1, 1};
    PresolvedProblem orig =
        problem_from_csr(Ax, Ai, Ap, m, n, nnz, lhs, rhs, lbs, ubs, c);

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, m, n, nnz, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);
    mu_assert("presolve reduces", run_presolver(presolver) == REDUCED);
    mu_assert("reduced dims",
              presolver->reduced_prob->m == 2 && presolver->reduced_prob->n == 3);

    // reduced Farkas certificate: y >= 0 on an rhs row, <= 0 on an lhs row
    double y[] = {-1.0, 1.0};
    double y_orig[] = {0, 0, 0};
    postsolve_primal_infeas_ray(presolver, y, y_orig);
    mu_assert("original primal ray certificate",
              is_primal_ray_certificate(&orig, y_orig, 1e-9));
    // y_0 makes A'y + z = 0 for the eliminated (free) column x0
    mu_assert("owner row multiplier", fabs(y_orig[0] - 1.0) < 1e-12);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* Dual infeasibility ray through a substitution: x0 = x1, minimise
   -x0 with x1 - x2 <= 10 and x >= 0 is unbounded along (1, 1, 1). */
static char *test_dton_dual_ray()
{
    // r0: x0 - x1 = 0, r1: x1 - x2 <= 10
    double Ax[] = {1, -1, 1, -1};
    int Ai[] = {0, 1, 1, 2};
    int Ap[] = {0, 2, 4};
    int m = 2, n = 3, nnz = 4;
    double lhs[] = {0, -INF};
    double rhs[] = {0, 10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {INF, INF, INF};
    double c[] = {-1, 0, 0};
    PresolvedProblem orig =
        problem_from_csr(Ax, Ai, Ap, m, n, nnz, lhs, rhs, lbs, ubs, c);

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, m, n, nnz, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);
    mu_assert("presolve reduces", run_presolver(presolver) == REDUCED);
    mu_assert("reduced dims",
              presolver->reduced_prob->m == 1 && presolver->reduced_prob->n == 2);

    double d[] = {1.0, 1.0};
    double d_orig[] = {0, 0, 0};
    postsolve_dual_infeas_ray(presolver, d, d_orig);
    mu_assert("original dual ray certificate",
              is_dual_ray_certificate(&orig, d_orig, 1e-9));
    mu_assert("eliminated component follows the survivor",
              fabs(d_orig[0] - 1.0) < 1e-12);

    PS_FREE(stgs);
    free_presolver(presolver);
    return 0;
}

/* ------------------------------------------------------------------------
   Transpose update (phase 5). The transpose is compact, so a target whose
   length does not grow is merged in place and one that grows is moved to the tail.
   Each test sets the workspace's rebuild_dirty_frac to 1e9 (small LPs
   always trip the default quarter-of-nnz guard) and ends with the debugger's
   comparison against a fresh transpose.
   ------------------------------------------------------------------------ */
static Presolver *dton_new_presolver(double *Ax, int *Ai, int *Ap, int m, int n,
                                     int nnz, double *lhs, double *rhs, double *lbs,
                                     double *ubs, double *c, Settings **stgs_out,
                                     double dirty_frac)
{
    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, m, n, nnz, lhs, rhs, lbs, ubs, c, stgs);
    if (presolver != NULL)
    {
        Work *work = presolver->prob->constraints->state->work;
        work->dton->rebuild_dirty_frac = dirty_frac;
    }
    *stgs_out = stgs;
    return presolver;
}

static bool dton_col_is(const Matrix *AT, int c, const int *rows, const double *vals,
                        int len)
{
    if (AT->p[c].end - AT->p[c].start != len)
    {
        return false;
    }
    for (int q = 0; q < len; ++q)
    {
        if (AT->i[AT->p[c].start + q] != rows[q] ||
            AT->x[AT->p[c].start + q] != vals[q])
        {
            return false;
        }
    }
    return true;
}

static bool dton_spans_ordered(const Matrix *AT)
{
    for (int c = 0; c < (int) AT->m; ++c)
    {
        if (AT->p[c + 1].start < AT->p[c].end)
        {
            return false;
        }
    }
    return true;
}

/* Same substitution with four fill-in rows: x0's column {r0, r5} grows to
   5 > 2 entries, so it moves to the tail reserved by the initial transpose. */
static char *test_dton_merge_relocate()
{
    // r0: x0 + 2 x1 = 1, r1..r4: x1 + x_{2..5} <= 5, r5: x0 + x6 <= 5
    double Ax[] = {1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 1, 3, 1, 4, 1, 5, 0, 6};
    int Ap[] = {0, 2, 4, 6, 8, 10, 12};
    double lhs[] = {1, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {1, 5, 5, 5, 5, 5};
    double lbs[] = {0, 0, 0, 0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10, 10, 10, 10};
    double c[] = {1, 1, 1, 1, 1, 1, 1};
    Settings *stgs;
    Presolver *ps =
        dton_new_presolver(Ax, Ai, Ap, 6, 7, 12, lhs, rhs, lbs, ubs, c, &stgs, 1e9);
    mu_assert("presolver allocation failed", ps != NULL);
    Constraints *constraints = ps->prob->constraints;
    const Matrix *AT = constraints->AT;
    DtonWorkspace *dton_work = constraints->state->work->dton;
    size_t old_alloc = AT->n_alloc;
    int old_start_x2 = AT->p[2].start;

    remove_dton_eq_rows(ps->prob);

    mu_assert("merge path taken", !dton_work->last_round_rebuilt);
    mu_assert("x1 span empty", AT->p[1].end == AT->p[1].start);
    mu_assert("x0 moved to the tail",
              AT->p[0].start >= dton_work->at.tail_base &&
                  dton_work->at.tail_base == AT->p[AT->m].start);
    mu_assert("tail reserved up front, no realloc",
              AT->n_alloc == old_alloc &&
                  (size_t) dton_work->at.tail_next <= AT->n_alloc &&
                  (size_t) dton_work->at.tail_base < AT->n_alloc);
    mu_assert("relocated slot fits exactly", dton_work->at.cap[0] == 5);
    mu_assert("neighbour untouched", AT->p[2].start == old_start_x2);
    int rows[] = {1, 2, 3, 4, 5};
    double vals[] = {-0.5, -0.5, -0.5, -0.5, 1};
    mu_assert("x0 column content", dton_col_is(AT, 0, rows, vals, 5));
    mu_assert("x0 size", constraints->state->col_sizes[0] == 5);
    mu_assert("spans now out of order", !dton_spans_ordered(AT));

    DEBUG(run_debugger(constraints, false));

    // a rebuild restores the ordered compact layout and invalidates the caps
    dton_rebuild_AT(ps->prob, dton_work);
    AT = constraints->AT;
    mu_assert("rebuild orders spans", dton_spans_ordered(AT));
    mu_assert("caps invalidated", !dton_work->at_valid);
    mu_assert("x0 content after rebuild", dton_col_is(AT, 0, rows, vals, 5));
    DEBUG(run_debugger(constraints, false));

    free_presolver(ps);
    PS_FREE(stgs);
    return 0;
}

/* Shrink and cancellation: x0's column {r0, r1, r2} loses the owner row r0
   and the r1 entry that cancels (2 x0 - 2 x0), keeping only r2. */
static char *test_dton_merge_shrink()
{
    // r0: x0 + 2 x1 = 1, r1: 0.5 x0 + x1 + x2 <= 5, r2: x0 + x3 <= 3
    double Ax[] = {1, 2, 0.5, 1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 1, 2, 0, 3};
    int Ap[] = {0, 2, 5, 7};
    double lhs[] = {1, -INF, -INF};
    double rhs[] = {1, 5, 3};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {1, 1, 1, 1};
    Settings *stgs;
    Presolver *ps =
        dton_new_presolver(Ax, Ai, Ap, 3, 4, 7, lhs, rhs, lbs, ubs, c, &stgs, 1e9);
    mu_assert("presolver allocation failed", ps != NULL);
    Constraints *constraints = ps->prob->constraints;
    const Matrix *AT = constraints->AT;
    DtonWorkspace *dton_work = constraints->state->work->dton;
    int old_start = AT->p[0].start;

    remove_dton_eq_rows(ps->prob);

    mu_assert("merge path taken", !dton_work->last_round_rebuilt);
    mu_assert("x0 stayed in place", AT->p[0].start == old_start);
    int rows[] = {2};
    double vals[] = {1};
    mu_assert("x0 column content", dton_col_is(AT, 0, rows, vals, 1));
    mu_assert("x0 size", constraints->state->col_sizes[0] == 1);
    mu_assert("x0 now a singleton column",
              iVec_contains(constraints->state->ston_cols, 0));
    mu_assert("r1 shrank to x2", constraints->state->row_sizes[1] == 1);
    DEBUG(run_debugger(constraints, false));

    free_presolver(ps);
    PS_FREE(stgs);
    return 0;
}

/* The relocation LP through the full presolver (all explorers) with a
   postsolve round trip: later explorers must cope with a relocated column. */
static char *test_dton_merge_full_presolve()
{
    // the relocation LP plus a dense row r6 so that no column is a singleton
    // (otherwise the singleton-column pass removes everything before the
    // doubleton code runs): r6: x0 + ... + x6 <= 20
    double Ax[] = {2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 1, 3, 1, 4, 1, 5, 0, 6, 0, 1, 2, 3, 4, 5, 6};
    int Ap[] = {0, 2, 4, 6, 8, 10, 12, 19};
    double lhs[] = {1, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {1, 5, 5, 5, 5, 5, 20};
    double lbs[] = {0, 0, 0, 0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10, 10, 10, 10};
    double c[] = {1, 2, 3, 4, 5, 6, 7};
    return dton_check_postsolve(Ax, Ai, Ap, 7, 7, 19, lhs, rhs, lbs, ubs, c, true,
                                1e9);
}

/* The chain LP through the rebuild fallback on every round. */
static char *test_dton_postsolve_chain_rebuild_path()
{
    double Ax[] = {1, 1, 2, 1, 2, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 0, 1, 2, 3};
    int Ap[] = {0, 2, 4, 8};
    double lhs[] = {1, 4, -INF};
    double rhs[] = {1, 4, 10};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10};
    double c[] = {1, 2, 3, 4};
    return dton_check_postsolve(Ax, Ai, Ap, 3, 4, 8, lhs, rhs, lbs, ubs, c, false,
                                0.0);
}

/* More claims in one round than the initial record capacity: 300 independent
   doubleton rows x_{2i} + 2 x_{2i+1} = 1 grow the claim, target and log arrays. */
static char *test_dton_record_growth()
{
    enum
    {
        N = 300
    };
    double Ax[2 * N];
    int Ai[2 * N];
    int Ap[N + 1];
    double lhs[N], rhs[N], lbs[2 * N], ubs[2 * N], c[2 * N];
    for (int i = 0; i < N; ++i)
    {
        Ax[2 * i] = 1;
        Ax[2 * i + 1] = 2;
        Ai[2 * i] = 2 * i;
        Ai[2 * i + 1] = 2 * i + 1;
        Ap[i] = 2 * i;
        lhs[i] = rhs[i] = 1;
    }
    Ap[N] = 2 * N;
    for (int j = 0; j < 2 * N; ++j)
    {
        lbs[j] = 0;
        ubs[j] = 10;
        c[j] = 1;
    }
    Settings *stgs;
    Presolver *ps = dton_new_presolver(Ax, Ai, Ap, N, 2 * N, 2 * N, lhs, rhs, lbs,
                                       ubs, c, &stgs, 1e9);
    mu_assert("presolver allocation failed", ps != NULL);
    Constraints *constraints = ps->prob->constraints;
    DtonWorkspace *dton_work = constraints->state->work->dton;
    mu_assert("initial capacity below the round", dton_work->substs.cap < N);

    remove_dton_eq_rows(ps->prob);

    mu_assert("records grew", dton_work->substs.cap >= N);
    for (int i = 0; i < N; ++i)
    {
        mu_assert("every row eliminated",
                  HAS_TAG(constraints->row_tags[i], R_TAG_INACTIVE));
        mu_assert("odd columns eliminated",
                  HAS_TAG(constraints->col_tags[2 * i + 1], C_TAG_INACTIVE));
        mu_assert("claim state reset", dton_work->substs.col_subst[2 * i + 1] == -1);
    }
    mu_assert("matrix empty", constraints->A->nnz == 0);
    DEBUG(run_debugger(constraints, false));

    free_presolver(ps);
    PS_FREE(stgs);
    return 0;
}

/* ------------------------------------------------------------------------
   Long row. r0 spans columns [0, L) minus column 5, and doubleton rows
   substitute sources sitting in r0, covering every placement of a target
   relative to its sources:
     A: x110 + x120 = 1     -> x110 = 1 - x120       (sparser column; target
                                                      x120 present in r0)
     B: x120 + 2 x130 = 1   -> x130 = (1 - x120) / 2 (larger coefficient;
        target x120 again, so it sits between its two sources)
     C: x5 + 2 x150 = 1     -> x150 = (1 - x5) / 2   (target x5 absent from
                                                      r0: inserted)
     D: x30 + 2 x_{L+1} = 1 -> x_{L+1} eliminated    (not in r0; target x30
                                                      present without a source)
   Row r_last: x120 + x_L <= 5 keeps x120 in a non-owner row.
   ------------------------------------------------------------------------ */
typedef struct
{
    int L, n, m, nnz;
    double *Ax, *lhs, *rhs, *lbs, *ubs, *c;
    int *Ai, *Ap;
} DtonLongLP;

static void dton_long_lp_build(DtonLongLP *lp, int L, double coef120)
{
    lp->L = L;
    lp->n = L + 2;
    lp->m = 1 + 4 + 1;
    lp->nnz = (L - 1) + 2 * 4 + 2;
    lp->Ax = (double *) calloc((size_t) lp->nnz, sizeof(double));
    lp->Ai = (int *) calloc((size_t) lp->nnz, sizeof(int));
    lp->Ap = (int *) calloc((size_t) lp->m + 1, sizeof(int));
    lp->lhs = (double *) calloc((size_t) lp->m, sizeof(double));
    lp->rhs = (double *) calloc((size_t) lp->m, sizeof(double));
    lp->lbs = (double *) calloc((size_t) lp->n, sizeof(double));
    lp->ubs = (double *) calloc((size_t) lp->n, sizeof(double));
    lp->c = (double *) calloc((size_t) lp->n, sizeof(double));
    int p = 0, r = 0;
    // r0: every column but 5, coefficient 1 except x120 (coef120)
    lp->Ap[r] = p;
    for (int j = 0; j < L; ++j)
    {
        if (j == 5) continue;
        lp->Ai[p] = j;
        lp->Ax[p] = (j == 120) ? coef120 : 1.0;
        p++;
    }
    lp->lhs[r] = -INF;
    lp->rhs[r] = 1000;
    r++;
    int dton[4][3] = {{110, 120, 1}, {120, 130, 2}, {5, 150, 2}, {30, L + 1, 2}};
    for (int d = 0; d < 4; ++d)
    {
        lp->Ap[r] = p;
        lp->Ai[p] = dton[d][0];
        lp->Ax[p] = 1.0;
        p++;
        lp->Ai[p] = dton[d][1];
        lp->Ax[p] = dton[d][2];
        p++;
        lp->lhs[r] = lp->rhs[r] = 1;
        r++;
    }
    // r_last: x120 + x_L <= 5
    lp->Ap[r] = p;
    lp->Ai[p] = 120;
    lp->Ax[p] = 1;
    p++;
    lp->Ai[p] = L;
    lp->Ax[p] = 1;
    p++;
    lp->lhs[r] = -INF;
    lp->rhs[r] = 5;
    r++;
    lp->Ap[r] = p;
    assert(r == lp->m && p == lp->nnz);
    for (int j = 0; j < lp->n; ++j)
    {
        lp->lbs[j] = 0;
        lp->ubs[j] = 10;
        lp->c[j] = 1 + (j % 3);
    }
}

static void dton_long_lp_free(DtonLongLP *lp)
{
    free(lp->Ax);
    free(lp->Ai);
    free(lp->Ap);
    free(lp->lhs);
    free(lp->rhs);
    free(lp->lbs);
    free(lp->ubs);
    free(lp->c);
}

/* Coefficient of column j in row q, or 0 when absent. */
static double dton_coeff(const Constraints *cs, int q, int j)
{
    const Matrix *A = cs->A;
    int pos = sorted_find(A->i + A->p[q].start, cs->state->row_sizes[q], j);
    return pos < 0 ? 0.0 : A->x[A->p[q].start + pos];
}

/* r0's x120 entry (coefficient 1.5) cancels against its two sources
   (x110 = 1 - x120 contributes -1, x130 = (1 - x120)/2 contributes -0.5): the
   target is dropped from r0 and the merge applies the delete tuple. */
static char *test_dton_cancellation_log()
{
    DtonLongLP lp;
    dton_long_lp_build(&lp, 200, 1.5);
    Settings *stgs;
    Presolver *ps =
        dton_new_presolver(lp.Ax, lp.Ai, lp.Ap, lp.m, lp.n, lp.nnz, lp.lhs, lp.rhs,
                           lp.lbs, lp.ubs, lp.c, &stgs, 1e9);
    mu_assert("presolver allocation failed", ps != NULL);
    remove_dton_eq_rows(ps->prob);
    const Constraints *cs = ps->prob->constraints;
    mu_assert("merge path taken", !cs->state->work->dton->last_round_rebuilt);
    mu_assert("x120 dropped from r0", dton_coeff(cs, 0, 120) == 0.0);
    // 199 - 3 eliminated (x110, x130, x150) + x5 inserted - x120 dropped
    mu_assert("r0 length", cs->state->row_sizes[0] == 196);
    mu_assert("x120 column is r_last only",
              cs->state->col_sizes[120] == 1 &&
                  cs->AT->i[cs->AT->p[120].start] == lp.m - 1);
    DEBUG(run_debugger(cs, false));
    free_presolver(ps);
    PS_FREE(stgs);
    dton_long_lp_free(&lp);
    return 0;
}

/* ------------------------------------------------------------------------
   Infeasible bound transfers. A transferred bound that contradicts the stay
   column's bounds means the row and the two columns' bounds are
   inconsistent: the eliminator must report INFEASIBLE instead of dropping
   the substituted column's bound with the row.
   ------------------------------------------------------------------------ */

/* One row: x0 + x1 = 1 with x0 in [0, 0.2] forces x1 >= 0.8, but x1 <= 0.5.
   The eliminator reports INFEASIBLE and leaves A untouched. */
static char *test_dton_infeasible_transfer()
{
    // r0: x0 + x1 = 1, r1: x0 + x2 <= 10, r2: x1 + x2 <= 10 (pads)
    double Ax[] = {1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 0, 2, 1, 2};
    int Ap[] = {0, 2, 4, 6};
    double lhs[] = {1, -INF, -INF};
    double rhs[] = {1, 10, 10};
    double lbs[] = {0, 0, 0};
    double ubs[] = {0.2, 0.5, 10};
    double c[] = {1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 3, 6, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);
    Constraints *constraints = presolver->prob->constraints;

    mu_assert("infeasible", remove_dton_eq_rows(presolver->prob) == INFEASIBLE);
    mu_assert("r0 still active", !HAS_TAG(constraints->row_tags[0], R_TAG_INACTIVE));
    mu_assert("x0 still active", !HAS_TAG(constraints->col_tags[0], C_TAG_INACTIVE));
    mu_assert("A untouched", constraints->A->nnz == 6);

    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* Two rows onto one column, each feasible on its own: x0 + x1 = 1 with
   x0 in [0, 0.4] gives x1 in [0.6, 1]; x3 + x1 = 1 with x3 in [0.5, 1] gives
   x1 in [0, 0.5]. Neither row is activity-infeasible, so only the transfer
   can see it. Through the presolver, with the eliminator alone and with the
   default explorers minus the singleton-column one (which happened to mask
   the wrong answer on this LP). */
static char *test_dton_infeasible_two_rows()
{
    // r0: x0 + x1 = 1, r1: x1 + x3 = 1, r2: x0 + x2 <= 10, r3: x2 + x3 <= 10
    double Ax[] = {1, 1, 1, 1, 1, 1, 1, 1};
    int Ai[] = {0, 1, 1, 3, 0, 2, 2, 3};
    int Ap[] = {0, 2, 4, 6, 8};
    double lhs[] = {1, 1, -INF, -INF};
    double rhs[] = {1, 1, 10, 10};
    double lbs[] = {0, 0, 0, 0.5};
    double ubs[] = {0.4, 1, 10, 1};
    double c[] = {1, 1, 1, 1};

    for (int variant = 0; variant < 2; ++variant)
    {
        Settings *stgs = default_settings();
        stgs->verbose = false;
        if (variant == 0)
        {
            set_settings_false(stgs);
            stgs->dton_eq = true;
        }
        else
        {
            stgs->ston_cols = false;
        }
        Presolver *presolver =
            new_presolver(Ax, Ai, Ap, 4, 4, 8, lhs, rhs, lbs, ubs, c, stgs);
        mu_assert("presolver allocation failed", presolver != NULL);
        mu_assert("infeasible", run_presolver(presolver) == INFEASIBLE);
        free_presolver(presolver);
        PS_FREE(stgs);
    }
    return 0;
}

/* The contradiction appears at the second link of a chain: x0 -> x1 gives
   x1 in [0.6, 1], then x1 -> x2 through 2 x1 + x2 = 4 gives x2 >= 2, but
   x2 <= 1. The deeper link's tightening stays, the round is abandoned. */
static char *test_dton_infeasible_chain()
{
    // r0: x0 + x1 = 1, r1: 2 x1 + x2 = 4, r2: x2 + x3 <= 10 (pad)
    double Ax[] = {1, 1, 2, 1, 1, 1};
    int Ai[] = {0, 1, 1, 2, 2, 3};
    int Ap[] = {0, 2, 4, 6};
    double lhs[] = {1, 4, -INF};
    double rhs[] = {1, 4, 10};
    double lbs[] = {0, 0, 0, 0};
    double ubs[] = {0.4, 10, 1, 10};
    double c[] = {1, 1, 1, 1};

    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->dton_eq = true;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 3, 4, 6, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);
    Constraints *constraints = presolver->prob->constraints;

    mu_assert("infeasible", remove_dton_eq_rows(presolver->prob) == INFEASIBLE);
    mu_assert("first link transferred",
              constraints->bounds[1].lb == 0.6 && constraints->bounds[1].ub == 1.0);
    mu_assert("x2 bounds untouched",
              constraints->bounds[2].lb == 0.0 && constraints->bounds[2].ub == 1.0);
    mu_assert("A untouched", constraints->A->nnz == 6);

    free_presolver(presolver);
    PS_FREE(stgs);
    return 0;
}

/* End-to-end: the full presolver with ALL default explorers plus the
   doubleton eliminator, on the substitution/chain LP. Only coarse outcomes are
   asserted (other explorers legitimately reshape the exact reduction), plus
   a postsolve round trip on a stationary reduced point. */
static char *test_dton_e2e_full_settings()
{
    // feasible: x1 = 0, x4 = 1, x0 = 4 satisfies every row (r0/r2 force
    // x0 >= 10/3, so r1's rhs must exceed 10)
    double Ax[] = {1, 2, 3, 1, 1, 2, 7, 1, 1, 1};
    int Ai[] = {0, 1, 0, 2, 3, 0, 1, 4, 1, 3};
    int Ap[] = {0, 2, 5, 8, 10};
    double lhs[] = {4, -INF, 9, -INF};
    double rhs[] = {4, 20, 9, 10};
    double lbs[] = {0, 0, 0, 0, 0};
    double ubs[] = {10, 10, 10, 10, 10};
    double c[] = {1, 1, 1, 1, 1};

    Settings *stgs = default_settings();
    stgs->verbose = false;
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, 4, 5, 10, lhs, rhs, lbs, ubs, c, stgs);
    mu_assert("presolver allocation failed", presolver != NULL);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("presolve must reduce", status == REDUCED);
    mu_assert("reduced problem exists", presolver->reduced_prob != NULL);
    mu_assert("rows reduced", presolver->reduced_prob->m < 4);
    mu_assert("cols reduced", presolver->reduced_prob->n < 5);
    mu_assert("nnz reduced", presolver->reduced_prob->nnz < 10);

    free_presolver(presolver);
    PS_FREE(stgs);

    return dton_check_postsolve(Ax, Ai, Ap, 4, 5, 10, lhs, rhs, lbs, ubs, c, true,
                                1e9);
}

static const char *all_tests_dton()
{
    mu_run_test(test_dton_workspace, counter_dton);
    mu_run_test(test_dton_choose_subst, counter_dton);
    mu_run_test(test_dton_claim_conflict, counter_dton);
    mu_run_test(test_dton_chain_depth2, counter_dton);
    mu_run_test(test_dton_cycle_break, counter_dton);
    mu_run_test(test_dton_pivot_large, counter_dton);
    mu_run_test(test_dton_eliminate_isolated, counter_dton);
    mu_run_test(test_dton_apply_substitution, counter_dton);
    mu_run_test(test_dton_cancellation, counter_dton);
    mu_run_test(test_dton_chain_empties_matrix, counter_dton);
    mu_run_test(test_dton_chain_bounds_live, counter_dton);
    mu_run_test(test_dton_chain_bounds_live_tags, counter_dton);
    mu_run_test(test_dton_postsolve_single_link, counter_dton);
    mu_run_test(test_dton_postsolve_chain, counter_dton);
    mu_run_test(test_dton_postsolve_tree, counter_dton);
    mu_run_test(test_dton_postsolve_two_rounds, counter_dton);
    mu_run_test(test_dton_postsolve_cycle, counter_dton);
    mu_run_test(test_dton_postsolve_cancellation, counter_dton);
    mu_run_test(test_dton_free_subst, counter_dton);
    mu_run_test(test_dton_primal_ray, counter_dton);
    mu_run_test(test_dton_dual_ray, counter_dton);
    mu_run_test(test_dton_merge_relocate, counter_dton);
    mu_run_test(test_dton_merge_shrink, counter_dton);
    mu_run_test(test_dton_merge_full_presolve, counter_dton);
    mu_run_test(test_dton_postsolve_chain_rebuild_path, counter_dton);
    mu_run_test(test_dton_record_growth, counter_dton);
    mu_run_test(test_dton_cancellation_log, counter_dton);
    mu_run_test(test_dton_infeasible_transfer, counter_dton);
    mu_run_test(test_dton_infeasible_two_rows, counter_dton);
    mu_run_test(test_dton_infeasible_chain, counter_dton);
    mu_run_test(test_dton_e2e_full_settings, counter_dton);
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
