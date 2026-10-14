#ifndef TEST_FOLD_H
#define TEST_FOLD_H

#include "Bounds.h"
#include "Constraints.h"
#include "Matrix.h"
#include "PSLP_API.h"
#include "PSLP_sol.h"
#include "PSLP_stats.h"
#include "PSLP_status.h"
#include "Postsolver.h"
#include "Problem.h"
#include "State.h"
#include "Tags.h"
#include "Workspace.h"
#include "fold.h"
#include "fold_apply.h"
#include "fold_internal.h"
#include "fold_postsolve.h"
#include "fold_ref.h"
#include "glbopts.h"
#include "kkt.h"
#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int counter_fold = 0;
static char fold_errbuf[256];

/* ---- small LP container --------------------------------------------------- */

typedef struct FoldTestLP
{
    int m, n, nnz;
    int *Ap, *Ai;
    double *Ax, *lhs, *rhs, *lb, *ub, *c;
} FoldTestLP;

static FoldTestLP flp_alloc(int m, int n, int nnz)
{
    FoldTestLP p;
    p.m = m;
    p.n = n;
    p.nnz = nnz;
    p.Ap = (int *) calloc((size_t) m + 1, sizeof(int));
    p.Ai = (int *) calloc((size_t) nnz + 1, sizeof(int));
    p.Ax = (double *) calloc((size_t) nnz + 1, sizeof(double));
    p.lhs = (double *) calloc((size_t) m + 1, sizeof(double));
    p.rhs = (double *) calloc((size_t) m + 1, sizeof(double));
    p.lb = (double *) calloc((size_t) n + 1, sizeof(double));
    p.ub = (double *) calloc((size_t) n + 1, sizeof(double));
    p.c = (double *) calloc((size_t) n + 1, sizeof(double));
    return p;
}

static void flp_free(FoldTestLP *p)
{
    free(p->Ap);
    free(p->Ai);
    free(p->Ax);
    free(p->lhs);
    free(p->rhs);
    free(p->lb);
    free(p->ub);
    free(p->c);
}

/* An LP with the given matrix and zero line data. */
static FoldTestLP flp_new(int m, int n, int nnz, const int *Ap, const int *Ai,
                          const double *Ax)
{
    FoldTestLP p = flp_alloc(m, n, nnz);
    memcpy(p.Ap, Ap, ((size_t) m + 1) * sizeof(int));
    memcpy(p.Ai, Ai, (size_t) nnz * sizeof(int));
    memcpy(p.Ax, Ax, (size_t) nnz * sizeof(double));
    return p;
}

static void flp_fill_rows(FoldTestLP *p, double lhs, double rhs)
{
    for (int i = 0; i < p->m; i++)
    {
        p->lhs[i] = lhs;
        p->rhs[i] = rhs;
    }
}

static void flp_fill_cols(FoldTestLP *p, double lb, double ub, double c)
{
    for (int j = 0; j < p->n; j++)
    {
        p->lb[j] = lb;
        p->ub[j] = ub;
        p->c[j] = c;
    }
}

/* A copy of P with room for more rows, columns and entries after P's. The
   caller writes the new rows and the final Ap. */
static FoldTestLP flp_extend(const FoldTestLP *P, int rows, int cols, int nnz)
{
    FoldTestLP p = flp_alloc(P->m + rows, P->n + cols, P->nnz + nnz);
    memcpy(p.Ap, P->Ap, (size_t) P->m * sizeof(int));
    memcpy(p.Ai, P->Ai, (size_t) P->nnz * sizeof(int));
    memcpy(p.Ax, P->Ax, (size_t) P->nnz * sizeof(double));
    memcpy(p.lhs, P->lhs, (size_t) P->m * sizeof(double));
    memcpy(p.rhs, P->rhs, (size_t) P->m * sizeof(double));
    memcpy(p.lb, P->lb, (size_t) P->n * sizeof(double));
    memcpy(p.ub, P->ub, (size_t) P->n * sizeof(double));
    memcpy(p.c, P->c, (size_t) P->n * sizeof(double));
    return p;
}

/* Block-diagonal k copies of P. */
static FoldTestLP flp_replicate(const FoldTestLP *P, int k)
{
    FoldTestLP r = flp_alloc(P->m * k, P->n * k, P->nnz * k);
    int pos = 0;
    for (int b = 0; b < k; b++)
    {
        for (int i = 0; i < P->m; i++)
        {
            r.Ap[b * P->m + i] = pos;
            for (int q = P->Ap[i]; q < P->Ap[i + 1]; q++)
            {
                r.Ai[pos] = b * P->n + P->Ai[q];
                r.Ax[pos] = P->Ax[q];
                pos++;
            }
            r.lhs[b * P->m + i] = P->lhs[i];
            r.rhs[b * P->m + i] = P->rhs[i];
        }
        for (int j = 0; j < P->n; j++)
        {
            r.lb[b * P->n + j] = P->lb[j];
            r.ub[b * P->n + j] = P->ub[j];
            r.c[b * P->n + j] = P->c[j];
        }
    }
    r.Ap[r.m] = pos;
    return r;
}

static FoldOptions fold_test_options(double min_reduction)
{
    FoldOptions options;
    options.min_reduction = min_reduction;
    options.work_budget = 100.0;
    options.max_rounds = 8192;
    return options;
}

/* Runs the fold on P. */
static FoldReason run_fold(const FoldTestLP *P, const FoldOptions *options,
                           FoldMaps **maps, FoldStats *stats)
{
    Matrix *A = matrix_new(P->Ax, P->Ai, P->Ap, (size_t) P->m, (size_t) P->n,
                           (size_t) P->nnz);
    int *work = (int *) calloc((size_t) P->n + 1, sizeof(int));
    Matrix *AT = transpose(A, work, 0);
    Bound *bounds = (Bound *) malloc(((size_t) P->n + 1) * sizeof(Bound));
    for (int j = 0; j < P->n; j++)
    {
        bounds[j].lb = P->lb[j];
        bounds[j].ub = P->ub[j];
    }
    FoldReason reason =
        fold_detect(A, AT, P->lhs, P->rhs, bounds, P->c, options, maps, stats);
    free(bounds);
    free(work);
    free_matrix(AT);
    free_matrix(A);
    return reason;
}

/* NULL if the fold's partition (with the gate off) equals the reference's. The
   maps are handed back through 'maps_out' when given, else freed. */
static const char *compare_with_ref(const FoldTestLP *P, FoldReason *reason_out,
                                    FoldStats *stats, FoldMaps **maps_out)
{
    FoldOptions options = fold_test_options(0.0);
    FoldMaps *maps = NULL;
    FoldReason reason = run_fold(P, &options, &maps, stats);
    *reason_out = reason;
    if (reason != FOLD_FOLDABLE && reason != FOLD_NOTHING)
    {
        fold_maps_free(maps);
        snprintf(fold_errbuf, sizeof(fold_errbuf),
                 "unexpected reason %s (%d rounds, work %.1f)",
                 fold_reason_string(reason), stats->rounds, stats->work);
        return fold_errbuf;
    }

    int *row_ref = (int *) malloc(((size_t) P->m + 1) * sizeof(int));
    int *col_ref = (int *) malloc(((size_t) P->n + 1) * sizeof(int));
    int rounds;
    const char *err = NULL;
    if (!fold_ref(P->m, P->n, P->Ap, P->Ai, P->Ax, P->lhs, P->rhs, P->lb, P->ub,
                  P->c, 100000, row_ref, col_ref, &rounds))
    {
        err = "reference did not converge";
    }
    /* A fold that removes no entries is not made: the reference may then only
       merge lines without entries. */
    if (!err && reason == FOLD_NOTHING)
    {
        int *col_count = (int *) calloc((size_t) P->n + 1, sizeof(int));
        for (int k = 0; k < P->nnz; k++)
        {
            col_count[P->Ai[k]]++;
        }
        for (int i = 0; i < P->m && !err; i++)
        {
            bool empty = P->Ap[i + 1] == P->Ap[i];
            if (row_ref[i] != i && !empty)
            {
                err = "reference merges a row with entries, fold has none";
            }
        }
        for (int j = 0; j < P->n && !err; j++)
        {
            if (col_ref[j] != j && col_count[j] > 0)
            {
                err = "reference merges a column with entries, fold has none";
            }
        }
        free(col_count);
    }
    if (!err && reason == FOLD_FOLDABLE)
    {
        for (int i = 0; i < P->m && !err; i++)
        {
            if (maps->row_rep[i] != row_ref[i])
            {
                snprintf(fold_errbuf, sizeof(fold_errbuf),
                         "row %d: fold class %d, reference %d", i, maps->row_rep[i],
                         row_ref[i]);
                err = fold_errbuf;
            }
        }
        for (int j = 0; j < P->n && !err; j++)
        {
            if (maps->col_rep[j] != col_ref[j])
            {
                snprintf(fold_errbuf, sizeof(fold_errbuf),
                         "col %d: fold class %d, reference %d", j, maps->col_rep[j],
                         col_ref[j]);
                err = fold_errbuf;
            }
        }
    }
    free(row_ref);
    free(col_ref);
    if (maps_out && !err)
    {
        *maps_out = maps;
    }
    else
    {
        fold_maps_free(maps);
    }
    return err;
}

/* Two rows, two columns, no symmetry inside the block. */
static FoldTestLP base_block(void)
{
    static const int Ap[3] = {0, 2, 4}, Ai[4] = {0, 1, 0, 1};
    static const double Ax[4] = {1, 2, 3, -1};
    FoldTestLP p = flp_new(2, 2, 4, Ap, Ai, Ax);
    p.lhs[0] = -INF;
    p.rhs[0] = 4;
    p.lhs[1] = -2;
    p.rhs[1] = INF;
    p.lb[0] = 0;
    p.ub[0] = 10;
    p.c[0] = 1;
    p.lb[1] = 0;
    p.ub[1] = 5;
    p.c[1] = 2;
    return p;
}

/* ---- tests ---------------------------------------------------------------- */

/* Rows x1 + x2 <= 1 and x1 + x2 <= 2 never share a class (different rhs); the
   two identical columns do. */
static const char *test_fold_sides(void)
{
    static const int Ap[3] = {0, 2, 4}, Ai[4] = {0, 1, 0, 1};
    static const double Ax[4] = {1, 1, 1, 1};
    FoldTestLP p = flp_new(2, 2, 4, Ap, Ai, Ax);
    flp_fill_rows(&p, -INF, 1);
    p.rhs[1] = 2;
    flp_fill_cols(&p, 0, INF, 1);
    FoldReason reason;
    FoldStats stats;
    FoldMaps *maps = NULL;
    const char *err = compare_with_ref(&p, &reason, &stats, &maps);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("expected foldable", reason == FOLD_FOLDABLE);
    mu_assert("rows should be singletons, columns one class",
              maps->n_rows_folded == 2 && maps->n_cols_folded == 1);
    fold_maps_free(maps);
    flp_free(&p);
    return 0;
}

/* Two identical blocks fold onto one. */
static const char *test_fold_symmetric_blocks(void)
{
    FoldTestLP b = base_block();
    FoldTestLP p = flp_replicate(&b, 2);
    FoldReason reason;
    FoldStats stats;
    FoldMaps *maps = NULL;
    const char *err = compare_with_ref(&p, &reason, &stats, &maps);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("expected foldable", reason == FOLD_FOLDABLE);
    mu_assert("expected 2 x 2 folded",
              maps->n_rows_folded == 2 && maps->n_cols_folded == 2);
    mu_assert("expected ids {0,1,0,1}",
              maps->row_rep[2] == 0 && maps->row_rep[3] == 1 &&
                  maps->col_rep[2] == 0 && maps->col_rep[3] == 1);
    fold_maps_free(maps);
    flp_free(&p);
    flp_free(&b);
    return 0;
}

/* The same blocks plus one row touching block 0 only: the asymmetry propagates
   and nothing folds. */
static const char *test_fold_asymmetric_blocks(void)
{
    FoldTestLP b = base_block();
    FoldTestLP r = flp_replicate(&b, 2);
    FoldTestLP p = flp_extend(&r, 1, 0, 1);
    p.Ap[r.m] = r.nnz;
    p.Ai[r.nnz] = 0;
    p.Ax[r.nnz] = 1.0;
    p.Ap[p.m] = r.nnz + 1;
    p.lhs[r.m] = -INF;
    p.rhs[r.m] = 3;
    FoldReason reason;
    FoldStats stats;
    const char *err = compare_with_ref(&p, &reason, &stats, NULL);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("expected nothing to fold", reason == FOLD_NOTHING);
    mu_assert("expected several rounds", stats.rounds >= 2);
    flp_free(&p);
    flp_free(&r);
    flp_free(&b);
    return 0;
}

/* Degenerate dimensions. */
static const char *test_fold_degenerate(void)
{
    static const int Ap[2] = {0, 1}, Ai[1] = {0};
    static const double Ax[1] = {2.0};
    FoldTestLP p = flp_new(1, 1, 1, Ap, Ai, Ax);
    p.lhs[0] = 1;
    p.rhs[0] = 1;
    p.lb[0] = 0;
    p.ub[0] = 3;
    p.c[0] = 1;
    FoldReason reason;
    FoldStats stats;
    const char *err = compare_with_ref(&p, &reason, &stats, NULL);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("expected nothing to fold", reason == FOLD_NOTHING);
    flp_free(&p);

    /* Two empty rows with equal sides and no entries: nothing to fold. */
    FoldTestLP q = flp_alloc(2, 1, 0);
    flp_fill_rows(&q, -INF, 1);
    flp_fill_cols(&q, 0, 1, 1);
    FoldOptions options = fold_test_options(0.0);
    FoldMaps *maps = NULL;
    reason = run_fold(&q, &options, &maps, &stats);
    mu_assert("expected nothing to fold for nnz == 0", reason == FOLD_NOTHING);
    mu_assert("maps empty", maps == NULL);
    flp_free(&q);

    /* No rows, or no columns: runs through and folds nothing. */
    FoldTestLP no_rows = flp_alloc(0, 2, 0);
    flp_fill_cols(&no_rows, 0, 1, 1);
    reason = run_fold(&no_rows, &options, &maps, &stats);
    mu_assert("expected nothing to fold without rows", reason == FOLD_NOTHING);
    flp_free(&no_rows);
    FoldTestLP no_cols = flp_alloc(2, 0, 0);
    flp_fill_rows(&no_cols, -INF, 1);
    reason = run_fold(&no_cols, &options, &maps, &stats);
    mu_assert("expected nothing to fold without columns", reason == FOLD_NOTHING);
    flp_free(&no_cols);
    return 0;
}

/* x_{t+1} - x_t <= 1 for t = 1..T-1 plus x_1 <= 0. The asymmetry travels one
   period per round. */
static FoldTestLP staircase(int T)
{
    FoldTestLP p = flp_alloc(T, T, 2 * (T - 1) + 1);
    int pos = 0;
    for (int idx = 0; idx < T - 1; idx++)
    {
        p.Ap[idx] = pos;
        p.Ai[pos] = idx;
        p.Ax[pos] = -1.0;
        pos++;
        p.Ai[pos] = idx + 1;
        p.Ax[pos] = 1.0;
        pos++;
        p.lhs[idx] = -INF;
        p.rhs[idx] = 1;
    }
    p.Ap[T - 1] = pos;
    p.Ai[pos] = 0;
    p.Ax[pos] = 1.0;
    pos++;
    p.lhs[T - 1] = -INF;
    p.rhs[T - 1] = 0;
    p.Ap[T] = pos;
    flp_fill_cols(&p, 0, 10, 1);
    return p;
}

/* The round limit and the work budget both stop the fold. */
static const char *test_fold_not_converged(void)
{
    FoldTestLP p = staircase(64);
    FoldOptions options = fold_test_options(0.0);
    options.max_rounds = 4;
    FoldMaps *maps = NULL;
    FoldStats stats;
    FoldReason reason = run_fold(&p, &options, &maps, &stats);
    mu_assert("expected the round limit", reason == FOLD_LIMIT);
    mu_assert("expected exactly 4 rounds", stats.rounds == 4);
    mu_assert("maps empty", maps == NULL);

    options = fold_test_options(0.0);
    options.work_budget = 0.5;
    reason = run_fold(&p, &options, &maps, &stats);
    mu_assert("expected the work budget", reason == FOLD_LIMIT);

    const char *err = compare_with_ref(&p, &reason, &stats, NULL);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("expected nothing to fold", reason == FOLD_NOTHING);
    mu_assert("expected ~T/2 rounds", stats.rounds >= 30);
    flp_free(&p);
    return 0;
}

/* Unique rows with duplicate columns, and the mirror. */
static const char *test_fold_one_sided(void)
{
    FoldTestLP p = flp_alloc(3, 2, 6);
    int pos = 0;
    for (int i = 0; i < 3; i++)
    {
        p.Ap[i] = pos;
        p.Ai[pos] = 0;
        p.Ax[pos] = 1.0;
        pos++;
        p.Ai[pos] = 1;
        p.Ax[pos] = 1.0;
        pos++;
        p.lhs[i] = -INF;
        p.rhs[i] = (double) i;
    }
    p.Ap[3] = pos;
    flp_fill_cols(&p, 0, INF, 1);
    FoldReason reason;
    FoldStats stats;
    FoldMaps *maps = NULL;
    const char *err = compare_with_ref(&p, &reason, &stats, &maps);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("expected foldable", reason == FOLD_FOLDABLE);
    mu_assert("expected 3 x 1",
              maps->n_rows_folded == 3 && maps->n_cols_folded == 1);
    fold_maps_free(maps);
    flp_free(&p);

    FoldTestLP q = flp_alloc(2, 3, 6);
    int qpos = 0;
    for (int i = 0; i < 2; i++)
    {
        q.Ap[i] = qpos;
        for (int j = 0; j < 3; j++)
        {
            q.Ai[qpos] = j;
            q.Ax[qpos] = 1.0;
            qpos++;
        }
    }
    q.Ap[2] = qpos;
    flp_fill_rows(&q, -INF, 5);
    flp_fill_cols(&q, 0, INF, 0);
    q.c[1] = 1;
    q.c[2] = 2;
    err = compare_with_ref(&q, &reason, &stats, &maps);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("expected foldable", reason == FOLD_FOLDABLE);
    mu_assert("expected 1 x 3",
              maps->n_rows_folded == 1 && maps->n_cols_folded == 3);
    fold_maps_free(maps);
    flp_free(&q);
    return 0;
}

/* Two identical blocks beside a large unique diagonal part: 8 removable
   entries out of 2008 is below the 1% gate. */
static const char *test_fold_gate(void)
{
    FoldTestLP b = base_block();
    FoldTestLP r = flp_replicate(&b, 2);
    int extra = 2000;
    FoldTestLP p = flp_extend(&r, extra, extra, extra);
    for (int ii = 0; ii < extra; ii++)
    {
        p.Ap[r.m + ii] = r.nnz + ii;
        p.Ai[r.nnz + ii] = r.n + ii;
        p.Ax[r.nnz + ii] = 1.0 + ii;
        p.lhs[r.m + ii] = -INF;
        p.rhs[r.m + ii] = ii;
        p.lb[r.n + ii] = 0;
        p.ub[r.n + ii] = INF;
        p.c[r.n + ii] = 1;
    }
    p.Ap[p.m] = p.nnz;
    FoldOptions options = fold_test_options(0.01);
    FoldMaps *maps = NULL;
    FoldStats stats;
    FoldReason reason = run_fold(&p, &options, &maps, &stats);
    mu_assert("expected nothing to fold below the gate", reason == FOLD_NOTHING);
    mu_assert("maps empty", maps == NULL);
    options.min_reduction = 0.0;
    reason = run_fold(&p, &options, &maps, &stats);
    mu_assert("expected foldable", reason == FOLD_FOLDABLE);
    fold_maps_free(maps);
    flp_free(&p);
    flp_free(&r);
    flp_free(&b);
    return 0;
}

/* ---- random instances ------------------------------------------------------ */

static uint64_t fold_rng_state = 0x9E3779B97F4A7C15ULL;

static unsigned fold_rnd(void)
{
    fold_rng_state =
        fold_rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (unsigned) (fold_rng_state >> 33);
}

/* A random block of m0 x n0 with about 2/3 density. */
static FoldTestLP random_block(int m0, int n0, int n_values, int n_sides)
{
    static const double vals[] = {1.0, 2.0, -1.0, 0.5, 3.0, -2.5};
    static const double sides[] = {0.0, 1.0, 5.0};
    FoldTestLP b = flp_alloc(m0, n0, m0 * n0);
    int pos = 0;
    for (int i = 0; i < m0; i++)
    {
        b.Ap[i] = pos;
        for (int j = 0; j < n0; j++)
        {
            if (fold_rnd() % 3 != 0)
            {
                b.Ai[pos] = j;
                b.Ax[pos] = vals[fold_rnd() % (unsigned) n_values];
                pos++;
            }
        }
        int kind = (int) (fold_rnd() % 3);
        double s = sides[fold_rnd() % (unsigned) n_sides];
        b.lhs[i] = (kind == 0) ? -INF : s;
        b.rhs[i] = (kind == 1) ? INF : s;
    }
    b.Ap[m0] = pos;
    b.nnz = pos;
    for (int j = 0; j < n0; j++)
    {
        b.lb[j] = (fold_rnd() % 2) ? 0.0 : -INF;
        b.ub[j] = (fold_rnd() % 2) ? INF : 10.0;
        b.c[j] = (double) (fold_rnd() % 3);
    }
    return b;
}

/* A random small block replicated 1-3 times, with one symmetry possibly broken.
   The draws are sequenced explicitly, since the order of evaluation of
   function arguments depends on the compiler. The perturbation cannot cancel
   an entry: explicit zeros would be dropped by the matrix but not by the
   reference. */
static FoldTestLP random_instance(void)
{
    int m0 = 2 + (int) (fold_rnd() % 6);
    int n0 = 2 + (int) (fold_rnd() % 6);
    FoldTestLP b = random_block(m0, n0, 6, 3);
    int k = 1 + (int) (fold_rnd() % 3);
    FoldTestLP p = flp_replicate(&b, k);
    if (p.nnz > 0 && (fold_rnd() % 2))
    {
        p.Ax[fold_rnd() % (unsigned) p.nnz] += 0.25;
    }
    flp_free(&b);
    return p;
}

static const char *test_fold_random_vs_reference(void)
{
    for (int trial = 0; trial < 60; trial++)
    {
        FoldTestLP p = random_instance();
        FoldReason reason;
        FoldStats stats;
        const char *err =
            (p.nnz > 0) ? compare_with_ref(&p, &reason, &stats, NULL) : NULL;
        flp_free(&p);
        if (err)
        {
            static char b2[320];
            snprintf(b2, sizeof(b2), "trial %d: %s", trial, err);
            mu_assert(b2, 0);
        }
    }
    return 0;
}

/* Large instances: a block replicated into one giant class, and a random
   sparse LP with many distinct sides. */
static const char *test_fold_large(void)
{
    FoldTestLP b = random_block(5, 6, 6, 3);
    FoldTestLP p = flp_replicate(&b, 30000);
    p.Ax[7] += 1.0;
    FoldReason reason;
    FoldStats stats;
    const char *err = compare_with_ref(&p, &reason, &stats, NULL);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("expected foldable", reason == FOLD_FOLDABLE);
    flp_free(&p);
    flp_free(&b);

    int m = 150000, n = 150000;
    FoldTestLP q = flp_alloc(m, n, 3 * m);
    int pos = 0;
    for (int i = 0; i < m; i++)
    {
        q.Ap[i] = pos;
        int j0 = (int) (fold_rnd() % (unsigned) (n - 2));
        for (int k = 0; k < 3; k++)
        {
            q.Ai[pos] = j0 + k;
            q.Ax[pos] = 1.0 + (double) (fold_rnd() % 2);
            pos++;
        }
        q.lhs[i] = -INF;
        q.rhs[i] = (double) (fold_rnd() % 50000);
    }
    q.Ap[m] = pos;
    flp_fill_cols(&q, 0, INF, 1);
    err = compare_with_ref(&q, &reason, &stats, NULL);
    mu_assert(err ? err : "", err == NULL);
    flp_free(&q);
    return 0;
}

/* ---- verification --------------------------------------------------------- */

/* Refines P through the workspace so the partition can be tampered with before
   verification. */
typedef struct FoldVerifyFixture
{
    Matrix *A;
    Matrix *AT;
    int *work;
    Bound *bounds;
    FoldWorkspace *workspace;
} FoldVerifyFixture;

static bool fixture_new(FoldVerifyFixture *fx, const FoldTestLP *P)
{
    memset(fx, 0, sizeof(*fx));
    fx->A = matrix_new(P->Ax, P->Ai, P->Ap, (size_t) P->m, (size_t) P->n,
                       (size_t) P->nnz);
    fx->work = (int *) calloc((size_t) P->n + 1, sizeof(int));
    fx->AT = transpose(fx->A, fx->work, 0);
    fx->bounds = (Bound *) malloc(((size_t) P->n + 1) * sizeof(Bound));
    for (int j = 0; j < P->n; j++)
    {
        fx->bounds[j].lb = P->lb[j];
        fx->bounds[j].ub = P->ub[j];
    }
    fx->workspace = fold_workspace_new(fx->A, fx->AT);
    if (!fx->workspace)
    {
        return false;
    }
    fold_start_colours(fx->workspace, P->m, P->n, P->lhs, P->rhs, fx->bounds, P->c);
    ColourRefinementOptions options = {100.0, 8192};
    ColourRefinementStats stats;
    return colour_refinement_run(
               fx->workspace->refine, fx->workspace->row_colours,
               fx->workspace->col_colours, &options, fx->workspace->row_classes,
               fx->workspace->col_classes, &stats) == COLOUR_REFINEMENT_CONVERGED;
}

static void fixture_free(FoldVerifyFixture *fx)
{
    fold_workspace_free(fx->workspace);
    free(fx->bounds);
    free(fx->work);
    free_matrix(fx->AT);
    free_matrix(fx->A);
}

/* The class ids, writable so that a test can tamper with the partition. */
static int *fixture_row_ids(const FoldVerifyFixture *fx)
{
    return fx->workspace->row_classes;
}

static int *fixture_col_ids(const FoldVerifyFixture *fx)
{
    return fx->workspace->col_classes;
}

/* Whether the fixture's current partition is equitable. */
static bool fixture_verify(FoldVerifyFixture *fx)
{
    return colour_refinement_verify(fx->workspace->refine,
                                    fx->workspace->row_classes,
                                    fx->workspace->col_classes);
}

/* Two identical blocks: the refined partition verifies; moving one row of
   block 0 into the class of the other row of block 0 (different entries) is
   rejected, and so is moving a column. */
static const char *test_fold_verify_entries(void)
{
    FoldTestLP b = base_block();
    FoldTestLP p = flp_replicate(&b, 2);
    FoldVerifyFixture fx;
    mu_assert("refinement should find the blocks", fixture_new(&fx, &p));
    mu_assert("refined partition should verify", fixture_verify(&fx));
    mu_assert("rows 0 and 2 should share a class",
              fixture_row_ids(&fx)[0] == fixture_row_ids(&fx)[2]);

    int saved = fixture_row_ids(&fx)[1];
    fixture_row_ids(&fx)[1] = fixture_row_ids(&fx)[0];
    mu_assert("rows with different entries must be rejected", !fixture_verify(&fx));
    fixture_row_ids(&fx)[1] = saved;
    mu_assert("restored partition should verify", fixture_verify(&fx));

    fixture_col_ids(&fx)[1] = fixture_col_ids(&fx)[0];
    mu_assert("columns with different entries must be rejected",
              !fixture_verify(&fx));
    fixture_free(&fx);
    flp_free(&p);
    flp_free(&b);
    return 0;
}

/* Rows x0 <= 1 and x1 <= 1 look alike but their columns are in different
   classes (different upper bounds), so merging the rows must be rejected on the
   class of the other end. */
static const char *test_fold_verify_other_end(void)
{
    static const int Ap[3] = {0, 1, 2}, Ai[2] = {0, 1};
    static const double Ax[2] = {1, 1};
    FoldTestLP p = flp_new(2, 2, 2, Ap, Ai, Ax);
    flp_fill_rows(&p, -INF, 1);
    flp_fill_cols(&p, 0, 10, 1);
    p.ub[1] = 11;
    FoldVerifyFixture fx;
    FoldOptions options = fold_test_options(0.0);
    FoldMaps *maps = NULL;
    FoldStats stats;
    mu_assert("expected nothing to fold",
              run_fold(&p, &options, &maps, &stats) == FOLD_NOTHING);
    fixture_new(&fx, &p);
    fixture_row_ids(&fx)[1] = fixture_row_ids(&fx)[0];
    mu_assert("rows touching different column classes must be rejected",
              !fixture_verify(&fx));
    fixture_free(&fx);
    flp_free(&p);
    return 0;
}

/* Two identical blocks, the second with -0.0 for the lower bound 0.0 of its
   columns: the fold passes the check. Changing the line data of a line that is
   not a representative makes the check fail. */
static const char *test_fold_check(void)
{
    FoldTestLP b = base_block();
    FoldTestLP p = flp_replicate(&b, 2);
    p.lb[2] = -0.0;
    p.lb[3] = -0.0;
    FoldOptions options = fold_test_options(0.0);
    FoldMaps *maps = NULL;
    FoldStats stats;
    mu_assert("expected foldable",
              run_fold(&p, &options, &maps, &stats) == FOLD_FOLDABLE);
    mu_assert("row 2 should be represented by row 0", maps->row_rep[2] == 0);
    mu_assert("column 3 should be represented by column 1", maps->col_rep[3] == 1);

    Bound bounds[4];
    for (int j = 0; j < 4; j++)
    {
        bounds[j].lb = p.lb[j];
        bounds[j].ub = p.ub[j];
    }
    mu_assert("the fold should pass the check",
              fold_check(maps, 4, 4, p.lhs, p.rhs, bounds, p.c));

    p.rhs[2] += 1.0;
    mu_assert("a row with a different rhs must fail the check",
              !fold_check(maps, 4, 4, p.lhs, p.rhs, bounds, p.c));
    p.rhs[2] -= 1.0;

    bounds[3].ub += 1.0;
    mu_assert("a column with a different ub must fail the check",
              !fold_check(maps, 4, 4, p.lhs, p.rhs, bounds, p.c));
    fold_maps_free(maps);
    flp_free(&p);
    flp_free(&b);
    return 0;
}

/* ---- apply ----------------------------------------------------------------- */

/* The folded problem of P, with the line data folded on copies of P's. */
typedef struct FoldApplyResult
{
    Matrix *A;
    Matrix *folded;
    double *lhs, *rhs, *c;
    RowTag *row_tags;
    ColTag *col_tags;
    Bound *bounds;
    FoldMaps *maps;
    FoldStats stats;
    FoldReason detect_reason;
    FoldReason reason;
} FoldApplyResult;

static void apply_result_free(FoldApplyResult *r)
{
    free_matrix(r->A);
    free_matrix(r->folded);
    free(r->lhs);
    free(r->rhs);
    free(r->c);
    free(r->row_tags);
    free(r->col_tags);
    free(r->bounds);
    fold_maps_free(r->maps);
}

static FoldApplyResult run_apply(const FoldTestLP *P)
{
    FoldApplyResult r;
    memset(&r, 0, sizeof(r));
    r.A = matrix_new(P->Ax, P->Ai, P->Ap, (size_t) P->m, (size_t) P->n,
                     (size_t) P->nnz);
    int *work = (int *) calloc((size_t) P->n + 1, sizeof(int));
    Matrix *AT = transpose(r.A, work, 0);
    r.lhs = (double *) malloc(((size_t) P->m + 1) * sizeof(double));
    r.rhs = (double *) malloc(((size_t) P->m + 1) * sizeof(double));
    r.c = (double *) malloc(((size_t) P->n + 1) * sizeof(double));
    r.row_tags = (RowTag *) calloc((size_t) P->m + 1, sizeof(RowTag));
    r.col_tags = (ColTag *) calloc((size_t) P->n + 1, sizeof(ColTag));
    r.bounds = (Bound *) malloc(((size_t) P->n + 1) * sizeof(Bound));
    memcpy(r.lhs, P->lhs, (size_t) P->m * sizeof(double));
    memcpy(r.rhs, P->rhs, (size_t) P->m * sizeof(double));
    memcpy(r.c, P->c, (size_t) P->n * sizeof(double));
    for (int i = 0; i < P->m; i++)
    {
        r.row_tags[i] = (RowTag) (i % 7);
    }
    for (int j = 0; j < P->n; j++)
    {
        r.bounds[j].lb = P->lb[j];
        r.bounds[j].ub = P->ub[j];
        r.col_tags[j] = (ColTag) (j % 5);
    }
    FoldOptions options = fold_test_options(0.0);
    r.reason = fold_detect(r.A, AT, r.lhs, r.rhs, r.bounds, r.c, &options, &r.maps,
                           &r.stats);
    r.detect_reason = r.reason;
    if (r.reason == FOLD_FOLDABLE)
    {
        r.reason = fold_apply(r.A, r.lhs, r.rhs, r.row_tags, r.bounds, r.col_tags,
                              r.c, r.maps, &r.folded);
    }
    free_matrix(AT);
    free(work);
    return r;
}

/* NULL if the folded problem equals E with costs scaled by 'cost_factor'. */
static const char *same_as_lp(const FoldApplyResult *r, const FoldTestLP *E,
                              double cost_factor)
{
    const Matrix *B = r->folded;
    if ((int) B->m != E->m || (int) B->n != E->n)
    {
        snprintf(fold_errbuf, sizeof(fold_errbuf),
                 "dims %zu x %zu, expected %d x %d", B->m, B->n, E->m, E->n);
        return fold_errbuf;
    }
    if ((int) B->nnz != E->nnz)
    {
        snprintf(fold_errbuf, sizeof(fold_errbuf), "nnz %zu, expected %d", B->nnz,
                 E->nnz);
        return fold_errbuf;
    }
    for (int i = 0; i < E->m; i++)
    {
        if (B->p[i].start != E->Ap[i] || B->p[i].end != E->Ap[i + 1])
        {
            snprintf(fold_errbuf, sizeof(fold_errbuf), "row %d range", i);
            return fold_errbuf;
        }
        for (int k = E->Ap[i]; k < E->Ap[i + 1]; k++)
        {
            if (B->i[k] != E->Ai[k] || B->x[k] != E->Ax[k])
            {
                snprintf(fold_errbuf, sizeof(fold_errbuf),
                         "row %d entry %d: (%d, %g), expected (%d, %g)", i, k,
                         B->i[k], B->x[k], E->Ai[k], E->Ax[k]);
                return fold_errbuf;
            }
        }
        if (r->lhs[i] != E->lhs[i] || r->rhs[i] != E->rhs[i])
        {
            return "row sides";
        }
    }
    for (int j = 0; j < E->n; j++)
    {
        if (r->bounds[j].lb != E->lb[j] || r->bounds[j].ub != E->ub[j])
        {
            return "bounds";
        }
        if (r->c[j] != cost_factor * E->c[j])
        {
            snprintf(fold_errbuf, sizeof(fold_errbuf), "cost %d: %g, expected %g", j,
                     r->c[j], cost_factor * E->c[j]);
            return fold_errbuf;
        }
    }
    return NULL;
}

/* Three copies of the block fold to the block with costs times three; tags
   travel with the representatives and keepers. */
static const char *test_fold_apply_blocks(void)
{
    FoldTestLP b = base_block();
    FoldTestLP p = flp_replicate(&b, 3);
    FoldApplyResult r = run_apply(&p);
    mu_assert("expected foldable", r.reason == FOLD_FOLDABLE);
    const char *err = same_as_lp(&r, &b, 3.0);
    mu_assert(err ? err : "", err == NULL);
    mu_assert("row tags of the representatives",
              r.row_tags[0] == 0 && r.row_tags[1] == 1);
    mu_assert("column tags of the keepers",
              r.col_tags[0] == 0 && r.col_tags[1] == 1);
    mu_assert("class sizes",
              r.maps->row_class_size[0] == 3 && r.maps->col_class_size[1] == 3);
    mu_assert("fold index of a member",
              r.maps->folded_row[3] == 1 && r.maps->folded_col[4] == 0);
    mu_assert("row sentinel", r.folded->p[2].start == 4 && r.folded->p[2].end == 4);
    apply_result_free(&r);
    flp_free(&p);
    flp_free(&b);
    return 0;
}

/* One row over three identical columns: the folded coefficient is the sum. */
static const char *test_fold_apply_column_sum(void)
{
    FoldTestLP p = flp_alloc(1, 3, 3);
    p.Ap[0] = 0;
    p.Ap[1] = 3;
    for (int j = 0; j < 3; j++)
    {
        p.Ai[j] = j;
        p.Ax[j] = 0.5;
    }
    flp_fill_rows(&p, -INF, 3);
    flp_fill_cols(&p, -1, 2, 0.25);
    FoldApplyResult r = run_apply(&p);
    mu_assert("expected foldable", r.reason == FOLD_FOLDABLE);
    mu_assert("1 x 1 with one entry",
              r.folded->m == 1 && r.folded->n == 1 && r.folded->nnz == 1);
    mu_assert("coefficient 1.5", r.folded->i[0] == 0 && r.folded->x[0] == 1.5);
    mu_assert("cost 0.75", r.c[0] == 0.75);
    mu_assert("bounds", r.bounds[0].lb == -1 && r.bounds[0].ub == 2);
    mu_assert("sides", r.lhs[0] == -INF && r.rhs[0] == 3);
    apply_result_free(&r);
    flp_free(&p);
    return 0;
}

/* [[1, -1], [-1, 1]] with equal data folds to one row and one column whose
   coefficient cancels exactly, so the entry is dropped. */
static const char *test_fold_apply_exact_zero(void)
{
    static const int Ap[3] = {0, 2, 4}, Ai[4] = {0, 1, 0, 1};
    static const double Ax[4] = {1, -1, -1, 1};
    FoldTestLP p = flp_new(2, 2, 4, Ap, Ai, Ax);
    flp_fill_cols(&p, 0, 1, 1);
    FoldApplyResult r = run_apply(&p);
    mu_assert("expected foldable", r.reason == FOLD_FOLDABLE);
    mu_assert("1 x 1 without entries",
              r.folded->m == 1 && r.folded->n == 1 && r.folded->nnz == 0);
    mu_assert("empty row range",
              r.folded->p[0].start == 0 && r.folded->p[0].end == 0);
    apply_result_free(&r);
    flp_free(&p);
    return 0;
}

/* Rows x1 + x2 and x0 + x1 fold to one row; columns 0 and 2 form folded
   column 0 and column 1 is folded column 1. The representative row touches
   folded column 1 before folded column 0, so its entries need sorting. */
static const char *test_fold_apply_unsorted_row(void)
{
    static const int Ap[3] = {0, 2, 4}, Ai[4] = {1, 2, 0, 1};
    static const double Ax[4] = {1, 1, 1, 1};
    FoldTestLP p = flp_new(2, 3, 4, Ap, Ai, Ax);
    flp_fill_rows(&p, -INF, 1);
    flp_fill_cols(&p, 0, 1, 1);
    p.c[1] = 2;
    FoldApplyResult r = run_apply(&p);
    mu_assert("expected foldable", r.reason == FOLD_FOLDABLE);
    mu_assert("1 x 2 with two entries",
              r.folded->m == 1 && r.folded->n == 2 && r.folded->nnz == 2);
    mu_assert("sorted columns", r.folded->i[r.folded->p[0].start] == 0 &&
                                    r.folded->i[r.folded->p[0].start + 1] == 1);
    mu_assert("coefficients", r.folded->x[r.folded->p[0].start] == 1 &&
                                  r.folded->x[r.folded->p[0].start + 1] == 1);
    apply_result_free(&r);
    flp_free(&p);
    return 0;
}

/* Random blocks without symmetry of their own, replicated: the folded problem
   is the block with the costs scaled by the number of copies. The draws are
   sequenced explicitly. */
static const char *test_fold_apply_random(void)
{
    int trials = 0;
    while (trials < 30)
    {
        int m0 = 2 + (int) (fold_rnd() % 6);
        int n0 = 2 + (int) (fold_rnd() % 6);
        FoldTestLP b = random_block(m0, n0, 6, 3);
        FoldOptions options = fold_test_options(0.0);
        FoldMaps *maps = NULL;
        FoldStats stats;
        FoldReason reason =
            b.nnz > 0 ? run_fold(&b, &options, &maps, &stats) : FOLD_FOLDABLE;
        fold_maps_free(maps);
        if (reason != FOLD_NOTHING)
        {
            flp_free(&b);
            continue;
        }
        int k = 2 + (int) (fold_rnd() % 3);
        FoldTestLP p = flp_replicate(&b, k);
        FoldApplyResult r = run_apply(&p);
        const char *err = NULL;
        if (r.reason == FOLD_FOLDABLE)
        {
            err = same_as_lp(&r, &b, (double) k);
        }
        else
        {
            snprintf(fold_errbuf, sizeof(fold_errbuf),
                     "expected foldable, got %s after %s (block %d x %d, %d nnz, %d "
                     "copies)",
                     fold_reason_string(r.reason),
                     fold_reason_string(r.detect_reason), b.m, b.n, b.nnz, k);
            err = fold_errbuf;
        }
        apply_result_free(&r);
        flp_free(&p);
        flp_free(&b);
        if (err)
        {
            static char b2[320];
            snprintf(b2, sizeof(b2), "trial %d: %s", trials, err);
            mu_assert(b2, 0);
        }
        trials++;
    }
    return 0;
}

/* The folded problem of r as an LP, for the KKT checks. */
static FoldTestLP folded_lp(const FoldApplyResult *r)
{
    const Matrix *B = r->folded;
    FoldTestLP p = flp_alloc((int) B->m, (int) B->n, (int) B->nnz);
    int pos = 0;
    for (int i = 0; i < p.m; i++)
    {
        p.Ap[i] = pos;
        for (int k = B->p[i].start; k < B->p[i].end; k++)
        {
            p.Ai[pos] = B->i[k];
            p.Ax[pos] = B->x[k];
            pos++;
        }
        p.lhs[i] = r->lhs[i];
        p.rhs[i] = r->rhs[i];
    }
    p.Ap[p.m] = pos;
    for (int j = 0; j < p.n; j++)
    {
        p.lb[j] = r->bounds[j].lb;
        p.ub[j] = r->bounds[j].ub;
        p.c[j] = r->c[j];
    }
    return p;
}

static PresolvedProblem lp_problem(FoldTestLP *p)
{
    return problem_from_csr(p->Ax, p->Ai, p->Ap, p->m, p->n, p->nnz, p->lhs, p->rhs,
                            p->lb, p->ub, p->c);
}

/* Three copies of  min -x0 - 2 x1  s.t.  x0 + x1 <= 1,  0 <= x <= 1  fold to
   x_f0 + x_f1 <= 1 with costs (-3, -6), whose KKT point x_f = (0, 1),
   y_f = -3, z_f = (0, -3) lifts to x = (0, 1), y = -1, z = (0, -1) per copy. */
static const char *test_fold_postsolve_lift(void)
{
    static const int Ap[2] = {0, 2}, Ai[2] = {0, 1};
    static const double Ax[2] = {1, 1};
    FoldTestLP b = flp_new(1, 2, 2, Ap, Ai, Ax);
    flp_fill_rows(&b, -INF, 1);
    flp_fill_cols(&b, 0, 1, -1);
    b.c[1] = -2;
    FoldTestLP p = flp_replicate(&b, 3);
    FoldApplyResult r = run_apply(&p);
    mu_assert("expected foldable", r.reason == FOLD_FOLDABLE);
    FoldTestLP f = folded_lp(&r);
    mu_assert("folded 1 x 2", f.m == 1 && f.n == 2 && f.nnz == 2);
    mu_assert("folded costs", f.c[0] == -3.0 && f.c[1] == -6.0);

    double x_f[2] = {0, 1}, y_f[1] = {-3}, z_f[2] = {0, -3};
    PresolvedProblem folded = lp_problem(&f);
    mu_assert("folded point is KKT", is_kkt_point(&folded, x_f, y_f, z_f, 1e-9));
    double x[6], y[3], z[6];
    fold_postsolve(r.maps, p.m, p.n, x_f, y_f, z_f, x, y, z);
    for (int j = 0; j < 6; j++)
    {
        mu_assert("x copied", x[j] == (double) (j % 2));
        mu_assert("z divided by the column class size",
                  z[j] == (j % 2 ? -1.0 : 0.0));
    }
    for (int i = 0; i < 3; i++)
    {
        mu_assert("y divided by the row class size", y[i] == -1.0);
    }
    PresolvedProblem original = lp_problem(&p);
    mu_assert("lifted point is KKT", is_kkt_point(&original, x, y, z, 1e-9));
    flp_free(&f);
    apply_result_free(&r);
    flp_free(&p);
    flp_free(&b);
    return 0;
}

/* Two equal rows x0 + x1 + x2 = 1.5 with 0 <= x <= 0.5 and costs -1 fold to
   3 x_f = 1.5 with cost -3. The KKT point x_f = 0.5, y_f = -0.5, z_f = -1.5
   lifts with y divided by 2 and z by 3. */
static const char *test_fold_postsolve_class_sizes(void)
{
    static const int Ap[3] = {0, 3, 6}, Ai[6] = {0, 1, 2, 0, 1, 2};
    static const double Ax[6] = {1, 1, 1, 1, 1, 1};
    FoldTestLP p = flp_new(2, 3, 6, Ap, Ai, Ax);
    flp_fill_rows(&p, 1.5, 1.5);
    flp_fill_cols(&p, 0, 0.5, -1);
    FoldApplyResult r = run_apply(&p);
    mu_assert("expected foldable", r.reason == FOLD_FOLDABLE);
    FoldTestLP f = folded_lp(&r);
    mu_assert("folded 1 x 1", f.m == 1 && f.n == 1 && f.nnz == 1);
    mu_assert("folded data", f.Ax[0] == 3.0 && f.c[0] == -3.0);

    double x_f[1] = {0.5}, y_f[1] = {-0.5}, z_f[1] = {-1.5};
    PresolvedProblem folded = lp_problem(&f);
    mu_assert("folded point is KKT", is_kkt_point(&folded, x_f, y_f, z_f, 1e-9));
    double x[3], y[2], z[3];
    fold_postsolve(r.maps, p.m, p.n, x_f, y_f, z_f, x, y, z);
    for (int j = 0; j < 3; j++)
    {
        mu_assert("x copied", x[j] == 0.5);
        mu_assert("z divided by the column class size", z[j] == -0.5);
    }
    for (int i = 0; i < 2; i++)
    {
        mu_assert("y divided by the row class size", y[i] == -0.25);
    }
    PresolvedProblem original = lp_problem(&p);
    mu_assert("lifted point is KKT", is_kkt_point(&original, x, y, z, 1e-9));
    flp_free(&f);
    apply_result_free(&r);
    flp_free(&p);
    return 0;
}

/* ---- through the presolver ----------------------------------------------- */

static Presolver *fold_presolver(const FoldTestLP *p, const Settings *stgs)
{
    return new_presolver(p->Ax, p->Ai, p->Ap, (size_t) p->m, (size_t) p->n,
                         (size_t) p->nnz, p->lhs, p->rhs, p->lb, p->ub, p->c, stgs);
}

/* Everything off but the fold. */
static Settings *fold_only_settings(void)
{
    Settings *stgs = default_settings();
    set_settings_false(stgs);
    stgs->lp_fold = true;
    return stgs;
}

/* The LP of test_fold_postsolve_lift: three copies of
   min -x0 - 2 x1  s.t.  x0 + x1 <= 1,  0 <= x <= 1. */
static FoldTestLP simplex_copies(void)
{
    static const int Ap[2] = {0, 2}, Ai[2] = {0, 1};
    static const double Ax[2] = {1, 1};
    FoldTestLP b = flp_new(1, 2, 2, Ap, Ai, Ax);
    flp_fill_rows(&b, -INF, 1);
    flp_fill_cols(&b, 0, 1, -1);
    b.c[1] = -2;
    FoldTestLP p = flp_replicate(&b, 3);
    flp_free(&b);
    return p;
}

static bool all_equal(const double *values, int len, double value)
{
    for (int idx = 0; idx < len; idx++)
    {
        if (values[idx] != value) return false;
    }
    return true;
}

/* The fold happens in new_presolver, which rebuilds the transpose, and
   postsolve lifts the point of the folded problem to the original. */
static const char *test_fold_presolver_lift(void)
{
    FoldTestLP p = simplex_copies();
    Settings *stgs = fold_only_settings();
    Presolver *presolver = fold_presolver(&p, stgs);
    mu_assert("presolver", presolver != NULL);
    const Constraints *constraints = presolver->prob->constraints;
    mu_assert("folded 1 x 2", constraints->A->m == 1 && constraints->A->n == 2 &&
                                  constraints->A->nnz == 2);
    mu_assert("transpose rebuilt", constraints->AT->m == 2 &&
                                       constraints->AT->n == 1 &&
                                       constraints->AT->nnz == 2);
    mu_assert("maps kept", constraints->state->work->fold_maps != NULL);
    mu_assert("presolve", run_presolver(presolver) != INFEASIBLE);
    const PresolvedProblem *reduced = presolver->reduced_prob;
    mu_assert("reduced 1 x 2",
              reduced->m == 1 && reduced->n == 2 && reduced->nnz == 2);
    mu_assert("folded data", reduced->c[0] == -3.0 && reduced->c[1] == -6.0 &&
                                 reduced->rhs[0] == 1.0);
    const PresolveStats *stats = presolver->stats;
    mu_assert("original dimensions", stats->n_rows_original == 3 &&
                                         stats->n_cols_original == 6 &&
                                         stats->nnz_original == 6);
    mu_assert("folded dimensions", stats->n_rows_folded == 1 &&
                                       stats->n_cols_folded == 2 &&
                                       stats->nnz_removed_lp_fold == 4);

    double x_f[2] = {0, 1}, y_f[1] = {-3}, z_f[2] = {0, -3};
    mu_assert("folded point is KKT", is_kkt_point(reduced, x_f, y_f, z_f, 1e-9));
    postsolve(presolver, x_f, y_f, z_f);
    const Solution *sol = presolver->sol;
    for (int j = 0; j < 6; j++)
    {
        mu_assert("x copied", sol->x[j] == (double) (j % 2));
        mu_assert("z divided by the column class size",
                  sol->z[j] == (j % 2 ? -1.0 : 0.0));
    }
    for (int i = 0; i < 3; i++)
    {
        mu_assert("y divided by the row class size", sol->y[i] == -1.0);
    }
    PresolvedProblem original = lp_problem(&p);
    mu_assert("lifted point is KKT",
              is_kkt_point(&original, sol->x, sol->y, sol->z, 1e-9));
    free_presolver(presolver);
    free_settings(stgs);
    flp_free(&p);
    return 0;
}

/* The forward map and the ray postsolves are not supported with a fold: they
   print an error and leave their outputs untouched. */
static const char *test_fold_refuses_map_to_reduced(void)
{
    FoldTestLP p = simplex_copies();
    Settings *stgs = fold_only_settings();
    Presolver *presolver = fold_presolver(&p, stgs);
    mu_assert("presolver", presolver && run_presolver(presolver) != INFEASIBLE);
    double x[6] = {0, 1, 0, 1, 0, 1}, y[3] = {-1, -1, -1};
    double x_red[2] = {DUMMY_VALUE, DUMMY_VALUE}, y_red[1] = {DUMMY_VALUE};
    double z_red[2] = {DUMMY_VALUE, DUMMY_VALUE};
    map_original_sol_to_reduced(presolver, x, y, x_red, y_red, z_red);
    mu_assert("x_red untouched", all_equal(x_red, 2, DUMMY_VALUE));
    mu_assert("y_red untouched", all_equal(y_red, 1, DUMMY_VALUE));
    mu_assert("z_red untouched", all_equal(z_red, 2, DUMMY_VALUE));
    free_presolver(presolver);
    free_settings(stgs);
    flp_free(&p);
    return 0;
}

static const char *test_fold_refuses_primal_ray(void)
{
    FoldTestLP p = simplex_copies();
    Settings *stgs = fold_only_settings();
    Presolver *presolver = fold_presolver(&p, stgs);
    mu_assert("presolver", presolver && run_presolver(presolver) != INFEASIBLE);
    double y_f[1] = {1};
    double y_orig[3] = {DUMMY_VALUE, DUMMY_VALUE, DUMMY_VALUE};
    postsolve_primal_infeas_ray(presolver, y_f, y_orig);
    mu_assert("y_orig untouched", all_equal(y_orig, 3, DUMMY_VALUE));
    free_presolver(presolver);
    free_settings(stgs);
    flp_free(&p);
    return 0;
}

static const char *test_fold_refuses_dual_ray(void)
{
    FoldTestLP p = simplex_copies();
    Settings *stgs = fold_only_settings();
    Presolver *presolver = fold_presolver(&p, stgs);
    mu_assert("presolver", presolver && run_presolver(presolver) != INFEASIBLE);
    double x_f[2] = {1, 1};
    double x_orig[6] = {DUMMY_VALUE, DUMMY_VALUE, DUMMY_VALUE,
                        DUMMY_VALUE, DUMMY_VALUE, DUMMY_VALUE};
    postsolve_dual_infeas_ray(presolver, x_f, x_orig);
    mu_assert("x_orig untouched", all_equal(x_orig, 6, DUMMY_VALUE));
    free_presolver(presolver);
    free_settings(stgs);
    flp_free(&p);
    return 0;
}

/* Three copies of  min -x0 - x1  s.t.  x0 + x1 >= 1,  0 <= x <= 1  with every
   reduction on: rows and columns fold to  2 x_f >= 1  with cost -6, the dual
   fix solves that, and the lift gives x = 1, y = 0, z = -1. */
static const char *test_fold_presolver_full_settings(void)
{
    static const int Ap[2] = {0, 2}, Ai[2] = {0, 1};
    static const double Ax[2] = {1, 1};
    FoldTestLP b = flp_new(1, 2, 2, Ap, Ai, Ax);
    flp_fill_rows(&b, 1, INF);
    flp_fill_cols(&b, 0, 1, -1);
    FoldTestLP p = flp_replicate(&b, 3);
    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->verbose = false;
    stgs->lp_fold = true;
    Presolver *presolver = fold_presolver(&p, stgs);
    mu_assert("presolver", presolver && run_presolver(presolver) != INFEASIBLE);
    const PresolveStats *stats = presolver->stats;
    mu_assert("folded 1 x 1",
              stats->n_rows_folded == 1 && stats->n_cols_folded == 1);
    mu_assert("solved by presolve",
              presolver->reduced_prob->m == 0 && presolver->reduced_prob->n == 0);

    double x_f[1], y_f[1], z_f[1]; /* the reduced problem is empty */
    postsolve(presolver, x_f, y_f, z_f);
    const Solution *sol = presolver->sol;
    mu_assert("x lifted", all_equal(sol->x, 6, 1.0));
    mu_assert("y lifted", all_equal(sol->y, 3, 0.0));
    mu_assert("z lifted", all_equal(sol->z, 6, -1.0));
    PresolvedProblem original = lp_problem(&p);
    mu_assert("lifted point is KKT",
              is_kkt_point(&original, sol->x, sol->y, sol->z, 1e-9));
    free_presolver(presolver);
    free_settings(stgs);
    flp_free(&p);
    flp_free(&b);
    return 0;
}

static const char *all_tests_fold(void)
{
    mu_run_test(test_fold_sides, counter_fold);
    mu_run_test(test_fold_symmetric_blocks, counter_fold);
    mu_run_test(test_fold_asymmetric_blocks, counter_fold);
    mu_run_test(test_fold_degenerate, counter_fold);
    mu_run_test(test_fold_not_converged, counter_fold);
    mu_run_test(test_fold_one_sided, counter_fold);
    mu_run_test(test_fold_gate, counter_fold);
    mu_run_test(test_fold_random_vs_reference, counter_fold);
    mu_run_test(test_fold_large, counter_fold);
    mu_run_test(test_fold_verify_entries, counter_fold);
    mu_run_test(test_fold_verify_other_end, counter_fold);
    mu_run_test(test_fold_check, counter_fold);
    mu_run_test(test_fold_apply_blocks, counter_fold);
    mu_run_test(test_fold_apply_column_sum, counter_fold);
    mu_run_test(test_fold_apply_exact_zero, counter_fold);
    mu_run_test(test_fold_apply_unsorted_row, counter_fold);
    mu_run_test(test_fold_apply_random, counter_fold);
    mu_run_test(test_fold_postsolve_lift, counter_fold);
    mu_run_test(test_fold_postsolve_class_sizes, counter_fold);
    mu_run_test(test_fold_presolver_lift, counter_fold);
    mu_run_test(test_fold_refuses_map_to_reduced, counter_fold);
    mu_run_test(test_fold_refuses_primal_ray, counter_fold);
    mu_run_test(test_fold_refuses_dual_ray, counter_fold);
    mu_run_test(test_fold_presolver_full_settings, counter_fold);
    return 0;
}

int test_fold(void)
{
    const char *result = all_tests_fold();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("fold: TEST FAILED!\n");
    }
    else
    {
        printf("fold: ALL TESTS PASSED\n");
    }
    printf("fold: Tests run: %d\n", counter_fold);
    return result == 0;
}

#endif // TEST_FOLD_H
