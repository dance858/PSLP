#ifndef TEST_ROW_DIFFERENCING_H
#define TEST_ROW_DIFFERENCING_H

#include "Matrix.h"
#include "Memory_wrapper.h"
#include "PSLP_API.h"
#include "PSLP_sol.h"
#include "PSLP_stats.h"
#include "RowDifferencing.h"
#include "kkt.h"
#include "minunit.h"
#include <math.h>
#include <string.h>

static int counter_row_differencing = 0;

/* A chain of running sums. Row k (k = 0, ..., K - 1) is
       sign_k (x_0 + ... + x_{(k+1)w-1} - s_k) = side_k
   over K w columns x followed by K columns s. With K = 6 and w = 10 the rows
   have 11, 21, ..., 61 entries, and a row minus the one before it has w + 2:
   rows 2 to 5 qualify (12 <= len / 2), row 1 does not (12 > 10). */
#define RD_K 6
#define RD_W 10
#define RD_N (RD_K * RD_W + RD_K)
#define RD_NNZ (RD_W * RD_K * (RD_K + 1) / 2 + RD_K)

typedef struct
{
    double Ax[RD_NNZ + 2000];
    int Ai[RD_NNZ + 2000];
    int Ap[RD_K + 40 + 1];
    double lhs[RD_K + 40], rhs[RD_K + 40];
    int m, nnz;
} RdChain;

static void rd_chain(RdChain *ch, const double *signs, const double *sides)
{
    int nnz = 0;
    for (int k = 0; k < RD_K; ++k)
    {
        ch->Ap[k] = nnz;
        for (int j = 0; j < (k + 1) * RD_W; ++j)
        {
            ch->Ai[nnz] = j;
            ch->Ax[nnz++] = signs[k];
        }
        ch->Ai[nnz] = RD_K * RD_W + k;
        ch->Ax[nnz++] = -signs[k];
        ch->lhs[k] = ch->rhs[k] = sides[k];
    }
    ch->Ap[RD_K] = nnz;
    ch->m = RD_K;
    ch->nnz = nnz;
}

static void rd_dense_row(const double *x, const int *idx, int start, int end,
                         double scale, double *dense)
{
    for (int q = start; q < end; ++q)
    {
        dense[idx[q]] += scale * x[q];
    }
}

/* True iff every row of A is the input row with its sides, or the input row
   minus alpha times its parent with the side to match, the rows are strictly
   increasing without zeros, and A'S'y = A'y for a test vector y. */
static bool rd_matches(const RdChain *ch, const double *sides_in, const Matrix *A,
                       const double *lhs, const double *rhs,
                       const RowDifferencing *rd)
{
    double dense[RD_N], y[RD_K + 40], ATy_new[RD_N] = {0}, ATy_old[RD_N] = {0};
    size_t k = 0;
    for (int r = 0; r < ch->m; ++r)
    {
        memset(dense, 0, sizeof(dense));
        rd_dense_row(ch->Ax, ch->Ai, ch->Ap[r], ch->Ap[r + 1], 1.0, dense);
        if (k < rd->rows->len && rd->rows->data[k] == r)
        {
            int p = rd->parents->data[k];
            double alpha = (double) rd->alphas->data[k];
            if (p >= r || (alpha != 1.0 && alpha != -1.0)) return false;
            rd_dense_row(ch->Ax, ch->Ai, ch->Ap[p], ch->Ap[p + 1], -alpha, dense);
            double side = sides_in[r] - alpha * sides_in[p];
            if (lhs[r] != side || rhs[r] != side) return false;
            ++k;
        }
        else if (lhs[r] != ch->lhs[r] || rhs[r] != ch->rhs[r])
        {
            return false;
        }
        for (int q = A->p[r].start; q < A->p[r].end; ++q)
        {
            if (A->x[q] == 0.0 || dense[A->i[q]] != A->x[q]) return false;
            if (q > A->p[r].start && A->i[q - 1] >= A->i[q]) return false;
            dense[A->i[q]] = 0.0;
        }
        for (int j = 0; j < RD_N; ++j)
        {
            if (dense[j] != 0.0) return false;
        }
    }
    if (k != rd->rows->len) return false;

    /* the duals: y of the differenced rows, mapped back */
    for (int r = 0; r < ch->m; ++r)
    {
        y[r] = 0.5 + r;
        rd_dense_row(A->x, A->i, A->p[r].start, A->p[r].end, y[r], ATy_new);
    }
    row_differencing_postsolve_duals(rd, y);
    for (int r = 0; r < ch->m; ++r)
    {
        rd_dense_row(ch->Ax, ch->Ai, ch->Ap[r], ch->Ap[r + 1], y[r], ATy_old);
    }
    for (int j = 0; j < RD_N; ++j)
    {
        if (fabs(ATy_new[j] - ATy_old[j]) > 1e-12) return false;
    }
    return true;
}

/* The chain collapses: rows 2 to 5 become the differences. */
static char *test_row_differencing_chain(void)
{
    static RdChain ch;
    double signs[RD_K] = {1, 1, 1, 1, 1, 1};
    double sides[RD_K] = {3, 3, 3, 3, 3, 3};
    rd_chain(&ch, signs, sides);
    {
        double lhs[RD_K], rhs[RD_K];
        memcpy(lhs, ch.lhs, sizeof(lhs));
        memcpy(rhs, ch.rhs, sizeof(rhs));
        RowDifferencing *rd = NULL;
        Matrix *A = row_differencing_matrix_new(
            ch.Ax, ch.Ai, ch.Ap, RD_K, RD_K, RD_N, (size_t) ch.nnz, lhs, rhs, &rd);
        mu_assert("the chain must be rewritten", A != NULL && rd != NULL);
        mu_assert("rows 2 to 5", rd->rows->len == 4 && rd->rows->data[0] == 2 &&
                                     rd->rows->data[3] == 5);
        mu_assert("parents", rd->parents->data[0] == 1 && rd->parents->data[3] == 4);
        mu_assert("entries", A->nnz == 11 + 21 + 4 * 12);
        mu_assert("removed", rd->nnz_removed == (size_t) ch.nnz - A->nnz);
        mu_assert("equal sides give a zero side", rhs[2] == 0.0 && rhs[5] == 0.0);
        mu_assert("matrix", rd_matches(&ch, sides, A, lhs, rhs, rd));
        free_matrix(A);
        row_differencing_free(rd);
    }

    /* explicit zeros in the input: the slack entries of rows 3 and 5. They
       are dropped from the differences and from the copied rows. */
    {
        double lhs[RD_K], rhs[RD_K];
        ch.Ax[ch.Ap[4] - 1] = 0.0;
        ch.Ax[ch.Ap[6] - 1] = 0.0;
        ch.Ax[ch.Ap[1] - 1] = 0.0;
        memcpy(lhs, ch.lhs, sizeof(lhs));
        memcpy(rhs, ch.rhs, sizeof(rhs));
        RowDifferencing *rd = NULL;
        Matrix *A = row_differencing_matrix_new(
            ch.Ax, ch.Ai, ch.Ap, RD_K, RD_K, RD_N, (size_t) ch.nnz, lhs, rhs, &rd);
        mu_assert("rewritten with zeros",
                  A != NULL && rd != NULL && rd->rows->len == 4);
        mu_assert("entries without zeros", A->nnz == 10 + 21 + 12 + 11 + 11 + 11);
        mu_assert("matrix without zeros", rd_matches(&ch, sides, A, lhs, rhs, rd));
        free_matrix(A);
        row_differencing_free(rd);
    }
    return 0;
}

/* alpha = -1 and the three exact cases of the side: equal up to alpha, a
   zero parent side, a zero own side. */
static char *test_row_differencing_signs_and_sides(void)
{
    static RdChain ch;
    double signs[RD_K] = {1, -1, 1, -1, -1, 1};
    double sides[RD_K] = {0, 0, 0, 5, 0, 7};
    rd_chain(&ch, signs, sides);
    /* row 2 (0 vs 0), row 3 (5 vs 0), row 4 (0 vs 5, alpha 1), row 5 (7 vs 0) */
    double lhs[RD_K], rhs[RD_K];
    memcpy(lhs, ch.lhs, sizeof(lhs));
    memcpy(rhs, ch.rhs, sizeof(rhs));
    RowDifferencing *rd = NULL;
    Matrix *A = row_differencing_matrix_new(ch.Ax, ch.Ai, ch.Ap, RD_K, RD_K, RD_N,
                                            (size_t) ch.nnz, lhs, rhs, &rd);
    mu_assert("rewritten", A != NULL && rd != NULL && rd->rows->len == 4);
    mu_assert("alphas", rd->alphas->data[0] == -1 && rd->alphas->data[1] == -1 &&
                            rd->alphas->data[2] == 1 && rd->alphas->data[3] == -1);
    mu_assert("sides",
              rhs[2] == 0.0 && rhs[3] == 5.0 && rhs[4] == -5.0 && rhs[5] == 7.0);
    mu_assert("matrix", rd_matches(&ch, sides, A, lhs, rhs, rd));
    free_matrix(A);
    row_differencing_free(rd);

    /* two nonzero sides that differ: row 5 is not a candidate and keeps its
       side */
    double sides2[RD_K] = {1, 1, 1, 1, 1, 3};
    double signs2[RD_K] = {1, 1, 1, 1, 1, 1};
    rd_chain(&ch, signs2, sides2);
    memcpy(lhs, ch.lhs, sizeof(lhs));
    memcpy(rhs, ch.rhs, sizeof(rhs));
    A = row_differencing_matrix_new(ch.Ax, ch.Ai, ch.Ap, RD_K, RD_K, RD_N,
                                    (size_t) ch.nnz, lhs, rhs, &rd);
    mu_assert("rows 2 to 4 only",
              A != NULL && rd->rows->len == 3 && rd->rows->data[2] == 4);
    mu_assert("side kept", rhs[5] == 3.0);
    mu_assert("matrix", rd_matches(&ch, sides2, A, lhs, rhs, rd));
    free_matrix(A);
    row_differencing_free(rd);
    return 0;
}

/* Nothing is rewritten, and the sides are left alone. */
static bool rd_rejected(RdChain *ch, int n_lin_rows)
{
    double lhs[RD_K + 40], rhs[RD_K + 40];
    memcpy(lhs, ch->lhs, sizeof(lhs));
    memcpy(rhs, ch->rhs, sizeof(rhs));
    RowDifferencing *rd = NULL;
    Matrix *A = row_differencing_matrix_new(ch->Ax, ch->Ai, ch->Ap, (size_t) ch->m,
                                            (size_t) n_lin_rows, RD_N,
                                            (size_t) ch->nnz, lhs, rhs, &rd);
    bool ok = A == NULL && rd == NULL &&
              memcmp(lhs, ch->lhs, (size_t) ch->m * sizeof(double)) == 0 &&
              memcmp(rhs, ch->rhs, (size_t) ch->m * sizeof(double)) == 0;
    free_matrix(A);
    row_differencing_free(rd);
    return ok;
}

static char *test_row_differencing_rejections(void)
{
    static RdChain ch;
    double signs[RD_K] = {1, 1, 1, 1, 1, 1};
    double sides[RD_K] = {3, 3, 3, 3, 3, 3};

    /* inequalities */
    rd_chain(&ch, signs, sides);
    for (int k = 0; k < RD_K; ++k) ch.lhs[k] = -INFINITY;
    mu_assert("inequalities", rd_rejected(&ch, RD_K));

    /* a shared entry that does not cancel in every row */
    rd_chain(&ch, signs, sides);
    for (int k = 0; k < RD_K; ++k) ch.Ax[ch.Ap[k] + 5] = 1.0 + 0.125 * k;
    mu_assert("no exact cancellation", rd_rejected(&ch, RD_K));

    /* leading coefficients that are not equal up to sign */
    rd_chain(&ch, signs, sides);
    for (int k = 0; k < RD_K; ++k) ch.Ax[ch.Ap[k]] = 1.0 + k;
    mu_assert("leading coefficients", rd_rejected(&ch, RD_K));

    /* cone rows are not linear rows */
    rd_chain(&ch, signs, sides);
    mu_assert("no linear rows", rd_rejected(&ch, 0));

    /* Other rows after the chain. With 4 rows of 50 entries 136 of 416
       entries go (33%) and the chain is rewritten; with 40 rows it is 136 of
       2216 (6%), below ROW_DIFFERENCING_MIN_GAIN, and nothing is. */
    rd_chain(&ch, signs, sides);
    for (int r = RD_K; r < RD_K + 40; ++r)
    {
        for (int j = 0; j < 50; ++j)
        {
            ch.Ai[ch.nnz] = j;
            ch.Ax[ch.nnz++] = 2.0 + r + j;
        }
        ch.Ap[r + 1] = ch.nnz;
        ch.lhs[r] = -INFINITY;
        ch.rhs[r] = 1.0;
        if (r == RD_K + 3)
        {
            double lhs[RD_K + 40], rhs[RD_K + 40];
            ch.m = r + 1;
            memcpy(lhs, ch.lhs, sizeof(lhs));
            memcpy(rhs, ch.rhs, sizeof(rhs));
            RowDifferencing *rd = NULL;
            Matrix *A = row_differencing_matrix_new(
                ch.Ax, ch.Ai, ch.Ap, (size_t) ch.m, (size_t) ch.m, RD_N,
                (size_t) ch.nnz, lhs, rhs, &rd);
            mu_assert("rows after the last pair are copied",
                      A != NULL && rd != NULL && rd->rows->len == 4 &&
                          rd_matches(&ch, ch.rhs, A, lhs, rhs, rd));
            free_matrix(A);
            row_differencing_free(rd);
        }
    }
    ch.m = RD_K + 40;
    mu_assert("gain below the minimum", rd_rejected(&ch, ch.m));
    return 0;
}

/* End to end. The x are fixed by their bounds and the slacks s are free with
   a cost, so presolve solves the problem and the duals of the original rows
   are determined: the slack of row k is in that row only, so y_k = -c_k.
   In the differenced rows a slack is in two rows, so the KKT test fails
   unless the duals are mapped back. */
static char *test_row_differencing_presolve_postsolve(void)
{
    static RdChain ch;
    double signs[RD_K] = {1, 1, 1, 1, 1, 1};
    double sides[RD_K] = {3, 3, 3, 3, 3, 3};
    rd_chain(&ch, signs, sides);
    double lbs[RD_N], ubs[RD_N], c[RD_N];
    for (int j = 0; j < RD_N; ++j)
    {
        bool is_x = j < RD_K * RD_W;
        lbs[j] = is_x ? 0.25 * (j % 7) : -INFINITY;
        ubs[j] = is_x ? lbs[j] : INFINITY;
        c[j] = is_x ? 1.0 : 1.0 + (j - RD_K * RD_W);
    }
    ProblemData orig = problem_from_csr(ch.Ax, ch.Ai, ch.Ap, RD_K, RD_N, ch.nnz,
                                        ch.lhs, ch.rhs, lbs, ubs, c);
    for (int on = 1; on >= 0; --on)
    {
        Settings *stgs = default_settings();
        stgs->verbose = false;
        stgs->row_differencing = on;
        Presolver *presolver =
            new_presolver(ch.Ax, ch.Ai, ch.Ap, RD_K, RD_N, (size_t) ch.nnz, ch.lhs,
                          ch.rhs, lbs, ubs, c, stgs);
        mu_assert("presolver", presolver != NULL);
        mu_assert("rows differenced",
                  presolver->stats->n_rows_differenced == (on ? 4 : 0));
        mu_assert("entries removed",
                  presolver->stats->nnz_removed_row_differencing ==
                      (size_t) (on ? ch.nnz - (11 + 21 + 4 * 12) : 0));
        run_presolver(presolver);
        mu_assert("presolve solves the problem",
                  presolver->stats->n_rows_reduced == 0 &&
                      presolver->stats->n_cols_reduced == 0);
        double none = 0.0;
        postsolve(presolver, &none, &none, &none);
        for (int k = 0; k < RD_K; ++k)
        {
            mu_assert("dual of an original row",
                      fabs(presolver->sol->y[k] + c[RD_K * RD_W + k]) < 1e-9);
        }
        mu_assert("postsolved point is not optimal",
                  is_kkt_point(&orig, presolver->sol->x, presolver->sol->y,
                               presolver->sol->z, 1e-9));
        free_presolver(presolver);
        PS_FREE(stgs);
    }
    return 0;
}

/* mapping duals to the differenced rows inverts the dual postsolve, also
   along a chain */
static char *test_row_differencing_map_duals(void)
{
    int rows[4] = {1, 2, 4, 5}, parents[4] = {0, 1, 3, 4};
    int alphas[4] = {1, -1, 1, 1};
    iVec rows_vec = {.data = rows, .len = 4, .capacity = 4};
    iVec parents_vec = {.data = parents, .len = 4, .capacity = 4};
    iVec alphas_vec = {.data = alphas, .len = 4, .capacity = 4};
    RowDifferencing rd = {
        .rows = &rows_vec, .parents = &parents_vec, .alphas = &alphas_vec};
    double y0[6] = {0.5, -2.0, 3.25, 1.0, 0.0, -7.5}, y[6];
    memcpy(y, y0, sizeof(y));
    row_differencing_map_duals(&rd, y);
    /* y = S^T y': y'_5 = y_5, y'_4 = y_4 + y'_5, y'_3 = y_3 + y'_4, ... */
    mu_assert("mapped duals", y[5] == -7.5 && y[4] == -7.5 && y[3] == -6.5 &&
                                  y[2] == 3.25 && y[1] == -5.25 && y[0] == -4.75);
    row_differencing_postsolve_duals(&rd, y);
    for (int k = 0; k < 6; ++k)
    {
        mu_assert("round trip", y[k] == y0[k]);
    }
    return 0;
}

static const char *all_tests_row_differencing(void)
{
    mu_run_test(test_row_differencing_chain, counter_row_differencing);
    mu_run_test(test_row_differencing_signs_and_sides, counter_row_differencing);
    mu_run_test(test_row_differencing_rejections, counter_row_differencing);
    mu_run_test(test_row_differencing_presolve_postsolve, counter_row_differencing);
    mu_run_test(test_row_differencing_map_duals, counter_row_differencing);
    return 0;
}

int test_row_differencing(void)
{
    const char *result = all_tests_row_differencing();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("row differencing: TEST FAILED!\n");
    }
    else
    {
        printf("row differencing: ALL TESTS PASSED\n");
    }
    printf("row differencing: Tests run: %d\n", counter_row_differencing);
    return result == 0;
}

#endif // TEST_ROW_DIFFERENCING_H
