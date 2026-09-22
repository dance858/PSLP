#include "kkt.h"
#include "Numerics.h"
#include "glbopts.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

PresolvedProblem problem_from_csr(double *Ax, int *Ai, int *Ap, int n_rows,
                                  int n_cols, int nnz, double *lhs, double *rhs,
                                  double *lbs, double *ubs, double *c)
{
    PresolvedProblem prob;
    prob.Ax = Ax;
    prob.Ai = Ai;
    prob.Ap = Ap;
    prob.m = (size_t) n_rows;
    prob.n = (size_t) n_cols;
    prob.nnz = (size_t) nnz;
    prob.lhs = lhs;
    prob.rhs = rhs;
    prob.c = c;
    prob.lbs = lbs;
    prob.ubs = ubs;
    prob.obj_offset = 0.0;
    return prob;
}

/* p(s; l, u) = sum_i (s_i < 0 ? s_i l_i : s_i > 0 ? s_i u_i : 0).
   A term with an infinite bound is skipped; for a dual feasible point the
   corresponding multiplier is zero anyway. */
static double support_term(double s, double l, double u)
{
    if (s < 0.0) return isinf(l) ? 0.0 : s * l;
    if (s > 0.0) return isinf(u) ? 0.0 : s * u;
    return 0.0;
}

/* Complementary-slackness projection: the part of the multiplier v that is
   admissible when the activity a sits in [l, u]. */
static double admissible_multiplier(double v, double a, double l, double u)
{
    if (!isinf(l) && !isinf(u) && fabs(l - u) <= FEAS_TOL) return v;
    if (!isinf(l) && fabs(a - l) <= FEAS_TOL) return v > 0.0 ? v : 0.0;
    if (!isinf(u) && fabs(a - u) <= FEAS_TOL) return v < 0.0 ? v : 0.0;
    return 0.0;
}

static double clip(double v, double l, double u)
{
    if (!isinf(l) && v < l) v = l;
    if (!isinf(u) && v > u) v = u;
    return v;
}

void kkt_residuals(const PresolvedProblem *prob, const double *x, const double *y,
                   const double *z, KKTResiduals *out)
{
    const size_t m = prob->m, n = prob->n;
    const int *Ap = prob->Ap, *Ai = prob->Ai;
    const double *Ax = prob->Ax, *c = prob->c;
    const double *lhs = prob->lhs, *rhs = prob->rhs;
    const double *lbs = prob->lbs, *ubs = prob->ubs;

    double *Ax_vec = (double *) calloc(m > 0 ? m : 1, sizeof(double));
    double *ATy = (double *) calloc(n > 0 ? n : 1, sizeof(double));
    if (!Ax_vec || !ATy)
    {
        free(Ax_vec);
        free(ATy);
        out->dual_res_abs = out->primal_res_abs = out->gap_abs = NAN;
        out->dual_res_rel = out->primal_res_rel = out->gap_rel = NAN;
        out->y_comp_slack = out->z_comp_slack = NAN;
        out->y_norm = out->z_norm = out->primal_obj = NAN;
        out->bound_viol = out->x_norm = NAN;
        return;
    }

    for (size_t i = 0; i < m; i++)
    {
        double acc = 0.0;
        for (int k = Ap[i]; k < Ap[i + 1]; k++)
        {
            acc += Ax[k] * x[Ai[k]];
            ATy[Ai[k]] += Ax[k] * y[i];
        }
        Ax_vec[i] = acc;
    }

    // objective, norms of the data, bound violation of x
    double cx = 0.0, c_norm2 = 0.0, sides_norm2 = 0.0;
    double x_norm2 = 0.0, bviol2 = 0.0;
    for (size_t j = 0; j < n; j++)
    {
        cx += c[j] * x[j];
        c_norm2 += c[j] * c[j];
        x_norm2 += x[j] * x[j];
        double proj = clip(x[j], lbs[j], ubs[j]);
        bviol2 += (x[j] - proj) * (x[j] - proj);
    }
    for (size_t i = 0; i < m; i++)
    {
        if (!isinf(lhs[i])) sides_norm2 += lhs[i] * lhs[i];
        if (!isinf(rhs[i])) sides_norm2 += rhs[i] * rhs[i];
    }

    // dual residual, z norm, bound support term, z complementary slackness
    double dual2 = 0.0, z_norm2 = 0.0, p_z = 0.0, zcs2 = 0.0;
    for (size_t j = 0; j < n; j++)
    {
        double r = ATy[j] + z[j] - c[j];
        dual2 += r * r;
        z_norm2 += z[j] * z[j];
        p_z += support_term(-z[j], lbs[j], ubs[j]);
        double zj = admissible_multiplier(z[j], x[j], lbs[j], ubs[j]);
        zcs2 += (z[j] - zj) * (z[j] - zj);
    }

    // primal residual, y norm, row support term, y complementary slackness
    double primal2 = 0.0, y_norm2 = 0.0, p_y = 0.0, ycs2 = 0.0;
    for (size_t i = 0; i < m; i++)
    {
        double a = Ax_vec[i];
        double proj = clip(a, lhs[i], rhs[i]);
        primal2 += (a - proj) * (a - proj);
        y_norm2 += y[i] * y[i];
        p_y += support_term(-y[i], lhs[i], rhs[i]);
        double yi = admissible_multiplier(y[i], a, lhs[i], rhs[i]);
        ycs2 += (y[i] - yi) * (y[i] - yi);
    }

    double p_total = p_y + p_z;
    out->primal_obj = cx + prob->obj_offset;
    out->dual_res_abs = sqrt(dual2);
    out->primal_res_abs = sqrt(primal2);
    out->gap_abs = fabs(cx + p_total);
    out->dual_res_rel = out->dual_res_abs / (1.0 + sqrt(c_norm2));
    out->primal_res_rel = out->primal_res_abs / (1.0 + sqrt(sides_norm2));
    out->gap_rel = out->gap_abs / (1.0 + fabs(p_total) + fabs(cx));
    out->y_comp_slack = sqrt(ycs2);
    out->z_comp_slack = sqrt(zcs2);
    out->y_norm = sqrt(y_norm2);
    out->z_norm = sqrt(z_norm2);
    out->bound_viol = sqrt(bviol2);
    out->x_norm = sqrt(x_norm2);

    free(Ax_vec);
    free(ATy);
}

bool is_kkt_point(const PresolvedProblem *prob, const double *x, const double *y,
                  const double *z, double tol)
{
    KKTResiduals r;
    kkt_residuals(prob, x, y, z, &r);

    const char *names[] = {
        "bound violation",           "primal residual",           "dual residual",
        "z complementary slackness", "y complementary slackness", "gap"};
    const double values[] = {r.bound_viol / (1.0 + r.x_norm),
                             r.primal_res_rel,
                             r.dual_res_rel,
                             r.z_comp_slack / (1.0 + r.z_norm),
                             r.y_comp_slack / (1.0 + r.y_norm),
                             r.gap_rel};

    for (size_t k = 0; k < sizeof(values) / sizeof(values[0]); ++k)
    {
        // written so that a NaN fails instead of passing
        if (!(values[k] <= tol))
        {
            printf("KKT: %s = %.3e > tol %.3e\n", names[k], values[k], tol);
            return false;
        }
    }
    return true;
}

bool is_primal_ray_certificate(const PresolvedProblem *prob, const double *y,
                               double tol)
{
    const size_t m = prob->m, n = prob->n;
    double *z = (double *) calloc(n > 0 ? n : 1, sizeof(double));
    if (!z)
    {
        printf("primal ray: allocation failed\n");
        return false;
    }

    // z = -A^T y
    for (size_t i = 0; i < m; ++i)
    {
        for (int p = prob->Ap[i]; p < prob->Ap[i + 1]; ++p)
        {
            z[prob->Ai[p]] -= prob->Ax[p] * y[i];
        }
    }

    double support = 0.0;
    bool ok = true;
    for (size_t i = 0; i < m && ok; ++i)
    {
        if (y[i] > tol)
        {
            if (IS_POS_INF(prob->rhs[i]))
            {
                printf("primal ray: y_%zu = %g > 0 but rhs is infinite\n", i, y[i]);
                ok = false;
            }
            else
                support += prob->rhs[i] * y[i];
        }
        else if (y[i] < -tol)
        {
            if (IS_NEG_INF(prob->lhs[i]))
            {
                printf("primal ray: y_%zu = %g < 0 but lhs is infinite\n", i, y[i]);
                ok = false;
            }
            else
                support += prob->lhs[i] * y[i];
        }
    }

    for (size_t j = 0; j < n && ok; ++j)
    {
        if (z[j] > tol)
        {
            if (IS_POS_INF(prob->ubs[j]))
            {
                printf("primal ray: z_%zu = %g > 0 but ub is infinite\n", j, z[j]);
                ok = false;
            }
            else
                support += prob->ubs[j] * z[j];
        }
        else if (z[j] < -tol)
        {
            if (IS_NEG_INF(prob->lbs[j]))
            {
                printf("primal ray: z_%zu = %g < 0 but lb is infinite\n", j, z[j]);
                ok = false;
            }
            else
                support += prob->lbs[j] * z[j];
        }
    }
    free(z);

    if (ok && support >= -tol)
    {
        printf("primal ray: support %g is not negative\n", support);
        ok = false;
    }
    return ok;
}

/* Recession-direction test of a single row or column: the activity a must be
   zero if both sides are finite, <= 0 if only the upper side is finite and
   >= 0 if only the lower side is finite. */
static bool is_recession_component(double a, double l, double u, double tol)
{
    if (!IS_NEG_INF(l) && !IS_POS_INF(u)) return fabs(a) <= tol;
    if (!IS_POS_INF(u)) return a <= tol;
    if (!IS_NEG_INF(l)) return a >= -tol;
    return true;
}

bool is_dual_ray_certificate(const PresolvedProblem *prob, const double *d,
                             double tol)
{
    for (size_t i = 0; i < prob->m; ++i)
    {
        double ad = 0.0;
        for (int p = prob->Ap[i]; p < prob->Ap[i + 1]; ++p)
        {
            ad += prob->Ax[p] * d[prob->Ai[p]];
        }
        if (!is_recession_component(ad, prob->lhs[i], prob->rhs[i], tol))
        {
            printf("dual ray: (A d)_%zu = %g is not a recession direction of "
                   "[%g, %g]\n",
                   i, ad, prob->lhs[i], prob->rhs[i]);
            return false;
        }
    }

    for (size_t j = 0; j < prob->n; ++j)
    {
        if (!is_recession_component(d[j], prob->lbs[j], prob->ubs[j], tol))
        {
            printf("dual ray: d_%zu = %g is not a recession direction of "
                   "[%g, %g]\n",
                   j, d[j], prob->lbs[j], prob->ubs[j]);
            return false;
        }
    }

    double cd = 0.0;
    for (size_t j = 0; j < prob->n; ++j)
    {
        cd += prob->c[j] * d[j];
    }
    if (cd >= -tol)
    {
        printf("dual ray: c'd = %g is not negative\n", cd);
        return false;
    }
    return true;
}
