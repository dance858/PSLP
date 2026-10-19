#ifndef KKT_H
#define KKT_H

#include "PSLP_API.h"
#include <stdbool.h>

/* Optimality and certificate checks for the LP

       min  c^T x   s.t.  lhs <= A x <= rhs,   lbs <= x <= ubs,

   given as a ProblemData (CSR matrix, sides, bounds, objective). The
   dual convention is  A^T y + z = c,  y_i >= 0 on an active lhs, y_i <= 0 on
   an active rhs, z_j >= 0 at a lower bound, z_j <= 0 at an upper bound, i.e.
   the convention of postsolve().

   The residuals of a primal-dual point (x, y, z) are

     dual_res_abs   = || A^T y + z - c - P x ||_2
     primal_res_abs = || A x - proj_[lhs, rhs](A x) ||_2
     gap_abs        = | x^T P x + c^T x + p(-y; lhs, rhs) + p(-z; lbs, ubs) |
                      where p(s; l, u) = sum_i (s_i < 0 ? s_i l_i : s_i > 0 ? s_i u_i
   : 0) y_comp_slack   = || y - y(x) ||_2,  y_i(x) = y_i on an equality row, max(y_i,
   0) if (A x)_i = lhs_i, min(y_i, 0) if (A x)_i = rhs_i, 0 otherwise z_comp_slack =
   || z - z(x) ||_2,  analogously with the variable bounds bound_viol     = || x -
   proj_[lbs, ubs](x) ||_2

     dual_res_rel   = dual_res_abs   / (1 + ||c||_2 + ||P x||_2)
     primal_res_rel = primal_res_abs / (1 + ||(lhs_fin, rhs_fin)||_2)
     gap_rel        = gap_abs / (1 + |p(-y; lhs, rhs) + p(-z; lbs, ubs)| +
                                 |x^T P x + c^T x|)

   P is the matrix of the quadratic term (1/2) x^T P x of the objective; for
   an LP (Pp == NULL or Pnnz == 0) the P terms vanish and the multipliers are
   z = c - A^T y. Equality tests in the complementary-slackness terms use
   FEAS_TOL. */

typedef struct
{
    double primal_obj; // (1/2) x^T P x + c^T x + obj_offset
    double dual_res_abs;
    double primal_res_abs;
    double gap_abs;
    double dual_res_rel;
    double primal_res_rel;
    double gap_rel;
    double y_comp_slack;
    double z_comp_slack;
    double y_norm;     // ||y||_2, for scaling the comp-slack thresholds
    double z_norm;     // ||z||_2
    double bound_viol; // || x - proj_[lbs, ubs](x) ||_2
    double x_norm;     // ||x||_2, for scaling the bound-violation threshold
} KKTResiduals;

/* Wraps CSR arrays as a ProblemData with obj_offset = 0 and P = 0. The
   argument order matches new_presolver so a test builds it from its local
   arrays. */
ProblemData problem_from_csr(double *Ax, int *Ai, int *Ap, int n_rows, int n_cols,
                             int nnz, double *lhs, double *rhs, double *lbs,
                             double *ubs, double *c);

/* Same with the quadratic objective matrix P (full symmetric CSR). */
ProblemData problem_from_csr_qp(double *Ax, int *Ai, int *Ap, int n_rows, int n_cols,
                                int nnz, double *lhs, double *rhs, double *lbs,
                                double *ubs, double *c, double *Px, int *Pi, int *Pp,
                                int Pnnz);

void kkt_residuals(const ProblemData *prob, const double *x, const double *y,
                   const double *z, KKTResiduals *out);

/* True iff bound_viol / (1 + ||x||), primal_res_rel, dual_res_rel,
   z_comp_slack / (1 + ||z||), y_comp_slack / (1 + ||y||) and gap_rel are all
   <= tol. Prints the first violated quantity. */
bool is_kkt_point(const ProblemData *prob, const double *x, const double *y,
                  const double *z, double tol);

/* True iff y is a Farkas certificate of primal infeasibility: with
   z = -A^T y, every positive multiplier pairs with a finite rhs/ub, every
   negative one with a finite lhs/lb, and the support
   sum_{y_i > 0} rhs_i y_i + sum_{y_i < 0} lhs_i y_i + sum_{z_j > 0} ub_j z_j +
   sum_{z_j < 0} lb_j z_j is < -tol. This is the sign convention of
   postsolve_primal_infeas_ray (y_i >= 0 on an rhs-active row), the opposite
   of postsolve(). Prints the reason on failure. */
bool is_primal_ray_certificate(const ProblemData *prob, const double *y, double tol);

/* True iff d is a certificate of dual infeasibility (unboundedness): a
   recession direction of the constraint set ((A d)_i = 0 on two-sided rows,
   <= 0 on rhs-only rows, >= 0 on lhs-only rows, analogously for the bounds)
   with P d = 0 and c^T d < -tol. Prints the reason on failure. */
bool is_dual_ray_certificate(const ProblemData *prob, const double *d, double tol);

#endif // KKT_H
