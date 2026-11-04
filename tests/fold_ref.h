#ifndef FOLD_REF_H
#define FOLD_REF_H

/* Plain reference of the fold refinement for the unit tests: start colouring
   by the bit patterns of (lhs, rhs) and (lb, ub, c), then synchronous rounds
   (rows, then columns) in which two lines stay together iff their sorted
   sequences of (class of the other end, value bits) are identical. Class id =
   lowest member index. Returns 1 if converged within max_rounds. */
int fold_ref(int m, int n, const int *Ap, const int *Ai, const double *Ax,
             const double *lhs, const double *rhs, const double *lb,
             const double *ub, const double *c, int max_rounds, int *row_class,
             int *col_class, int *rounds_out);

#endif // FOLD_REF_H
