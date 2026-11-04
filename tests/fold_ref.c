#include "fold_ref.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct RefPair
{
    uint64_t cls;
    uint64_t bits;
} RefPair;

/* qsort has no user pointer in C99. */
static const uint64_t *ref_key0, *ref_key1, *ref_key2;
static const int *ref_ptr;
static const RefPair *ref_pairs;
static const int *ref_cls;

static uint64_t ref_bits(double value)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits == 0x8000000000000000ULL ? 0 : bits;
}

static int cmp_u64(uint64_t a, uint64_t b)
{
    return (a < b) ? -1 : (a > b) ? 1 : 0;
}

static int cmp_start(const void *pa, const void *pb)
{
    int a = *(const int *) pa, b = *(const int *) pb;
    int r = cmp_u64(ref_key0[a], ref_key0[b]);
    if (r == 0) r = cmp_u64(ref_key1[a], ref_key1[b]);
    if (r == 0 && ref_key2 != NULL) r = cmp_u64(ref_key2[a], ref_key2[b]);
    if (r == 0) r = (a < b) ? -1 : (a > b) ? 1 : 0;
    return r;
}

static int cmp_pair(const void *pa, const void *pb)
{
    const RefPair *a = (const RefPair *) pa, *b = (const RefPair *) pb;
    int r = cmp_u64(a->cls, b->cls);
    if (r == 0) r = cmp_u64(a->bits, b->bits);
    return r;
}

/* Lines compare by (old class, length, sorted pairs), index tie-break. */
static int cmp_line(const void *pa, const void *pb)
{
    int a = *(const int *) pa, b = *(const int *) pb;
    if (ref_cls[a] != ref_cls[b])
    {
        return (ref_cls[a] < ref_cls[b]) ? -1 : 1;
    }
    int la = ref_ptr[a + 1] - ref_ptr[a], lb = ref_ptr[b + 1] - ref_ptr[b];
    if (la != lb)
    {
        return (la < lb) ? -1 : 1;
    }
    for (int k = 0; k < la; k++)
    {
        int r = cmp_pair(&ref_pairs[ref_ptr[a] + k], &ref_pairs[ref_ptr[b] + k]);
        if (r != 0)
        {
            return r;
        }
    }
    return (a < b) ? -1 : (a > b) ? 1 : 0;
}

/* idx sorted by cmp_line (is_line_cmp) or cmp_start; consecutive lines that
   compare equal apart from the tie-break form a group whose id is its first
   member. Returns 1 if any class id changed. */
static int assign_groups(int *idx, int n_lines, int *cls, int is_line_cmp)
{
    int changed = 0;
    int first = idx[0];
    for (int i = 0; i < n_lines; i++)
    {
        int line = idx[i];
        if (i > 0)
        {
            int prev = idx[i - 1];
            int same;
            if (is_line_cmp)
            {
                same = (ref_cls[prev] == ref_cls[line]) &&
                       (ref_ptr[prev + 1] - ref_ptr[prev] ==
                        ref_ptr[line + 1] - ref_ptr[line]);
                for (int k = 0; same && k < ref_ptr[line + 1] - ref_ptr[line]; k++)
                {
                    if (cmp_pair(&ref_pairs[ref_ptr[prev] + k],
                                 &ref_pairs[ref_ptr[line] + k]) != 0)
                    {
                        same = 0;
                    }
                }
            }
            else
            {
                same = (ref_key0[prev] == ref_key0[line]) &&
                       (ref_key1[prev] == ref_key1[line]) &&
                       (ref_key2 == NULL || ref_key2[prev] == ref_key2[line]);
            }
            if (!same)
            {
                first = line;
            }
        }
        if (cls[line] != first)
        {
            changed = 1;
        }
        cls[line] = first;
    }
    return changed;
}

static int all_singletons(int n_lines, const int *cls)
{
    for (int i = 0; i < n_lines; i++)
    {
        if (cls[i] != i)
        {
            return 0;
        }
    }
    return 1;
}

static void start_colouring(int n_lines, const uint64_t *k0, const uint64_t *k1,
                            const uint64_t *k2, int *idx, int *cls)
{
    for (int i = 0; i < n_lines; i++)
    {
        idx[i] = i;
        cls[i] = -1;
    }
    ref_key0 = k0;
    ref_key1 = k1;
    ref_key2 = k2;
    qsort(idx, (size_t) n_lines, sizeof(int), cmp_start);
    assign_groups(idx, n_lines, cls, 0);
}

static int half_round(int n_lines, const int *ptr, const int *idx_arr,
                      const double *val, const int *other_cls, int *cls,
                      int *old_cls, RefPair *pairs, int *order)
{
    memcpy(old_cls, cls, (size_t) n_lines * sizeof(int));
    for (int i = 0; i < n_lines; i++)
    {
        for (int k = ptr[i]; k < ptr[i + 1]; k++)
        {
            pairs[k].cls = (uint64_t) other_cls[idx_arr[k]];
            pairs[k].bits = ref_bits(val[k]);
        }
        qsort(pairs + ptr[i], (size_t) (ptr[i + 1] - ptr[i]), sizeof(RefPair),
              cmp_pair);
        order[i] = i;
    }
    ref_ptr = ptr;
    ref_pairs = pairs;
    ref_cls = old_cls;
    qsort(order, (size_t) n_lines, sizeof(int), cmp_line);
    return assign_groups(order, n_lines, cls, 1);
}

int fold_ref(int m, int n, const int *Ap, const int *Ai, const double *Ax,
             const double *lhs, const double *rhs, const double *lb,
             const double *ub, const double *c, int max_rounds, int *row_class,
             int *col_class, int *rounds_out)
{
    int nnz = Ap[m];
    int converged = 0;
    int rounds = 0;
    int mn = (m > n) ? m : n;

    int *ATp = (int *) calloc((size_t) n + 1, sizeof(int));
    int *ATi = (int *) malloc(((size_t) nnz + 1) * sizeof(int));
    double *ATx = (double *) malloc(((size_t) nnz + 1) * sizeof(double));
    int *fill = (int *) calloc((size_t) n + 1, sizeof(int));
    for (int k = 0; k < nnz; k++)
    {
        ATp[Ai[k] + 1]++;
    }
    for (int j = 0; j < n; j++)
    {
        ATp[j + 1] += ATp[j];
    }
    for (int i = 0; i < m; i++)
    {
        for (int k = Ap[i]; k < Ap[i + 1]; k++)
        {
            int j = Ai[k];
            int pos = ATp[j] + fill[j]++;
            ATi[pos] = i;
            ATx[pos] = Ax[k];
        }
    }
    free(fill);

    uint64_t *k0 = (uint64_t *) malloc((size_t) mn * sizeof(uint64_t));
    uint64_t *k1 = (uint64_t *) malloc((size_t) mn * sizeof(uint64_t));
    uint64_t *k2 = (uint64_t *) malloc((size_t) mn * sizeof(uint64_t));
    int *order = (int *) malloc((size_t) mn * sizeof(int));
    int *old_cls = (int *) malloc((size_t) mn * sizeof(int));
    RefPair *pairs = (RefPair *) malloc(((size_t) nnz + 1) * sizeof(RefPair));

    for (int i = 0; i < m; i++)
    {
        k0[i] = ref_bits(lhs[i]);
        k1[i] = ref_bits(rhs[i]);
    }
    start_colouring(m, k0, k1, NULL, order, row_class);
    for (int j = 0; j < n; j++)
    {
        k0[j] = ref_bits(lb[j]);
        k1[j] = ref_bits(ub[j]);
        k2[j] = ref_bits(c[j]);
    }
    start_colouring(n, k0, k1, k2, order, col_class);

    if (nnz > 0)
    {
        while (!(all_singletons(m, row_class) && all_singletons(n, col_class)))
        {
            if (rounds >= max_rounds)
            {
                break;
            }
            rounds++;
            int changed = 0;
            changed |= half_round(m, Ap, Ai, Ax, col_class, row_class, old_cls,
                                  pairs, order);
            changed |= half_round(n, ATp, ATi, ATx, row_class, col_class, old_cls,
                                  pairs, order);
            if (!changed)
            {
                converged = 1;
                break;
            }
        }
        if (all_singletons(m, row_class) && all_singletons(n, col_class))
        {
            converged = 1;
        }
    }
    else
    {
        converged = 1;
    }

    free(ATp);
    free(ATi);
    free(ATx);
    free(k0);
    free(k1);
    free(k2);
    free(order);
    free(old_cls);
    free(pairs);
    *rounds_out = rounds;
    return converged;
}
