#ifndef PSLP_PSLP_INFS_H
#define PSLP_PSLP_INFS_H

#ifdef __cplusplus
extern "C"
{
#endif

#define PSLP_INF 1e20
#define IS_POS_INF(x) ((x) >= PSLP_INF)
#define IS_NEG_INF(x) ((x) <= -PSLP_INF)
#define IS_ABS_INF(x) (IS_POS_INF(x) || IS_NEG_INF(x))

#ifdef __cplusplus
}
#endif

#endif /* PSLP_PSLP_INFS_H */