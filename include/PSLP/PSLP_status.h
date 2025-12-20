#ifndef PSLP_STATUS_H
#define PSLP_STATUS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

    typedef uint8_t PresolveStatus;

    enum PresolveStatus_
    {
        UNCHANGED = 0,
        REDUCED = 1 << 0,
        INFEASIBLE = 1 << 1,
        UNBNDORINFEAS = INFEASIBLE | (1 << 2),
    };

#ifdef __cplusplus
}
#endif

#endif // PSLP_STATUS_H
