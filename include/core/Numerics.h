/*
 * Copyright 2025 Daniel Cederberg
 *
 * This file is part of the PSLP project (LP Presolver).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CORE_NUMERICS_H
#define CORE_NUMERICS_H

#include <math.h>

#define FEAS_TOL 1e-6
#define BOUND_MARGINAL 0.5 * FEAS_TOL
#define ZERO_TOL 1e-10
#define ZERO_TOL_DUAL_POSTSOLVE 1e-6
#define HUGE_VAL_PS 1e7

// used for model clean up according to Gurobi paper
#define CLEAN1 1e-3
#define CLEAN2 1e-10

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#define IS_INTEGRAL(x) ((x) == (int) (x))

#define IS_EQUAL_FEAS_TOL(x, y) (ABS((x) - (y)) <= FEAS_TOL)
#define IS_ZERO_FEAS_TOL(x) (ABS(x) <= FEAS_TOL)
#define IS_GT_FEAS_TOL(x, y) ((x) - (y) >= FEAS_TOL)
#define IS_LT_FEAS_TOL(x, y) ((y) - (x) >= FEAS_TOL)
#define IS_HUGE(x) (ABS(x) >= HUGE_VAL_PS)

#define SET_ZERO_IF_SMALL_DUAL_POSTSOLVE(x)                                         \
    do                                                                              \
    {                                                                               \
        if (fabs(x) < ZERO_TOL_DUAL_POSTSOLVE) x = 0.0;                             \
    } while (0)

#endif // CORE_NUMERICS_H
