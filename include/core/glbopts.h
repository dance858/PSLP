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

#ifndef GLB_H_GUARD
#define GLB_H_GUARD

#include "PSLP_inf.h"

#define INF PSLP_INF
#define IS_POS_INF(x) ((x) >= PSLP_INF)
#define IS_NEG_INF(x) ((x) <= -PSLP_INF)
#define IS_ABS_INF(x) (IS_POS_INF(x) || IS_NEG_INF(x))

#ifdef TESTING
#define EXTRA_ROW_SPACE 2
#define EXTRA_MEMORY_RATIO 1
#else
#define EXTRA_ROW_SPACE 4
// this value on EXTRA_MEMORY_RATIO matters for neos-5093327-huahum
#define EXTRA_MEMORY_RATIO 2
#endif

#define RETURN_IF_INFEASIBLE(x)                                                     \
    if ((x) == INFEASIBLE) return INFEASIBLE
#define RETURN_IF_UNBNDORINFEAS(x)                                                  \
    if ((x) == UNBNDORINFEAS) return UNBNDORINFEAS
#define RETURN_IF_NOT_UNCHANGED(x)                                                  \
    if ((x) != UNCHANGED) return x

#define SIZE_INACTIVE_ROW -1
#define SIZE_INACTIVE_COL -1
#define MAX_RATIO_PIVOT 1e3
#define PSLP_presolve_VERSION "0.0.1"

#endif
