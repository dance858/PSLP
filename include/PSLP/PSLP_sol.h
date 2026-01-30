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

#ifndef PSLP_SOL_H
#define PSLP_SOL_H

#ifdef __cplusplus
#include <cstddef> // size_t
extern "C"
{
#else
#include <stddef.h> // size_t
#endif

    typedef struct Solution
    {
        double *x;
        double *y;
        double *z;
        size_t dim_x;
        size_t dim_y;
    } Solution;

#ifdef __cplusplus
}
#endif

#endif // PSLP_SOL_H
