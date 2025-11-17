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

#ifndef DEBUG_MACROS_H
#define DEBUG_MACROS_H

#include "Tags.h"
#include <stdbool.h>
#include <stddef.h>

#ifdef NDEBUG
#define HUGE_BOUND_IS_OK
#define HUGE_BOUND_IS_NOT_OK
#define HUGE_BOUND_PARAM
#define HUGE_BOUND_ARG
#else
#define HUGE_BOUND_IS_OK , true
#define HUGE_BOUND_IS_NOT_OK , false
#define HUGE_BOUND_PARAM , bool huge_bound_ok
#define HUGE_BOUND_ARG , huge_bound_ok
#endif

// OBS: The macros should be wrapped inside a DEBUG when use

#ifndef NDEBUG
#define DEBUG(code)                                                                 \
    do                                                                              \
    {                                                                               \
        code;                                                                       \
    } while (0)
#else
#define DEBUG(code)                                                                 \
    do                                                                              \
    {                                                                               \
    } while (0)
#endif

static inline bool ARRAYS_EQUAL_INT(int *arr1, int *arr2, int size)
{
    for (int i = 0; i < size; i++)
    {
        if (arr1[i] != arr2[i]) return false;
    }
    return true;
}

static inline bool ARRAYS_EQUAL_DOUBLE(double *arr1, double *arr2, int size)
{
    for (int i = 0; i < size; i++)
    {
        if (arr1[i] != arr2[i]) return false;
    }
    return true;
}

static inline bool ARRAYS_EQUAL_ROWTAG(RowTag *arr1, const RowTag *arr2, int size)
{
    for (int i = 0; i < size; i++)
    {
        if (arr1[i] != arr2[i]) return false;
    }
    return true;
}

static inline bool ARRAYS_EQUAL_COLTAG(ColTag *arr1, const ColTag *arr2, int size)
{
    for (int i = 0; i < size; i++)
    {
        if (arr1[i] != arr2[i]) return false;
    }
    return true;
}

#define ASSERT_NO_INACTIVE_ROWS(rows, row_tags, len)                                \
    do                                                                              \
    {                                                                               \
        for (int i = 0; i < (len); i++)                                             \
        {                                                                           \
            assert(!HAS_TAG(row_tags[rows[i]], R_TAG_INACTIVE));                    \
        }                                                                           \
    } while (0)

#endif // DEBUG_MACROS_H