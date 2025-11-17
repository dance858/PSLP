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

#ifndef MEMORY_WRAPPER_H
#define MEMORY_WRAPPER_H

#include <stdio.h>
#include <stdlib.h>

static inline void *ps_malloc(int n, size_t size)
{
    return malloc(n * size);
}

static inline void *ps_calloc(int n, size_t size)
{
    return calloc(n, size);
}

#define PS_FREE(p)                                                                  \
    do                                                                              \
    {                                                                               \
        free(p);                                                                    \
        p = NULL;                                                                   \
    } while (0)

static inline void *ps_realloc(void *p, int n, size_t size)
{
    return realloc(p, n * size);
}

#define RETURN_PTR_IF_NULL(ptr, ret_val)                                            \
    do                                                                              \
    {                                                                               \
        if ((ptr) == NULL)                                                          \
        {                                                                           \
            return (ret_val);                                                       \
        }                                                                           \
    } while (0)

#define RETURN_IF_NULL(ptr)                                                         \
    do                                                                              \
    {                                                                               \
        if ((ptr) == NULL) return;                                                  \
    } while (0)

#endif // MEMORY_WRAPPER_H
