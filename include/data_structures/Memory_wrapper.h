/*
 * Copyright 2025-2026 Daniel Cederberg
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

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Allocations of at least this many bytes ask the kernel for transparent huge
   pages (Linux, madvise). First touch of a large fresh block is then one
   fault per 2 MB instead of one per 4 KB. */
#define PS_HUGE_PAGE_MIN_BYTES ((size_t) 4 << 20)

/* No-op where there is nothing to ask for. Never fails. */
void ps_hint_huge_pages(void *p, size_t bytes);

static inline void *ps_malloc(size_t n, size_t size)
{
    void *p = malloc(n * size);
    if (n * size >= PS_HUGE_PAGE_MIN_BYTES) ps_hint_huge_pages(p, n * size);
    return p;
}

static inline void *ps_calloc(size_t n, size_t size)
{
    void *p = calloc(n, size);
    if (n * size >= PS_HUGE_PAGE_MIN_BYTES) ps_hint_huge_pages(p, n * size);
    return p;
}

#define PS_FREE(p)                                                                  \
    do                                                                              \
    {                                                                               \
        free(p);                                                                    \
        p = NULL;                                                                   \
    } while (0)

static inline void *ps_realloc(void *p, size_t n, size_t size)
{
    void *q = realloc(p, n * size);
    if (n * size >= PS_HUGE_PAGE_MIN_BYTES) ps_hint_huge_pages(q, n * size);
    return q;
}

/* Reallocates the array whose pointer lives at 'pp' (the address of the array
   pointer) to n elements of size bytes. On failure the old block is kept and
   false is returned. */
static inline bool ps_grow(void *pp, size_t n, size_t size)
{
    void *old, *grown;
    memcpy(&old, pp, sizeof(old));
    grown = ps_realloc(old, n, size);
    if (grown == NULL)
    {
        return false;
    }
    memcpy(pp, &grown, sizeof(grown));
    return true;
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
