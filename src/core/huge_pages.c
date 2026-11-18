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

/* Transparent huge pages for large allocations (see Memory_wrapper.h). */
#if defined(__linux__)
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE /* madvise */
#endif
#include <stdint.h>
#include <sys/mman.h>
#include <unistd.h>
#endif
#include <stddef.h>

void ps_hint_huge_pages(void *p, size_t bytes);

void ps_hint_huge_pages(void *p, size_t bytes)
{
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    if (p == NULL) return;
    /* madvise wants whole pages: advise the pages that lie inside the block */
    const uintptr_t page = (uintptr_t) sysconf(_SC_PAGESIZE);
    const uintptr_t begin = ((uintptr_t) p + page - 1) & ~(page - 1);
    const uintptr_t end = ((uintptr_t) p + bytes) & ~(page - 1);
    if (end > begin)
    {
        (void) madvise((void *) begin, (size_t) (end - begin), MADV_HUGEPAGE);
    }
#else
    (void) p;
    (void) bytes;
#endif
}
