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

#ifndef CORE_BINARY_SEARCH_H
#define CORE_BINARY_SEARCH_H

#define BSEARCH_LINEAR_THRESHOLD 8

/* Find the index of 'target' in a sorted array 'arr' of length 'len'.
   Returns the relative index (0-based) if found, or -1 if not found.
   Falls back to linear scan for short arrays. */
static inline int sorted_find(const int *arr, int len, int target)
{
    if (len <= BSEARCH_LINEAR_THRESHOLD)
    {
        for (int i = 0; i < len; ++i)
        {
            if (arr[i] == target) return i;
            if (arr[i] > target) return -1;
        }
        return -1;
    }

    int lo = 0, hi = len - 1;
    while (lo <= hi)
    {
        int mid = lo + (hi - lo) / 2;
        if (arr[mid] == target)
        {
            return mid;
        }
        if (arr[mid] < target)
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid - 1;
        }
    }
    return -1;
}

/* Find the first index where arr[i] >= target in a sorted array of length 'len'.
   Returns 'len' if all elements are less than target.
   Falls back to linear scan for short arrays. */
static inline int sorted_lower_bound(const int *arr, int len, int target)
{
    if (len <= BSEARCH_LINEAR_THRESHOLD)
    {
        for (int i = 0; i < len; ++i)
        {
            if (arr[i] >= target)
            {
                return i;
            }
        }
        return len;
    }

    int lo = 0, hi = len;
    while (lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if (arr[mid] < target)
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid;
        }
    }
    return lo;
}

#endif // CORE_BINARY_SEARCH_H
