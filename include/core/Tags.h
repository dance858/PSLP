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

#ifndef CORE_TAGS_H
#define CORE_TAGS_H
#include <stddef.h>
#include <stdint.h>

// Parallel cols assumes these are 8 bits. For the macros below
// RowTag and ColTag should be the same too.
typedef uint8_t RowTag;
typedef uint8_t ColTag;

enum RowTag
{
    R_TAG_NONE = 0,
    R_TAG_LHS_INF = 1 << 0,
    R_TAG_RHS_INF = 1 << 1,
    R_TAG_EQ = 1 << 2,
    R_TAG_INACTIVE = 1 << 3,
    R_TAG_INFEAS = 1 << 4
};

enum ColTag
{
    C_TAG_NONE = 0,
    C_TAG_LB_INF = 1 << 0,
    C_TAG_LB_HUGE = 1 << 1,
    C_TAG_UB_INF = 1 << 2,
    C_TAG_UB_HUGE = 1 << 3,
    C_TAG_FIXED = 1 << 4,
    C_TAG_SUBSTITUTED = 1 << 5,
    // for parallel cols to work correctly it is important
    // that C_TAG_INACTIVE contains R_TAG_INACTIVE
    C_TAG_INACTIVE = (C_TAG_FIXED | C_TAG_SUBSTITUTED)
};

RowTag *new_rowtags(double *lhs, double *rhs, size_t n_rows);

// #define UPDATE_TAG(tag, new_tag) (tag |= new_tag)
// #define REMOVE_TAG(tag, old_tag) (tag &= ~old_tag)
// #define HAS_TAG(tag, check_tag) (tag & check_tag)
// #define RESET_TAG(tag, new_tag) (tag = new_tag)
// #define HAS_STATUS(status, check_status) (status & check_status)
// #define HAS_INF_TAG(tag) HAS_TAG((tag), (C_TAG_LB_INF | C_TAG_UB_INF))

#define UPDATE_TAG(tag, new_tag) ((tag) |= (RowTag) (new_tag))
#define REMOVE_TAG(tag, old_tag) ((tag) &= (RowTag) ~(RowTag) (old_tag))
#define HAS_TAG(tag, check_tag) (((tag) & (RowTag) (check_tag)) != 0)
#define RESET_TAG(tag, new_tag) ((tag) = (RowTag) (new_tag))
#define HAS_STATUS(status, check_status) (((status) & (check_status)) != 0)
#define HAS_INF_TAG(tag) HAS_TAG((tag), (C_TAG_LB_INF | C_TAG_UB_INF))

#endif // CORE_TAGS_H
