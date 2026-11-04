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

/* Exact check of a partition. Every line is compared with the lowest member
   of its class: same number of entries and the same sorted pairs (class of
   the other end, value bits). */

#include "Memory_wrapper.h"
#include "colour_refinement_internal.h"
#include "hash.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* One entry of a line as seen by the check. */
typedef struct EntryPair
{
    uint64_t other_class;
    uint64_t bits;
} EntryPair;

static int compare_pairs(const void *a, const void *b)
{
    const EntryPair *x = (const EntryPair *) a;
    const EntryPair *y = (const EntryPair *) b;
    if (x->other_class != y->other_class)
    {
        return x->other_class < y->other_class ? -1 : 1;
    }
    if (x->bits != y->bits)
    {
        return x->bits < y->bits ? -1 : 1;
    }
    return 0;
}

/* Writes the sorted pairs of a line and returns their number. */
static int line_pairs(const RefineSide *side, const int *other_classes, int line,
                      EntryPair *pairs)
{
    int count = 0;
    for (int idx = side->p[line].start; idx < side->p[line].end; ++idx)
    {
        EntryPair pair = {(uint64_t) other_classes[side->i[idx]],
                          value_bits(side->x[idx])};
        pairs[count++] = pair;
    }
    qsort(pairs, (size_t) count, sizeof(EntryPair), compare_pairs);
    return count;
}

/* Checks the classes of one side against the classes of the other side. */
static ColourRefinementResult verify_side(const RefineSide *side, const int *classes,
                                          const int *other_classes)
{
    int n_lines = side->n_lines;
    int max_length = 0;
    for (int line = 0; line < n_lines; ++line)
    {
        int length = side->p[line].end - side->p[line].start;
        max_length = length > max_length ? length : max_length;
    }
    int *lowest = (int *) ps_malloc((size_t) n_lines + 1, sizeof(int));
    EntryPair *lowest_pairs =
        (EntryPair *) ps_malloc((size_t) max_length + 1, sizeof(EntryPair));
    EntryPair *pairs =
        (EntryPair *) ps_malloc((size_t) max_length + 1, sizeof(EntryPair));
    ColourRefinementResult result = COLOUR_REFINEMENT_ALLOCATION_FAILED;
    if (!lowest || !lowest_pairs || !pairs)
    {
        goto done;
    }

    for (int line = 0; line < n_lines; ++line)
    {
        lowest[line] = -1;
    }
    for (int line = 0; line < n_lines; ++line)
    {
        assert(classes[line] >= 0 && classes[line] < n_lines);
        if (lowest[classes[line]] < 0)
        {
            lowest[classes[line]] = line;
        }
    }

    result = COLOUR_REFINEMENT_EQUITABLE;
    for (int line = 0; line < n_lines; ++line)
    {
        int first = lowest[classes[line]];
        if (first == line)
        {
            continue;
        }
        int count = line_pairs(side, other_classes, line, pairs);
        int first_count = line_pairs(side, other_classes, first, lowest_pairs);
        if (count != first_count ||
            memcmp(pairs, lowest_pairs, (size_t) count * sizeof(EntryPair)) != 0)
        {
            result = COLOUR_REFINEMENT_NOT_EQUITABLE;
            break;
        }
    }

done:
    PS_FREE(lowest);
    PS_FREE(lowest_pairs);
    PS_FREE(pairs);
    return result;
}

ColourRefinementResult colour_refinement_verify(const ColourRefinement *refine,
                                                const int *row_classes,
                                                const int *col_classes)
{
    ColourRefinementResult result =
        verify_side(refine->rows, row_classes, col_classes);
    if (result != COLOUR_REFINEMENT_EQUITABLE)
    {
        return result;
    }
    return verify_side(refine->cols, col_classes, row_classes);
}
