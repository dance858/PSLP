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

/* Colour refinement.

   A round refines the rows against the classes of the columns, then the
   columns against the classes of the rows. To refine one side, every line gets
   a signature: its pairs (class of the other end, value bits) in sorted order.
   The lines are sorted by (class, number of entries, signature) and every run
   of equal keys becomes a class. */

#include "colour_refinement.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include <stdlib.h>
#include <string.h>

typedef struct RefineSide
{
    int n_lines;
    const RowRange *ranges;
    const int *index; /* The line on the other side of each entry. */
    const double *values;
    int n_entries;
    uint64_t *colour; /* Start colours, written by the caller. */
    uint32_t *class_id;
    int n_classes;
    int duplicate_entries;
} RefineSide;

/* The sort key of a line. */
typedef struct RefineKey
{
    uint64_t key; /* The start colour, later the class id. */
    int count;    /* Entries of the line. */
    int line;
    const uint64_t *signature; /* count sorted pairs, NULL in the start pass. */
} RefineKey;

struct ColourRefinement
{
    RefineSide rows;
    RefineSide cols;
    RefineKey *keys;     /* Of the larger side. */
    uint64_t *signature; /* Of the side with more entries. */
    uint32_t *row_copy;  /* Class ids refined by verify. */
    uint32_t *col_copy;
};

/* Bit pattern of a value with -0.0 read as 0.0. */
static uint64_t value_bits(double value)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits == 0x8000000000000000ULL ? 0 : bits;
}

static int compare_pair(const void *a, const void *b)
{
    return memcmp(a, b, 2 * sizeof(uint64_t));
}

static int compare_key(const void *a, const void *b)
{
    const RefineKey *x = (const RefineKey *) a;
    const RefineKey *y = (const RefineKey *) b;
    if (x->key != y->key)
    {
        return x->key < y->key ? -1 : 1;
    }
    if (x->count != y->count)
    {
        return x->count < y->count ? -1 : 1;
    }
    if (!x->signature)
    {
        return 0;
    }
    return memcmp(x->signature, y->signature,
                  (size_t) x->count * 2 * sizeof(uint64_t));
}

static int line_count(const RefineSide *side, int line)
{
    return side->ranges[line].end - side->ranges[line].start;
}

/* Sorts the keys and gives every run of equal keys a class id. Returns the
   number of distinct keys the lines came with. */
static int assign_classes(RefineSide *side, RefineKey *keys)
{
    qsort(keys, (size_t) side->n_lines, sizeof(RefineKey), compare_key);
    int n_old = 0;
    int n_classes = 0;
    int duplicate_entries = 0;
    for (int ii = 0; ii < side->n_lines; ii++)
    {
        if (ii == 0 || keys[ii - 1].key != keys[ii].key)
        {
            n_old++;
        }
        if (ii == 0 || compare_key(&keys[ii - 1], &keys[ii]) != 0)
        {
            n_classes++;
        }
        else
        {
            duplicate_entries += keys[ii].count;
        }
        side->class_id[keys[ii].line] = (uint32_t) (n_classes - 1);
    }
    side->n_classes = n_classes;
    side->duplicate_entries = duplicate_entries;
    return n_old;
}

/* Groups the lines by start colour and number of entries. */
static void group_side(RefineSide *side, RefineKey *keys)
{
    for (int line = 0; line < side->n_lines; line++)
    {
        RefineKey key = {side->colour[line], line_count(side, line), line, NULL};
        keys[line] = key;
    }
    assign_classes(side, keys);
}

/* Splits the classes of a side by the classes of the other side. Returns the
   number of classes before the split. */
static int refine_side(RefineSide *side, const uint32_t *other_class,
                       RefineKey *keys, uint64_t *signature)
{
    uint64_t *pair = signature;
    for (int line = 0; line < side->n_lines; line++)
    {
        int count = line_count(side, line);
        RefineKey key = {side->class_id[line], count, line, pair};
        keys[line] = key;
        for (int idx = side->ranges[line].start; idx < side->ranges[line].end; idx++)
        {
            *pair++ = other_class[side->index[idx]];
            *pair++ = value_bits(side->values[idx]);
        }
        qsort(pair - 2 * count, (size_t) count, 2 * sizeof(uint64_t), compare_pair);
    }
    return assign_classes(side, keys);
}

ColourRefinementResult colour_refinement_run(ColourRefinement *refine,
                                             const ColourRefinementOptions *options,
                                             ColourRefinementStats *stats)
{
    RefineSide *rows = &refine->rows;
    RefineSide *cols = &refine->cols;
    ColourRefinementResult result = COLOUR_REFINEMENT_NOT_CONVERGED;
    int rounds = 0;
    double work = 0;
    group_side(rows, refine->keys);
    group_side(cols, refine->keys);
    while (rounds < options->max_rounds && work <= options->work_budget)
    {
        int n_row_classes = rows->n_classes;
        int n_col_classes = cols->n_classes;
        refine_side(rows, cols->class_id, refine->keys, refine->signature);
        refine_side(cols, rows->class_id, refine->keys, refine->signature);
        rounds++;
        /* A round touches the entries of both sides. */
        work += 2.0;
        if (rows->n_classes == n_row_classes && cols->n_classes == n_col_classes)
        {
            result = COLOUR_REFINEMENT_CONVERGED;
            break;
        }
    }
    stats->rounds = rounds;
    stats->work = work;
    return result;
}

bool colour_refinement_verify(ColourRefinement *refine)
{
    /* Copies of both sides are refined against the original classes. The
       partition is equitable when no class splits, whatever the ids. */
    RefineSide rows = refine->rows;
    RefineSide cols = refine->cols;
    rows.class_id = refine->row_copy;
    cols.class_id = refine->col_copy;
    memcpy(rows.class_id, refine->rows.class_id,
           (size_t) rows.n_lines * sizeof(uint32_t));
    memcpy(cols.class_id, refine->cols.class_id,
           (size_t) cols.n_lines * sizeof(uint32_t));
    int n_row_classes =
        refine_side(&rows, refine->cols.class_id, refine->keys, refine->signature);
    int n_col_classes =
        refine_side(&cols, refine->rows.class_id, refine->keys, refine->signature);
    return rows.n_classes == n_row_classes && cols.n_classes == n_col_classes;
}

const uint32_t *colour_refinement_row_classes(const ColourRefinement *refine)
{
    return refine->rows.class_id;
}

const uint32_t *colour_refinement_col_classes(const ColourRefinement *refine)
{
    return refine->cols.class_id;
}

int colour_refinement_n_row_classes(const ColourRefinement *refine)
{
    return refine->rows.n_classes;
}

int colour_refinement_n_col_classes(const ColourRefinement *refine)
{
    return refine->cols.n_classes;
}

uint64_t *colour_refinement_row_colours(ColourRefinement *refine)
{
    return refine->rows.colour;
}

uint64_t *colour_refinement_col_colours(ColourRefinement *refine)
{
    return refine->cols.colour;
}

int colour_refinement_duplicate_entries(const ColourRefinement *refine)
{
    return refine->rows.duplicate_entries + refine->cols.duplicate_entries;
}

void colour_refinement_lowest_member_ids(const ColourRefinement *refine, bool rows,
                                         int *scratch, int *ids)
{
    const RefineSide *side = rows ? &refine->rows : &refine->cols;
    for (int class_index = 0; class_index < side->n_classes; class_index++)
    {
        scratch[class_index] = -1;
    }
    for (int line = 0; line < side->n_lines; line++)
    {
        uint32_t class_index = side->class_id[line];
        if (scratch[class_index] < 0)
        {
            scratch[class_index] = line;
        }
        ids[line] = scratch[class_index];
    }
}

double colour_refinement_memory_estimate(size_t n_rows, size_t n_cols, size_t nnz)
{
    double lines = (double) n_rows + (double) n_cols;
    double larger = (double) (n_rows > n_cols ? n_rows : n_cols);
    /* Colour, class id and its copy per line; a key per line of the larger
       side; a signature pair per entry. */
    return (sizeof(uint64_t) + 2 * sizeof(uint32_t)) * lines +
           sizeof(RefineKey) * larger + 2 * sizeof(uint64_t) * (double) nnz;
}

static bool side_new(RefineSide *side, const Matrix *M)
{
    side->n_lines = (int) M->m;
    side->ranges = M->p;
    side->index = M->i;
    side->values = M->x;
    side->n_entries = 0;
    for (int line = 0; line < side->n_lines; line++)
    {
        side->n_entries += line_count(side, line);
    }
    side->colour =
        (uint64_t *) ps_calloc((size_t) side->n_lines + 1, sizeof(uint64_t));
    side->class_id =
        (uint32_t *) ps_calloc((size_t) side->n_lines + 1, sizeof(uint32_t));
    side->n_classes = 0;
    side->duplicate_entries = 0;
    return side->colour && side->class_id;
}

ColourRefinement *colour_refinement_new(const Matrix *A, const Matrix *AT)
{
    ColourRefinement *refine =
        (ColourRefinement *) ps_calloc(1, sizeof(ColourRefinement));
    if (!refine)
    {
        return NULL;
    }
    if (!side_new(&refine->rows, A) || !side_new(&refine->cols, AT))
    {
        colour_refinement_free(refine);
        return NULL;
    }
    int n_lines = (int) (A->m > A->n ? A->m : A->n);
    int n_entries = refine->rows.n_entries > refine->cols.n_entries
                        ? refine->rows.n_entries
                        : refine->cols.n_entries;
    refine->keys = (RefineKey *) ps_malloc((size_t) n_lines + 1, sizeof(RefineKey));
    refine->signature =
        (uint64_t *) ps_malloc(2 * (size_t) n_entries + 1, sizeof(uint64_t));
    refine->row_copy = (uint32_t *) ps_malloc(A->m + 1, sizeof(uint32_t));
    refine->col_copy = (uint32_t *) ps_malloc(A->n + 1, sizeof(uint32_t));
    if (!refine->keys || !refine->signature || !refine->row_copy ||
        !refine->col_copy)
    {
        colour_refinement_free(refine);
        return NULL;
    }
    return refine;
}

void colour_refinement_free(ColourRefinement *refine)
{
    if (!refine)
    {
        return;
    }
    PS_FREE(refine->rows.colour);
    PS_FREE(refine->rows.class_id);
    PS_FREE(refine->cols.colour);
    PS_FREE(refine->cols.class_id);
    PS_FREE(refine->keys);
    PS_FREE(refine->signature);
    PS_FREE(refine->row_copy);
    PS_FREE(refine->col_copy);
    PS_FREE(refine);
}
