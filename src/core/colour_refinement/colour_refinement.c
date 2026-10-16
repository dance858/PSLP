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

   A line is a row or a column; a side is all rows or all columns, read from A or
   from its transpose.

   Every line has a class id and a signature: the sum over its entries of a
   hash of (class of the other end, value bits). Every class stores the
   signature its members had when it was last split. A half-round refines one
   side against the current classes of the other side: it recomputes the
   signature of every line in a class with more than one line, and the lines
   whose signature differs from that of their class are the candidates. Within
   a class, the candidates are grouped by signature. If the class has lines
   that are not candidates, those keep the class id and every group gets a new
   class; otherwise the largest group keeps the id. The half-rounds alternate
   between the sides until one changes no class.

   The candidates are partitioned by class id range into buckets, so all
   candidates of a class are in one bucket. The buckets are grouped
   independently, on the threads of the pool, each with a hash table
   (colour_refinement_group.c). */

#include "Memory_wrapper.h"
#include "colour_refinement_internal.h"
#include "group_table.h"
#include "hash.h"
#include <assert.h>
#include <string.h>

#define LINES_PER_CHUNK 4096
#define LINES_PER_BUCKET (1 << 14)

typedef enum RegroupMode
{
    REGROUP_START, /* Every line is a candidate with its start colour. */
    REGROUP_DENSE  /* The lines of non-singleton classes, signatures recomputed. */
} RegroupMode;

/* Context of the passes of one half-round. */
typedef struct RegroupPass
{
    RefineSide *side;
    const int *other_class;
    RegroupMode mode;
    RefineWorkspace *work;
    int shift; /* Class c is in bucket c >> shift. */
    int failed;
} RegroupPass;

/* The most buckets a side with 'n_lines' lines is split into: four per
   LINES_PER_BUCKET lines. */
static inline int max_bucket_count(int n_lines)
{
    return 4 * (n_lines / LINES_PER_BUCKET + 1);
}

/* The bucket that contains class 'id'. */
static inline int bucket_of(const RegroupPass *pass, int id)
{
    /* buckets are blocks of 2^shift consecutive ids */
    return id >> pass->shift;
}

/* Signature contribution of an entry with label 'bits' whose other end is in
   class 'other_class'. */
static inline uint64_t entry_hash(int other_class, uint64_t bits)
{
    return mix(bits * 0xff51afd7ed558ccdULL +
               (uint64_t) other_class * 0x9e3779b97f4a7c15ULL);
}

static inline int line_size(const RefineSide *side, int line)
{
    return side->p[line].end - side->p[line].start;
}

/* A line is a candidate if its class has more than one member and its new
   signature is not equal to the signature of the class. */
static inline bool is_candidate(const RefineSide *side, int line)
{
    const RefineClass *class_data = side->classes + side->class_id[line];
    return class_data->size > 1 && side->signature[line] != class_data->signature;
}

/* The signature of a line is the sum over its entries of the hash of (class of
   the other end, value bits). */
static uint64_t line_signature(const RefineSide *side, const int *other_class,
                               int line)
{
    uint64_t signature = 0;
    for (int idx = side->p[line].start; idx < side->p[line].end; ++idx)
    {
        signature += entry_hash(other_class[side->i[idx]], value_bits(side->x[idx]));
    }
    return signature;
}

/* The bucket counts of a chunk, one per bucket. */
static int *chunk_counts(const RefineWorkspace *work, int chunk)
{
    return work->bucket_counts + (size_t) chunk * (size_t) work->n_buckets;
}

/* Counts the candidates of a chunk of lines per bucket. In a dense half-round
   it first computes the signatures of the lines of non-singleton classes. */
static void count_chunk(void *context, int chunk, int thread)
{
    (void) thread;
    RegroupPass *pass = (RegroupPass *) context;
    RefineSide *side = pass->side;
    RefineWorkspace *work = pass->work;
    int *counts = chunk_counts(work, chunk);
    int start = work->chunk_bounds[chunk];
    int end = work->chunk_bounds[chunk + 1];

    if (pass->mode == REGROUP_START)
    {
        /* Every line is in class 0 and is a candidate. */
        assert(side->n_classes <= 1);
        counts[0] = end - start;
        return;
    }

    int entries = 0;
    for (int line = start; line < end; ++line)
    {
        /* which class does the current row (column) belong to */
        int id = side->class_id[line];

        /* a class with just one member can't be split */
        if (side->classes[id].size == 1)
        {
            continue;
        }

        /* compute the signature of the line given the classes of the other side */
        side->signature[line] = line_signature(side, pass->other_class, line);
        entries += line_size(side, line);

        /* count the candidate in the bucket of its class */
        if (is_candidate(side, line))
        {
            counts[bucket_of(pass, id)]++;
        }
    }

    work->chunk_entries[chunk] = entries;
}

/* Writes the candidates of a chunk of lines to their buckets. */
static void scatter_chunk(void *context, int chunk, int thread)
{
    (void) thread;
    RegroupPass *pass = (RegroupPass *) context;
    RefineSide *side = pass->side;
    RefineWorkspace *work = pass->work;
    int *positions = chunk_counts(work, chunk);
    int start = work->chunk_bounds[chunk];
    int end = work->chunk_bounds[chunk + 1];

    for (int line = start; line < end; ++line)
    {
        if (pass->mode == REGROUP_DENSE && !is_candidate(side, line))
        {
            continue;
        }
        int id = side->class_id[line];
        RefineCandidate candidate = {side->signature[line], id, line};
        work->candidates[positions[bucket_of(pass, id)]++] = candidate;
    }
}

/* Collects the candidates of the side into buckets. */
static void collect_candidates(RegroupPass *pass, ThreadPool *pool)
{
    RefineSide *side = pass->side;
    RefineWorkspace *work = pass->work;

    // --------------------------------------------------------------------------
    //                   split the lines into chunks
    // --------------------------------------------------------------------------
    int n_chunks =
        thread_pool_chunk_count(work->n_threads, side->n_lines, LINES_PER_CHUNK);
    thread_pool_chunk_bounds(side->n_lines, n_chunks, work->chunk_bounds);

    // ---------------------------------------------------------------------------
    //   choose the buckets so all candidates of a class are in one bucket
    // ----------------------------------------------------------------------------
    int bucket_goal = max_bucket_count(side->n_lines);
    int shift = 0;
    while ((side->n_classes >> shift) + 1 > bucket_goal)
    {
        shift++;
    }
    pass->shift = shift;
    work->n_buckets = (side->n_classes >> shift) + 1;
    assert(work->n_buckets <= work->max_buckets && n_chunks <= work->max_chunks);

    memset(work->bucket_counts, 0,
           (size_t) n_chunks * (size_t) work->n_buckets * sizeof(int));
    memset(work->chunk_entries, 0, (size_t) n_chunks * sizeof(int));

    /* run count_chunk on chunks in parallel */
    thread_pool_run(pool, count_chunk, pass, n_chunks);

    int n_candidates = thread_pool_bucket_offsets(
        work->bucket_counts, n_chunks, work->n_buckets, work->bucket_offsets);
    thread_pool_run(pool, scatter_chunk, pass, n_chunks);

    work->work_done += n_candidates;
    for (int chunk = 0; chunk < n_chunks; ++chunk)
    {
        work->work_done += work->chunk_entries[chunk];
    }
    for (int ii = 0; ii < work->n_buckets; ++ii)
    {
        RefineBucket *bucket = work->buckets + ii;
        memset(bucket, 0, sizeof(*bucket));
        bucket->start = work->bucket_offsets[ii];
        bucket->end = work->bucket_offsets[ii + 1];
    }
}

static void group_task(void *context, int bucket, int thread)
{
    RegroupPass *pass = (RegroupPass *) context;
    if (!refine_group_bucket(pass->side, pass->work, pass->work->buckets + bucket,
                             thread))
    {
        ps_atomic_store_int(&pass->failed, 1);
    }
}

static void apply_task(void *context, int bucket, int thread)
{
    (void) thread;
    RegroupPass *pass = (RegroupPass *) context;
    refine_apply_bucket(pass->side, pass->work, pass->work->buckets + bucket);
}

/* Refines one side. Returns false when a table cannot grow. */
static bool regroup(ColourRefinement *refine, RefineSide *side,
                    const int *other_class, RegroupMode mode)
{
    RefineWorkspace *work = refine->work;
    RegroupPass pass = {side, other_class, mode, work, 0, 0};

    collect_candidates(&pass, refine->pool);
    thread_pool_run(refine->pool, group_task, &pass, work->n_buckets);
    if (pass.failed)
    {
        return false;
    }

    /* The new classes of the buckets get consecutive ids. */
    int next_class = side->n_classes;
    int n_changed = 0;
    for (int ii = 0; ii < work->n_buckets; ++ii)
    {
        RefineBucket *bucket = work->buckets + ii;
        bucket->first_new_class = next_class;
        next_class += bucket->n_new_classes;
        n_changed += bucket->n_changed_lines;
    }
    thread_pool_run(refine->pool, apply_task, &pass, work->n_buckets);
    side->n_classes = next_class;
    side->n_changed = n_changed;
    return true;
}

/* Splits the lines of a side by their start colours. */
static bool start_side(ColourRefinement *refine, RefineSide *side,
                       const uint64_t *colours)
{
    /* put every line in class 0 */
    memset(side->class_id, 0, (size_t) side->n_lines * sizeof(int));
    side->classes[0].signature = 0;
    side->classes[0].size = side->n_lines;
    side->n_classes = 1;

    /* assign start colours as initial signature */
    memcpy(side->signature, colours, (size_t) side->n_lines * sizeof(uint64_t));

    return regroup(refine, side, NULL, REGROUP_START);
}

/* Sum of the entries of the lines that are not the lowest member of their
   class. 'seen' is scratch with n_classes entries. */
static int duplicate_entries(const RefineSide *side, int *seen)
{
    memset(seen, 0, (size_t) side->n_classes * sizeof(int));
    int entries = 0;
    for (int line = 0; line < side->n_lines; ++line)
    {
        if (seen[side->class_id[line]])
        {
            entries += line_size(side, line);
        }
        seen[side->class_id[line]] = 1;
    }
    return entries;
}

ColourRefinementResult colour_refinement_run(ColourRefinement *refine,
                                             const uint64_t *row_colours,
                                             const uint64_t *col_colours,
                                             const ColourRefinementOptions *options,
                                             int *row_classes, int *col_classes,
                                             ColourRefinementStats *stats)
{
    RefineSide *rows = refine->rows;
    RefineSide *cols = refine->cols;
    RefineWorkspace *work = refine->work;
    double budget = options->work_budget * (double) refine->nnz;
    ColourRefinementResult result = COLOUR_REFINEMENT_ALLOCATION_FAILED;
    int half_rounds = 0;
    work->work_done = 0;

    if (!start_side(refine, rows, row_colours) ||
        !start_side(refine, cols, col_colours))
    {
        goto done;
    }

    /* The first round refines both sides. */
    if (!regroup(refine, rows, cols->class_id, REGROUP_DENSE) ||
        !regroup(refine, cols, rows->class_id, REGROUP_DENSE))
    {
        goto done;
    }
    half_rounds = 2;

    /* 'from' was refined last and 'to' is refined against it. */
    RefineSide *from = cols;
    RefineSide *to = rows;
    while (true)
    {
        if (from->n_changed == 0)
        {
            result = COLOUR_REFINEMENT_CONVERGED;
            break;
        }
        if (half_rounds / 2 >= options->max_rounds ||
            (double) work->work_done > budget)
        {
            result = COLOUR_REFINEMENT_NOT_CONVERGED;
            break;
        }
        if (!regroup(refine, to, from->class_id, REGROUP_DENSE))
        {
            goto done;
        }
        half_rounds++;
        RefineSide *swap = from;
        from = to;
        to = swap;
    }

    memcpy(row_classes, rows->class_id, (size_t) rows->n_lines * sizeof(int));
    memcpy(col_classes, cols->class_id, (size_t) cols->n_lines * sizeof(int));
    stats->n_row_classes = rows->n_classes;
    stats->n_col_classes = cols->n_classes;
    stats->duplicate_entries = duplicate_entries(rows, work->class_seen) +
                               duplicate_entries(cols, work->class_seen);

done:
    stats->rounds = (half_rounds + 1) / 2;
    stats->work = (double) work->work_done / (double) refine->nnz;
    return result;
}

static void side_free(RefineSide *side)
{
    if (!side)
    {
        return;
    }
    PS_FREE(side->class_id);
    PS_FREE(side->signature);
    PS_FREE(side->classes);
    PS_FREE(side);
}

static RefineSide *side_new(const Matrix *M)
{
    RefineSide *side = (RefineSide *) ps_calloc(1, sizeof(RefineSide));
    if (!side)
    {
        return NULL;
    }
    size_t size = M->m + 1;
    side->n_lines = (int) M->m;
    side->p = M->p;
    side->i = M->i;
    side->x = M->x;
    side->class_id = (int *) ps_malloc(size, sizeof(int));
    side->signature = (uint64_t *) ps_malloc(size, sizeof(uint64_t));
    side->classes = (RefineClass *) ps_malloc(size, sizeof(RefineClass));
    if (!side->class_id || !side->signature || !side->classes)
    {
        side_free(side);
        return NULL;
    }
    return side;
}

static void workspace_free(RefineWorkspace *work)
{
    if (!work)
    {
        return;
    }
    if (work->tables)
    {
        for (int thread = 0; thread < work->n_threads; ++thread)
        {
            group_table_free(work->tables + thread);
        }
    }
    PS_FREE(work->candidates);
    PS_FREE(work->candidate_group);
    PS_FREE(work->group_first);
    PS_FREE(work->group_size);
    PS_FREE(work->group_class);
    PS_FREE(work->class_seen);
    PS_FREE(work->class_largest_group);
    PS_FREE(work->tables);
    PS_FREE(work->buckets);
    PS_FREE(work->chunk_bounds);
    PS_FREE(work->chunk_entries);
    PS_FREE(work->bucket_counts);
    PS_FREE(work->bucket_offsets);
    PS_FREE(work);
}

static RefineWorkspace *workspace_new(int max_lines, int n_threads)
{
    RefineWorkspace *work =
        (RefineWorkspace *) ps_calloc(1, sizeof(RefineWorkspace));
    if (!work)
    {
        return NULL;
    }
    size_t lines = (size_t) max_lines + 1;
    work->n_threads = n_threads;
    work->candidates = (RefineCandidate *) ps_malloc(lines, sizeof(RefineCandidate));
    work->candidate_group = (int *) ps_malloc(lines, sizeof(int));
    work->group_first = (int *) ps_malloc(lines, sizeof(int));
    work->group_size = (int *) ps_malloc(lines, sizeof(int));
    work->group_class = (int *) ps_malloc(lines, sizeof(int));
    work->class_seen = (int *) ps_malloc(lines, sizeof(int));
    work->class_largest_group = (int *) ps_malloc(lines, sizeof(int));
    work->tables = (GroupTable *) ps_calloc((size_t) n_threads, sizeof(GroupTable));

    work->max_buckets = max_bucket_count(max_lines);
    work->buckets =
        (RefineBucket *) ps_malloc((size_t) work->max_buckets, sizeof(RefineBucket));
    work->max_chunks =
        thread_pool_chunk_count(n_threads, max_lines, LINES_PER_CHUNK);
    work->chunk_bounds =
        (int *) ps_malloc((size_t) work->max_chunks + 1, sizeof(int));
    work->chunk_entries = (int *) ps_malloc((size_t) work->max_chunks, sizeof(int));
    work->bucket_counts = (int *) ps_malloc(
        (size_t) work->max_chunks * (size_t) work->max_buckets, sizeof(int));
    work->bucket_offsets =
        (int *) ps_malloc((size_t) work->max_buckets + 1, sizeof(int));

    if (!work->candidates || !work->candidate_group || !work->group_first ||
        !work->group_size || !work->group_class || !work->class_seen ||
        !work->class_largest_group || !work->tables || !work->buckets ||
        !work->chunk_bounds || !work->chunk_entries || !work->bucket_counts ||
        !work->bucket_offsets)
    {
        workspace_free(work);
        return NULL;
    }
    return work;
}

ColourRefinement *colour_refinement_new(const Matrix *A, const Matrix *AT,
                                        ThreadPool *pool)
{
    assert(A->nnz > 0);
    ColourRefinement *refine =
        (ColourRefinement *) ps_calloc(1, sizeof(ColourRefinement));
    if (!refine)
    {
        return NULL;
    }
    refine->pool = pool;
    refine->nnz = (int) A->nnz;
    int max_lines = (int) (A->m > A->n ? A->m : A->n);
    refine->rows = side_new(A);
    refine->cols = side_new(AT);
    refine->work = workspace_new(max_lines, thread_pool_n_threads(pool));
    if (!refine->rows || !refine->cols || !refine->work)
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
    side_free(refine->rows);
    side_free(refine->cols);
    workspace_free(refine->work);
    PS_FREE(refine);
}
