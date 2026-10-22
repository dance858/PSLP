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

#ifndef COLOUR_REFINEMENT_INTERNAL_H
#define COLOUR_REFINEMENT_INTERNAL_H

/* State of the colour refinement, shared by colour_refinement.c and
   colour_refinement_verify.c. */

#include "Matrix.h"
#include "colour_refinement.h"
#include "hash.h"
#include "thread_pool.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct GroupTable GroupTable;

/* Capacity, in lines, of the buffer of a thread of the parallel push. A thread
   collects the lines it reaches in its buffer and appends them to the frontier
   in one go when the buffer is full and when the push ends. */
#define PUSH_BUFFER_SIZE 1024

/* For each class, we store the signature its members had when it was last split
   and its number of members. */
typedef struct RefineClass
{
    uint64_t signature;
    int size;
} RefineClass;

/* A line is a candidate when its current class has more than one line and its new
   signature differs from the stored signature of its current class. */
typedef struct RefineCandidate
{
    uint64_t signature;
    int class_id;
    int line;
} RefineCandidate;

/* The rows or the columns. The signature of a line is the sum over its
   entries of a hash of (class of the other end, value bits). */

/* The side of the matrix, either the rows or the columns.

    * n_lines:     number of lines (rows or columns)
    * p, i, x:     pointer representation of A or AT
    * class_id:    class id of every line
    * signature:   latest signature of every line (not updated if the line is
                   alone in its class). Bit 0 marks a line the push reached,
                   and a candidate between the count and the scatter pass
                   (colour_refinement.c).
    * classes:     signature and size of every class, indexed by class id
    * n_classes:   number of classes, the ids are 0, ..., n_classes - 1
    * frontier:    lines whose signature the last push changed
    * changed:     lines that changed class in the last half-round
    * changed_old_class: class of each changed line before the half-round
    * changed_entries:   entries of the changed lines
    * nonsingleton_entries: entries of the lines in classes with more than
                            one line
*/
typedef struct RefineSide
{
    int n_lines;
    const RowRange *p;
    const int *i;
    const double *x;
    int *class_id;        /* size: n_lines */
    uint64_t *signature;  /* size: n_lines */
    RefineClass *classes; /* size: n_lines */
    int n_classes;
    int *frontier; /* size: n_lines */
    int n_frontier;
    int *changed;           /* size: n_lines */
    int *changed_old_class; /* size: n_lines */
    int n_changed;
    int changed_entries;
    int nonsingleton_entries;
} RefineSide;

/* The candidates of a range of consecutive class ids, at positions [start, end)
   of the candidate array. All candidates of a class are in one bucket. One
   thread groups it. */
typedef struct RefineBucket
{
    int start;
    int end;
    int n_groups;
    int n_new_classes;
    int n_changed_lines;
    int n_changed_entries;
    int singleton_entries; /* Entries of the lines left alone in a class. */
    int first_new_class;   /* Id of its first new class. */
    int change_start;      /* Position of its first changed line. */
} RefineBucket;

typedef struct RefineWorkspace
{
    /* Chunk bookkeeping (chunk c covers the lines chunk_bounds[c], ...,
       chunk_bounds[c + 1] - 1. */
    int n_threads;
    int max_chunks;
    int *chunk_bounds;  /* size: max_chunks + 1 */
    int *chunk_entries; /* Entries read by each chunk. */

    /* Buckets bookkeeping */
    RefineBucket *buckets;
    int max_buckets;
    int n_buckets;
    int *bucket_counts;          /* size: max_chunks * max_buckets */
    int *bucket_offsets;         /* size: max_buckets + 1 */
    RefineCandidate *candidates; /* size: max lines of a side */

    /* Group bookkeeping within a bucket. A group is a set of candidates of a bucket
       that have have the same old class and the same new signature.

          * candidate_group[ii] is the group of candidate ii
          * group_first[ii] is the position of the first candidate of group ii.
          * group_size[group] is the number of candidates in a group.
          * group_class[group] is the class the lines of a group get, NEW_CLASS
             until its new class exists.
          * tables[thread] is the hash table for the groups of a bucket, used by
            thread 'thread'.
    */
    int *candidate_group;
    int *group_first;
    int *group_size;
    int *group_class;
    GroupTable *tables;

    /* Scratch indexed by class. */
    int *class_seen;          /* Scratch of the statistics. */
    int *class_largest_group; /* Largest group of the class in its bucket. */

    /* Lines reached by the push, buffered per thread (PUSH_BUFFER_SIZE each). */
    int **push_buffers;
    int *push_counts;

    /* for keeping track of the work budget */
    int64_t work_done;
} RefineWorkspace;

/* Contribution of one entry to the signature of a line. An entry joins a row
   and a column. 'bits' is its label and 'other_class' the class of the line
   at its other end: the column of the entry when the signature of a row is
   computed, the row when the signature of a column is computed. Even, so bit
   0 of a signature is free. */
static inline uint64_t entry_hash(int other_class, uint64_t bits)
{
    return mix(bits * 0xff51afd7ed558ccdULL +
               (uint64_t) other_class * 0x9e3779b97f4a7c15ULL)
           << 1;
}

/* We use bit 0 of a stored signature as a mark. Between half-rounds it is set
   on the lines the push reached, the frontier. Within a half-round the count
   pass leaves it set on the candidates and the scatter pass clears it. */
#define SIGNATURE_MARK ((uint64_t) 1)

static inline bool has_mark(uint64_t signature)
{
    return (signature & SIGNATURE_MARK) != 0;
}

static inline uint64_t without_mark(uint64_t signature)
{
    return signature & ~SIGNATURE_MARK;
}

static inline int line_size(const RefineSide *side, int line)
{
    return side->p[line].end - side->p[line].start;
}

/* Groups the candidates of a bucket and chooses the keepers of its classes
   (colour_refinement_group.c). False when the table of 'thread' cannot
   grow. */
bool refine_group_bucket(RefineSide *side, RefineWorkspace *work,
                         RefineBucket *bucket, int thread);

/* Creates the new classes of a bucket, from its first new id on, and moves
   their lines into them. */
void refine_apply_bucket(RefineSide *side, RefineWorkspace *work,
                         const RefineBucket *bucket);

/* The lines of 'from' that changed class push the change of their id into the
   signatures of 'to' (colour_refinement_push.c). The lines reached become the
   frontier of 'to'. */
void refine_push(const RefineSide *from, RefineSide *to, RefineWorkspace *work,
                 ThreadPool *pool);

struct ColourRefinement
{
    RefineSide *rows;
    RefineSide *cols;
    RefineWorkspace *work;
    ThreadPool *pool;
    int nnz;
};

#endif // COLOUR_REFINEMENT_INTERNAL_H
