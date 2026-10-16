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
#include "thread_pool.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct GroupTable GroupTable;

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
                   alone in its class)
    * classes:     signature and size of every class, indexed by class id
    * n_classes:   number of classes, the ids are 0, ..., n_classes - 1
    * n_changed:   number of lines that changed class in the last half-round
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
    int n_changed;
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
    int first_new_class; /* Id of its first new class. */
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

    /* for keeping track of the work budget */
    int64_t work_done;
} RefineWorkspace;

/* Groups the candidates of a bucket and chooses the keepers of its classes
   (colour_refinement_group.c). False when the table of 'thread' cannot
   grow. */
bool refine_group_bucket(RefineSide *side, RefineWorkspace *work,
                         RefineBucket *bucket, int thread);

/* Creates the new classes of a bucket, from its first new id on, and moves
   their lines into them. */
void refine_apply_bucket(RefineSide *side, RefineWorkspace *work,
                         const RefineBucket *bucket);

struct ColourRefinement
{
    RefineSide *rows;
    RefineSide *cols;
    RefineWorkspace *work;
    ThreadPool *pool;
    int nnz;
};

#endif // COLOUR_REFINEMENT_INTERNAL_H
