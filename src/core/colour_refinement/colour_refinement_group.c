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

/* Grouping of one bucket of candidates.

   A bucket is a range of consecutive class ids together with the candidates
   of those classes, so all candidates of a class are in one bucket.

   The candidates of a bucket are grouped by (class, signature) with the hash
   table of the thread. In every class with candidates, one group keeps the
   class id: the lines that are not candidates if there are any, else the
   largest group. Every other group becomes a new class.

   A bucket writes only to its own classes and to its own range of the
   candidate and group arrays, so the pool groups buckets in parallel. */

#include "colour_refinement_internal.h"
#include "group_table.h"

/* The group_class of a group that becomes a new class, until it gets its id. */
#define NEW_CLASS -1

/* Takes the candidates of one bucket and sorts them into groups. Candidates
   with the same class and the same new signature end up in the same group. It
   does this in one pass with a hash table.*/
static int group_bucket(RefineWorkspace *work, const RefineBucket *bucket,
                        int thread)
{
    GroupTable *table = work->tables + thread;
    if (!group_table_reset(table, bucket->end - bucket->start))
    {
        return -1;
    }

    /* The groups of the bucket are numbered 0, 1, ... and described in its own
       range of the group arrays. */
    const RefineCandidate *candidates = work->candidates;
    int *group_first = work->group_first + bucket->start;
    int *group_size = work->group_size + bucket->start;
    int *candidate_group = work->candidate_group;
    int end = bucket->end;
    int n_groups = 0;

    for (int ii = bucket->start; ii < end; ++ii)
    {
        size_t slot =
            group_table_find(table, &candidates[ii], candidates, group_first);

        /* The first candidate of a group creates it. */
        if (group_table_is_empty(table, slot))
        {
            /* fill hash table with the new group */
            group_table_set(table, slot, n_groups);

            /* the first candidate of the new group */
            group_first[n_groups] = ii;
            group_size[n_groups] = 0;
            n_groups++;
        }

        int group = group_table_get(table, slot);
        group_size[group]++;

        /* record the group of candidate ii, numbered within the bucket */
        candidate_group[ii] = group;
    }
    return n_groups;
}

/* Decides for every class of a bucket which group, if any, keeps the class id
   and takes the lines of the other groups out of the class. Counts the new
   classes, the lines that change class and their entries, and the entries of
   the lines left alone in a class. The lines of a group have equal entry
   counts, and so have all lines of a class once the side has been refined. */
static void choose_keepers(RefineSide *side, RefineWorkspace *work,
                           RefineBucket *bucket)
{
    const RefineCandidate *candidates = work->candidates;
    const int *group_first = work->group_first;
    const int *group_size = work->group_size;
    int *largest = work->class_largest_group;
    int first_group = bucket->start;
    int end_group = first_group + bucket->n_groups;

    /* reset the largest group of every class that appears in the bucket */
    for (int group = first_group; group < end_group; ++group)
    {
        largest[candidates[group_first[group]].class_id] = -1;
    }

    /* Take every group out of its class and find the largest group of each
       class. */
    for (int group = first_group; group < end_group; ++group)
    {
        /* the class this group belonged to */
        int id = candidates[group_first[group]].class_id;

        /* remove the group's lines from the class */
        side->classes[id].size -= group_size[group];

        /* update the largest group of the class if this group is larger */
        if (largest[id] < 0 || group_size[group] > group_size[largest[id]])
        {
            largest[id] = group;
        }
    }

    /* A class left empty had only candidates, so its largest group goes back
       and keeps the id. Every other group becomes a new class. */
    int *group_class = work->group_class;
    int n_new_classes = 0;
    int n_changed_lines = 0;
    int n_changed_entries = 0;
    int singleton_entries = 0;
    for (int group = first_group; group < end_group; ++group)
    {
        const RefineCandidate *first = candidates + group_first[group];
        RefineClass *class_data = side->classes + first->class_id;
        int line_entries = line_size(side, first->line);
        if (largest[first->class_id] == group && class_data->size == 0)
        {
            class_data->size = group_size[group];
            class_data->signature = first->signature;
            group_class[group] = first->class_id;
        }
        else
        {
            group_class[group] = NEW_CLASS;
            n_new_classes++;
            n_changed_lines += group_size[group];
            n_changed_entries += group_size[group] * line_entries;
        }

        /* the group is alone in its class */
        if (group_size[group] == 1)
        {
            singleton_entries += line_entries;
        }

        /* one line that is not a candidate is left in the class, counted once
           through the largest group */
        if (largest[first->class_id] == group && class_data->size == 1 &&
            group_class[group] == NEW_CLASS)
        {
            singleton_entries += line_entries;
        }
    }
    bucket->n_new_classes += n_new_classes;
    bucket->n_changed_lines += n_changed_lines;
    bucket->n_changed_entries += n_changed_entries;
    bucket->singleton_entries += singleton_entries;
}

bool refine_group_bucket(RefineSide *side, RefineWorkspace *work,
                         RefineBucket *bucket, int thread)
{
    /* assign the candidates in 'bucket' to different groups */
    bucket->n_groups = group_bucket(work, bucket, thread);

    /* false when the table cannot grow */
    if (bucket->n_groups < 0)
    {
        return false;
    }

    /* decide, for every group in the bucket, whether its lines stay in their current
       class or move to a new class */
    choose_keepers(side, work, bucket);

    return true;
}

/* create the bucket's new classes and move every cand into its group's class. */
void refine_apply_bucket(RefineSide *side, RefineWorkspace *work,
                         const RefineBucket *bucket)
{
    RefineClass *classes = side->classes;
    const RefineCandidate *candidates = work->candidates;
    const int *group_first = work->group_first;
    const int *candidate_group = work->candidate_group;
    int *group_class = work->group_class;
    int *class_id = side->class_id;
    int *changed = side->changed;
    int *changed_old_class = side->changed_old_class;
    int next_class = bucket->first_new_class;
    int first_group = bucket->start;
    int end_group = first_group + bucket->n_groups;
    int end = bucket->end;

    /* create new classes for the groups that become new classes */
    for (int group = first_group; group < end_group; ++group)
    {
        if (group_class[group] == NEW_CLASS)
        {
            classes[next_class].signature = candidates[group_first[group]].signature;
            classes[next_class].size = work->group_size[group];
            group_class[group] = next_class++;
        }
    }

    /* assign each candidate to its group's class and record the lines that
       change class */
    int position = bucket->change_start;
    for (int ii = bucket->start; ii < end; ++ii)
    {
        int group = first_group + candidate_group[ii];
        int line = candidates[ii].line;
        if (group_class[group] != candidates[ii].class_id)
        {
            changed[position] = line;
            changed_old_class[position++] = candidates[ii].class_id;
        }
        class_id[line] = group_class[group];
    }
}
