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

#ifndef GROUP_TABLE_H
#define GROUP_TABLE_H

/* A very basic open-addressing hash table with linear probing. It maps the key
   (class, signature) of a candidate to a group number. A slot holds a group
   number, or -1 when it is empty; the key of a group is read from its first
   candidate, so the table stores no keys. */

#include "colour_refinement_internal.h"
#include "hash.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

struct GroupTable
{
    int *slots;      /* A group number, or -1 when the slot is empty. */
    size_t capacity; /* Allocated slots. */
    size_t mask;     /* Slots in use minus one, a power of two minus one. */
};

/* Empties the table and gives it at least twice as many slots in use as
   'n_keys', a power of two and at least 16, so a search stops after a few
   slots and 'mask' reduces a hash to a slot. The table grows when it has too
   few slots and is reused otherwise. False when it cannot grow. */
bool group_table_reset(GroupTable *table, int n_keys);

void group_table_free(GroupTable *table);

/* Hash of the key of a candidate. */
static inline uint64_t group_table_hash(const RefineCandidate *key)
{
    return mix(key->signature ^ ((uint64_t) key->class_id * 0xd6e8feb86659fd93ULL));
}

/* Whether two candidates have the same key. */
static inline bool group_table_same_key(const RefineCandidate *a,
                                        const RefineCandidate *b)
{
    return a->signature == b->signature && a->class_id == b->class_id;
}

/* The slot of the group of 'key', or the empty slot where its group goes. The
   first candidate of group g is candidates[group_first[g]]. Starts at the slot
   of the hash of the key and moves on past the slots of other groups (linear
   probing). */
static inline size_t group_table_find(const GroupTable *table,
                                      const RefineCandidate *key,
                                      const RefineCandidate *candidates,
                                      const int *group_first)
{
    /* slot where the search for the candidate's group starts */
    size_t slot = group_table_hash(key) & table->mask;

    /* move past the slots of other groups until the key's group or an empty slot */
    while (table->slots[slot] >= 0 &&
           !group_table_same_key(&candidates[group_first[table->slots[slot]]], key))
    {
        slot = (slot + 1) & table->mask;
    }
    return slot;
}

static inline bool group_table_is_empty(const GroupTable *table, size_t slot)
{
    return table->slots[slot] < 0;
}

static inline int group_table_get(const GroupTable *table, size_t slot)
{
    return table->slots[slot];
}

static inline void group_table_set(GroupTable *table, size_t slot, int group)
{
    table->slots[slot] = group;
}

#endif // GROUP_TABLE_H
