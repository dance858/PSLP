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

#include "group_table.h"
#include "Memory_wrapper.h"
#include <string.h>

bool group_table_reset(GroupTable *table, int n_keys)
{
    size_t capacity = 16;
    while (capacity < 2 * (size_t) n_keys)
    {
        capacity <<= 1;
    }
    if (capacity > table->capacity)
    {
        if (!ps_grow(&table->slots, capacity, sizeof(int)))
        {
            return false;
        }
        table->capacity = capacity;
    }
    memset(table->slots, 0xff, capacity * sizeof(int));
    table->mask = capacity - 1;
    return true;
}

void group_table_free(GroupTable *table)
{
    PS_FREE(table->slots);
    table->capacity = 0;
}
