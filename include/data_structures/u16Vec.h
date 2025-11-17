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

#include <stdint.h>

// This macro defines a vector of uint16_t.
DEFINE_VECTOR(uint16_t, u16)

__attribute__((unused)) static void print_u16Vec(const u16Vec *vec)
{
    for (int i = 0; i < vec->len; ++i)
    {
        printf("%d ", vec->data[i]);
    }
    printf("\n");
}
