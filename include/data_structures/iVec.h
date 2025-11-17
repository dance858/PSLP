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

#ifndef IVEC_H
#define IVEC_H

#include <stdio.h>

#include "Vec_macros.h"

// This macro defines a vector of integers.
DEFINE_VECTOR(int, i)

__attribute__((unused)) static void print_ivec(iVec *vec)
{
    for (size_t i = 0; i < vec->len; ++i)
    {
        printf("%d ", vec->data[i]);
    }
    printf("\n");
}

#endif // IVEC_H
