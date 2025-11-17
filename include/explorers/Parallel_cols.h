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

#ifndef PARALLEL_COLS_H
#define PARALLEL_COLS_H

#include "PresolveStatus.h"

// forward declaration
struct Problem;

/* This function finds parallel columns and merges them. It always
   returns UNCHANGED. When this function finishes, the problem is
   up to date, and no extra synchronization is needed.
*/
PresolveStatus remove_parallel_cols(struct Problem *prob);

#endif // PARALLEL_COLS_H
