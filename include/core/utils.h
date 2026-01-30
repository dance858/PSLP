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

#ifndef CORE_UTILS_H
#define CORE_UTILS_H
#include "Tags.h"
#include "iVec.h"

/* Functions for shrinking (ie. removing data belonging to inactive
   cols/rows) */
void dPtr_shrink(double *ptr, const int *map, size_t len);
void iPtr_shrink(int *ptr, const int *map, size_t len);
void rowTagPtr_shrink(RowTag *ptr, const int *map, size_t len);
void colTagPtr_shrink(ColTag *ptr, const int *map, size_t len);
void shrink_idx_vector(iVec *vec, const int *map);

/* Other utility functions */
double get_max_abs(const double *vals, size_t len);

#endif // CORE_UTILS_H
