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

#ifndef FOLD_APPLY_H
#define FOLD_APPLY_H

#include "Bounds.h"
#include "Matrix.h"
#include "Tags.h"
#include "fold.h"

/* Builds the folded problem of 'maps' into a new matrix: the rows are the
   representatives with their coefficients summed over each column class, lhs,
   rhs and row tags those of the representatives, bounds and column tags those
   of the keepers, and c the keeper's cost times the class size, all compacted
   in place. The caller rebuilds the transpose. Returns FOLD_FOLDABLE with
   *folded set, or FOLD_FAILED when an allocation fails, in which case *folded
   is NULL and nothing is touched. */
FoldReason fold_apply(const Matrix *A, double *lhs, double *rhs, RowTag *row_tags,
                      Bound *bounds, ColTag *col_tags, double *c,
                      const FoldMaps *maps, Matrix **folded);

/* Finds and applies the fold. Returns true when the problem was folded: *A is
   then the folded matrix (the given one freed), the line data is compacted in
   place and *maps holds the maps; the caller rebuilds AT. Otherwise nothing is
   touched and *maps is NULL. 'stats' may be NULL; 'pool' runs the
   refinement. */
bool fold_problem(Matrix **A, const Matrix *AT, double *lhs, double *rhs,
                  RowTag *row_tags, Bound *bounds, ColTag *col_tags, double *c,
                  const FoldOptions *options, struct ThreadPool *pool,
                  FoldMaps **maps, FoldStats *stats);

#endif // FOLD_APPLY_H
