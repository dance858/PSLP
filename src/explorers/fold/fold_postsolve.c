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

#include "fold_postsolve.h"

void fold_postsolve(const FoldMaps *maps, int n_rows, int n_cols,
                    const double *x_folded, const double *y_folded,
                    const double *z_folded, double *x, double *y, double *z)
{
    for (int col = 0; col < n_cols; ++col)
    {
        int folded_col = maps->col_fold_index[col];
        x[col] = x_folded[folded_col];
        z[col] = z_folded[folded_col] / maps->col_class_size[folded_col];
    }
    for (int row = 0; row < n_rows; ++row)
    {
        int folded_row = maps->row_fold_index[row];
        y[row] = y_folded[folded_row] / maps->row_class_size[folded_row];
    }
}
