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

#include "Matrix.h"
#include "Tags.h"

typedef struct RowView
{
    double *vals;
    int *cols;
    int *len;
    RowRange *range;
    double *lhs;
    double *rhs;
    RowTag *tag;
    int i;
} RowView;

typedef struct ConstRowView
{
    const double *vals;
    const int *cols;
    const int *len;
    const RowRange *range;
    const double *lhs;
    const double *rhs;
    const RowTag *tag;
    int i;
} ConstRowView;

typedef struct ColView
{
    double *vals;
    int *rows;
    int *len;
    RowRange *range;
    double *lb;
    double *ub;
    ColTag *tag;
    int k;
} ColView;

typedef struct ConstColView
{
    const double *vals;
    const int *rows;
    const int *len;
    const RowRange *range;
    const double *lb;
    const double *ub;
    const ColTag *tag;
    int k;
} ConstColView;

static inline RowView new_rowview(double *vals, int *cols, int *len,
                                  RowRange *row_range, double *lhs, double *rhs,
                                  RowTag *tag, int i)
{
    RowView view;
    view.vals = vals;
    view.cols = cols;
    view.len = len;
    view.range = row_range;
    view.lhs = lhs;
    view.rhs = rhs;
    view.tag = tag;
    view.i = i;
    return view;
}

static inline ConstRowView new_const_rowview(const double *vals, const int *cols,
                                             const int *len,
                                             const RowRange *row_range,
                                             const double *lhs, const double *rhs,
                                             const RowTag *tag, int i)
{
    ConstRowView view;
    view.vals = vals;
    view.cols = cols;
    view.len = len;
    view.range = row_range;
    view.lhs = lhs;
    view.rhs = rhs;
    view.tag = tag;
    view.i = i;
    return view;
}

static inline ColView new_colview(double *vals, int *rows, int *len,
                                  RowRange *col_range, double *lb, double *ub,
                                  ColTag *tag, int k)
{
    ColView view;
    view.vals = vals;
    view.rows = rows;
    view.len = len;
    view.range = col_range;
    view.lb = lb;
    view.ub = ub;
    view.tag = tag;
    view.k = k;
    return view;
}

static inline ConstColView new_const_colview(const double *vals, const int *rows,
                                             const int *len,
                                             const RowRange *col_range,
                                             const double *lb, const double *ub,
                                             const ColTag *tag, int k)
{
    ConstColView view;
    view.vals = vals;
    view.rows = rows;
    view.len = len;
    view.range = col_range;
    view.lb = lb;
    view.ub = ub;
    view.tag = tag;
    view.k = k;
    return view;
}