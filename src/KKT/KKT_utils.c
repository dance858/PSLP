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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#include "KKT_utils.h"
#include <math.h>

void csr_matvec(const Matrix *A, const double *x,
                double *result)
{
    for (size_t k = 0; k < A->m; ++k)
    {
        double sum = 0.0;
        for (int j = A->p[k].start; j < A->p[k].end; ++j)
        {
            sum += A->x[j] * x[A->i[j]];
        }
        result[k] = sum;
    }
}

double norm2(const double *v, size_t len)
{
    double sum = 0.0;
    for (size_t i = 0; i < len; ++i)
    {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}
