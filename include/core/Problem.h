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

#ifndef CORE_PROBLEM_H
#define CORE_PROBLEM_H

#include <stdbool.h>
struct Constraints;

typedef struct
{
    double *c;
    double offset;
} Objective;

typedef struct Problem
{
    struct Constraints *constraints;
    Objective *obj;
} Problem;

/* Constructor and destructor */
Problem *new_problem(struct Constraints *constraints, Objective *obj);
void free_problem(Problem *problem);

/* Cleans the problem in the sense that all data associated with
   inactive rows and inactive columns is removed. Rows and columns
   are also re-indexed to take the removal of rows/columns into
   account. */
void problem_clean(Problem *problem, bool remove_all);

/* Constructor and destructor */
Objective *objective_new(double *c);
void free_objective(Objective *obj);

/* Updates the offset when a variable is fixed */
void fix_var_in_obj(Objective *obj, int col, double value);

/* Substitutes variable 'k' from the objective using row 'i'
   (assumed to be an equality). Row 'i' is specified by
   (vals, cols, len, rhs) */
void sub_var_in_obj(Objective *obj, const double *vals, const int *cols, int len,
                    int k, double aik, double rhs);

#endif // CORE_PROBLEM_H