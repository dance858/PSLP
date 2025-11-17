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

#include "Problem.h"
#include "Constraints.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "State.h"
#include "Workspace.h"
#include "utils.h"

Problem *new_problem(Constraints *constraints, Objective *obj)
{
    Problem *problem = (Problem *) ps_malloc(1, sizeof(Problem));
    RETURN_PTR_IF_NULL(problem, NULL);

    problem->constraints = constraints;
    problem->obj = obj;

    return problem;
}

void free_problem(Problem *problem)
{
    RETURN_IF_NULL(problem);
    RETURN_IF_NULL(problem->constraints);
    State *data = problem->constraints->state;

    if (data)
    {
        free_work(data->work);
        PS_FREE(data->col_locks);
        PS_FREE(data->activities);
        PS_FREE(data->row_sizes);
        PS_FREE(data->col_sizes);
        free_state(data);
    }

    free_constraints(problem->constraints);
    free_objective(problem->obj);
    PS_FREE(problem);
}

Objective *objective_new(double *c)
{
    Objective *obj = (Objective *) ps_malloc(1, sizeof(Objective));
    RETURN_PTR_IF_NULL(obj, NULL);
    obj->c = c;
    obj->offset = 0.0;

    return obj;
}

void free_objective(Objective *obj)
{
    RETURN_IF_NULL(obj);
    PS_FREE(obj->c);
    PS_FREE(obj);
}

void objective_shrink(double *c, int *map, int len)
{
    dPtr_shrink(c, map, len);
}
void fix_var_in_obj(Objective *obj, int col, double value)
{
    obj->offset += obj->c[col] * value;
}

void problem_clean(Problem *prob, bool remove_all)
{
    Mapping *maps = prob->constraints->state->work->mappings;
    int n_cols_old = prob->constraints->n;
    int n_rows_old = prob->constraints->m;
    constraints_clean(prob->constraints, maps, remove_all);
    clean_state(prob->constraints->state, maps, n_rows_old, n_cols_old);
    objective_shrink(prob->obj->c, maps->cols, n_cols_old);
}

void sub_var_in_obj(Objective *obj, const double *vals, const int *cols, int len,
                    int k, double aik, double rhs)
{
    if (obj->c[k] == 0.0)
    {
        return;
    }

    double ratio = obj->c[k] / aik;
    for (int i = 0; i < len; ++i)
    {
        obj->c[cols[i]] -= ratio * vals[i];
    }

    obj->offset += rhs * ratio;
}
