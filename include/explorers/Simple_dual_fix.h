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

#ifndef SIMPLE_DUAL_FIX_H
#define SIMPLE_DUAL_FIX_H

#include "PSLP_status.h"

// forward declaration
struct Problem;

/* Uses locks and objective coefficients to fix variables at one of their
   bounds. Note that a variable may be fixed to infinity, making the
   constraints it appears in redundant. See Gurobi paper Section 4.4.

   When this function finishes, the problem is up to date, and no extra
   synchronization is needed. */
PresolveStatus simple_dual_fix(struct Problem *prob);

#endif // SIMPLE_DUAL_FIX_H
