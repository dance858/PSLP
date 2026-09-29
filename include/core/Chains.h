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

#ifndef CORE_CHAINS_H
#define CORE_CHAINS_H

/* depth of a node dropped to break a cycle */
#define CHAINS_DROPPED (-2)

/* Consider a graph with n nodes. Node i points to succ[i], or to nothing when
   succ[i] == -1, so following the pointers from a node traces a chain.
   Several nodes can point to the same successor, so their chains merge and
   share the rest of the path. A cycle is broken by dropping its node with the
   highest drop_priority. Nodes pointing at a dropped node treat it as
   pointing to nothing.

   Computes depth[i], the number of steps from node i to the last node of its
   chain. A dropped node gets depth CHAINS_DROPPED. The dropped nodes are
   stored in 'dropped' and their number is returned. 'order' stores the nodes
   that are not dropped, ordered so that every node comes before the node it
   points to. 'order' and 'dropped' need room for n entries. 'stamp' and
   'path' are scratch of n ints. */
int compute_chain_depths(int n, const int *succ, const int *drop_priority,
                         int *depth, int *order, int *dropped, int *stamp,
                         int *path);

#endif // CORE_CHAINS_H
