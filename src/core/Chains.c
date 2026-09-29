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

#include "Chains.h"

int compute_chain_depths(int n, const int *succ, const int *drop_priority,
                         int *depth, int *order, int *dropped, int *stamp, int *path)
{
    int n_dropped = 0;
    int walk = 0;

    for (int i = 0; i < n; ++i)
    {
        depth[i] = -1;
        stamp[i] = 0;
    }

    for (int start = 0; start < n; ++start)
    {
        // an unresolved start is walked again after a cycle was broken on its
        // chain (the dropped node may be the start itself, which ends the loop)
        while (depth[start] == -1)
        {
            // follow the pointers until nothing, a resolved node, a dropped
            // node, or a node stamped by this same walk (a cycle)
            int wid = ++walk;
            int len = 0;
            int node = start;
            int base = -1;
            int cycle_at = -1;
            for (;;)
            {
                stamp[node] = wid;
                path[len++] = node;
                int t = succ[node];

                if (t < 0 || depth[t] == CHAINS_DROPPED)
                {
                    break; // the chain ends here
                }
                if (depth[t] >= 0)
                {
                    base = depth[t];
                    break;
                }
                if (stamp[t] == wid)
                {
                    cycle_at = t;
                    break;
                }
                node = t;
            }

            if (cycle_at >= 0)
            {
                // the cycle is the path from cycle_at on: drop its node with
                // the highest drop priority
                int ci = 0;
                while (path[ci] != cycle_at)
                {
                    ci++;
                }
                int drop = path[ci];
                for (int q = ci + 1; q < len; ++q)
                {
                    if (drop_priority[path[q]] > drop_priority[drop])
                    {
                        drop = path[q];
                    }
                }
                depth[drop] = CHAINS_DROPPED;
                dropped[n_dropped++] = drop;
                continue;
            }

            // the path gets its depths from the end
            for (int q = len - 1; q >= 0; --q)
            {
                depth[path[q]] = ++base;
            }
        }
    }

    // order: stable counting sort of the resolved nodes by descending depth
    // (stamp[] is free now; a depth is below n)
    int maxdepth = -1;
    for (int i = 0; i < n; ++i)
    {
        if (depth[i] > maxdepth)
        {
            maxdepth = depth[i];
        }
    }
    int *count = stamp;
    for (int d = 0; d <= maxdepth; ++d)
    {
        count[d] = 0;
    }
    for (int i = 0; i < n; ++i)
    {
        if (depth[i] >= 0)
        {
            count[depth[i]]++;
        }
    }
    int run = 0;
    for (int d = maxdepth; d >= 0; --d)
    {
        int c = count[d];
        count[d] = run;
        run += c;
    }
    for (int i = 0; i < n; ++i)
    {
        if (depth[i] >= 0)
        {
            order[count[depth[i]]++] = i;
        }
    }
    return n_dropped;
}
