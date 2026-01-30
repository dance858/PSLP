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

#ifndef TIMER_H
#define TIMER_H

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
typedef struct
{
    LARGE_INTEGER start, end;
} Timer;
static inline void clock_gettime_monotonic(Timer *timer, int is_start)
{
    if (is_start)
    {
        QueryPerformanceCounter(&timer->start);
    }
    else
    {
        QueryPerformanceCounter(&timer->end);
    }
}
static inline double get_elapsed_seconds(const Timer *timer)
{
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return (double) (timer->end.QuadPart - timer->start.QuadPart) /
           (double) freq.QuadPart;
}
#define clock_gettime(CLOCK, PTR)                                                   \
    clock_gettime_monotonic((Timer *) (PTR), ((PTR) == &((Timer *) (PTR))->start))
#define GET_ELAPSED_SECONDS(timer) get_elapsed_seconds(&(timer))
#else
#include <time.h>
typedef struct
{
    struct timespec start, end;
} Timer;
// #define GET_ELAPSED_SECONDS(timer)
//     (((timer).end.tv_sec - (timer).start.tv_sec) +
//      ((double) ((timer).end.tv_nsec - (timer).start.tv_nsec) * 1e-9))

#define GET_ELAPSED_SECONDS(timer)                                                  \
    ((double) ((timer).end.tv_sec - (timer).start.tv_sec) +                         \
     (double) ((timer).end.tv_nsec - (timer).start.tv_nsec) / 1e9)
#endif

#define RUN_AND_TIME(func, timer, time_variable, result_var, ...)                   \
    do                                                                              \
    {                                                                               \
        clock_gettime(CLOCK_MONOTONIC, &timer.start);                               \
        (result_var) = func(__VA_ARGS__);                                           \
        clock_gettime(CLOCK_MONOTONIC, &timer.end);                                 \
        (time_variable) += GET_ELAPSED_SECONDS(timer);                              \
    } while (0)

#endif // TIMER_H
