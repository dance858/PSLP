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

#ifndef HASH_H
#define HASH_H

#include <stdint.h>
#include <string.h>

/* Bit pattern of a value with -0.0 read as 0.0. */
static inline uint64_t value_bits(double value)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits == 0x8000000000000000ULL ? 0 : bits;
}

/* The splitmix64 finaliser. */
static inline uint64_t mix(uint64_t x)
{
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

/* Hash of two values. */
static inline uint64_t hash_two_values(double a, double b)
{
    return mix(value_bits(a) ^ mix(value_bits(b)));
}

/* Hash of three values. */
static inline uint64_t hash_three_values(double a, double b, double c)
{
    return mix(value_bits(a) ^ hash_two_values(b, c));
}

#endif // HASH_H
