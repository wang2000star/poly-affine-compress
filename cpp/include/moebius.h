#ifndef MOEBIUS_H
#define MOEBIUS_H

#include <cstdint>

// In-place packed Möbius transform (truth table ↔ ANF, self-inverse).
void moebius_packed(uint64_t* data, int n);
void moebius_packed_mt(uint64_t* data, int n, int n_threads);

// In-place packed upward Möbius transform (superset-sum, self-inverse).
// Applied to ANF coefficients, gives ANF of g(z) = f(z ⊕ 1).
// This is its own inverse (involution over GF(2)).
void moebius_upward_packed(uint64_t* data, int n);

#endif
