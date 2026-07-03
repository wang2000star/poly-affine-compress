#ifndef MOEBIUS_H
#define MOEBIUS_H

#include <cstdint>

// In-place packed Möbius transform (truth table ↔ ANF, self-inverse).
void moebius_packed(uint64_t* data, int n);
void moebius_packed_mt(uint64_t* data, int n, int n_threads);

#endif
