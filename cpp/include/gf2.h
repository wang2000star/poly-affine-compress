#ifndef GF2_H
#define GF2_H

#include <cstdint>

// Rank of m×n matrix over GF(2), rows packed as uint32_t[n] (n ≤ 32).
int gf2_rank(const uint32_t* rows, int m, int n);

// Invert n×n matrix M over GF(2). Returns false if singular.
bool gf2_invert(const uint32_t* M, int n, uint32_t* M_inv);

#endif
