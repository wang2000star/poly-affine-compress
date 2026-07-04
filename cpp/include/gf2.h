#ifndef GF2_H
#define GF2_H

#include <cstdint>

// Rank of m×n matrix over GF(2), rows packed as uint32_t[n] (n ≤ 32).
int gf2_rank(const uint32_t* rows, int m, int n);

// Invert n×n matrix M over GF(2). Returns false if singular.
bool gf2_invert(const uint32_t* M, int n, uint32_t* M_inv);

// Right-inverse P (n×m) of m×n matrix M (n ≤ 32, m ≤ n). Returns false if rank < m.
// P is stored as n uint32_t values, where P[j] is an m-bit mask: P[j][k] = coefficient of z_k in x_j.
bool gf2_right_inverse(const uint32_t* M, int m, int n, uint32_t* P);

// Kernel basis of m×n matrix M. Returns dimension (n - rank).
// kernel[i] (for i in 0..dim-1) is an n-bit mask: kernel[i][j] = 1 if basis vector i has bit j set.
// pivot_cols (output): pivot column indices. pivot_cols_size = rank.
// Workspace: M_work needs m entries. Can be NULL.
int gf2_kernel_basis(const uint32_t* M, int m, int n, int* pivot_cols,
                     uint32_t* kernel, uint32_t* M_work);

#endif
