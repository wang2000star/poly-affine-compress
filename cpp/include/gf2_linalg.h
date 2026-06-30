#ifndef ANF_GF2_LINALG_H
#define ANF_GF2_LINALG_H

#include <cstdint>
#include <vector>

// GF(2) matrix: each uint64_t row packs up to 64 columns (bit i = column i)
using GF2Matrix = std::vector<uint64_t>;
using GF2Vector = uint64_t;

// Rank of m×n matrix (n ≤ 64)
int gf2_rank(const GF2Matrix& M, int n);

// Inverse of n×n matrix, returns empty vector on failure
GF2Matrix gf2_invert(const GF2Matrix& M, int n);

// Extend m×n matrix (m < n, full row rank) to n×n invertible
GF2Matrix gf2_extend_to_invertible(const GF2Matrix& M, int m, int n);

// Extend m×n matrix (m > n, full column rank) to m×m invertible by adding columns
GF2Matrix gf2_extend_columns_to_invertible(const GF2Matrix& M, int m, int n);

// Multiply GF(2) matrix (n×m) by vector (m-bit): x_j = Σ_i M[j][i]·v_i
// Each M[j] is an n-bit mask, v is an m-bit vector
uint64_t gf2_mat_vec_mul(const GF2Matrix& M, uint64_t v);

// Print matrix for debugging
void gf2_print_matrix(const GF2Matrix& M, int rows, int cols);

#endif // ANF_GF2_LINALG_H
