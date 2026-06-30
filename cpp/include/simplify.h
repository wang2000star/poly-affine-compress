#ifndef ANF_SIMPLIFY_H
#define ANF_SIMPLIFY_H

#include "sparse_anf.h"
#include <vector>
#include <cstdint>

// Result of simplification: g(z) = f(Mz ⊕ b) with g having m variables
struct SimplifyResult {
    SparseANF g;
    std::vector<uint64_t> M;  // m×n matrix, each row is n-bit mask
    uint64_t b = 0;           // m-bit vector
};

// ---- Greedy XOR Merge ----

// Simplify f by repeatedly merging x_i → x_i⊕x_j when it reduces T.
// M tracks accumulated transformation as m×n matrix.
SimplifyResult greedy_merge_simplify(const SparseANF& f,
                                     int max_iter = 100,
                                     bool verbose = false);

// ---- Complement Search ----

// For n ≤ 16: exhaustive Gray code over all 2ⁿ complement patterns.
// For n > 16: greedy bit flipping.
SimplifyResult simplify_by_complement(const SparseANF& f,
                                      bool verbose = false);

// ---- Random M,b search ----

// Try random M,b candidates, keeping the best result
SimplifyResult search_random(const SparseANF& f,
                             int max_m = 8,
                             int n_trials = 200,
                             uint64_t seed = 1,
                             bool verbose = false);

// ---- Combined Pipeline ----

// Complement → Greedy merge → Random search → Done
SimplifyResult simplify(const SparseANF& f, bool verbose = false);

#endif // ANF_SIMPLIFY_H
