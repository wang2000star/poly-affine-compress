#ifndef ANF_VECTOR_INT_POLY_H
#define ANF_VECTOR_INT_POLY_H

#include "int_poly.h"
#include <vector>
#include <cstdint>

// Vector integer polynomial: k components over n variables.
class VectorIntPoly {
public:
    std::vector<IntPoly> components;
    int n = 0;
    int k = 0;

    VectorIntPoly() = default;
    VectorIntPoly(const std::vector<IntPoly>& components);

    int union_T() const;
    int T() const { return union_T(); }

    // Apply z = Mx + b to all components
    VectorIntPoly substitute_affine(const std::vector<std::vector<int64_t>>& M,
                                     const std::vector<int64_t>& b) const;

    // Union framework: Z = Z1 ∪ Z2 where Z1 = X, Z2 = Mx+b
    VectorIntPoly substitute_affine_union(const std::vector<std::vector<int64_t>>& M,
                                           const std::vector<int64_t>& b) const;

    // Union of variables used across all components
    std::vector<int> variables_used() const;
};

// ---- Simplification ----

struct VectorIntSimplifyResult {
    VectorIntPoly g;
    std::vector<std::vector<int64_t>> M;  // m×n
    std::vector<int64_t> b;              // m
};

// Greedy merge for vector integer polynomial
VectorIntSimplifyResult vector_greedy_merge_int(const VectorIntPoly& vec,
                                                 int max_iter = 50,
                                                 bool verbose = false);

// Random search with union for vector integer polynomial
VectorIntSimplifyResult vector_search_random_int(const VectorIntPoly& vec,
                                                   int max_m, int n_trials,
                                                   uint64_t seed = 1,
                                                   bool verbose = false);

// Combined pipeline: greedy merge + random search with union
VectorIntSimplifyResult vector_simplify_int(const VectorIntPoly& vec,
                                            bool verbose = false);

#endif // ANF_VECTOR_INT_POLY_H
