#ifndef ANF_VECTOR_ANF_H
#define ANF_VECTOR_ANF_H

#include "sparse_anf.h"
#include "simplify.h"
#include <vector>
#include <cstdint>

// Vector Boolean function: k components over n variables.
// Goal: find M,b minimizing |Union_i supp(g_i)|.
class VectorANF {
public:
    std::vector<SparseANF> components;
    int n = 0;
    int k = 0;

    VectorANF() = default;
    VectorANF(const std::vector<SparseANF>& components);

    int union_T() const;
    int T() const { return union_T(); }

    // Apply z = Mx ⊕ b to all components
    VectorANF substitute_affine(const std::vector<uint64_t>& M, uint64_t b) const;

    // Union framework: Z = Z1 ∪ Z2 where Z1 = X, Z2 = Mx⊕b
    VectorANF substitute_affine_union(const std::vector<uint64_t>& M,
                                      uint64_t b) const;

    // Replace x_var with Σ coeffs[j]*x_j (Boolean ring: OR for masks)
    VectorANF substitute_linear(int var, const std::vector<uint8_t>& coeffs) const;

    // Gradients of each component: returns k×n list
    std::vector<std::vector<SparseANF>> gradient() const;

    // Union of variables used across all components
    std::vector<int> variables_used() const;
};

// ---- Simplification ----

struct VectorSimplifyResult {
    VectorANF g;
    std::vector<uint64_t> M;  // m×n
    uint64_t b = 0;           // m-bit
};

// Greedy XOR merge for vector Boolean function
VectorSimplifyResult vector_greedy_merge(const VectorANF& vec,
                                         int max_iter = 50,
                                         bool verbose = false);

// Complement search for vector Boolean function
VectorSimplifyResult vector_simplify_by_complement(const VectorANF& vec,
                                                    bool verbose = false);

// Combined pipeline: complement + greedy merge
VectorSimplifyResult vector_simplify(const VectorANF& vec, bool verbose = false);

// Per-component simplification (strategy 2)
struct VectorPerCompResult {
    std::vector<SimplifyResult> per_component;
    VectorANF joint_vec;
    std::vector<uint64_t> M_joint;
    uint64_t b_joint;
};
VectorPerCompResult vector_simplify_per_component(const VectorANF& vec,
                                                   bool verbose = false);

#endif // ANF_VECTOR_ANF_H
