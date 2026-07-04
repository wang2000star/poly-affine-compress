#ifndef SPARSE_ANF_H
#define SPARSE_ANF_H

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>

// SparseANF: dict-based {mask: coeff} representation.
// mask is a uint64_t bitmask: bit j = 1 means variable z_j appears.
// Multiplication in Boolean ring: x_i^2 = x_i → masks combine with OR.
using AnfMap = std::unordered_map<uint64_t, int>;

class SparseANF {
    AnfMap terms_;
    int n_;  // number of variables

public:
    SparseANF() : n_(0) {}
    explicit SparseANF(int n) : n_(n) {}
    SparseANF(const AnfMap& terms, int n);

    // Query
    int T() const { return (int)terms_.size(); }
    int T_nonlinear() const;
    int degree() const;
    int n() const { return n_; }
    bool is_zero() const { return terms_.empty(); }
    bool is_one() const;
    int eval_mask(uint64_t x) const;
    const AnfMap& terms() const { return terms_; }
    std::string to_string(const std::vector<std::string>& var_names) const;

    // Mutators (in-place)
    void not_op();                              // g = 1 ⊕ g
    void xor_op(const SparseANF& other);        // g = g ⊕ other
    void xor_const();                           // g = g ⊕ 1  (same as not_op but without the name confusion)
    void and_op(const SparseANF& other);        // Boolean ring AND (OR for masks)
    void shift(int amount);                     // shift all masks by `amount` bits
    void complement_var(int j);                 // substitute z_j = z_j ⊕ 1
    void substitute_zi_eq_zi_xor_zj(int i, int j);  // substitute z_i = z_i ⊕ z_j

    // Substitute z = M_rows * x ⊕ b, where M_rows is (m_new × n_) matrix
    // Each entry of M_rows is a bitmask of old variables.
    SparseANF substitute_affine(const std::vector<uint32_t>& M_rows,
                                uint64_t b, int new_n) const;

    // Extract linear terms (degree=1), return L vector length n_
    // Removes them from this ANF.
    std::vector<int> extract_linear_terms();

    // Complement search: greedy bit-flipping O(T) per flip
    // Returns new (g, M_rows, b) where M_rows corresponds to complement pattern
    // The caller composes M = M_comp @ M_old, b = M_comp @ b_old ⊕ b_comp
    void complement_search_greedy(std::vector<uint32_t>& out_M_rows,
                                  uint64_t& out_b);

private:
    AnfMap _expand_monomial(uint64_t mask,
                            const std::vector<uint32_t>& M_rows,
                            uint64_t b, int new_n) const;
};

#endif // SPARSE_ANF_H
