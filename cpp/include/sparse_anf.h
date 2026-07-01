#ifndef ANF_SPARSE_ANF_H
#define ANF_SPARSE_ANF_H

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>
#include <random>

// Low-bit mask: bits ones. Safe for bits=64 (returns UINT64_MAX).
static inline uint64_t low_mask(int bits) {
    return bits >= 64 ? UINT64_MAX : ((uint64_t)1 << bits) - 1;
}

// Boolean function ANF over F₂.
// Masks are uint64_t: bit i = 1 means x_i appears in the monomial.
// n ≤ 64.  Boolean ring: x_i² = x_i (OR for combining monomials).
class SparseANF {
public:
    std::unordered_map<uint64_t, uint8_t> terms;  // {mask: coeff}, coeff ∈ {0,1}
    int n = 0;  // number of variables

    SparseANF() = default;
    SparseANF(const std::unordered_map<uint64_t, uint8_t>& terms, int n);
    SparseANF(std::unordered_map<uint64_t, uint8_t>&& terms, int n);

    // Constructors
    static SparseANF from_mask_set(const std::vector<uint64_t>& masks, int n);
    static SparseANF random_cubic(int n, double density = 0.15, uint64_t seed = 0);

    int T() const { return (int)terms.size(); }
    int degree() const;

    // Evaluate f(x) for x given as raw bitmask
    uint8_t eval_mask(uint64_t x) const;

    // Substitute z = Mx ⊕ b, compute g(z) = f(x).
    // M: m×n matrix (vector of m uint64_t, each is n-bit row mask)
    // b: m-bit vector (b[j] = (b>>j)&1)
    // Returns g(z) in m variables.
    SparseANF substitute_affine(const std::vector<uint64_t>& M, uint64_t b) const;

    // Union framework: Z = Z1 ∪ Z2 where Z1 = X (n vars), Z2 = Mx⊕b (m vars).
    // Builds M_ext = [M; I_n] ((m+n)×n), b_ext = [b; 0], calls substitute_affine.
    // Result g has (m+n) vars; m ≤ n recommended (Z1 provides n extra vars).
    SparseANF substitute_affine_union(const std::vector<uint64_t>& M,
                                      uint64_t b) const;

    // Variables used in any monomial
    std::vector<int> variables_used() const;

    // Partial derivative ∂f/∂x_var (Boolean derivative: f(x) ⊕ f(x⊕e_var))
    SparseANF partial_deriv(int var) const;

    // Gradient: list of partial derivatives
    std::vector<SparseANF> gradient() const;

    // Expand f(N(z ⊕ c)) where N is n×m, c is n-bit vector
    SparseANF expand_with(const std::vector<uint64_t>& N, uint64_t c, int m) const;

    // Debug: print terms
    void print() const;
    std::string to_string() const;

private:
    // Expand a single monomial through linear forms
    std::unordered_map<uint64_t, uint8_t> expand_monomial(
        uint64_t x_mask, const std::vector<uint64_t>& N, uint64_t c, int m) const;

    // Direct substitution for structured M (each row has ≤1 bit)
    SparseANF substitute_affine_structured(const std::vector<uint64_t>& M,
                                           uint64_t b) const;
};

// Apply merge x_i → x_i⊕x_j, return new terms
std::unordered_map<uint64_t, uint8_t> apply_xor_merge(
    const std::unordered_map<uint64_t, uint8_t>& terms, int i, int j);

// Verify g(Mx⊕b) == f(x) for random tests
bool verify_substitution(const SparseANF& f, const SparseANF& g,
                         const std::vector<uint64_t>& M, uint64_t b,
                         int n_tests = 50);

#endif // ANF_SPARSE_ANF_H
