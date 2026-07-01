#ifndef ANF_INT_POLY_H
#define ANF_INT_POLY_H

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>

// Exponent vector: e[i] = degree of x_i
using ExpVector = std::vector<int>;

struct ExpHash {
    std::size_t operator()(const ExpVector& v) const {
        std::size_t h = 0;
        for (int e : v) {
            h ^= std::hash<int>{}(e) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

// Sparse multivariate integer polynomial over Z.
// terms: exponent vector -> coefficient (non-zero).
class IntPoly {
public:
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    int n = 0;  // number of variables

    IntPoly() = default;
    IntPoly(const std::unordered_map<ExpVector, int64_t, ExpHash>& terms, int n);
    IntPoly(std::unordered_map<ExpVector, int64_t, ExpHash>&& terms, int n);

    int T() const { return (int)terms.size(); }
    int degree() const;
    std::vector<int> deg_by_var() const;

    // Arithmetic
    IntPoly operator+(const IntPoly& other) const;
    IntPoly operator-(const IntPoly& other) const;
    IntPoly operator*(const IntPoly& other) const;
    IntPoly operator*(int64_t scalar) const;
    friend IntPoly operator*(int64_t scalar, const IntPoly& poly) {
        return poly * scalar;
    }

    // Calculus
    IntPoly partial_deriv(int var) const;
    std::vector<IntPoly> gradient() const;

    // Which variables actually appear
    std::vector<int> variables_used() const;

    // Replace x_var with Σ coeffs[j] * x_j (no constant term)
    IntPoly substitute_linear(int var, const std::vector<int64_t>& coeffs) const;

    // --- Affine substitution (z = Mx + b direction, matches SparseANF convention) ---
    // M is m×n, b is m-vector.
    // Computes g(z) in m variables such that f(x) = g(Mx + b).
    // Supports m ≤ n (right inverse) and m > n with structured or full-rank M.
    IntPoly substitute_affine(const std::vector<std::vector<int64_t>>& M,
                              const std::vector<int64_t>& b) const;

    // Union framework: Z = Z1 ∪ Z2 where Z1 = X (n vars), Z2 = Mx+b (m vars).
    // Builds M_ext = [M; I_n] ((m+n)×n), b_ext = [b; 0], calls substitute_affine.
    IntPoly substitute_affine_union(const std::vector<std::vector<int64_t>>& M,
                                    const std::vector<int64_t>& b) const;

    // --- Forward expansion (x = Nz + c direction, always works) ---
    // N is n×m, c is n-vector.
    // Computes g(z) = f(Nz + c) for ANY m (including m > n).
    // This is the "forward" substitution that never needs inversion.
    IntPoly expand_affine(const std::vector<std::vector<int64_t>>& N,
                          const std::vector<int64_t>& c) const;

    // Evaluate at given integer values
    int64_t eval(const std::vector<int64_t>& values) const;

    // Build polynomial Σ coeffs[j] * x_j in n variables
    static IntPoly linear_form_poly(const std::vector<int64_t>& coeffs, int n);

    std::string to_string() const;
    void print() const;

private:
    // Build polynomial (Σ row[j] * z_j + const_val) in m variables
    static IntPoly single_affine(const std::vector<int64_t>& row,
                                 int64_t const_val, int m);

    // Structured substitute_affine for M with ≤1 non-zero per row
    IntPoly substitute_affine_structured(
        const std::vector<std::vector<int64_t>>& M,
        const std::vector<int64_t>& b) const;
};

// ---- Simplification ----

struct IntSimplifyResult {
    IntPoly g;
    std::vector<std::vector<int64_t>> M;  // m×n
    std::vector<int64_t> b;              // m
};

// Try merge: replace x_i with x_i - k*x_j in f (inverse of basis merge)
IntPoly try_merge_int(const IntPoly& f, int i, int j, int64_t k);

// Greedy pairwise merge for integer polynomials
IntSimplifyResult greedy_merge_simplify_int(const IntPoly& f,
                                            int max_iter = 50,
                                            bool verbose = false);

// Gradient-guided simplification
IntSimplifyResult simplify_by_gradient_int(const IntPoly& f,
                                           bool verbose = false);

// Random M,b search for IntPoly (structured + small general matrices)
IntSimplifyResult search_random_int(const IntPoly& f, int max_m, int n_trials,
                                     uint64_t seed = 1, bool verbose = false);

// Combined pipeline: gradient-guided + exhaustive merge + random search
IntSimplifyResult simplify_int(const IntPoly& f, bool verbose = false);

// Verify f(x) == g(Mx + b) for random integer inputs
int verify_int_poly(const IntPoly& f, const IntPoly& g,
                    const std::vector<std::vector<int64_t>>& M,
                    const std::vector<int64_t>& b,
                    int n_tests = 100);

#endif // ANF_INT_POLY_H
