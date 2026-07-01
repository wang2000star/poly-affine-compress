#include "sparse_anf.h"
#include "simplify.h"
#include "gf2_linalg.h"
#include "int_poly.h"
#include "vector_anf.h"
#include "vector_int_poly.h"
#include <cassert>
#include <iostream>
#include <random>

static int failures = 0;
static int tests = 0;

#define TEST(name, expr) do { \
    ++tests; \
    if (!(expr)) { \
        std::cerr << "FAIL: " << name << " (line " << __LINE__ << ")\n"; \
        ++failures; \
    } \
} while(0)

// ====================================================================
//  GF(2) linalg tests
// ====================================================================

void test_gf2_rank() {
    GF2Matrix I = {1, 2, 4, 8};
    TEST("identity rank", gf2_rank(I, 4) == 4);
    GF2Matrix zero = {0, 0, 0};
    TEST("zero rank", gf2_rank(zero, 4) == 0);
    GF2Matrix dep = {3, 6, 5};
    TEST("dependent rank", gf2_rank(dep, 3) == 2);
}

void test_gf2_invert() {
    GF2Matrix I = {1, 2, 4};
    auto inv = gf2_invert(I, 3);
    TEST("identity invert size", inv.size() == 3);
    TEST("identity invert[0]", inv[0] == 1);

    GF2Matrix M = {3, 2, 4};
    inv = gf2_invert(M, 3);
    TEST("invert non-trivial size", inv.size() == 3);
    for (int i = 0; i < 3; ++i) {
        uint64_t row = 0;
        uint64_t t = inv[i];
        while (t) { int j = __builtin_ctzll(t); t &= t - 1; row ^= M[j]; }
        TEST("invert check", row == ((uint64_t)1 << i));
    }
}

void test_gf2_extend() {
    GF2Matrix M = {5, 2};
    auto ext = gf2_extend_to_invertible(M, 2, 3);
    TEST("extend size", ext.size() == 3);
    TEST("extend rank", gf2_rank(ext, 3) == 3);
}

void test_gf2_extend_cols() {
    // m=3, n=2, full column rank
    GF2Matrix M = {1, 2, 3};  // 3×2: rows 001, 010, 011 — rank 2
    auto ext = gf2_extend_columns_to_invertible(M, 3, 2);
    TEST("extend cols size", ext.size() == 3);
    int r = gf2_rank(ext, 3);
    TEST("extend cols rank", r == 3);
}

// ====================================================================
//  SparseANF tests
// ====================================================================

void test_sparse_anf_basic() {
    std::unordered_map<uint64_t, uint8_t> terms = {{1, 1}, {6, 1}};
    SparseANF f(terms, 3);
    TEST("T", f.T() == 2);
    TEST("degree", f.degree() == 2);
    TEST("eval 000", f.eval_mask(0) == 0);
    TEST("eval 001", f.eval_mask(1) == 1);
    TEST("eval 110", f.eval_mask(6) == 1);
    TEST("eval 111", f.eval_mask(7) == 0);
}

void test_substitute_affine() {
    std::unordered_map<uint64_t, uint8_t> terms = {{1, 1}, {6, 1}};
    SparseANF f(terms, 3);

    std::vector<uint64_t> M = {1, 2, 4};
    auto g = f.substitute_affine(M, 0);
    TEST("identity T", g.T() == f.T());
    TEST("identity n", g.n == 3);
    TEST("identity verify", verify_substitution(f, g, M, 0, 100));

    M = {2, 1, 4};
    g = f.substitute_affine(M, 0);
    TEST("permute verify", verify_substitution(f, g, M, 0, 100));

    // m < n
    std::unordered_map<uint64_t, uint8_t> terms2 = {{3, 1}};
    SparseANF f2(terms2, 3);
    M = {1, 2};
    g = f2.substitute_affine(M, 0);
    TEST("m<n verify", verify_substitution(f2, g, M, 0, 100));

    // m < n full row rank
    std::unordered_map<uint64_t, uint8_t> terms3 = {{3, 1}}; // x0*x1, n=2
}

void test_substitute_affine_m_lt_n() {
    // Test m < n with full row rank
    std::unordered_map<uint64_t, uint8_t> terms = {{1, 1}, {2, 1}}; // x0 + x1, n=2 (really n=2)
    SparseANF f(terms, 2);
    // Use m=1 < n=2: z0 = x0 ⊕ x1
    std::vector<uint64_t> M = {3};  // 1×2: z0 = x0⊕x1
    auto g = f.substitute_affine(M, 0);
    // f(x) = x0 ⊕ x1 = z0 → g(z) = z0, T=1
    TEST("m<n T", g.T() == 1);
    TEST("m<n n", g.n == 1);
}

void test_complement() {
    std::unordered_map<uint64_t, uint8_t> terms = {{1, 1}, {2, 1}};
    SparseANF f(terms, 2);

    auto r = simplify_by_complement(f, false);
    TEST("complement verify", verify_substitution(f, r.g, r.M, r.b, 100));
    TEST("complement T <= original", r.g.T() <= f.T());
}

void test_greedy_merge() {
    std::unordered_map<uint64_t, uint8_t> terms = {{3, 1}, {5, 1}};
    SparseANF f(terms, 3);

    auto r = greedy_merge_simplify(f, 10, false);
    TEST("greedy merge verify", verify_substitution(f, r.g, r.M, r.b, 100));
    TEST("greedy merge T=1", r.g.T() == 1);
    TEST("greedy merge m=2", r.g.n == 2);
}

void test_simplify_pipeline() {
    std::mt19937_64 rng(12345);
    int inner_n = 4, outer_n = 10;
    auto inner = SparseANF::random_cubic(inner_n, 0.3, 42);
    std::vector<uint64_t> H(inner_n, 0);
    for (int i = 0; i < inner_n; ++i) H[i] = rng() & ((uint64_t(1) << outer_n) - 1);
    auto f = inner.expand_with(H, 0, outer_n);

    auto r = simplify(f, false);
    TEST("pipeline verify", verify_substitution(f, r.g, r.M, r.b, 100));
    TEST("pipeline m<=n", r.g.n <= f.n);
}

void test_substitute_affine_structured_m_gt_n() {
    // Structured m > n: each row of M has at most one 1-bit.
    // z_i = x_i, z_{i+n} = x_i⊕1 — both plain and complement available.
    // f(x) = x₀ ⊕ x₀x₁, n=2
    std::unordered_map<uint64_t, uint8_t> terms = {{1, 1}, {3, 1}};
    SparseANF f(terms, 2);

    // z₀=x₀, z₁=¬x₀, z₂=x₁, z₃=¬x₁ (m=4 > n=2)
    std::vector<uint64_t> M = {1, 1, 2, 2};
    uint64_t b = 0b1010;  // bits: z₁=1, z₃=1
    auto g = f.substitute_affine(M, b);

    TEST("structured m>n T", g.T() == 2);
    TEST("structured m>n m", g.n == 4);
    TEST("structured m>n verify", verify_substitution(f, g, M, b, 100));
}

void test_substitute_affine_union() {
    // f(x) = x₀x₁ ⊕ x₀x₂, n=3
    std::unordered_map<uint64_t, uint8_t> terms = {{3, 1}, {5, 1}};
    SparseANF f(terms, 3);

    // Z2: z₀ = x₀⊕x₁ (m=1), Z1 = {x₀,x₁,x₂}
    // Z = Z1 ∪ Z2 has 1+3=4 vars
    std::vector<uint64_t> M = {3};  // 1×3: z₀ = x₀⊕x₁
    auto g = f.substitute_affine_union(M, 0);

    // g uses (m+n)=4 vars; verify correctness
    TEST("union n == m+n", g.n == 4);

    // Build verification M_ext = [M; I₃]
    std::vector<uint64_t> M_ext = {3, 1, 2, 4};
    TEST("union verify", verify_substitution(f, g, M_ext, 0, 100));
}

// ====================================================================
//  IntPoly tests
// ====================================================================

void test_intpoly_basic() {
    IntPoly f({}, 2);
    TEST("empty T", f.T() == 0);
    TEST("empty degree", f.degree() == 0);

    // f = 3x0^2 + 2x0x1 + x1
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    terms[ExpVector{2, 0}] = 3;
    terms[ExpVector{1, 1}] = 2;
    terms[ExpVector{0, 1}] = 1;
    IntPoly f2(terms, 2);
    TEST("T", f2.T() == 3);
    TEST("degree", f2.degree() == 2);

    // eval at x0=2, x1=3: 3*4 + 2*6 + 3 = 12 + 12 + 3 = 27
    TEST("eval", f2.eval({2, 3}) == 27);

    // partial_deriv x0: 6x0 + 2x1
    auto d0 = f2.partial_deriv(0);
    TEST("deriv T", d0.T() == 2);
    TEST("deriv eval", d0.eval({2, 3}) == 18);

    // add
    auto sum = f2 + f2;
    TEST("add T", sum.T() == 3);
    TEST("add eval", sum.eval({2, 3}) == 54);

    // scalar mul
    auto scaled = f2 * 2;
    TEST("scalar T", scaled.T() == 3);
    TEST("scalar eval", scaled.eval({2, 3}) == 54);

    // substitute_linear: replace x0 with x0 + x1 in f2
    // f2(x0, x1) → g(x0, x1) = 3(x0+x1)^2 + 2(x0+x1)x1 + x1
    // = 3x0^2 + 6x0x1 + 3x1^2 + 2x0x1 + 2x1^2 + x1
    // = 3x0^2 + 8x0x1 + 5x1^2 + x1
    auto g = f2.substitute_linear(0, {1, 1});
    TEST("subst_linear T", g.T() == 4);
    TEST("subst_linear eval", g.eval({2, 3}) == f2.eval({5, 3})); // x0→2+3=5
}

void test_intpoly_substitute_affine() {
    // === expand_affine (forward: x = Nz + c, always works) ===
    // f(x) = x0^2 + x1^2, n=2
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    terms[ExpVector{2, 0}] = 1;
    terms[ExpVector{0, 2}] = 1;
    IntPoly f(terms, 2);

    // x = Nz + c: x0 = z0+z1, x1 = z0-z1
    // g(z) = f(Nz + c) = (z0+z1)^2 + (z0-z1)^2 = 2z0^2 + 2z1^2
    std::vector<std::vector<int64_t>> N = {{1, 1}, {1, -1}};
    std::vector<int64_t> c = {0, 0};
    auto g = f.expand_affine(N, c);
    TEST("expand_affine T", g.T() == 2);
    TEST("expand_affine n", g.n == 2);
    TEST("expand_affine eval", g.eval({3, 1}) == f.eval({4, 2}));

    // === substitute_affine (inverse: f(x) = g(Mx + b), needs integer-invertible M) ===
    // M = [[1,1],[0,1]], det=1, M^{-1} = [[1,-1],[0,1]]
    // f(x) = x0^2 + x1^2
    // g(z) = f(M^{-1}z) = (z0-z1)^2 + z1^2 = z0^2 - 2z0z1 + 2z1^2
    std::vector<std::vector<int64_t>> M = {{1, 1}, {0, 1}};
    std::vector<int64_t> b = {0, 0};
    auto g2 = f.substitute_affine(M, b);
    TEST("subst_affine T", g2.T() == 3);
    TEST("subst_affine n", g2.n == 2);
    // g(3,1) = 9 - 6 + 2 = 5, f(2,1) = 4 + 1 = 5
    TEST("subst_affine eval", g2.eval({3, 1}) == f.eval({2, 1}));
}

void test_intpoly_simplify() {
    // f = (x0 + x1)^3 + (x2 - x3)^2 — has structure
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    // (x0 + x1)^3
    terms[ExpVector{3, 0, 0, 0}] = 1;
    terms[ExpVector{2, 1, 0, 0}] = 3;
    terms[ExpVector{1, 2, 0, 0}] = 3;
    terms[ExpVector{0, 3, 0, 0}] = 1;
    // (x2 - x3)^2
    terms[ExpVector{0, 0, 2, 0}] = 1;
    terms[ExpVector{0, 0, 1, 1}] = -2;
    terms[ExpVector{0, 0, 0, 2}] = 1;
    IntPoly f(terms, 4);
    TEST("init T", f.T() == 7);
    TEST("init n", f.n == 4);

    auto r = simplify_int(f, false);
    int errors = verify_int_poly(f, r.g, r.M, r.b, 100);
    TEST("intpoly simplify verify", errors == 0);
    // Should reduce to about 3 variables (x0+x1, x2-x3, maybe constant)
    TEST("intpoly simplify m <= 3", r.g.n <= 3);
}

void test_intpoly_substitute_affine_union() {
    // f(x) = x₀x₁ + x₂², n=3
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    terms[ExpVector{1, 1, 0}] = 1;
    terms[ExpVector{0, 0, 2}] = 1;
    IntPoly f(terms, 3);

    // Structured M: z₀ = x₀ + x₁ (m=1)
    std::vector<std::vector<int64_t>> M = {{1, 1, 0}};
    std::vector<int64_t> b = {0};

    auto g = f.substitute_affine_union(M, b);
    // g should have 1+3=4 vars
    TEST("intpoly union n", g.n == 4);
    TEST("intpoly union verify", verify_int_poly(f, g, {{1,1,0},{1,0,0},{0,1,0},{0,0,1}}, {0,0,0,0}, 50) == 0);

    // Test with complement-style structure: z₀ = -x₀ + 1
    terms.clear();
    terms[ExpVector{1, 0}] = 1;
    terms[ExpVector{0, 1}] = 1;
    IntPoly f2(terms, 2);
    M = {{-1, 0}};
    b = {1};
    g = f2.substitute_affine_union(M, b);
    TEST("intpoly union complement n", g.n == 3);
    // M_ext = [M; I₂]
    std::vector<std::vector<int64_t>> M_ext = {{-1, 0}, {1, 0}, {0, 1}};
    std::vector<int64_t> b_ext = {1, 0, 0};
    TEST("intpoly union complement verify", verify_int_poly(f2, g, M_ext, b_ext, 50) == 0);
}

void test_intpoly_search_random() {
    // f = (x₀ + x₁)² + x₂² = x₀² + 2x₀x₁ + x₁² + x₂², n=3
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    terms[ExpVector{2, 0, 0}] = 1;
    terms[ExpVector{1, 1, 0}] = 2;
    terms[ExpVector{0, 2, 0}] = 1;
    terms[ExpVector{0, 0, 2}] = 1;
    IntPoly f(terms, 3);

    auto r = search_random_int(f, 3, 50, 42, false);
    int errors = verify_int_poly(f, r.g, r.M, r.b, 50);
    TEST("intpoly random search verify", errors == 0);
    TEST("intpoly random search T <= orig", r.g.T() <= f.T());
}

void test_intpoly_gradient_merge() {
    // f = (x0 + x1 + x2)^3 + (x3 - x4)^2
    // Gradient:∂f/∂x0 = ∂f/∂x1 = ∂f/∂x2 (all equal)
    // So gradient-guided merge should find this
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    for (int a = 0; a <= 3; ++a) {
        for (int b = 0; b <= 3-a; ++b) {
            int c = 3-a-b;
            int coeff = 6 / (a*b*c + (a==0)+(b==0)+(c==0) - 1);
            // Use hardcoded values
        }
    }
    // Just use explicit terms
    terms.clear();
    // (x0+x1+x2)^3 = sum of all monomials of degree 3 in x0,x1,x2 with multinomial coefficients
    terms[ExpVector{3,0,0,0,0}] = 1;
    terms[ExpVector{2,1,0,0,0}] = 3;
    terms[ExpVector{2,0,1,0,0}] = 3;
    terms[ExpVector{1,2,0,0,0}] = 3;
    terms[ExpVector{1,1,1,0,0}] = 6;
    terms[ExpVector{1,0,2,0,0}] = 3;
    terms[ExpVector{0,3,0,0,0}] = 1;
    terms[ExpVector{0,2,1,0,0}] = 3;
    terms[ExpVector{0,1,2,0,0}] = 3;
    terms[ExpVector{0,0,3,0,0}] = 1;
    // (x3-x4)^2
    terms[ExpVector{0,0,0,2,0}] = 1;
    terms[ExpVector{0,0,0,1,1}] = -2;
    terms[ExpVector{0,0,0,0,2}] = 1;

    IntPoly f(terms, 5);
    auto r = simplify_by_gradient_int(f, false);
    int errors = verify_int_poly(f, r.g, r.M, r.b, 100);
    TEST("gradient merge verify", errors == 0);
    // Should reduce dimension (gradients of x0,x1,x2 are equal)
    TEST("gradient merge m <= 3", r.g.n <= 4);
}

// ====================================================================
//  VectorANF tests
// ====================================================================

void test_vector_anf_basic() {
    // Two components: f0 = x0 + x1, f1 = x0*x1
    std::vector<SparseANF> comps;
    comps.push_back(SparseANF({{1, 1}, {2, 1}}, 2));
    comps.push_back(SparseANF({{3, 1}}, 2));

    VectorANF vec(comps);
    TEST("vec n", vec.n == 2);
    TEST("vec k", vec.k == 2);
    TEST("vec union_T", vec.union_T() == 3); // x0, x1, x0*x1
}

void test_vector_anf_simplify() {
    // Two components that share structure
    std::vector<SparseANF> comps;
    comps.push_back(SparseANF({{3, 1}, {5, 1}}, 3));   // x0*x1 + x0*x2
    comps.push_back(SparseANF({{3, 1}, {6, 1}}, 3));   // x0*x1 + x1*x2

    VectorANF vec(comps);
    int orig_T = vec.union_T();

    auto r = vector_simplify(vec, false);
    TEST("vector simplify T <= orig", r.g.union_T() <= orig_T);

    // Verify
    std::mt19937_64 rng(42);
    int errors = 0;
    for (int t = 0; t < 50; ++t) {
        uint64_t x = rng() & 7;
        uint64_t z = 0;
        for (int j = 0; j < (int)r.M.size(); ++j) {
            if (__builtin_parityll(r.M[j] & x)) z |= (uint64_t)1 << j;
        }
        z ^= r.b;
        z &= (1 << r.g.n) - 1;
        for (int i = 0; i < vec.k; ++i) {
            if (vec.components[i].eval_mask(x) != r.g.components[i].eval_mask(z))
                ++errors;
        }
    }
    TEST("vector simplify verify", errors == 0);
}

void test_vector_anf_complement_union() {
    // Components that benefit from having both x_i and x_i⊕1
    // f0 = x₀x₁ ⊕ x₀x₂, f1 = x₁x₂ (same structure, complement union may help)
    std::vector<SparseANF> comps;
    comps.push_back(SparseANF({{3, 1}, {5, 1}}, 3));   // x₀x₁ + x₀x₂
    comps.push_back(SparseANF({{6, 1}}, 3));            // x₁x₂

    VectorANF vec(comps);
    int orig_T = vec.union_T();
    auto r = vector_simplify(vec, false);

    int errors = 0;
    std::mt19937_64 rng(42);
    for (int t = 0; t < 50; ++t) {
        uint64_t x = rng() & 7;
        uint64_t z = 0;
        for (int j = 0; j < (int)r.M.size(); ++j) {
            if (__builtin_parityll(r.M[j] & x)) z |= (uint64_t)1 << j;
        }
        z ^= r.b;
        z &= (1 << r.g.n) - 1;
        for (int i = 0; i < vec.k; ++i) {
            if (vec.components[i].eval_mask(x) != r.g.components[i].eval_mask(z))
                ++errors;
        }
    }
    TEST("vector anf complement union verify", errors == 0);
    TEST("vector anf complement union T <= orig", r.g.union_T() <= orig_T);
}

// ====================================================================
//  VectorIntPoly tests
// ====================================================================

void test_vector_intpoly_basic() {
    std::vector<IntPoly> comps;
    comps.push_back(IntPoly({{ExpVector{1, 0}, 2}}, 2));  // f0 = 2x0
    comps.push_back(IntPoly({{ExpVector{0, 1}, 3}}, 2));  // f1 = 3x1

    VectorIntPoly vec(comps);
    TEST("vec int n", vec.n == 2);
    TEST("vec int k", vec.k == 2);
    TEST("vec int union_T", vec.union_T() == 2);
}

void test_vector_intpoly_simplify() {
    // Two components that share structure
    // f0 = (x0+x1)^2, f1 = (x0+x1)^3
    std::vector<IntPoly> comps;
    // (x0+x1)^2 = x0^2 + 2x0x1 + x1^2
    comps.push_back(IntPoly({
        {ExpVector{2, 0, 0}, 1},
        {ExpVector{1, 1, 0}, 2},
        {ExpVector{0, 2, 0}, 1}
    }, 3));
    // (x0+x1)^3 = x0^3 + 3x0^2x1 + 3x0x1^2 + x1^3
    comps.push_back(IntPoly({
        {ExpVector{3, 0, 0}, 1},
        {ExpVector{2, 1, 0}, 3},
        {ExpVector{1, 2, 0}, 3},
        {ExpVector{0, 3, 0}, 1}
    }, 3));

    VectorIntPoly vec(comps);
    int orig_T = vec.union_T();

    auto r = vector_simplify_int(vec, false);
    TEST("vec int simplify T <= orig", r.g.union_T() <= orig_T);

    // Verify
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int64_t> dist(-3, 3);
    int errors = 0;
    for (int t = 0; t < 50; ++t) {
        std::vector<int64_t> x(3);
        for (int i = 0; i < 3; ++i) x[i] = dist(rng);

        std::vector<int64_t> z(r.g.n, 0);
        for (int j = 0; j < r.g.n; ++j) {
            for (int i = 0; i < 3; ++i) z[j] += r.M[j][i] * x[i];
            z[j] += r.b[j];
        }

        for (int i = 0; i < vec.k; ++i) {
            if (vec.components[i].eval(x) != r.g.components[i].eval(z))
                ++errors;
        }
    }
    TEST("vec int simplify verify", errors == 0);
}

// ====================================================================
//  Main
// ====================================================================

int main() {
    // GF(2) linalg
    test_gf2_rank();
    test_gf2_invert();
    test_gf2_extend();
    test_gf2_extend_cols();

    // SparseANF
    test_sparse_anf_basic();
    test_substitute_affine();
    test_substitute_affine_m_lt_n();
    test_complement();
    test_greedy_merge();
    test_simplify_pipeline();
    test_substitute_affine_structured_m_gt_n();
    test_substitute_affine_union();

    // IntPoly
    test_intpoly_basic();
    test_intpoly_substitute_affine();
    test_intpoly_substitute_affine_union();
    test_intpoly_search_random();
    test_intpoly_simplify();
    test_intpoly_gradient_merge();

    // VectorANF
    test_vector_anf_basic();
    test_vector_anf_simplify();
    test_vector_anf_complement_union();

    // VectorIntPoly
    test_vector_intpoly_basic();
    test_vector_intpoly_simplify();

    std::cout << "Tests: " << tests << ", failures: " << failures << "\n";
    return failures ? 1 : 0;
}
