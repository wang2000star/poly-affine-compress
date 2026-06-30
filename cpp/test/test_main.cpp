#include "sparse_anf.h"
#include "simplify.h"
#include "gf2_linalg.h"
#include <cassert>
#include <iostream>
#include <random>

static int failures = 0;
static int tests = 0;

#define TEST(name, expr) do { \
    ++tests; \
    if (!(expr)) { \
        std::cerr << "FAIL: " << name << "\n"; \
        ++failures; \
    } \
} while(0)

void test_gf2_rank() {
    // Rank of identity
    GF2Matrix I = {1, 2, 4, 8};
    TEST("identity rank", gf2_rank(I, 4) == 4);

    // Rank of zero rows
    GF2Matrix zero = {0, 0, 0};
    TEST("zero rank", gf2_rank(zero, 4) == 0);

    // Rank of dependent rows
    GF2Matrix dep = {3, 6, 5};  // 011, 110, 101 -> rank 2 in 3 cols
    TEST("dependent rank", gf2_rank(dep, 3) == 2);
}

void test_gf2_invert() {
    // Invert identity
    GF2Matrix I = {1, 2, 4};
    auto inv = gf2_invert(I, 3);
    TEST("identity invert size", inv.size() == 3);
    TEST("identity invert[0]", inv[0] == 1);

    // Invert a known matrix
    GF2Matrix M = {3, 2, 4};  // rows: 011, 010, 100
    inv = gf2_invert(M, 3);
    TEST("invert non-trivial size", inv.size() == 3);

    // Check M^{-1} * M = I
    for (int i = 0; i < 3; ++i) {
        uint64_t row = 0;
        uint64_t t = inv[i];
        while (t) {
            int j = __builtin_ctzll(t);
            t &= t - 1;
            row ^= M[j];
        }
        TEST("invert check row " + std::to_string(i), row == ((uint64_t)1 << i));
    }
}

void test_gf2_extend() {
    GF2Matrix M = {5, 2};  // m=2, n=3: 101, 010
    auto ext = gf2_extend_to_invertible(M, 2, 3);
    TEST("extend size", ext.size() == 3);
    int r = gf2_rank(ext, 3);
    TEST("extend rank", r == 3);
}

void test_sparse_anf_basic() {
    // f = x0 + x1*x2
    std::unordered_map<uint64_t, uint8_t> terms = {{1, 1}, {6, 1}};
    SparseANF f(terms, 3);
    TEST("T", f.T() == 2);
    TEST("degree", f.degree() == 2);
    TEST("eval 000", f.eval_mask(0) == 0);
    TEST("eval 001", f.eval_mask(1) == 1);   // x0=1
    TEST("eval 110", f.eval_mask(6) == 1);  // x1=1,x2=1
    TEST("eval 111", f.eval_mask(7) == 0);  // x0+x1*x2 = 1+1*1 = 0
}

void test_substitute_affine() {
    // f = x0 + x1*x2
    std::unordered_map<uint64_t, uint8_t> terms = {{1, 1}, {6, 1}};
    SparseANF f(terms, 3);

    // Identity: z_i = x_i (M = I, b = 0)
    std::vector<uint64_t> M = {1, 2, 4};
    auto g = f.substitute_affine(M, 0);
    TEST("identity substitute T", g.T() == f.T());
    TEST("identity substitute n", g.n == 3);
    bool ok = verify_substitution(f, g, M, 0, 100);
    TEST("identity verify", ok);

    // Permutation: z0 = x1, z1 = x0, z2 = x2
    M = {2, 1, 4};
    g = f.substitute_affine(M, 0);
    ok = verify_substitution(f, g, M, 0, 100);
    TEST("permute verify", ok);

    // m < n: f = x0*x1 (depends on x0, x1 only), n=3, M = {1, 2}
    // z0 = x0, z1 = x1, g(z0, z1) = z0*z1
    std::unordered_map<uint64_t, uint8_t> terms2 = {{3, 1}};
    SparseANF f2(terms2, 3);
    M = {1, 2};
    g = f2.substitute_affine(M, 0);
    ok = verify_substitution(f2, g, M, 0, 100);
    TEST("m<n verify", ok);

    // m < n with non-trivial M: f = x0*x2 (depends on x2, but M maps z0 = x0, z1 = x2)
    // So M = {1, 4}, m=2, g(z0, z1) = z0*z1
    std::unordered_map<uint64_t, uint8_t> terms3 = {{5, 1}};  // x0*x2: bits 0 and 2
    SparseANF f3(terms3, 3);
    M = {1, 4};
    g = f3.substitute_affine(M, 0);
    ok = verify_substitution(f3, g, M, 0, 100);
    TEST("m<n non-trivial M verify", ok);
}

void test_complement() {
    // f = x0 + x1
    std::unordered_map<uint64_t, uint8_t> terms = {{1, 1}, {2, 1}};
    SparseANF f(terms, 2);

    auto r = simplify_by_complement(f, false);
    bool ok = verify_substitution(f, r.g, r.M, r.b, 100);
    TEST("complement simple verify", ok);
    // T should be at most 2 (no bigger than original)
    TEST("complement T <= original", r.g.T() <= f.T());
}

void test_greedy_merge() {
    // f = x0*x1 + x0*x2  ->  merge x2->x2⊕x1 to get x0*x1 + x0*(x1⊕x2) = x0*x1 + x0*x1 + x0*x2 = x0*x2
    // Actually let's test the known reduction
    std::unordered_map<uint64_t, uint8_t> terms = {{3, 1}, {5, 1}};
    SparseANF f(terms, 3);

    auto r = greedy_merge_simplify(f, 10, false);
    bool ok = verify_substitution(f, r.g, r.M, r.b, 100);
    TEST("greedy merge simple verify", ok);
    // Should reduce to 1 term
    TEST("greedy merge reduced T", r.g.T() == 1);
    TEST("greedy merge reduced m", r.g.n == 2);
}

void test_simplify_pipeline() {
    // LinStr-style: f = inner.expand_with(H, 0, outer_n)
    std::mt19937_64 rng(12345);
    int inner_n = 4, outer_n = 10;
    auto inner = SparseANF::random_cubic(inner_n, 0.3, 42);
    std::vector<uint64_t> H(inner_n, 0);
    for (int i = 0; i < inner_n; ++i) {
        H[i] = rng() & ((uint64_t(1) << outer_n) - 1);
    }
    auto f = inner.expand_with(H, 0, outer_n);

    auto r = simplify(f, false);
    bool ok = verify_substitution(f, r.g, r.M, r.b, 100);
    TEST("pipeline linstr verify", ok);
    TEST("pipeline reduced m", r.g.n <= f.n);
}

int main() {
    test_gf2_rank();
    test_gf2_invert();
    test_gf2_extend();
    test_sparse_anf_basic();
    test_substitute_affine();
    test_complement();
    test_greedy_merge();
    test_simplify_pipeline();

    std::cout << "Tests: " << tests << ", failures: " << failures << "\n";
    return failures ? 1 : 0;
}
