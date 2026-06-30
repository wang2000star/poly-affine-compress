#include "sparse_anf.h"
#include "simplify.h"
#include <iostream>
#include <iomanip>
#include <random>

static void run_test(const SparseANF& f, const std::string& label, bool verbose) {
    auto result = simplify(f, verbose);
    bool ok = verify_substitution(f, result.g, result.M, result.b, 100);
    std::cout << std::left << std::setw(30) << label
              << "  T: " << std::setw(4) << f.T() << " -> " << std::setw(4) << result.g.T()
              << "  m: " << f.n << " -> " << result.g.n
              << "  verify: " << (ok ? "OK" : "FAIL")
              << "\n";
}

int main() {
    // Test 1: x0*x1 + x0*x2  (expect T: 2 -> 1, m: 3 -> 2)
    {
        std::unordered_map<uint64_t, uint8_t> terms = {{3, 1}, {5, 1}};
        SparseANF f(terms, 3);
        run_test(f, "x0*x1 + x0*x2", false);
    }

    // Test 2: x0*x1 + x0*x2 + x1*x2
    {
        std::unordered_map<uint64_t, uint8_t> terms = {{3, 1}, {5, 1}, {6, 1}};
        SparseANF f(terms, 3);
        run_test(f, "x0*x1 + x0*x2 + x1*x2", false);

        // Debug: check intermediate results
        auto r1 = simplify_by_complement(f, false);
        std::cout << "  After complement: T=" << r1.g.T() << " m=" << r1.g.n << "\n";
        auto r2 = greedy_merge_simplify(r1.g, 100, false);
        // Compose M
        int m2 = (int)r2.M.size();
        std::vector<uint64_t> M_acc = r1.M;
        std::vector<uint64_t> M_new;
        for (int row = 0; row < m2; ++row) {
            uint64_t combined = 0;
            uint64_t r2_row = r2.M[row];
            uint64_t t = r2_row;
            while (t) {
                int i = __builtin_ctzll(t);
                t &= t - 1;
                combined ^= M_acc[i];
            }
            M_new.push_back(combined);
        }
        uint64_t b_new = 0;
        for (int j = 0; j < m2; ++j) {
            if (__builtin_parityll(r2.M[j] & r1.b)) {
                b_new |= (uint64_t)1 << j;
            }
        }
        b_new ^= r2.b;
        std::cout << "  After merge: T=" << r2.g.T() << " m=" << r2.g.n << "\n";

        // Now call search_random with verbose
        auto r3 = search_random(r2.g, 3, 200, 1, true);
        if (r3.g.T() < r2.g.T()) {
            std::cout << "  Random search found improvement: T=" << r3.g.T() << " m=" << r3.M.size() << "\n";
            std::cout << "    M[0]=" << r3.M[0] << "\n";
            if (r3.M.size() > 1) std::cout << "    M[1]=" << r3.M[1] << "\n";
            if (r3.M.size() > 2) std::cout << "    M[2]=" << r3.M[2] << "\n";
            bool ok = verify_substitution(r2.g, r3.g, r3.M, r3.b, 100);
            std::cout << "    verify: " << (ok ? "OK" : "FAIL") << "\n";
        } else {
            std::cout << "  No improvement from random search\n";
        }
    }

    // Test 3: Random cubic n=8
    {
        auto f = SparseANF::random_cubic(8, 0.15, 42);
        run_test(f, "Random cubic n=8", true);
    }

    // Test 4: Random cubic n=12
    {
        auto f = SparseANF::random_cubic(12, 0.12, 123);
        run_test(f, "Random cubic n=12", false);
    }

    // Test 5: Hidden linear structure n=10, inner=4
    {
        std::mt19937_64 rng(999);
        int inner_n = 4, outer_n = 10;
        auto inner = SparseANF::random_cubic(inner_n, 0.3, 42);
        std::vector<uint64_t> H(inner_n, 0);
        for (int i = 0; i < inner_n; ++i) {
            H[i] = rng() & ((uint64_t(1) << outer_n) - 1);
        }
        auto f = inner.expand_with(H, 0, outer_n);
        run_test(f, "Hidden structure n=10", false);
    }

    // Test 6: Hidden linear structure n=16, inner=5
    {
        std::mt19937_64 rng(7777);
        int inner_n = 5, outer_n = 16;
        auto inner = SparseANF::random_cubic(inner_n, 0.3, 42);
        std::vector<uint64_t> H(inner_n, 0);
        for (int i = 0; i < inner_n; ++i) {
            H[i] = rng() & ((uint64_t(1) << outer_n) - 1);
        }
        auto f = inner.expand_with(H, 0, outer_n);
        run_test(f, "Hidden structure n=16", false);
    }

    // Test 7: LinStr n=24, inner=5
    {
        std::mt19937_64 rng(12345);
        int inner_n = 5, outer_n = 24;
        // G = random cubic with inner_n variables
        auto inner = SparseANF::random_cubic(inner_n, 0.3, 42);
        // H = random linear injective map: inner_n × outer_n
        std::vector<uint64_t> H(inner_n, 0);
        for (int i = 0; i < inner_n; ++i) {
            H[i] = rng() & ((uint64_t(1) << outer_n) - 1);
        }
        auto f = inner.expand_with(H, 0, outer_n);
        run_test(f, "LinStr n=24 inner=5", false);
    }

    std::cout << "\nAll tests done.\n";
    return 0;
}
