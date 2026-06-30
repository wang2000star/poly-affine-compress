#include "sparse_anf.h"
#include "simplify.h"
#include <iostream>
#include <chrono>
#include <random>

struct BenchResult {
    std::string label;
    int T_before;
    int T_after;
    int n_before;
    int n_after;
    double time_ms;
    bool verified;
};

static BenchResult bench_one(const std::string& label, const SparseANF& f) {
    auto t0 = std::chrono::steady_clock::now();
    auto result = simplify(f, false);
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    bool ok = verify_substitution(f, result.g, result.M, result.b, 50);
    return {label, f.T(), result.g.T(), f.n, result.g.n, ms, ok};
}

int main() {
    std::cout << "=== ANF Simplify Benchmarks ===\n\n";

    // 1) Random cubics
    for (int n : {8, 12, 16, 20}) {
        auto f = SparseANF::random_cubic(n, 0.15, 42);
        auto r = bench_one("Random cubic n=" + std::to_string(n), f);
        printf("  %-30s T: %4d -> %-4d  m: %d -> %-2d  %8.1f ms  %s\n",
               r.label.c_str(), r.T_before, r.T_after,
               r.n_before, r.n_after, r.time_ms,
               r.verified ? "OK" : "FAIL");
    }

    // 2) Hidden linear structure
    for (int outer_n : {10, 16, 24}) {
        int inner_n = 5;
        std::mt19937_64 rng(outer_n * 1000);
        auto inner = SparseANF::random_cubic(inner_n, 0.3, 42);
        std::vector<uint64_t> H(inner_n, 0);
        for (int i = 0; i < inner_n; ++i) {
            H[i] = rng() & ((uint64_t(1) << outer_n) - 1);
        }
        auto f = inner.expand_with(H, 0, outer_n);
        auto r = bench_one("LinStr n=" + std::to_string(outer_n), f);
        printf("  %-30s T: %4d -> %-4d  m: %d -> %-2d  %8.1f ms  %s\n",
               r.label.c_str(), r.T_before, r.T_after,
               r.n_before, r.n_after, r.time_ms,
               r.verified ? "OK" : "FAIL");
    }

    // 3) Large random cubic (greedy merge only)
    for (int n : {32, 48, 64}) {
        auto f = SparseANF::random_cubic(n, 0.02, 123);
        auto t0 = std::chrono::steady_clock::now();
        auto result = greedy_merge_simplify(f, 100, false);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        bool ok = verify_substitution(f, result.g, result.M, result.b, 50);
        printf("  %-30s T: %4d -> %-4d  m: %d -> %-2d  %8.1f ms  %s\n",
               ("Greedy n=" + std::to_string(n)).c_str(),
               f.T(), result.g.T(), f.n, result.g.n, ms,
               ok ? "OK" : "FAIL");
    }

    std::cout << "\nDone.\n";
    return 0;
}
