#include "simplify.h"
#include "gf2_linalg.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <random>

// Identity matrix of size n
static std::vector<uint64_t> identity_matrix(int n) {
    std::vector<uint64_t> I(n, 0);
    for (int i = 0; i < n; ++i) I[i] = (uint64_t)1 << i;
    return I;
}

// ====================================================================
//  Internal helpers
// ====================================================================

// Drop unused variables from f, compacting M and b
static SimplifyResult drop_unused(const SparseANF& f,
                                   const std::vector<uint64_t>& M,
                                   uint64_t b) {
    auto used = f.variables_used();
    std::sort(used.begin(), used.end());

    if ((int)used.size() == f.n) {
        return {f, M, b};
    }

    // Compact polynomial
    std::unordered_map<uint64_t, uint8_t> new_terms;
    for (auto& [mask, v] : f.terms) {
        uint64_t new_mask = 0;
        for (int new_idx = 0; new_idx < (int)used.size(); ++new_idx) {
            if ((mask >> used[new_idx]) & 1) {
                new_mask |= (uint64_t)1 << new_idx;
            }
        }
        new_terms[new_mask] ^= v;
    }
    std::unordered_map<uint64_t, uint8_t> clean;
    for (auto& [m, v] : new_terms) {
        if (v & 1) clean[m] = 1;
    }

    // Compact M and b
    std::vector<uint64_t> new_M;
    uint64_t new_b = 0;
    for (int new_idx = 0; new_idx < (int)used.size(); ++new_idx) {
        int old_idx = used[new_idx];
        new_M.push_back(M[old_idx]);
        if ((b >> old_idx) & 1) new_b |= (uint64_t)1 << new_idx;
    }

    return {SparseANF(clean, (int)used.size()), new_M, new_b};
}

// Score of merge x_i → x_i⊕x_j: ΔT = T(after) - T(before)
static int score_merge(const std::unordered_map<uint64_t, uint8_t>& terms,
                       int i, int j) {
    std::unordered_map<uint64_t, uint8_t> toggles;
    for (auto& [m, v] : terms) {
        if ((m >> i) & 1) {
            toggles[m] ^= v;
            uint64_t rest = m ^ (uint64_t(1) << i);
            uint64_t new_m = ((m >> j) & 1) ? rest : (rest | (uint64_t(1) << j));
            toggles[new_m] ^= v;
        }
    }

    int n_changed = 0;
    int n_toggles = 0;
    int n_overlap = 0;

    for (auto& [m, v] : terms) {
        if ((m >> i) & 1) ++n_changed;
    }
    for (auto& [m, c] : toggles) {
        if (c & 1) {
            ++n_toggles;
            if (terms.count(m) && !((m >> i) & 1)) {
                ++n_overlap;
            }
        }
    }
    return n_toggles - 2 * n_overlap - n_changed;
}

// ====================================================================
//  Greedy XOR merge
// ====================================================================

SimplifyResult greedy_merge_simplify(const SparseANF& f,
                                     int max_iter, bool verbose) {
    std::unordered_map<uint64_t, uint8_t> terms = f.terms;
    int m = f.n;
    std::vector<uint64_t> M(m, 0);
    for (int i = 0; i < m; ++i) M[i] = (uint64_t)1 << i;
    uint64_t b = 0;
    int orig_T = f.T();

    if (verbose) {
        std::cout << "\nGreedy merge: n=" << m << ", T₀=" << orig_T << "\n";
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        if (terms.size() <= 1) break;

        // Find active variables
        uint64_t active_mask = 0;
        for (auto& [mask, v] : terms) active_mask |= mask;
        std::vector<int> active_list;
        uint64_t t = active_mask;
        while (t) {
            int i = __builtin_ctzll(t);
            t &= t - 1;
            active_list.push_back(i);
        }
        if (active_list.size() <= 1) break;

        // Score all pairs
        int best_delta = 0;
        int best_i = -1, best_j = -1;

        for (int i : active_list) {
            // Precompute terms containing i
            std::vector<uint64_t> ti;
            for (auto& [mask, v] : terms) {
                if ((mask >> i) & 1) ti.push_back(mask);
            }
            int n_changed_i = (int)ti.size();

            for (int j : active_list) {
                if (i == j) continue;

                // Toggle computation
                std::unordered_map<uint64_t, uint8_t> toggle_counts;
                for (auto m : ti) {
                    toggle_counts[m] ^= 1;
                    uint64_t rest = m ^ (uint64_t(1) << i);
                    uint64_t new_m = ((m >> j) & 1) ? rest : (rest | (uint64_t(1) << j));
                    toggle_counts[new_m] ^= 1;
                }

                int n_toggles = 0;
                for (auto& [m, c] : toggle_counts) {
                    if (c & 1) ++n_toggles;
                }
                int n_overlap = 0;
                for (auto& [m, c] : toggle_counts) {
                    if ((c & 1) && terms.count(m) && !((m >> i) & 1)) {
                        ++n_overlap;
                    }
                }

                int delta = n_toggles - 2 * n_overlap - n_changed_i;
                if (delta < best_delta) {
                    best_delta = delta;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        if (best_delta >= 0) break;

        // Apply merge
        M[best_i] ^= M[best_j];
        terms = apply_xor_merge(terms, best_i, best_j);

        // Drop unused variables
        SparseANF cur(terms, m);
        auto result = drop_unused(cur, M, b);
        terms = result.g.terms;
        m = result.g.n;
        M = result.M;
        b = result.b;

        if (verbose) {
            double pct = (orig_T - (int)terms.size()) * 100.0 / orig_T;
            std::cout << "  iter " << iter << ": T=" << terms.size()
                      << "/" << orig_T << " (" << pct << "%↓) m=" << m << "\n";
        }
    }

    return {SparseANF(terms, m), M, b};
}

// ====================================================================
//  Complement search
// ====================================================================

SimplifyResult simplify_by_complement(const SparseANF& f, bool verbose) {
    int n = f.n;
    std::vector<uint64_t> M(n, 0);
    for (int i = 0; i < n; ++i) M[i] = (uint64_t)1 << i;
    uint64_t b = 0;

    if (n <= 0) return {f, M, b};

    if (n > 14) {
        // Greedy per-variable
        if (verbose) {
            std::cout << "Greedy complement (n=" << n << "):\n";
        }
        SparseANF cur = f;
        uint64_t cur_b = 0;
        for (int var = 0; var < n; ++var) {
            uint64_t test_b = (uint64_t)1 << var;
            auto g = cur.substitute_affine(M, test_b);
            if (g.T() < cur.T()) {
                cur = g;
                cur_b ^= (uint64_t)1 << var;
            }
        }
        auto result = drop_unused(cur, M, cur_b);
        return {result.g, result.M, result.b};
    }

    // Exhaustive Gray code O(2ⁿ·T(f))
    if (verbose) {
        std::cout << "Complement search (Gray, 2^" << n << "=" << (n >= 64 ? 0 : (uint64_t)1 << n) << "):\n";
    }

    int best_T = f.T();
    uint64_t best_b = 0;
    auto cur_terms = f.terms;
    uint64_t cur_b = 0;

    // Gray code iteration — in-place map updates
    for (int step = 1; step < (1 << n); ++step) {
        int gray = step ^ (step >> 1);
        int prev_gray = (step - 1) ^ ((step - 1) >> 1);
        int diff = gray ^ prev_gray;
        int bit = __builtin_ctz(diff);

        // Collect monomials containing x_bit
        std::vector<uint64_t> affected;
        for (auto& [m, v] : cur_terms) {
            if (v & 1 && ((m >> bit) & 1)) {
                affected.push_back(m);
            }
        }

        // Complement: each monomial m → m ⊕ (m\{x_i}). Only rest toggles.
        for (uint64_t m : affected) {
            uint64_t rest = m ^ (uint64_t(1) << bit);
            cur_terms[rest] ^= 1;
        }

        // Remove zero-coefficient entries
        for (auto it = cur_terms.begin(); it != cur_terms.end(); ) {
            if (it->second & 1) {
                it->second = 1;
                ++it;
            } else {
                it = cur_terms.erase(it);
            }
        }

        cur_b ^= (uint64_t(1) << bit);

        int T = (int)cur_terms.size();
        if (T < best_T) {
            best_T = T;
            best_b = cur_b;
        }
    }

    // Recompute best g from best_b
    SparseANF best_g;
    if (best_T < f.T()) {
        best_g = f.substitute_affine(M, best_b);
    } else {
        best_g = f;
        best_b = 0;
    }

    auto result = drop_unused(best_g, M, best_b);
    if (verbose) {
        double pct = (f.T() - best_T) * 100.0 / f.T();
        std::cout << "  best: T=" << best_T << "/" << f.T()
                  << " (" << pct << "%↓)\n";
    }
    return {result.g, result.M, result.b};
}

// ====================================================================
//  Random M,b search
// ====================================================================

SimplifyResult search_random(const SparseANF& f, int max_m, int n_trials,
                              uint64_t seed, bool verbose) {
    int n = f.n;
    int best_T = f.T();
    SimplifyResult best = {f, identity_matrix(n), 0};

    std::mt19937_64 rng(seed);
    for (int trial = 0; trial < n_trials; ++trial) {
        int m = (rng() % max_m) + 1;
        std::vector<uint64_t> M(m, 0);
        for (int j = 0; j < m; ++j) {
            M[j] = rng() & low_mask(n);
        }
        uint64_t b = rng() & low_mask(m);

        // Skip rank-deficient M (must be full row rank for m≤n, full column rank for m>n)
        if (gf2_rank(M, n) < std::min(m, n)) continue;

        auto g = f.substitute_affine(M, b);

        // Decide whether to use Z = Z1 ∪ Z2
        // Start with the plain version as default
        auto g_best = g;
        std::vector<uint64_t> M_best = M;
        uint64_t b_best = b;
        bool use_union = false;

        // Try union only if it improves T
        auto g_u = f.substitute_affine_union(M, b);
        if (g_u.T() < g.T()) {
            M_best = M;
            for (int i = 0; i < f.n; ++i)
                M_best.push_back((uint64_t)1 << i);
            b_best = b;
            g_best = g_u;
            use_union = true;
        }

        if (g_best.T() < best_T) {
            // Verify when dimension changes
            if (g_best.n != f.n && !verify_substitution(f, g_best, M_best, b_best, 50)) {
                continue;
            }
            best_T = g_best.T();
            best = {g_best, M_best, b_best};
            if (verbose) {
                double pct = (f.T() - best_T) * 100.0 / f.T();
                std::cout << "  trial " << trial << ": T=" << best_T
                          << "/" << f.T() << " (" << pct << "%↓)"
                          << " m=" << (int)M_best.size()
                          << (use_union ? " [union]" : "") << "\n";
            }
        }
    }
    return best;
}

// ====================================================================
//  Combined pipeline
// ====================================================================

SimplifyResult simplify(const SparseANF& f, bool verbose) {
    int n = f.n;
    int orig_T = f.T();

    if (verbose) {
        std::cout << "\n=== Simplify: n=" << n << ", T₀=" << orig_T << " ===\n";
    }

    // Phase 1: complement
    auto r1 = simplify_by_complement(f, verbose);
    auto g = r1.g;
    auto M_acc = r1.M;
    auto b_acc = r1.b;

    // Phase 2: greedy merge
    auto r2 = greedy_merge_simplify(g, 100, verbose);
    // Compose: M = r2.M @ M_acc (with uint64_t row ops)
    int m2 = (int)r2.M.size();
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
    // b_new = r2.M · b_acc ⊕ r2.b
    uint64_t b_new = 0;
    for (int j = 0; j < m2; ++j) {
        if (__builtin_parityll(r2.M[j] & b_acc)) {
            b_new |= (uint64_t)1 << j;
        }
    }
    b_new ^= r2.b;

    g = r2.g;
    M_acc = M_new;
    b_acc = b_new;

    if (verbose) {
        double pct = (orig_T - g.T()) * 100.0 / orig_T;
        std::cout << "  After merge: T=" << g.T() << "/" << orig_T
                  << " (" << pct << "%↓), m=" << g.n << "\n";
    }

    // Phase 3: random M,b search
    if (g.n > 0 && g.T() > 1) {
        int max_random_m = std::max(g.n, 1);
        auto r3 = search_random(g, max_random_m, 200, 1, false);
        if (r3.g.T() < g.T()) {
            // Compose: M = r3.M @ M_acc
            int m3 = (int)r3.M.size();
            std::vector<uint64_t> M3_new;
            for (int row = 0; row < m3; ++row) {
                uint64_t combined = 0;
                uint64_t r3_row = r3.M[row];
                uint64_t t = r3_row;
                while (t) {
                    int i = __builtin_ctzll(t);
                    t &= t - 1;
                    combined ^= M_acc[i];
                }
                M3_new.push_back(combined);
            }
            uint64_t b3_new = 0;
            for (int j = 0; j < m3; ++j) {
                if (__builtin_parityll(r3.M[j] & b_acc)) {
                    b3_new |= (uint64_t)1 << j;
                }
            }
            b3_new ^= r3.b;

            g = r3.g;
            M_acc = M3_new;
            b_acc = b3_new;

            if (verbose) {
                double pct = (orig_T - g.T()) * 100.0 / orig_T;
                std::cout << "  After random: T=" << g.T() << "/" << orig_T
                          << " (" << pct << "%↓), m=" << g.n << "\n";
            }
        }
    }

    if (verbose) {
        double pct = (orig_T - g.T()) * 100.0 / orig_T;
        std::cout << "  Final: T=" << g.T() << "/" << orig_T
                  << " (" << pct << "%↓), m=" << g.n << "\n";
    }

    return {g, M_acc, b_acc};
}
