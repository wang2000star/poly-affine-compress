#include "vector_int_poly.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>

// ====================================================================
//  VectorIntPoly implementation
// ====================================================================

VectorIntPoly::VectorIntPoly(const std::vector<IntPoly>& comps) {
    assert(!comps.empty());
    n = comps[0].n;
    k = (int)comps.size();
    for (auto& f : comps) {
        assert(f.n == n);
        components.push_back(f);
    }
}

int VectorIntPoly::union_T() const {
    std::unordered_map<ExpVector, int64_t, ExpHash> seen;
    for (auto& f : components) {
        for (auto& [exp, c] : f.terms) {
            seen[exp] = 1;
        }
    }
    return (int)seen.size();
}

VectorIntPoly VectorIntPoly::substitute_affine(
    const std::vector<std::vector<int64_t>>& M,
    const std::vector<int64_t>& b) const {
    std::vector<IntPoly> new_comps;
    for (auto& f : components) {
        new_comps.push_back(f.substitute_affine(M, b));
    }
    return VectorIntPoly(new_comps);
}

VectorIntPoly VectorIntPoly::substitute_affine_union(
    const std::vector<std::vector<int64_t>>& M,
    const std::vector<int64_t>& b) const {
    std::vector<IntPoly> new_comps;
    for (auto& f : components) {
        new_comps.push_back(f.substitute_affine_union(M, b));
    }
    return VectorIntPoly(new_comps);
}

std::vector<int> VectorIntPoly::variables_used() const {
    std::vector<bool> used_flag(n, false);
    for (auto& f : components) {
        for (auto& [exp, c] : f.terms) {
            for (int i = 0; i < n && i < (int)exp.size(); ++i) {
                if (exp[i] > 0) used_flag[i] = true;
            }
        }
    }
    std::vector<int> used;
    for (int i = 0; i < n; ++i) {
        if (used_flag[i]) used.push_back(i);
    }
    return used;
}

// ====================================================================
//  Simplification
// ====================================================================

static std::vector<std::vector<int64_t>> identity_int_mat(int n) {
    std::vector<std::vector<int64_t>> I(n, std::vector<int64_t>(n, 0));
    for (int i = 0; i < n; ++i) I[i][i] = 1;
    return I;
}

static VectorIntSimplifyResult vector_drop_unused_int(
    const VectorIntPoly& vec,
    const std::vector<std::vector<int64_t>>& M,
    const std::vector<int64_t>& b) {
    auto used = vec.variables_used();
    if ((int)used.size() == vec.n) return {vec, M, b};

    std::sort(used.begin(), used.end());

    std::vector<IntPoly> new_comps;
    for (auto& f : vec.components) {
        std::unordered_map<ExpVector, int64_t, ExpHash> new_terms;
        for (auto& [exp, c] : f.terms) {
            ExpVector new_exp;
            for (int idx : used)
                new_exp.push_back(idx < (int)exp.size() ? exp[idx] : 0);
            new_terms[std::move(new_exp)] = c;
        }
        new_comps.push_back(IntPoly(std::move(new_terms), (int)used.size()));
    }

    std::vector<std::vector<int64_t>> new_M;
    std::vector<int64_t> new_b;
    for (int idx : used) {
        if (idx < (int)M.size()) new_M.push_back(M[idx]);
        new_b.push_back(idx < (int)b.size() ? b[idx] : 0);
    }

    return {VectorIntPoly(new_comps), new_M, new_b};
}

static int vector_try_merge_int_T(const VectorIntPoly& vec, int i, int j, int64_t k) {
    std::vector<int64_t> coeffs(vec.n, 0);
    coeffs[i] = 1;
    coeffs[j] = -k;

    int total = 0;
    for (auto& f : vec.components) {
        auto g = f.substitute_linear(i, coeffs);
        total += g.T();
    }
    return total;
}

VectorIntSimplifyResult vector_greedy_merge_int(
    const VectorIntPoly& vec, int max_iter, bool verbose) {
    if (vec.components.empty() || vec.union_T() == 0) {
        return {vec, identity_int_mat(vec.n), std::vector<int64_t>(vec.n, 0)};
    }

    VectorIntPoly cur = vec;
    auto M = identity_int_mat(vec.n);
    std::vector<int64_t> b(vec.n, 0);
    int orig_T = cur.union_T();

    if (verbose) {
        std::cout << "  Vector greedy merge (int): n=" << vec.n
                  << ", k=" << vec.k << ", T₀=" << orig_T << "\n";
    }

    std::vector<int64_t> k_values = {1, -1, 2, -2};

    for (int iter = 0; iter < max_iter; ++iter) {
        auto active = cur.variables_used();
        if ((int)active.size() <= 1) break;

        int best_T = cur.union_T();
        int best_i = -1, best_j = -1;
        int64_t best_k = 0;

        // Use sum of T_i as proxy (faster than union_T for trial)
        int best_sum = 0;

        for (int i : active) {
            for (int j : active) {
                if (i == j) continue;
                for (int64_t k : k_values) {
                    int sum_T = vector_try_merge_int_T(cur, i, j, k);
                    if (best_i < 0 || sum_T < best_sum) {
                        best_sum = sum_T;
                        best_i = i; best_j = j; best_k = k;
                    }
                }
            }
        }

        if (best_i < 0) break;

        // Apply merge to all components
        std::vector<int64_t> coeffs(cur.n, 0);
        coeffs[best_i] = 1;
        coeffs[best_j] = -best_k;

        std::vector<IntPoly> new_comps;
        for (auto& f : cur.components) {
            new_comps.push_back(f.substitute_linear(best_i, coeffs));
        }
        cur = VectorIntPoly(new_comps);
        for (int col = 0; col < (int)M[best_i].size(); ++col)
            M[best_i][col] += best_k * M[best_j][col];
        b[best_i] += best_k * b[best_j];

        auto result = vector_drop_unused_int(cur, M, b);
        cur = result.g; M = result.M; b = result.b;

        if (verbose) {
            double pct = (orig_T - cur.union_T()) * 100.0 / orig_T;
            std::cout << "    iter " << iter << ": x" << best_i << "→x" << best_i
                      << "+" << best_k << "x" << best_j
                      << "  T=" << cur.union_T() << "/" << orig_T
                      << " (" << pct << "%↓)  m=" << cur.n << "\n";
        }

        if (cur.union_T() <= 1) break;
    }

    return {cur, M, b};
}

// ====================================================================
//  Random search with union
// ====================================================================

VectorIntSimplifyResult vector_search_random_int(const VectorIntPoly& vec,
                                                   int max_m, int n_trials,
                                                   uint64_t seed, bool verbose) {
    if (vec.components.empty() || vec.union_T() <= 1)
        return {vec, identity_int_mat(vec.n), std::vector<int64_t>(vec.n, 0)};

    int best_T = vec.union_T();
    VectorIntPoly best_vec = vec;
    auto best_M = identity_int_mat(vec.n);
    std::vector<int64_t> best_b(vec.n, 0);

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> entry_dist(-1, 1);
    std::uniform_int_distribution<int> m_dist(1, max_m);

    for (int trial = 0; trial < n_trials; ++trial) {
        int m = m_dist(rng);
        std::vector<std::vector<int64_t>> M(m, std::vector<int64_t>(vec.n, 0));
        std::vector<int64_t> b(m, 0);

        // Structured: each row at most one non-zero entry ±1
        bool ok = true;
        std::vector<bool> used_var(vec.n, false);
        for (int row = 0; row < m; ++row) {
            int var = rng() % vec.n;
            if (used_var[var]) { ok = false; break; }
            used_var[var] = true;
            M[row][var] = (rng() % 2) ? 1LL : -1LL;
            if (rng() % 3 == 0) b[row] = entry_dist(rng);
        }
        if (!ok) continue;

        // Try each component with union
        std::vector<IntPoly> new_comps;
        int total_m = m + vec.n;
        for (auto& f : vec.components) {
            auto g = f.substitute_affine_union(M, b);
            new_comps.push_back(g);
        }
        VectorIntPoly candidate(new_comps);
        int T = candidate.union_T();

        if (T < best_T) {
            // Build M_ext = [M; I_n]
            std::vector<std::vector<int64_t>> M_ext = M;
            for (int i = 0; i < vec.n; ++i) {
                std::vector<int64_t> row(vec.n, 0);
                row[i] = 1;
                M_ext.push_back(row);
            }
            std::vector<int64_t> b_ext = b;
            b_ext.resize(total_m, 0);

            best_T = T;
            best_vec = candidate;
            best_M = M_ext;
            best_b = b_ext;

            if (verbose) {
                double pct = (vec.union_T() - T) * 100.0 / vec.union_T();
                std::cout << "  trial " << trial << ": T=" << T
                          << "/" << vec.union_T() << " (" << pct << "%↓)"
                          << " m=" << total_m << " [union]\n";
            }
        }
    }
    return {best_vec, best_M, best_b};
}

// ====================================================================
//  Combined pipeline
// ====================================================================

VectorIntSimplifyResult vector_simplify_int(const VectorIntPoly& vec, bool verbose) {
    if (vec.components.empty()) {
        return {vec, identity_int_mat(vec.n), std::vector<int64_t>(vec.n, 0)};
    }

    int n = vec.n;
    auto M_acc = identity_int_mat(n);
    std::vector<int64_t> b_acc(n, 0);
    VectorIntPoly cur = vec;

    if (verbose) {
        std::cout << "Vector simplify (int): n=" << n << ", k=" << vec.k
                  << ", T_union=" << vec.union_T() << "\n";
    }

    // Phase 1: Greedy merge
    auto r1 = vector_greedy_merge_int(cur, 50, verbose);
    if (r1.g.union_T() < cur.union_T()) {
        int m_new = (int)r1.M.size();
        int m_cur = (int)M_acc.size();

        std::vector<std::vector<int64_t>> M_new(m_new, std::vector<int64_t>(n, 0));
        for (int j = 0; j < m_new; ++j) {
            for (int k = 0; k < m_cur; ++k) {
                if (r1.M[j][k] == 0) continue;
                for (int i = 0; i < n; ++i)
                    M_new[j][i] += r1.M[j][k] * M_acc[k][i];
            }
        }

        std::vector<int64_t> b_new(m_new, 0);
        for (int j = 0; j < m_new; ++j) {
            for (int k = 0; k < m_cur; ++k)
                b_new[j] += r1.M[j][k] * b_acc[k];
            b_new[j] += r1.b[j];
        }

        cur = r1.g;
        M_acc = M_new;
        b_acc = b_new;

        if (verbose)
            std::cout << "  After merge: T=" << cur.union_T()
                      << ", m=" << cur.n << "\n";
    }

    // Phase 2: Random search with union
    if (cur.union_T() > 1 && cur.n > 0) {
        int max_random_m = std::max(std::min(cur.n, 6), 1);
        auto r2 = vector_search_random_int(cur, max_random_m, 50, 2, verbose);
        if (r2.g.union_T() < cur.union_T()) {
            int m_new = (int)r2.M.size();
            int m_cur = (int)M_acc.size();

            std::vector<std::vector<int64_t>> M_new(m_new, std::vector<int64_t>(n, 0));
            for (int j = 0; j < m_new; ++j) {
                for (int k = 0; k < m_cur; ++k) {
                    if (r2.M[j][k] == 0) continue;
                    for (int i = 0; i < n; ++i)
                        M_new[j][i] += r2.M[j][k] * M_acc[k][i];
                }
            }

            std::vector<int64_t> b_new(m_new, 0);
            for (int j = 0; j < m_new; ++j) {
                for (int k = 0; k < m_cur; ++k)
                    b_new[j] += r2.M[j][k] * b_acc[k];
                b_new[j] += r2.b[j];
            }

            cur = r2.g;
            M_acc = M_new;
            b_acc = b_new;

            if (verbose)
                std::cout << "  After random: T=" << cur.union_T()
                          << ", m=" << cur.n << "\n";
        }
    }

    return {cur, M_acc, b_acc};
}
