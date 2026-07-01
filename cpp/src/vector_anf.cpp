#include "vector_anf.h"
#include "simplify.h"
#include "gf2_linalg.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>

// ====================================================================
//  VectorANF implementation
// ====================================================================

VectorANF::VectorANF(const std::vector<SparseANF>& comps) {
    assert(!comps.empty());
    n = comps[0].n;
    k = (int)comps.size();
    for (auto& f : comps) {
        assert(f.n == n);
        components.push_back(f);
    }
}

int VectorANF::union_T() const {
    uint64_t seen_masks = 0;
    // Use a small set; since n ≤ 64 we track via bit operations
    // But monomials can overlap — need proper set
    std::unordered_map<uint64_t, uint8_t> all;
    for (auto& f : components) {
        for (auto& [m, v] : f.terms) {
            all[m] = 1;
        }
    }
    return (int)all.size();
}

VectorANF VectorANF::substitute_affine(const std::vector<uint64_t>& M, uint64_t b) const {
    std::vector<SparseANF> new_comps;
    for (auto& f : components) {
        new_comps.push_back(f.substitute_affine(M, b));
    }
    return VectorANF(new_comps);
}

VectorANF VectorANF::substitute_affine_union(const std::vector<uint64_t>& M,
                                              uint64_t b) const {
    std::vector<SparseANF> new_comps;
    for (auto& f : components) {
        new_comps.push_back(f.substitute_affine_union(M, b));
    }
    return VectorANF(new_comps);
}

VectorANF VectorANF::substitute_linear(int var, const std::vector<uint8_t>& coeffs) const {
    std::vector<SparseANF> new_comps;
    for (auto& f : components) {
        std::unordered_map<uint64_t, uint8_t> new_terms;
        for (auto& [mask, v] : f.terms) {
            if (!((mask >> var) & 1)) {
                new_terms[mask] ^= v;
                continue;
            }
            uint64_t rest = mask ^ ((uint64_t)1 << var);
            for (int j = 0; j < (int)coeffs.size() && j < n; ++j) {
                if (coeffs[j]) {
                    if (j == var) {
                        new_terms[mask] ^= v;
                    } else {
                        uint64_t new_m = (mask >> j) & 1 ? rest : (rest | ((uint64_t)1 << j));
                        new_terms[new_m] ^= v;
                    }
                }
            }
        }
        std::unordered_map<uint64_t, uint8_t> clean;
        for (auto& [m, val] : new_terms) {
            if (val & 1) clean[m] = 1;
        }
        new_comps.push_back(SparseANF(clean, n));
    }
    return VectorANF(new_comps);
}

std::vector<std::vector<SparseANF>> VectorANF::gradient() const {
    std::vector<std::vector<SparseANF>> grads;
    for (auto& f : components) {
        grads.push_back(f.gradient());
    }
    return grads;
}

std::vector<int> VectorANF::variables_used() const {
    uint64_t seen = 0;
    for (auto& f : components) {
        for (auto& [mask, v] : f.terms) {
            seen |= mask;
        }
    }
    std::vector<int> used;
    uint64_t t = seen;
    while (t) {
        int i = __builtin_ctzll(t);
        t &= t - 1;
        used.push_back(i);
    }
    return used;
}

// ====================================================================
//  Simplification helpers
// ====================================================================

static std::vector<uint64_t> identity_vec(int n) {
    std::vector<uint64_t> I(n, 0);
    for (int i = 0; i < n; ++i) I[i] = (uint64_t)1 << i;
    return I;
}

static VectorSimplifyResult vector_drop_unused(const VectorANF& vec,
                                                const std::vector<uint64_t>& M,
                                                uint64_t b) {
    auto used = vec.variables_used();
    std::sort(used.begin(), used.end());
    if ((int)used.size() == vec.n) return {vec, M, b};

    // Compact each component
    std::vector<SparseANF> new_comps;
    for (auto& f : vec.components) {
        std::unordered_map<uint64_t, uint8_t> new_terms;
        for (auto& [mask, v] : f.terms) {
            uint64_t new_mask = 0;
            for (int new_idx = 0; new_idx < (int)used.size(); ++new_idx) {
                if ((mask >> used[new_idx]) & 1)
                    new_mask |= (uint64_t)1 << new_idx;
            }
            new_terms[new_mask] ^= v;
        }
        std::unordered_map<uint64_t, uint8_t> clean;
        for (auto& [m, val] : new_terms) {
            if (val & 1) clean[m] = 1;
        }
        new_comps.push_back(SparseANF(clean, (int)used.size()));
    }

    // Compact M and b
    std::vector<uint64_t> new_M;
    uint64_t new_b = 0;
    for (int new_idx = 0; new_idx < (int)used.size(); ++new_idx) {
        int old_idx = used[new_idx];
        new_M.push_back(M[old_idx]);
        if ((b >> old_idx) & 1) new_b |= (uint64_t)1 << new_idx;
    }

    return {VectorANF(new_comps), new_M, new_b};
}

// ====================================================================
//  Greedy XOR merge
// ====================================================================

static int vector_try_merge_T(const VectorANF& vec, int i, int j) {
    std::vector<uint8_t> coeffs(vec.n, 0);
    coeffs[i] = 1;
    coeffs[j] = 1;
    auto g = vec.substitute_linear(i, coeffs);
    return g.union_T();
}

VectorSimplifyResult vector_greedy_merge(const VectorANF& vec,
                                          int max_iter, bool verbose) {
    if (vec.union_T() == 0) {
        return {vec, identity_vec(vec.n), 0};
    }

    VectorANF cur = vec;
    auto M = identity_vec(vec.n);
    uint64_t b = 0;
    int orig_T = cur.union_T();

    if (verbose) {
        std::cout << "  Vector greedy merge: n=" << vec.n << ", k=" << vec.k
                  << ", T₀=" << orig_T << "\n";
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        auto active = cur.variables_used();
        if ((int)active.size() <= 1) break;

        int best_T = cur.union_T();
        int best_i = -1, best_j = -1;
        VectorANF best_vec;

        for (int i : active) {
            for (int j : active) {
                if (i == j) continue;
                int T = vector_try_merge_T(cur, i, j);
                if (T < best_T) {
                    best_T = T;
                    best_i = i;
                    best_j = j;
                    std::vector<uint8_t> coeffs(cur.n, 0);
                    coeffs[i] = 1;
                    coeffs[j] = 1;
                    best_vec = cur.substitute_linear(i, coeffs);
                }
            }
        }

        if (best_i < 0) break;

        M[best_i] ^= M[best_j];
        cur = best_vec;

        auto result = vector_drop_unused(cur, M, b);
        cur = result.g;
        M = result.M;
        b = result.b;

        if (verbose) {
            double pct = (orig_T - cur.union_T()) * 100.0 / orig_T;
            std::cout << "    iter " << iter << ": x" << best_i << "→x" << best_i
                      << "⊕x" << best_j << "  T=" << cur.union_T()
                      << "/" << orig_T << " (" << pct << "%↓)  m=" << cur.n << "\n";
        }

        if (cur.union_T() <= 1) break;
    }

    return {cur, M, b};
}

// ====================================================================
//  Complement search
// ====================================================================

VectorSimplifyResult vector_simplify_by_complement(const VectorANF& vec,
                                                    bool verbose) {
    int n = vec.n;

    if (n > 16) {
        // Greedy per-variable
        if (verbose) std::cout << "  Greedy complement (n=" << n << "):\n";
        VectorANF cur = vec;
        uint64_t cur_b = 0;
        auto M = identity_vec(n);

        for (int var = 0; var < n; ++var) {
            uint64_t test_b = (uint64_t)1 << var;
            auto g = cur.substitute_affine(M, test_b);
            if (g.union_T() < cur.union_T()) {
                cur = g;
                cur_b ^= (uint64_t)1 << var;
            }
        }
        auto result = vector_drop_unused(cur, M, cur_b);
        return {result.g, result.M, result.b};
    }

    if (verbose) {
        std::cout << "  Complement search (Gray, 2^" << n << "):\n";
    }

    // Gray code over all 2^n complement patterns
    VectorANF cur = vec;
    uint64_t cur_b = 0;
    int best_T = vec.union_T();
    uint64_t best_b = 0;
    VectorANF best_vec = vec;

    for (int step = 1; step < (1 << n); ++step) {
        int gray = step ^ (step >> 1);
        int prev = (step - 1) ^ ((step - 1) >> 1);
        int diff = gray ^ prev;
        int bit = __builtin_ctz(diff);

        // Flip bit: complement x_bit in all components
        std::vector<SparseANF> new_comps;
        for (auto& f : cur.components) {
            std::unordered_map<uint64_t, uint8_t> new_terms;
            for (auto& [mask, v] : f.terms) {
                new_terms[mask] ^= v;
                if ((mask >> bit) & 1) {
                    uint64_t rest = mask ^ ((uint64_t)1 << bit);
                    new_terms[rest] ^= v;
                }
            }
            std::unordered_map<uint64_t, uint8_t> clean;
            for (auto& [m, val] : new_terms) {
                if (val & 1) clean[m] = 1;
            }
            new_comps.push_back(SparseANF(clean, f.n));
        }
        cur = VectorANF(new_comps);
        cur_b ^= (uint64_t)1 << bit;

        int T = cur.union_T();
        if (T < best_T) {
            best_T = T;
            best_vec = cur;
            best_b = cur_b;
        }
    }

    auto result = vector_drop_unused(best_vec, identity_vec(n), best_b);
    return {result.g, result.M, result.b};
}

// ====================================================================
//  Combined pipeline
// ====================================================================

VectorSimplifyResult vector_simplify(const VectorANF& vec, bool verbose) {
    int n = vec.n;
    uint64_t b_acc = 0;
    auto M_acc = identity_vec(n);
    VectorANF cur = vec;

    if (verbose) {
        std::cout << "Vector simplify: n=" << n << ", k=" << vec.k
                  << ", T_union=" << vec.union_T() << "\n";
    }

    // Phase 1: complement
    auto r1 = vector_simplify_by_complement(cur, verbose);
    if (r1.g.union_T() < cur.union_T()) {
        // Compose: M = M1 ∘ M_acc (GF2)
        int m1 = (int)r1.M.size();
        if (m1 == n) {
            // M1 is n×n identity (complement only changes b)
            M_acc = r1.M;
            b_acc ^= r1.b;
        } else {
            std::vector<uint64_t> M_new;
            for (int row = 0; row < m1; ++row) {
                uint64_t combined = 0;
                uint64_t t = r1.M[row];
                while (t) {
                    int i = __builtin_ctzll(t);
                    t &= t - 1;
                    combined ^= M_acc[i];
                }
                M_new.push_back(combined);
            }
            M_acc = M_new;
            b_acc ^= r1.b;
        }
        cur = r1.g;
        if (verbose) std::cout << "  After complement: T=" << cur.union_T() << "\n";
    }

    // Phase 2: greedy XOR merge
    auto r2 = vector_greedy_merge(cur, 50, verbose);
    if (r2.g.union_T() < cur.union_T()) {
        // Compose M
        int m2 = (int)r2.M.size();
        int m_cur = (int)M_acc.size();
        std::vector<uint64_t> M_new;
        for (int row = 0; row < m2; ++row) {
            uint64_t combined = 0;
            uint64_t t = r2.M[row];
            while (t) {
                int i = __builtin_ctzll(t);
                t &= t - 1;
                combined ^= M_acc[i];
            }
            M_new.push_back(combined);
        }
        uint64_t b_new = 0;
        for (int j = 0; j < m2; ++j) {
            if (__builtin_parityll(r2.M[j] & b_acc))
                b_new |= (uint64_t)1 << j;
        }
        b_new ^= r2.b;

        cur = r2.g;
        M_acc = M_new;
        b_acc = b_new;

        if (verbose) {
            std::cout << "  After merge: T=" << cur.union_T()
                      << ", m=" << cur.n << "\n";
        }
    }

    // Phase 3: complement union — keep both x_i AND x_i⊕1
    if (cur.union_T() > 1 && cur.n > 0) {
        uint64_t ones = ((uint64_t)1 << cur.n) - 1;
        // Build M_comp = I_{cur.n}, b_comp = all ones
        std::vector<uint64_t> M_comp(cur.n, 0);
        for (int i = 0; i < cur.n; ++i) M_comp[i] = (uint64_t)1 << i;
        auto union_vec = cur.substitute_affine_union(M_comp, ones);
        if (union_vec.union_T() < cur.union_T()) {
            // Compose: M_new = [M_comp; I_n] ∘ M_acc
            int m_comp = cur.n;  // M_comp rows
            int total = m_comp + cur.n;
            std::vector<uint64_t> M_total;
            for (int row = 0; row < total; ++row) {
                uint64_t row_mask = 0;
                uint64_t t = row < m_comp ? M_comp[row] : ((uint64_t)1 << (row - m_comp));
                while (t) {
                    int i = __builtin_ctzll(t);
                    t &= t - 1;
                    row_mask ^= M_acc[i];
                }
                M_total.push_back(row_mask);
            }
            uint64_t b_total = 0;
            for (int j = 0; j < m_comp; ++j) {
                if (__builtin_parityll(M_comp[j] & b_acc))
                    b_total |= (uint64_t)1 << j;
            }
            b_total ^= ones;
            for (int j = m_comp; j < total; ++j) {
                if (__builtin_parityll(((uint64_t)1 << (j - m_comp)) & b_acc))
                    b_total |= (uint64_t)1 << j;
            }

            cur = union_vec;
            M_acc = M_total;
            b_acc = b_total;

            if (verbose) {
                std::cout << "  After complement union: T=" << cur.union_T()
                          << ", m=" << cur.n << "\n";
            }
        }
    }

    return {cur, M_acc, b_acc};
}

// ====================================================================
//  Per-component simplification
// ====================================================================

VectorPerCompResult vector_simplify_per_component(const VectorANF& vec,
                                                   bool verbose) {
    std::vector<SimplifyResult> results;
    for (int idx = 0; idx < vec.k; ++idx) {
        if (verbose) {
            std::cout << "  Component " << idx << ": n=" << vec.components[idx].n
                      << ", T=" << vec.components[idx].T() << "\n";
        }
        auto f = vec.components[idx];
        // Complement + greedy merge
        auto r1 = simplify_by_complement(f, false);
        auto r2 = greedy_merge_simplify(r1.g, 50, false);

        // Compose
        int m2 = (int)r2.M.size();
        std::vector<uint64_t> M_acc;
        for (int row = 0; row < m2; ++row) {
            uint64_t combined = 0;
            uint64_t t = r2.M[row];
            while (t) {
                int i = __builtin_ctzll(t);
                t &= t - 1;
                combined ^= r1.M[i];
            }
            M_acc.push_back(combined);
        }
        uint64_t b_acc = 0;
        for (int j = 0; j < m2; ++j) {
            if (__builtin_parityll(r2.M[j] & r1.b))
                b_acc |= (uint64_t)1 << j;
        }
        b_acc ^= r2.b;

        results.push_back(SimplifyResult{r2.g, M_acc, b_acc});
        if (verbose) {
            std::cout << "    → T=" << r2.g.T() << ", m=" << r2.g.n << "\n";
        }
    }

    // Build shared M with distinct rows from all M_i
    // Since n ≤ 64, represent rows as uint64_t
    std::unordered_map<uint64_t, int> row_to_idx;
    std::vector<uint64_t> all_rows;
    std::vector<std::vector<int>> comp_maps;

    for (auto& r : results) {
        std::vector<int> cmap;
        for (int ri = 0; ri < (int)r.M.size(); ++ri) {
            uint64_t row = r.M[ri];
            auto it = row_to_idx.find(row);
            if (it == row_to_idx.end()) {
                int idx = (int)all_rows.size();
                row_to_idx[row] = idx;
                all_rows.push_back(row);
                cmap.push_back(idx);
            } else {
                cmap.push_back(it->second);
            }
        }
        comp_maps.push_back(cmap);
    }

    int m_total = (int)all_rows.size();

    // Joint b
    uint64_t b_joint = 0;
    for (int idx = 0; idx < (int)results.size(); ++idx) {
        auto& cmap = comp_maps[idx];
        uint64_t b_i = results[idx].b;
        for (int local = 0; local < (int)cmap.size(); ++local) {
            if ((b_i >> local) & 1)
                b_joint |= (uint64_t)1 << cmap[local];
        }
    }

    // Remap each g_i to shared variable space
    std::vector<SparseANF> joint_comps;
    for (int idx = 0; idx < (int)results.size(); ++idx) {
        auto& g = results[idx].g;
        auto& cmap = comp_maps[idx];
        std::unordered_map<uint64_t, uint8_t> remapped;
        for (auto& [mask, v] : g.terms) {
            uint64_t new_mask = 0;
            uint64_t t = mask;
            while (t) {
                int local = __builtin_ctzll(t);
                t &= t - 1;
                int global = cmap[local];
                new_mask |= (uint64_t)1 << global;
            }
            remapped[new_mask] = v;
        }
        joint_comps.push_back(SparseANF(remapped, m_total));
    }

    VectorANF joint_vec(joint_comps);
    return {results, joint_vec, all_rows, b_joint};
}
