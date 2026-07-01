#include "sparse_anf.h"
#include "gf2_linalg.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <sstream>

SparseANF::SparseANF(const std::unordered_map<uint64_t, uint8_t>& t, int n)
    : n(n) {
    for (auto& [m, v] : t) {
        if (v & 1) terms[m] = 1;
    }
}

SparseANF::SparseANF(std::unordered_map<uint64_t, uint8_t>&& t, int n)
    : n(n) {
    for (auto& [m, v] : t) {
        if (v & 1) terms[m] = 1;
    }
}

SparseANF SparseANF::from_mask_set(const std::vector<uint64_t>& masks, int n) {
    std::unordered_map<uint64_t, uint8_t> t;
    for (auto m : masks) t[m] = 1;
    return SparseANF(t, n);
}

SparseANF SparseANF::random_cubic(int n, double density, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::unordered_map<uint64_t, uint8_t> terms;

    // degree-1 terms
    for (int i = 0; i < n; ++i) {
        if (dist(rng) < density * 3) {
            terms[uint64_t(1) << i] = 1;
        }
    }
    // degree-2 terms
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (dist(rng) < density) {
                terms[(uint64_t(1) << i) | (uint64_t(1) << j)] = 1;
            }
        }
    }
    // degree-3 terms
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                if (dist(rng) < density) {
                    terms[(uint64_t(1) << i) | (uint64_t(1) << j) | (uint64_t(1) << k)] = 1;
                }
            }
        }
    }
    if (terms.empty()) terms[1] = 1;
    return SparseANF(terms, n);
}

int SparseANF::degree() const {
    int d = 0;
    for (auto& [m, v] : terms) {
        d = std::max(d, __builtin_popcountll(m));
    }
    return d;
}

uint8_t SparseANF::eval_mask(uint64_t x) const {
    uint8_t r = 0;
    for (auto& [m, v] : terms) {
        if ((m & x) == m) r ^= v;
    }
    return r;
}

// Expand f(N(z ⊕ c)): N is n×m (rows = x-vars, cols = z-vars)
SparseANF SparseANF::expand_with(const std::vector<uint64_t>& N, uint64_t c, int m) const {
    std::unordered_map<uint64_t, uint8_t> result;
    for (auto& [x_mask, coeff] : terms) {
        auto expanded = expand_monomial(x_mask, N, c, m);
        for (auto& [z_mask, v] : expanded) {
            result[z_mask] ^= v;
        }
    }
    // Remove zero coeffs
    std::unordered_map<uint64_t, uint8_t> clean;
    for (auto& [m, v] : result) {
        if (v) clean[m] = 1;
    }
    return SparseANF(clean, m);
}

// Expand monomial x^s through x_i = N_i·z ⊕ c_i
// N is n×m matrix, each N_i (row i) is m-bit mask of z variables
// c_i is bit i of c
std::unordered_map<uint64_t, uint8_t> SparseANF::expand_monomial(
    uint64_t x_mask, const std::vector<uint64_t>& N, uint64_t c, int m) const {
    // For each variable x_i in the monomial:
    // x_i = (Σ_{j: N[i][j]=1} z_j) ⊕ c_i
    // This gives a linear form: either L_i(z) or L_i(z) ⊕ 1
    std::vector<std::pair<uint64_t, uint8_t>> forms;  // {mask, const}

    uint64_t t = x_mask;
    while (t) {
        int i = __builtin_ctzll(t);
        t &= t - 1;
        uint64_t row_mask = N[i];  // which z_j appear in x_i's linear form
        uint8_t const_bit = (c >> i) & 1;
        forms.push_back({row_mask, const_bit});
    }

    if (forms.empty()) {
        return {{0, 1}};  // constant 1 monomial
    }

    // Multiply: (L₁⊕c₁)·(L₂⊕c₂)·...·(L_k⊕c_k)
    // Each (L_i⊕c_i) = L_i when c_i=0, or L_i⊕1 when c_i=1
    // In ANF: (L_i⊕1) = L_i ⊕ 1 (XOR of monomials)
    std::unordered_map<uint64_t, uint8_t> result;
    result[0] = 1;  // start with monomial 1

    for (auto& [mask, const_bit] : forms) {
        std::unordered_map<uint64_t, uint8_t> new_result;

        if (const_bit == 0) {
            // (L_i) * current = (z_{b1}⊕z_{b2}⊕...) * current
            // Distribute over each bit: L_i(z)·g = Σ_j z_j·g
            for (auto& [z_mask, v] : result) {
                uint64_t t = mask;
                while (t) {
                    int bit = __builtin_ctzll(t);
                    t &= t - 1;
                    new_result[z_mask | ((uint64_t)1 << bit)] ^= v;
                }
            }
        } else {
            // (L_i ⊕ 1) * current = L_i*current ⊕ 1*current
            for (auto& [z_mask, v] : result) {
                // 1 * current
                new_result[z_mask] ^= v;
                // L_i * current: distribute over each bit
                uint64_t t = mask;
                while (t) {
                    int bit = __builtin_ctzll(t);
                    t &= t - 1;
                    new_result[z_mask | ((uint64_t)1 << bit)] ^= v;
                }
            }
        }
        result = std::move(new_result);
    }

    // Remove zero coeffs and return
    std::unordered_map<uint64_t, uint8_t> clean;
    for (auto& [m, v] : result) {
        if (v & 1) clean[m] = 1;
    }
    return clean;
}

// Structured substitute_affine: each z_j = x_{i(j)} or z_j = x_{i(j)}⊕1
// (each row of M has at most one 1-bit). Direct substitution without inversion.
SparseANF SparseANF::substitute_affine_structured(const std::vector<uint64_t>& M,
                                                   uint64_t b) const {
    int m = (int)M.size();
    // Map each x_i to a preferred z_j (plain over complement)
    std::vector<int> z_of_var(n, -1);
    for (int j = 0; j < m; ++j) {
        if (M[j] == 0) continue;
        int i = __builtin_ctzll(M[j]);
        bool is_comp = (b >> j) & 1;
        if (z_of_var[i] < 0) {
            z_of_var[i] = j;
        } else if (((b >> z_of_var[i]) & 1) && !is_comp) {
            z_of_var[i] = j;  // prefer plain over complement
        }
    }

    std::unordered_map<uint64_t, uint8_t> g_terms;
    for (auto& [mask, v] : terms) {
        uint64_t t = mask;
        uint64_t g_base = 0;
        std::vector<int> comp_rows;  // rows where z_j = x_i⊕1
        bool ok = true;
        while (t) {
            int i = __builtin_ctzll(t);
            t &= t - 1;
            int j = z_of_var[i];
            if (j < 0) { ok = false; break; }
            g_base |= (uint64_t)1 << j;
            if ((b >> j) & 1) comp_rows.push_back(j);
        }
        if (!ok) continue;

        if (comp_rows.empty()) {
            g_terms[g_base] ^= v;
        } else {
            // Expand: Π (z⊕1) = Σ_{subset} Π_{in subset} z
            int k = (int)comp_rows.size();
            for (int sub = 0; sub < (1 << k); ++sub) {
                uint64_t gs = g_base;
                for (int p = 0; p < k; ++p) {
                    if (!((sub >> p) & 1))
                        gs &= ~((uint64_t)1 << comp_rows[p]);
                }
                g_terms[gs] ^= v;
            }
        }
    }
    std::unordered_map<uint64_t, uint8_t> clean;
    for (auto& [tm, tv] : g_terms) {
        if (tv & 1) clean[tm] = 1;
    }
    return SparseANF(clean, m);
}

SparseANF SparseANF::substitute_affine(const std::vector<uint64_t>& M, uint64_t b) const {
    // M is m×n: M[j] is n-bit mask (which x_i appear in z_j)
    // z_j = Σ_i M[j][i]·x_i ⊕ b_j
    int m = (int)M.size();

    if (m == 0) {
        uint8_t const_val = b ? eval_mask(0) : 0;
        // Can't call eval_mask(b) since b is the z-space constant...
        // Actually with m=0, f(x) = g() where g has 0 vars = constant.
        // g = f(0) or f(x) for any x? With m=0, z=() so g() = f(x) for ALL x.
        // This means f must be constant. Let's check: pick x=0 and x=1.
        uint8_t v0 = eval_mask(0);
        uint8_t v1 = eval_mask(b ? 1 : 0);  // hack for m=0
        if (v0 == v1) {
            return SparseANF(v0 ? std::unordered_map<uint64_t, uint8_t>{{0, 1}} : std::unordered_map<uint64_t, uint8_t>{}, 0);
        }
        // inconsistent — return empty
        return SparseANF({}, 0);
    }

    int rank = gf2_rank(M, n);

    if (m == n && rank == n) {
        // Invertible: x = M^{-1}(z ⊕ b)
        GF2Matrix N = gf2_invert(M, n);
        uint64_t c = 0;
        for (int i = 0; i < n; ++i) {
            // c = N·b
            if (__builtin_parityll(N[i] & b)) {
                c |= (uint64_t)1 << i;
            }
        }
        return expand_with(N, c, m);
    }

    if (rank < std::min(m, n)) {
        // Rank-deficient case: stub
        return SparseANF({}, m);
    }

    if (m > n) {
        // Check if M is "structured": each row has at most one 1-bit.
        // This covers z_j = x_{i(j)} or z_j = x_{i(j)}⊕1, no linear combinations.
        // In this case we can do direct substitution without matrix inversion.
        bool structured = true;
        for (auto& row : M) {
            if (row && (row & (row - 1))) { structured = false; break; }
        }
        if (structured) {
            return substitute_affine_structured(M, b);
        }
        if (rank == n) {
            // Full column rank: left inverse N·M = I_n
            GF2Matrix M_ext = gf2_extend_columns_to_invertible(M, m, n);
            GF2Matrix M_inv = gf2_invert(M_ext, m);
            if (M_inv.empty()) return SparseANF({}, m);
            GF2Matrix N(M_inv.begin(), M_inv.begin() + n);
            uint64_t c = gf2_mat_vec_mul(N, b);
            return expand_with(N, c, m);
        }
        // rank < n: not supported (information loss)
        return SparseANF({}, m);
    }

    // m < n, full row rank: extend to invertible, transform, restrict
    GF2Matrix M_aug = gf2_extend_to_invertible(M, m, n);
    uint64_t b_aug = b;  // b_missing = 0 for extended rows

    // Invert M_aug to get N
    GF2Matrix N_aug = gf2_invert(M_aug, n);

    // Compute c = N_aug · b_aug
    uint64_t c_aug = 0;
    for (int i = 0; i < n; ++i) {
        if (__builtin_parityll(N_aug[i] & b_aug)) {
            c_aug |= (uint64_t)1 << i;
        }
    }

    // Expand: f' = f(N_aug(z_aug ⊕ c_aug))
    SparseANF f_prime = expand_with(N_aug, c_aug, n);

    // Restrict to first m variables
    uint64_t z_mask = low_mask(m);
    std::unordered_map<uint64_t, uint8_t> g_terms;
    for (auto& [zm, v] : f_prime.terms) {
        uint64_t g_mask = zm & z_mask;
        g_terms[g_mask] ^= v;
    }

    std::unordered_map<uint64_t, uint8_t> clean;
    for (auto& [m, v] : g_terms) {
        if (v) clean[m] = 1;
    }
    return SparseANF(clean, m);
}

// Z = Z1 ∪ Z2: Z1 = X (n vars), Z2 = Mx⊕b (m vars, m ≤ n recommended)
// Uses left inverse path directly to avoid structured shortcut.
SparseANF SparseANF::substitute_affine_union(const std::vector<uint64_t>& M,
                                              uint64_t b) const {
    int m = (int)M.size();
    int total = m + n;
    // Build M_ext = [M; I_n] ((m+n)×n), b_ext = [b; 0] ((m+n)-bit)
    GF2Matrix M_ext = M;
    for (int i = 0; i < n; ++i)
        M_ext.push_back((uint64_t)1 << i);
    uint64_t b_ext = b;
    // Left inverse: N·M_ext = I_n
    auto M_sq = gf2_extend_columns_to_invertible(M_ext, total, n);
    auto M_inv = gf2_invert(M_sq, total);
    if (M_inv.empty()) return SparseANF({}, total);
    GF2Matrix N(M_inv.begin(), M_inv.begin() + n);
    uint64_t c = gf2_mat_vec_mul(N, b_ext);
    return expand_with(N, c, total);
}

std::vector<int> SparseANF::variables_used() const {
    std::vector<int> used;
    uint64_t seen = 0;
    for (auto& [mask, v] : terms) {
        seen |= mask;
    }
    uint64_t t = seen;
    while (t) {
        int i = __builtin_ctzll(t);
        t &= t - 1;
        used.push_back(i);
    }
    return used;
}

SparseANF SparseANF::partial_deriv(int var) const {
    std::unordered_map<uint64_t, uint8_t> result;
    for (auto& [mask, v] : terms) {
        if ((mask >> var) & 1) {
            uint64_t rest = mask ^ ((uint64_t)1 << var);
            result[rest] ^= v;
        }
        // variable i not in monomial → derivative contribution is 0
    }
    std::unordered_map<uint64_t, uint8_t> clean;
    for (auto& [m, val] : result) {
        if (val & 1) clean[m] = 1;
    }
    return SparseANF(clean, n);
}

std::vector<SparseANF> SparseANF::gradient() const {
    std::vector<SparseANF> grads(n);
    for (int i = 0; i < n; ++i) grads[i] = partial_deriv(i);
    return grads;
}

void SparseANF::print() const {
    std::cout << to_string() << "\n";
}

std::string SparseANF::to_string() const {
    std::ostringstream os;
    os << "SparseANF(n=" << n << ", T=" << T() << "): ";
    bool first = true;
    for (auto& [m, v] : terms) {
        if (!first) os << " ⊕ ";
        first = false;
        if (m == 0) {
            os << "1";
        } else {
            uint64_t t = m;
            while (t) {
                int i = __builtin_ctzll(t);
                t &= t - 1;
                os << "x" << i;
            }
        }
    }
    return os.str();
}

// Apply x_i → x_i⊕x_j
std::unordered_map<uint64_t, uint8_t> apply_xor_merge(
    const std::unordered_map<uint64_t, uint8_t>& terms, int i, int j) {
    std::unordered_map<uint64_t, uint8_t> new_terms;
    for (auto& [m, v] : terms) {
        if ((m >> i) & 1) {
            new_terms[m] ^= v;  // toggle off
            uint64_t rest = m ^ (uint64_t(1) << i);
            uint64_t new_m = ((m >> j) & 1) ? rest : (rest | (uint64_t(1) << j));
            new_terms[new_m] ^= v;  // toggle on
        } else {
            new_terms[m] ^= v;
        }
    }
    std::unordered_map<uint64_t, uint8_t> clean;
    for (auto& [m, v] : new_terms) {
        if (v & 1) clean[m] = 1;
    }
    return clean;
}

bool verify_substitution(const SparseANF& f, const SparseANF& g,
                         const std::vector<uint64_t>& M, uint64_t b,
                         int n_tests) {
    int m = (int)M.size();
    std::mt19937_64 rng(12345);
    for (int t = 0; t < n_tests; ++t) {
        uint64_t x = rng() & low_mask(f.n);
        // Compute z = Mx ⊕ b
        uint64_t z = 0;
        for (int j = 0; j < m; ++j) {
            if (__builtin_parityll(M[j] & x)) {
                z |= (uint64_t)1 << j;
            }
        }
        z ^= b;
        z &= low_mask(m);

        if (f.eval_mask(x) != g.eval_mask(z)) {
            return false;
        }
    }
    return true;
}
