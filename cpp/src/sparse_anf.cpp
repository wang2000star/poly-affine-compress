#include "sparse_anf.h"
#include <algorithm>
#include <sstream>
#include <cassert>
#include <cstdio>

// ====================================================================
//  Constructor
// ====================================================================

SparseANF::SparseANF(const AnfMap& terms, int n) : n_(n) {
    for (auto& [m, v] : terms) {
        if (v & 1) terms_[m] = 1;
    }
}

// ====================================================================
//  Query
// ====================================================================

int SparseANF::T_nonlinear() const {
    int cnt = 0;
    for (auto& [m, v] : terms_) {
        if (v && __builtin_popcountll(m) >= 2) cnt++;
    }
    return cnt;
}

int SparseANF::degree() const {
    int d = 0;
    for (auto& [m, v] : terms_) {
        if (v) {
            int pc = __builtin_popcountll(m);
            if (pc > d) d = pc;
        }
    }
    return d;
}

bool SparseANF::is_one() const {
    auto it = terms_.find(0);
    return it != terms_.end() && it->second == 1 && terms_.size() == 1;
}

int SparseANF::eval_mask(uint64_t x) const {
    int r = 0;
    for (auto& [m, v] : terms_) {
        if (v && (m & x) == m) r ^= 1;
    }
    return r;
}

std::string SparseANF::to_string(const std::vector<std::string>& var_names) const {
    if (terms_.empty()) return "0";
    std::vector<std::pair<uint64_t, int>> sorted(terms_.begin(), terms_.end());
    std::sort(sorted.begin(), sorted.end());
    std::vector<std::string> parts;
    for (auto& [m, v] : sorted) {
        if (!v) continue;
        if (m == 0) { parts.emplace_back("1"); continue; }
        std::string term;
        bool first = true;
        for (int j = 0; j < 64; j++) {
            if ((m >> j) & 1) {
                if (!first) term += " * ";
                term += (j < (int)var_names.size()) ? var_names[j] : ("z_" + std::to_string(j));
                first = false;
            }
        }
        parts.push_back(term);
    }
    if (parts.empty()) return "0";
    std::string r = parts[0];
    for (size_t i = 1; i < parts.size(); i++) r += " + " + parts[i];
    return r;
}

// ====================================================================
//  Mutators
// ====================================================================

void SparseANF::xor_op(const SparseANF& other) {
    if (this == &other) { terms_.clear(); return; }
    for (auto& [m, v] : other.terms_) {
        if (!v) continue;
        auto it = terms_.find(m);
        if (it != terms_.end()) {
            it->second ^= v;
            if (it->second == 0) terms_.erase(it);
        } else {
            terms_[m] = v & 1;
        }
    }
    if (other.n_ > n_) n_ = other.n_;
}

void SparseANF::xor_const() {
    auto it = terms_.find(0);
    if (it != terms_.end()) {
        it->second ^= 1;
        if (it->second == 0) terms_.erase(it);
    } else {
        terms_[0] = 1;
    }
}

void SparseANF::and_op(const SparseANF& other) {
    if (is_zero() || other.is_zero()) { terms_.clear(); return; }
    if (is_one()) { *this = other; return; }
    if (other.is_one()) return;

    AnfMap result;
    for (auto& [m1, v1] : terms_) {
        if (!v1) continue;
        for (auto& [m2, v2] : other.terms_) {
            if (!v2) continue;
            uint64_t m = m1 | m2;  // Boolean ring: OR
            result[m] ^= 1;
        }
    }
    for (auto it = result.begin(); it != result.end(); ) {
        if (it->second == 0) it = result.erase(it); else ++it;
    }
    terms_ = std::move(result);
    if (other.n_ > n_) n_ = other.n_;
}

void SparseANF::shift(int amount) {
    if (amount == 0) return;
    AnfMap shifted;
    for (auto& [m, v] : terms_) if (v) shifted[m << amount] = 1;
    terms_ = std::move(shifted);
    n_ += amount;
}

void SparseANF::complement_var(int j) {
    // z_j = z_j ⊕ 1 → for each monomial with j, toggle version without j
    uint64_t bit = 1ULL << j;
    // Collect toggles first
    AnfMap toggles;
    for (auto& [m, v] : terms_) {
        if (v && (m & bit)) toggles[m ^ bit] ^= 1;
    }
    for (auto& [m, v] : toggles) {
        if (!v) continue;
        auto it = terms_.find(m);
        if (it != terms_.end()) { it->second ^= 1; if (it->second == 0) terms_.erase(it); }
        else terms_[m] = 1;
    }
}

void SparseANF::substitute_zi_eq_zi_xor_zj(int i, int j) {
    // z_i = z_i ⊕ z_j → for each monomial with i, toggle (without-i, with-j)
    uint64_t bit_i = 1ULL << i;
    uint64_t bit_j = 1ULL << j;
    AnfMap toggles;
    for (auto& [m, v] : terms_) {
        if (v && (m & bit_i)) {
            uint64_t new_mask = (m ^ bit_i) | bit_j;
            toggles[new_mask] ^= 1;
        }
    }
    for (auto& [m, v] : toggles) {
        if (!v) continue;
        auto it = terms_.find(m);
        if (it != terms_.end()) { it->second ^= 1; if (it->second == 0) terms_.erase(it); }
        else terms_[m] = 1;
    }
}

// ====================================================================
//  extract_linear_terms
// ====================================================================

std::vector<int> SparseANF::extract_linear_terms() {
    std::vector<int> L(n_, 0);
    for (auto it = terms_.begin(); it != terms_.end(); ) {
        uint64_t m = it->first;
        if (m != 0 && __builtin_popcountll(m) == 1) {
            int j = __builtin_ctzll(m);
            L[j] ^= it->second;
            it = terms_.erase(it);
        } else {
            ++it;
        }
    }
    return L;
}

// ====================================================================
//  substitute_affine
// ====================================================================

SparseANF SparseANF::substitute_affine(
    const std::vector<uint32_t>& M_rows,
    uint64_t b, int new_n) const
{
    // For each monomial in g, expand through z = M_rows * x ⊕ b
    AnfMap result;
    for (auto& [mask, coeff] : terms_) {
        if (!coeff) continue;
        AnfMap expanded = _expand_monomial(mask, M_rows, b, new_n);
        for (auto& [m, v] : expanded) {
            if (v) result[m] ^= 1;
        }
    }
    for (auto it = result.begin(); it != result.end(); ) {
        if (it->second == 0) it = result.erase(it); else ++it;
    }
    return SparseANF(result, new_n);
}

AnfMap SparseANF::_expand_monomial(
    uint64_t mask, const std::vector<uint32_t>& M_rows,
    uint64_t b, int new_n) const
{
    // Expand one monomial through z_j = <M_rows[j], x> ⊕ b_j
    // Boolean ring: z_j² = z_j → masks combine with OR
    AnfMap result;
    result[0] = 1;

    for (int j = 0; j < (int)M_rows.size(); j++) {
        if (!((mask >> j) & 1)) continue;

        uint32_t row = M_rows[j];
        bool has_b = (b >> j) & 1;
        AnfMap next;

        for (auto& [z_mask, val] : result) {
            if (!val) continue;
            if (has_b)
                next[z_mask] ^= 1;  // constant contribution from b_j
            // Each set bit in row → a variable
            uint32_t r = row;
            while (r) {
                int k = __builtin_ctz(r);
                r &= r - 1;
                next[z_mask | (1ULL << k)] ^= 1;
            }
        }
        result = std::move(next);
    }
    return result;
}

// ====================================================================
//  Complement search: greedy
// ====================================================================

void SparseANF::complement_search_greedy(
    std::vector<uint32_t>& out_M_rows,
    uint64_t& out_b)
{
    out_M_rows.clear();
    out_b = 0;
    int n = n_;
    for (int j = 0; j < n; j++)
        out_M_rows.push_back(1U << j);

    int best_T = T_nonlinear();

    // Iterate multiple passes until no improvement
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i < n; i++) {
            // Try complementing z_i
            SparseANF saved = *this;
            complement_var(i);
            int cur_T = T_nonlinear();
            if (cur_T < best_T) {
                best_T = cur_T;
                out_b ^= (1ULL << i);
                improved = true;
                // keep the change
            } else {
                // revert
                *this = std::move(saved);
            }
        }
    }
}
