#include "search.h"
#include "affine.h"
#include "anf.h"
#include "gf2.h"
#include "io.h"
#include "moebius.h"
#include "walsh.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <thread>
#include <csignal>
#include <cstring>
#include <filesystem>

// ============================================================
//  CandidateGenerator
// ============================================================

struct CandidateGenerator {
    int n;
    int max_m_search;
    const std::vector<WalshInfo>& walsh;
    std::mt19937_64 rng;

    CandidateGenerator(int n, int max_m, const std::vector<WalshInfo>& w)
        : n(n), max_m_search(max_m), walsh(w), rng(42) {}

    // Single-row candidates from Walsh (each bit with high |W_i|)
    std::vector<std::pair<uint32_t, uint32_t>> gen_walsh_single_rows(int top_k = 40) {
        struct BitScore { int bit; double score; };
        std::vector<BitScore> scores(n);
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (auto& w : walsh)
                s += (double)w.walsh_mag[i] / (1 << (n - 1));
            scores[i] = {i, s};
        }
        std::sort(scores.begin(), scores.end(),
            [](auto& a, auto& b) { return a.score > b.score; });

        std::vector<std::pair<uint32_t, uint32_t>> candidates;
        int k = std::min(top_k, n);
        for (int j = 0; j < k; j++) {
            uint32_t row = 1u << scores[j].bit;
            candidates.emplace_back(row, 0u);
            candidates.emplace_back(row, 1u);
        }
        return candidates;
    }

    // Multi-row candidates by combining top Walsh bits
    std::vector<std::pair<std::vector<uint32_t>, uint32_t>> gen_multi_row(int max_rows = 12) {
        struct BitScore { int bit; double score; };
        std::vector<BitScore> scores(n);
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (auto& w : walsh)
                s += (double)w.walsh_mag[i] / (1 << (n - 1));
            scores[i] = {i, s};
        }
        std::sort(scores.begin(), scores.end(),
            [](auto& a, auto& b) { return a.score > b.score; });

        std::vector<std::pair<std::vector<uint32_t>, uint32_t>> candidates;

        for (int m = 1; m <= std::min(max_rows, n); m++) {
            std::vector<uint32_t> rows;
            for (int j = 0; j < m; j++)
                rows.push_back(1u << scores[j].bit);
            candidates.emplace_back(rows, 0u);
        }

        int top_n = std::min(8, n);
        for (int i = 0; i < top_n; i++) {
            for (int j = i + 1; j < top_n; j++) {
                uint32_t row = (1u << scores[i].bit) | (1u << scores[j].bit);
                candidates.emplace_back(std::vector<uint32_t>{row}, 0u);
                candidates.emplace_back(std::vector<uint32_t>{row}, 1u);
            }
        }
        return candidates;
    }

    // Random M,b candidates
    std::vector<std::pair<std::vector<uint32_t>, uint32_t>> gen_random(
        int n_candidates, int max_m)
    {
        std::vector<std::pair<std::vector<uint32_t>, uint32_t>> candidates;
        for (int c = 0; c < n_candidates; c++) {
            int m = (rng() % max_m) + 1;
            std::vector<uint32_t> rows(m, 0);
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    if (rng() & 1) rows[j] |= (1u << i);
                }
            }
            uint64_t b = 0;
            for (int j = 0; j < m; j++) {
                if (rng() & 1) b |= (1ULL << j);
            }
            candidates.emplace_back(rows, b);
        }
        return candidates;
    }
};

// ============================================================
//  Dependency set computation (theory-guided)
//
//  For each output, determine which input variables actually
//  affect it (i.e., e_j ∉ Aut(f_i)).
//
//  n ≤ 20: exact via truth table scan
//  n > 20: random sampling (fallback)
// ============================================================

std::vector<uint32_t> compute_dependency_set(
    const TruthTable& tt, int n, int n_samples)
{
    int k = tt.n_outputs;
    std::vector<uint32_t> dep(k, 0);

    if (n <= 20) {
        // Exact: scan truth table
        int64_t n_words = tt.n_words;
        for (int j = 0; j < n; j++) {
            uint32_t bit = (1u << j);
            uint32_t word_offset = (j < 6) ? 0 : (1u << (j - 6));
            int bit_in_word = j & 63;
            uint64_t bit_mask = (1ULL << bit_in_word);

            for (int oi = 0; oi < k; oi++) {
                bool depends = false;
                for (int64_t w = 0; w < n_words && !depends; w++) {
                    uint64_t val = tt.tt[oi][w];
                    uint64_t shifted;
                    if (j < 6) {
                        // Compare f(x) with f(x^(1<<j)) via paired bits within same word
                        int shift = 1 << j;
                        uint64_t mask0 = 0;
                        for (int b = 0; b < 64; b += 2 * shift)
                            mask0 |= ((1ULL << shift) - 1) << b;
                        uint64_t val0 = val & mask0;
                        uint64_t val1 = (val >> shift) & mask0;
                        if (val0 != val1) { depends = true; dep[oi] |= bit; continue; }
                        continue;
                    } else {
                        // Bit flip across words
                        int64_t w_pair = w ^ word_offset;
                        if (w_pair >= n_words) continue;
                        shifted = tt.tt[oi][w_pair];
                    }
                    if (val != shifted) {
                        depends = true;
                        dep[oi] |= bit;
                    }
                }
            }
        }
    } else {
        // n > 20: random sampling
        for (int j = 0; j < n; j++) {
            uint32_t bit = (1u << j);
            for (int oi = 0; oi < k; oi++) {
                for (int t = 0; t < n_samples; t++) {
                    // Random x: use deterministic hashing
                    uint64_t x = (uint64_t)t * 0x9e3779b97f4a7c15ULL;
                    x = (x >> 32) ^ x;
                    if (n < 32) x &= ((1ULL << n) - 1);
                    uint64_t x_flip = x ^ (1ULL << j);
                    // Read from truth table (only valid for n≤25)
                    // For n>25, would need circuit evaluation
                    uint64_t f_x = (tt.tt[oi][x >> 6] >> (x & 63)) & 1;
                    uint64_t f_xf = (tt.tt[oi][x_flip >> 6] >> (x_flip & 63)) & 1;
                    if (f_x != f_xf) {
                        dep[oi] |= bit;
                        break;
                    }
                }
            }
        }
    }
    return dep;
}

// ============================================================
//  Aut(F) computation (n ≤ 20 via Walsh power spectrum)
//
//  For each output:
//    1. Compute Walsh-Hadamard transform
//    2. Square magnitudes (power spectrum)
//    3. Inverse Walsh → autocorrelation r(Δ)
//    4. r(Δ) = 1 ± ε → Δ ∈ Aut(f_i)
//    5. Intersect across all outputs → Aut(F)
//  Output: basis vectors of Aut(F) as bitmasks
// ============================================================

std::vector<uint64_t> compute_aut_basis(
    const TruthTable& tt, int n, int n_threads)
{
    int k = tt.n_outputs;
    int64_t n_pts = int64_t(1) << n;
    int64_t n_words = tt.n_words;

    if (n > 20) {
        // For n > 20, Aut(F) basis can't be computed via Walsh easily
        // Return empty basis (= dim 0, no invariance)
        return {};
    }

    // Allocate working arrays
    std::vector<double> walsh_re(n_pts, 0.0);
    std::vector<double> walsh_im(n_pts, 0.0);
    std::vector<double> power(n_pts, 0.0);
    std::vector<double> autocorr(n_pts, 0.0);

    // For each output, compute autocorrelation
    std::vector<uint64_t> aut_set(n_pts, ~0ULL); // start as all-1

    for (int oi = 0; oi < k; oi++) {
        // Initialize real part with (-1)^f(x)
        for (int64_t i = 0; i < n_pts; i++) {
            uint64_t bit = (tt.tt[oi][i >> 6] >> (i & 63)) & 1;
            walsh_re[i] = bit ? -1.0 : 1.0;
            walsh_im[i] = 0.0;
        }

        // Walsh-Hadamard transform (in-place)
        for (int step = 1; step < n_pts; step <<= 1) {
            for (int64_t i = 0; i < n_pts; i += step * 2) {
                for (int64_t j = 0; j < step; j++) {
                    double u_re = walsh_re[i + j];
                    double u_im = walsh_im[i + j];
                    double v_re = walsh_re[i + j + step];
                    double v_im = walsh_im[i + j + step];
                    walsh_re[i + j] = u_re + v_re;
                    walsh_im[i + j] = u_im + v_im;
                    walsh_re[i + j + step] = u_re - v_re;
                    walsh_im[i + j + step] = u_im - v_im;
                }
            }
        }

        // Power spectrum S[α] = |Ŵ(α)|²
        for (int64_t i = 0; i < n_pts; i++) {
            power[i] = walsh_re[i] * walsh_re[i] + walsh_im[i] * walsh_im[i];
        }

        // Inverse Walsh on power spectrum → autocorrelation r(Δ)
        // r(Δ) = (1/2ⁿ) · Σ_α S[α] · (-1)^{α·Δ}
        // Using Walsh on power gives 2ⁿ·r(Δ) (Walsh is self-inverse up to scale)
        for (int64_t i = 0; i < n_pts; i++) {
            walsh_re[i] = power[i];
            walsh_im[i] = 0.0;
        }
        for (int step = 1; step < n_pts; step <<= 1) {
            for (int64_t i = 0; i < n_pts; i += step * 2) {
                for (int64_t j = 0; j < step; j++) {
                    double u_re = walsh_re[i + j];
                    double u_im = walsh_im[i + j];
                    double v_re = walsh_re[i + j + step];
                    double v_im = walsh_im[i + j + step];
                    walsh_re[i + j] = u_re + v_re;
                    walsh_im[i + j] = u_im + v_im;
                    walsh_re[i + j + step] = u_re - v_re;
                    walsh_im[i + j + step] = u_im - v_im;
                }
            }
        }

        // Normalize: r(Δ) = power_walsh[Δ] / 2ⁿ
        // And find Δ with r(Δ) close to 1
        uint64_t aut_i = 1; // Δ=0 always included
        for (int64_t d = 1; d < n_pts; d++) {
            double r = walsh_re[d] / (double)n_pts;
            // r(Δ) = 1 means perfect autocorrelation (= Δ ∈ Aut(f))
            if (r > 0.9999) {
                aut_i |= (1ULL << d);
            }
        }
        // Intersect with other outputs
        aut_set[oi] = aut_i;
    }

    // Intersection across all outputs
    uint64_t aut_all = aut_set[0];
    for (int oi = 1; oi < k; oi++) {
        aut_all &= aut_set[oi];
    }

    // Extract basis vectors from aut_all
    std::vector<uint64_t> basis;
    if (aut_all == 0) return basis;

    // aut_all is a bitmask: bit d = 1 means d ∈ Aut(F)
    // We need a basis of the subspace.
    // Convert from vector set to basis via Gaussian elimination.
    std::vector<uint64_t> vectors;
    uint64_t mask = aut_all;
    while (mask) {
        int d = __builtin_ctzll(mask);
        mask &= mask - 1;
        vectors.push_back((uint64_t)d);
    }

    if (vectors.empty()) return basis;

    // Gaussian elimination to find basis
    int dim = n;
    std::vector<uint64_t> basis_mat(dim, 0);

    for (uint64_t v : vectors) {
        uint64_t vec = v;
        for (int i = 0; i < dim; i++) {
            if (!((vec >> i) & 1)) continue;
            if (basis_mat[i] == 0) {
                basis_mat[i] = vec;
                break;
            }
            vec ^= basis_mat[i];
        }
    }

    for (int i = 0; i < dim; i++) {
        if (basis_mat[i] != 0)
            basis.push_back(basis_mat[i]);
    }

    return basis;
}

// ============================================================
//  Progressive M construction (theory-guided)
//
//  Start from identity (m=n, ker={0}, always consistent).
//  For each candidate row w from the 3n pool:
//    1. Check exclusivity: rank(M ∪ [w]) > rank(M)
//    2. Find Δ: M·Δ = 0, w·Δ = 1
//    3. Δ ∈ Aut(F)? → skip
//    4. Evaluate M ∪ [w]; keep if T improves
//  Then hill climb from best result.
// ============================================================

MbCandidate progressive_m_search(
    const TruthTable& tt,
    const std::vector<uint64_t>& aut_basis,
    int n, int max_m, int n_threads)
{
    int k = tt.n_outputs;
    int effective_max = (max_m <= 0) ? n : std::min(max_m, n);

    // Start from identity (m=n)
    std::vector<uint32_t> curr_M(n);
    uint64_t curr_b = 0;
    for (int i = 0; i < n; i++)
        curr_M[i] = (1u << i);

    MbCandidate best = evaluate_Mb(tt, curr_M.data(), curr_b, n, n_threads);

    if (effective_max <= n) {
        // Can't increase m, just return identity baseline
        // Hill climb from identity is already done in run_search
        return best;
    }

    // Build pool rows
    int m_target = std::min(effective_max, 2 * n); // max m target
    if (m_target <= n) return best;

    // Generate complement + XOR pair candidate rows
    struct PoolRow { uint32_t row; uint64_t b_bit; double score; };
    std::vector<PoolRow> pool;

    // Complement rows: e_i with b=1
    for (int i = 0; i < n; i++)
        pool.push_back({(1u << i), 1ULL, 0.0});

    // XOR pair rows: e_i ⊕ e_j with b=0
    for (int i = 0; i < std::min(n, 8); i++) {
        for (int j = i + 1; j < std::min(n, 8); j++) {
            pool.push_back({(1u << i) | (1u << j), 0ULL, 0.0});
        }
    }

    // Progressive construction
    for (auto& pr : pool) {
        if ((int)curr_M.size() >= m_target) break;

        uint32_t w = pr.row;
        uint64_t wb = pr.b_bit;

        // Step 1: exclusivity check - is w linearly independent of curr_M?
        int r = (int)curr_M.size();
        uint32_t test_M[33];
        for (int i = 0; i < r; i++) test_M[i] = curr_M[i];
        test_M[r] = w;
        if (gf2_rank(test_M, r + 1, n) == gf2_rank(curr_M.data(), r, n))
            continue; // w is in row space of curr_M

        // Step 2: Solve M·Δ = 0, w·Δ = 1
        // Build augmented matrix [M; w] and find nullspace basis, then a Δ
        // For simplicity, use GF(2) rank test + random search for Δ
        uint64_t delta = 0;
        bool found_delta = false;

        // Try single-bit Δs (e_j) as candidates
        for (int j = 0; j < n && !found_delta; j++) {
            uint64_t d = (1ULL << j);
            // Check M·d = 0
            bool m_dot_d_zero = true;
            for (int ri = 0; ri < r && m_dot_d_zero; ri++) {
                if (__builtin_popcount(curr_M[ri] & (uint32_t)d) & 1)
                    m_dot_d_zero = false;
            }
            if (!m_dot_d_zero) continue;
            // Check w·d = 1
            if ((__builtin_popcount(w & (uint32_t)d) & 1)) {
                delta = d;
                found_delta = true;
            }
        }

        if (!found_delta) continue;

        // Step 3: Check Δ ∈ span(Aut(F))?
        // Use Gaussian elimination: reduce delta by basis vectors
        uint64_t d_check = delta;
        for (uint64_t aut_vec : aut_basis) {
            if (d_check == 0) break;
            uint64_t bv = aut_vec;
            if (bv == 0) continue;
            int top = 63 - __builtin_clzll(bv);
            if ((d_check >> top) & 1)
                d_check ^= bv;
        }
        if (d_check == 0) continue; // delta ∈ Aut(F), skip

        // Step 4: Evaluate candidate
        std::vector<uint32_t> new_M = curr_M;
        new_M.push_back(w);
        uint64_t new_b = curr_b;
        if (wb) new_b |= (1ULL << r); // set b bit for new row

        MbCandidate cand = evaluate_Mb(tt, new_M.data(), new_b,
                                        (int)new_M.size(), n_threads);
        if (cand.union_T < best.union_T && cand.union_T != INT64_MAX) {
            best = cand;
            curr_M = new_M;
            curr_b = new_b;
            std::cout << "    progressive: m=" << new_M.size()
                      << " T_union=" << cand.union_T << "\n";
        }
    }

    return best;
}

// Detect complement pairs in M,b: pairs of single-bit rows with same bit but opposite b.
// Returns number of pairs found, fills pair_masks with (1<<row_a)|(1<<row_b) for each.
static int detect_complement_pairs(const uint32_t* M_rows, uint64_t b, int m,
                                    uint64_t* pair_masks_out, int max_pairs)
{
    int np = 0;
    for (int i = 0; i < m && np < max_pairs; i++) {
        uint32_t ri = M_rows[i];
        if (ri == 0) continue;
        int bit_i = __builtin_ctz(ri);
        if (ri != (1u << bit_i)) continue;
        for (int j = i + 1; j < m && np < max_pairs; j++) {
            uint32_t rj = M_rows[j];
            if (rj == 0) continue;
            int bit_j = __builtin_ctz(rj);
            if (rj != (1u << bit_j)) continue;
            if (bit_i == bit_j && ((b >> i) & 1) != ((b >> j) & 1))
                pair_masks_out[np++] = (1ULL << i) | (1ULL << j);
        }
    }
    return np;
}

// Evaluate M,b with complement-pair filtering (if any pairs detected).
static MbCandidate evaluate_pool_rowset(
    const TruthTable& tt,
    const uint32_t* M_rows, uint64_t b, int m,
    int n, int n_threads)
{
    // Detect complement pairs
    uint64_t pair_masks[16];
    int np = detect_complement_pairs(M_rows, b, m, pair_masks, 16);

    MbCandidate cand = evaluate_Mb(tt, M_rows, b, m, n_threads, np > 0);
    if (cand.union_T >= INT64_MAX) return cand;

    if (np == 0) return cand;  // no filtering needed

    // Filter: Möbius → remove pair monomials → inverse Möbius
    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    std::vector<uint64_t> union_anf(n_words_g, 0);
    int64_t total_T = 0;

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        auto& raw = cand.g_tt_raw[oi];
        if (raw.empty()) { cand.per_output_T[oi] = 0; continue; }
        moebius_packed(raw.data(), m);
        filter_pairs(raw.data(), m, pair_masks, np);
        int64_t T_filt = count_T(raw.data(), m);
        cand.per_output_T[oi] = T_filt;
        total_T += T_filt;
        for (int64_t w = 0; w < n_words_g; w++)
            union_anf[w] |= raw[w];
        moebius_packed(raw.data(), m);
    }
    cand.total_T = total_T;
    cand.union_T = count_T(union_anf.data(), m);
    return cand;
}

// Build the 3n-row pool and score each row by Walsh magnitude.
// Row types:
//   0..n-1:   identity, b=0
//   n..2n-1:  complement, b=1
//   2n..3n-1: XOR of two top-Walsh bits (each row has ≥2 ones)
// Returns a vector of {M_row, b_row, score}.
static void build_row_pool(int n,
    const std::vector<WalshInfo>& walsh,
    std::vector<std::pair<uint32_t, uint32_t>>& pool_rows,
    std::vector<double>& scores,
    uint32_t dep_union = 0)
{
    pool_rows.clear();
    scores.clear();

    // Score each input bit by Walsh magnitude
    struct BitScore { int bit; double score; };
    std::vector<BitScore> bit_scores(n);
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (auto& w : walsh)
            s += (double)w.walsh_mag[i] / (1 << (n - 1));
        bit_scores[i] = {i, s};
    }
    std::sort(bit_scores.begin(), bit_scores.end(),
        [](auto& a, auto& b) { return a.score > b.score; });

    // Group 1: identity rows (b=0) and complement rows (b=1)
    // Skip bits that are don't-care (not in dep_union)
    for (int i = 0; i < n; i++) {
        uint32_t row = (1u << i);
        if (dep_union != 0 && !(dep_union & row)) continue;
        pool_rows.emplace_back(row, 0u);
        scores.push_back(bit_scores[i].score);
        pool_rows.emplace_back(row, 1u);
        scores.push_back(bit_scores[i].score);
    }

    // Group 2: XOR pairs of top Walsh bits (≥2 ones per row)
    int top_n = std::min(n, 8);
    for (int a = 0; a < top_n; a++) {
        for (int b_idx = a + 1; b_idx < top_n; b_idx++) {
            uint32_t row = (1u << bit_scores[a].bit) | (1u << bit_scores[b_idx].bit);
            if (dep_union != 0 && (row & dep_union) == 0) continue;
            double score = bit_scores[a].score + bit_scores[b_idx].score;
            pool_rows.emplace_back(row, 0u);
            scores.push_back(score);
            pool_rows.emplace_back(row, 1u);
            scores.push_back(score);
        }
    }
}

// Generate candidates from the 3n-row pool.
// Strategy: start from identity transform, try swapping/replacing rows.
static void add_row_pool_candidates(
    const TruthTable& tt,
    const std::vector<WalshInfo>& walsh,
    int n, int n_threads, int max_candidates,
    std::vector<MbCandidate>& results,
    uint32_t dep_union = 0)
{
    if (walsh.empty() || max_candidates <= 0) return;
    int pool_size = 2 * n + 2 * std::min(n, 8) * (std::min(n, 8) - 1) / 2;
    if (pool_size <= 0) return;

    std::vector<std::pair<uint32_t, uint32_t>> pool_rows;
    std::vector<double> scores;
    build_row_pool(n, walsh, pool_rows, scores, dep_union);

    // Index each bit's identity row and complement row in the pool.
    // Pool layout: first n = identity(b=0), next n = identity(b=1), then XOR-pair rows.
    // So for bit i: id_row = i, comp_row = i + n.
    std::vector<int> id_idx(n), comp_idx(n);
    for (int i = 0; i < n; i++) {
        id_idx[i] = i;
        comp_idx[i] = i + n;
    }

    int generated = 0;
    auto try_Mb = [&](const std::vector<int>& selected, const char* label) {
        if (generated >= max_candidates) return;
        int m = (int)selected.size();
        uint32_t M[32] = {0};
        uint64_t bv = 0;
        for (int j = 0; j < m; j++) {
            M[j] = pool_rows[selected[j]].first;
            if (pool_rows[selected[j]].second) bv |= (1ULL << j);
        }
        auto cand = evaluate_pool_rowset(tt, M, bv, m, n, n_threads);
        if (cand.union_T < INT64_MAX) {
            results.push_back(cand);
            generated++;
            std::cout << "      pool " << label << " m=" << m
                      << " T_union=" << cand.union_T << "\n";
        }
    };

    // 1) Identity transform (baseline)
    {
        std::vector<int> sel(n);
        for (int i = 0; i < n; i++) sel[i] = id_idx[i];
        try_Mb(sel, "identity");
    }

    // 2) Full complement set: identity + complement rows (m=2n)
    if (generated < max_candidates) {
        std::vector<int> sel;
        for (int i = 0; i < n; i++) {
            sel.push_back(id_idx[i]);
            sel.push_back(comp_idx[i]);
        }
        try_Mb(sel, "full-complement");
    }

    // 3) Try swapping each identity row with its complement (one at a time)
    //    M has identity rows, but bit i is swapped to complement (b=1).
    if (generated < max_candidates) {
        for (int i = 0; i < n && generated < max_candidates; i++) {
            std::vector<int> sel;
            for (int j = 0; j < n; j++)
                sel.push_back(j == i ? comp_idx[j] : id_idx[j]);
            try_Mb(sel, "swap-comp");
        }
    }

    // 4) Try swapping pairs of identity rows with complements
    if (generated < max_candidates) {
        for (int i = 0; i < n && generated < max_candidates; i++) {
            for (int j = i + 1; j < n && generated < max_candidates; j++) {
                std::vector<int> sel;
                for (int k = 0; k < n; k++)
                    sel.push_back((k == i || k == j) ? comp_idx[k] : id_idx[k]);
                try_Mb(sel, "swap2-comp");
            }
        }
    }

    // 5) Greedy: start from best complement mix, try flipping more bits
    if (generated < max_candidates) {
        // Find best combination so far by scanning results
        int64_t best_union = INT64_MAX;
        std::vector<int> best_config(n, 0); // 0 = identity, 1 = complement
        for (auto& r : results) {
            if (r.m == n && r.union_T < best_union) {
                best_union = r.union_T;
            }
        }
        // Greedy: start from identity, flip each bit to complement if it helps
        std::vector<int> config(n, 0);
        bool improved = true;
        while (improved && generated < max_candidates) {
            improved = false;
            for (int i = 0; i < n && generated < max_candidates; i++) {
                if (config[i]) continue; // already flipped
                config[i] = 1; // try flip
                std::vector<int> sel;
                for (int j = 0; j < n; j++)
                    sel.push_back(config[j] ? comp_idx[j] : id_idx[j]);
                uint32_t M[32] = {0};
                uint64_t bv = 0;
                for (int j = 0; j < n; j++) {
                    M[j] = pool_rows[sel[j]].first;
                    if (pool_rows[sel[j]].second) bv |= (1ULL << j);
                }
                auto cand = evaluate_pool_rowset(tt, M, bv, n, n, n_threads);
                if (cand.union_T < best_union) {
                    best_union = cand.union_T;
                    best_config = config;
                    improved = true;
                    results.push_back(cand);
                    generated++;
                    std::cout << "      pool greedy-flip bit " << i
                              << " T_union=" << cand.union_T << "\n";
                } else {
                    config[i] = 0; // revert
                }
            }
        }
    }

    // 6) Try adding XOR-pair rows on top of identity (m > n)
    if (generated < max_candidates) {
        int n_xor_pairs = std::min(12, (int)pool_rows.size() - 2 * n);
        for (int xi = 0; xi < n_xor_pairs && generated < max_candidates; xi++) {
            int pool_idx = 2 * n + xi; // XOR-pair rows start at 2n
            std::vector<int> sel;
            for (int j = 0; j < n; j++) sel.push_back(id_idx[j]);
            sel.push_back(pool_idx); // add XOR-pair row
            try_Mb(sel, "add-xor");
        }
    }

    std::cout << "    Pool candidates: " << generated << "\n";
}

// ============================================================
//  Hill climbing
// ============================================================

// ============================================================
//  Hill climbing
// ============================================================

static MbCandidate evaluate_with_pair_filter(
    const TruthTable& tt,
    const uint32_t* M_rows, uint64_t b, int m, int n_threads,
    const uint64_t* pair_masks, int n_pairs)
{
    MbCandidate cand = evaluate_Mb(tt, M_rows, b, m, n_threads, true);
    if (cand.union_T >= INT64_MAX) return cand;

    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    std::vector<uint64_t> union_anf(n_words_g, 0);
    int64_t total_T = 0;

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        auto& raw = cand.g_tt_raw[oi];
        if (raw.empty()) { cand.per_output_T[oi] = 0; continue; }
        moebius_packed(raw.data(), m);
        filter_pairs(raw.data(), m, pair_masks, n_pairs);
        int64_t T_filt = count_T(raw.data(), m);
        cand.per_output_T[oi] = T_filt;
        total_T += T_filt;
        for (int64_t w = 0; w < n_words_g; w++)
            union_anf[w] |= raw[w];
        moebius_packed(raw.data(), m);
    }
    cand.total_T = total_T;
    cand.union_T = count_T(union_anf.data(), m);
    return cand;
}

// Hill climb from a complement candidate: modify non-paired identity rows
// but keep complement rows and their paired identity rows fixed.
static MbCandidate hill_climb_complement(
    const TruthTable& tt, const MbCandidate& start,
    int n, int n_threads,
    const uint64_t* pair_masks, int n_pairs)
{
    MbCandidate best = start;
    int m = start.m;
    int k = m - n;

    // Identify which identity rows are paired with complement rows (must stay fixed)
    bool row_locked[32] = {false};
    for (int j = 0; j < k; j++) {
        // Extract the complement bit from pair_masks[j]:
        // pair_masks[j] = (1<<comp_bit) | (1<<(n+j))
        int comp_bit = __builtin_ctzll(pair_masks[j]);
        row_locked[comp_bit] = true;
    }

    auto is_better = [](const MbCandidate& a, const MbCandidate& b) -> bool {
        if (a.union_T != b.union_T) return a.union_T < b.union_T;
        return a.total_T < b.total_T;
    };

    MbCandidate cur = best;
    bool improved;
    int pass = 0;
    const int MAX_PASSES = 50;

    do {
        improved = false;
        pass++;

        // Flip bits in non-locked identity rows only
        for (int r = 0; r < n && !improved; r++) {
            if (row_locked[r]) continue;
            for (int c = 0; c < n && !improved; c++) {
                cur.M_rows[r] ^= (1u << c);
                MbCandidate cand = evaluate_with_pair_filter(
                    tt, cur.M_rows.data(), cur.b, m, n_threads, pair_masks, n_pairs);
                if (cand.union_T != INT64_MAX && is_better(cand, best))
                    best = cand;
                if (cand.union_T != INT64_MAX && is_better(cand, cur)) {
                    cur = cand;
                    improved = true;
                } else {
                    cur.M_rows[r] ^= (1u << c);
                }
            }
        }

        // Flip b-bit in non-locked rows only
        for (int r = 0; r < n && !improved; r++) {
            if (row_locked[r]) continue;
            cur.b ^= (1ULL << r);
            MbCandidate cand = evaluate_with_pair_filter(
                tt, cur.M_rows.data(), cur.b, m, n_threads, pair_masks, n_pairs);
            if (cand.union_T != INT64_MAX && is_better(cand, best))
                best = cand;
            if (cand.union_T != INT64_MAX && is_better(cand, cur)) {
                cur = cand;
                improved = true;
            } else {
                cur.b ^= (1ULL << r);
            }
        }
    } while (improved && pass < MAX_PASSES);

    return best;
}

// Old hill_climb_complement kept for reference, redirect to new version
// (The new version above replaces the old one.)

MbCandidate hill_climb(const TruthTable& tt, const MbCandidate& start,
                       int n, int n_threads)
{
    MbCandidate best = start;
    int m = start.m;

    auto eval = [&](uint32_t* M, uint64_t b) -> MbCandidate {
        return evaluate_Mb(tt, M, b, m, n_threads);
    };

    // Compare: union_T primary, total_T tiebreaker
    auto is_better = [](const MbCandidate& a, const MbCandidate& b) -> bool {
        if (a.union_T != b.union_T) return a.union_T < b.union_T;
        return a.total_T < b.total_T;
    };

    MbCandidate cur = best;
    bool improved;
    int pass = 0;
    const int MAX_PASSES = 50;

    do {
        improved = false;
        pass++;

        for (int r = 0; r < m && !improved; r++) {
            for (int c = 0; c < n && !improved; c++) {
                cur.M_rows[r] ^= (1u << c);
                MbCandidate cand = eval(cur.M_rows.data(), cur.b);
                if (cand.union_T != INT64_MAX && is_better(cand, best))
                    best = cand;
                if (cand.union_T != INT64_MAX && is_better(cand, cur)) {
                    cur = cand;
                    improved = true;
                } else {
                    cur.M_rows[r] ^= (1u << c);
                }
            }
        }

        for (int r = 0; r < m && !improved; r++) {
            cur.b ^= (1ULL << r);
            MbCandidate cand = eval(cur.M_rows.data(), cur.b);
            if (cand.union_T != INT64_MAX && is_better(cand, best))
                best = cand;
            if (cand.union_T != INT64_MAX && is_better(cand, cur)) {
                cur = cand;
                improved = true;
            } else {
                cur.b ^= (1ULL << r);
            }
        }
    } while (improved && pass < MAX_PASSES);

    return best;
}

// ============================================================
// ============================================================
//  d1c: Product-based ANF simplification via exhaustive b-search
// ============================================================

int64_t count_T_nl_packed(const uint64_t* anf, int n) {
    int64_t words = int64_t(1) << (n - 6);
    uint64_t deg1_mask = 0;
    for (int b = 0; b < 64; b++)
        if (__builtin_popcountll(b) == 1) deg1_mask |= (1ULL << b);
    int64_t total = 0;
    for (int64_t w = 0; w < words; w++) {
        uint64_t val = anf[w];
        if (val == 0) continue;
        int pw = __builtin_popcountll(w);
        int cnt = __builtin_popcountll(val);
        int deg1 = 0;
        if (pw == 0)
            deg1 += __builtin_popcountll(val & deg1_mask);
        else if (pw == 1)
            if (val & 1) deg1++;
        total += (cnt - deg1);
    }
    if (anf[0] & 1) total--;
    return total;
}

// Apply upward Möbius pass for ONE bit i on packed ANF.
// For all indices with bit i = 0: data[idx] ^= data[idx | (1<<i)]
void upward_pass_bit(uint64_t* data, int n, int i) {
    int64_t words = int64_t(1) << (n - 6);
    if (i >= 6) {
        int64_t step_words = int64_t(1) << (i - 6);
        for (int64_t w = 0; w + step_words < words; w++)
            if (!((w / step_words) & 1))
                data[w] ^= data[w + step_words];
    } else {
        int b = 1 << i;
        uint64_t mask = (1ULL << b) - 1;
        for (int64_t w = 0; w < words; w++) {
            uint64_t v = data[w];
            uint64_t res = 0;
            for (int g = 0; g < 64; g += 2 * b) {
                uint64_t lo = (v >> g) & mask;
                uint64_t hi = (v >> (g + b)) & mask;
                lo ^= hi;
                res |= lo << g;
                res |= hi << (g + b);
            }
            data[w] = res;
        }
    }
}

// ============================================================
//  d1c pair-substitution: extend building blocks from xj+1 to xj+xk+1
// ============================================================

// Apply substitution xi ← xi⊕xj to ANF in-place.
// For each monomial containing i:
//   if also contains j: keep m, also add m\{i}   (cancellation by XOR)
//   if doesn't contain j: keep m, also add m\{i}∪{j}  (variable replacement)
static void substitute_pair_anf(uint64_t* anf, int n, int i, int j) {
    int64_t total = (n < 6) ? (1 << n) : (int64_t(1) << n);
    for (int64_t m = 0; m < total; m++) {
        if (!(m & (1ULL << i))) continue;
        int64_t w = m >> 6;
        uint64_t bit = 1ULL << (m & 63);
        if (!(anf[w] & bit)) continue;
        int64_t m2 = (m & (1ULL << j)) ? (m ^ (1ULL << i))
                                       : (m ^ (1ULL << i) ^ (1ULL << j));
        anf[m2 >> 6] ^= (1ULL << (m2 & 63));
    }
}

// Compute support mask from packed ANF (shared helper)
static uint64_t compute_support_from_anf(const uint64_t* anf, int n) {
    int64_t n_words = (n < 6) ? 1 : (int64_t(1) << (n - 6));
    uint64_t support = 0, lower = 0;
    for (int64_t w = 0; w < n_words; w++) {
        uint64_t val = anf[w];
        if (val == 0) continue;
        support |= ((uint64_t)w << 6);
        for (uint64_t v = val; v; v &= v - 1)
            lower |= (uint64_t)__builtin_ctzll(v);
    }
    support |= lower;
    if (n < 64) support &= (1ULL << n) - 1;
    return support;
}

// ============================================================
//  Exhaustive search for best b (existing, unchanged)
// ============================================================
// g(z) = f(z⊕b).  Tries all 2^t b-vectors via selective upward Möbius on f's ANF.
// Returns T_nl(g_best), sets best_b.  Returns -1 if t too large.
int64_t exhaustive_search_best_b(const uint64_t* f_anf, int n,
                                        uint64_t support_mask,
                                        int max_exhaustive_t,
                                        uint64_t& best_b)
{
    // Extract support variable indices
    std::vector<int> svars;
    for (int i = 0; i < n; i++)
        if (support_mask & (1ULL << i)) svars.push_back(i);
    int t = (int)svars.size();
    if (t > max_exhaustive_t) return -1;

    int64_t words = int64_t(1) << (n - 6);
    std::vector<uint64_t> buf(words);
    uint64_t n_masks = 1ULL << t;

    int64_t best_T = count_T_nl_packed(f_anf, n);
    best_b = 0;

    for (uint64_t mask = 0; mask < n_masks; mask++) {
        uint64_t b = 0;
        for (int j = 0; j < t; j++)
            if (mask & (1ULL << j)) b |= (1ULL << svars[j]);

        // Copy ANF and apply selective upward Möbius for bits where b_i=1
        memcpy(buf.data(), f_anf, words * sizeof(uint64_t));
        for (int i = 0; i < n; i++)
            if (b & (1ULL << i)) upward_pass_bit(buf.data(), n, i);

        int64_t T = count_T_nl_packed(buf.data(), n);
        if (T < best_T) { best_T = T; best_b = b; }
    }
    return best_T;
}

// Compute support mask + T_nl for each output (handles tt in-place: Möbius→analyze→restore)
void compute_output_info(TruthTable& tt, int n, std::vector<OutputInfo>& info) {
    int n_out = tt.n_outputs;
    info.resize(n_out);
    int64_t words = int64_t(1) << (n - 6);

    for (int oi = 0; oi < n_out; oi++) {
        moebius_packed(tt.tt[oi].data(), n);
        uint64_t* anf = tt.tt[oi].data();

        // Support set
        uint64_t support = 0, lower = 0;
        for (int64_t w = 0; w < words; w++) {
            uint64_t val = anf[w];
            if (val == 0) continue;
            support |= ((uint64_t)w << 6);
            if (val == UINT64_MAX) { lower |= 0x3F; }
            else { for (uint64_t v = val; v; v &= v - 1)
                       lower |= (uint64_t)__builtin_ctzll(v); }
        }
        support |= lower;
        if (n < 64) support &= (1ULL << n) - 1;

        info[oi].support_mask = support;
        info[oi].t = __builtin_popcountll(support);
        info[oi].T_nl = count_T_nl_packed(anf, n);

        moebius_packed(anf, n);  // restore truth table
    }
}

// ============================================================
//  Signal handling — graceful interrupt
// ============================================================

static volatile sig_atomic_t g_interrupted = 0;

static void on_signal(int) { g_interrupted = 1; }

// ============================================================
//  run_search — full pipeline (Phase 1..6)
// ============================================================

void run_search(const TruthTable& tt, const Circuit& circ,
                const std::vector<int>& output_indices,
                const SearchParams& params)
{
    // Register signal handlers for graceful Ctrl+C handling
    std::signal(SIGINT, on_signal);
    std::signal(SIGTERM, on_signal);

    int n = tt.n;
    auto t0 = std::chrono::steady_clock::now();
    std::vector<MbCandidate> results;

    // Phase 2: Raw ANF baseline
    std::cout << "\nPhase 2: Raw ANF baseline...\n";
    TruthTable tt_copy = tt;  // mutable copy for ANF transform
    RawANFInfo raw = compute_raw_anf_info(tt_copy);
    std::cout << "  Sum T = " << raw.sum_T << "\n";
    std::cout << "  Max deg = " << raw.overall_max_deg << "\n";
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        std::cout << "    " << circ.outputs[output_indices[oi]]
                  << ": T=" << raw.per_output_T[oi]
                  << ", m=" << raw.per_output_deg[oi] << "\n";
    }

    // Save raw ANF if requested
    if (!params.anf_out_dir.empty()) {
        namespace fs = std::filesystem;
        fs::create_directories(params.anf_out_dir);
        std::string base = params.anf_out_dir + "/" + params.inst_name + "_raw";
        save_raw_anf(tt_copy, circ, output_indices, base + ".poly");
        save_raw_T(tt_copy, circ, output_indices, base + "_stats.txt");
    }

    // Phase 2c: d1c — product-based ANF simplification via exhaustive b-search
    // For n ≤ 20 and t ≤ max_exhaustive_t, try all 2^t b-vectors and pick the best.
    {
        auto t_d1c0 = std::chrono::steady_clock::now();
        int max_exhaustive_t = 20;
        bool do_d1c = (n <= 20);
        std::cout << "\nPhase 2c: d1c product-based simplification"
                  << (do_d1c ? "" : " (n>20, skipped)") << "...\n";

        if (do_d1c) {
            // Get ANF + support info for all outputs (handles Moebius→restore)
            std::vector<OutputInfo> out_info;
            compute_output_info(tt_copy, n, out_info);

            // For each output with small t, run exhaustive b-search
            int n_searched = 0, n_improved = 0;
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                auto& oi_i = out_info[oi];
                if (oi_i.t > max_exhaustive_t) continue;
                n_searched++;

                // Apply Moebius to get ANF for exhaustive search
                moebius_packed(tt_copy.tt[oi].data(), n);
                uint64_t* anf = tt_copy.tt[oi].data();

                uint64_t best_b = 0;
                int64_t best_T = exhaustive_search_best_b(anf, n,
                    oi_i.support_mask, max_exhaustive_t, best_b);

                moebius_packed(anf, n);  // restore truth table

                if (best_T >= 0 && best_T < oi_i.T_nl) {
                    n_improved++;
                    std::cout << "    [d1c] output " << oi
                              << " (t=" << oi_i.t << "): T_nl "
                              << oi_i.T_nl << " → " << best_T
                              << " (b=0x" << std::hex << best_b << std::dec << ")\n";

                    // Build M=I, b=best_b candidate and evaluate on ALL outputs
                    uint32_t M[32] = {0};
                    for (int i = 0; i < n && i < 32; i++)
                        M[i] = (1u << i);

                    auto cand = evaluate_Mb(tt_copy, M, best_b, n, params.n_threads);
                    if (cand.union_T < INT64_MAX) {
                        results.push_back(cand);
                        std::cout << "      → T_union=" << cand.union_T << "\n";
                    }
                }
            }

            // Pair-substitution extension: for n≤16, also try xi ⊕ xj building blocks
            static const int MAX_PAIR_D1C_N = 16;
            static const int MAX_PAIR_TRIALS = 200;
            if (n <= MAX_PAIR_D1C_N) {
                int64_t n_words = (n < 6) ? 1 : (int64_t(1) << (n - 6));
                for (int oi = 0; oi < tt.n_outputs; oi++) {
                    auto& oi_i = out_info[oi];
                    if (oi_i.t == 0 || oi_i.t > max_exhaustive_t) continue;

                    moebius_packed(tt_copy.tt[oi].data(), n);
                    uint64_t* orig_anf = tt_copy.tt[oi].data();
                    int64_t baseline_T = oi_i.T_nl;

                    int n_pair_trials = 0;
                    for (int si_idx = 0; si_idx < n && n_pair_trials < MAX_PAIR_TRIALS; si_idx++) {
                        if (!(oi_i.support_mask & (1ULL << si_idx))) continue;
                        for (int sj = 0; sj < n && n_pair_trials < MAX_PAIR_TRIALS; sj++) {
                            if (sj == si_idx) continue;
                            n_pair_trials++;

                            // Copy ANF and apply pair substitution
                            std::vector<uint64_t> buf(n_words);
                            memcpy(buf.data(), orig_anf, n_words * sizeof(uint64_t));
                            substitute_pair_anf(buf.data(), n, si_idx, sj);

                            // Recompute support for pair-substituted ANF
                            uint64_t pair_support = compute_support_from_anf(buf.data(), n);
                            int pair_t = __builtin_popcountll(pair_support);
                            if (pair_t == 0 || pair_t > max_exhaustive_t) continue;

                            uint64_t pair_b = 0;
                            int64_t pair_T = exhaustive_search_best_b(
                                buf.data(), n, pair_support, max_exhaustive_t, pair_b);

                            if (pair_T >= 0 && pair_T < baseline_T) {
                                baseline_T = pair_T;
                                uint32_t M[32] = {0};
                                for (int k = 0; k < n && k < 32; k++)
                                    M[k] = (1u << k);
                                M[si_idx] = (1u << si_idx) | (1u << sj);

                                auto cand = evaluate_Mb(tt_copy, M, pair_b, n, params.n_threads);
                                if (cand.union_T < INT64_MAX) {
                                    results.push_back(cand);
                                    std::cout << "    [d1c-pair] output " << oi
                                              << ": x" << si_idx << "←x" << si_idx << "⊕x" << sj
                                              << " T_union=" << cand.union_T << "\n";
                                }
                            }
                        }
                    }
                    moebius_packed(orig_anf, n);  // restore truth table
                }
            }

            if (n_searched == 0)
                std::cout << "    (all outputs have t > " << max_exhaustive_t << ")\n";
            else
                std::cout << "    Searched " << n_searched << " output(s), "
                          << n_improved << " improved\n";
        }
        auto t_d1c1 = std::chrono::steady_clock::now();
        std::cout << "    d1c time: "
                  << std::chrono::duration<double>(t_d1c1 - t_d1c0).count() << " s\n";
    }

    // Phase 2b: Theory-guided analysis (Aut(F) basis + dependency set)
    std::vector<uint64_t> aut_basis;
    uint32_t dep_union = (1u << n) - 1; // default: all bits are relevant
    if (n <= 20) {
        auto t_aut0 = std::chrono::steady_clock::now();
        std::cout << "\nPhase 2b: Computing Aut(F) basis...\n";
        aut_basis = compute_aut_basis(tt, n, params.n_threads);
        std::cout << "  dim(Aut(F)) = " << aut_basis.size() << "\n";
        if (aut_basis.empty()) {
            std::cout << "  (no non-trivial invariance)\n";
        } else {
            std::cout << "  basis vectors:";
            for (auto& v : aut_basis) std::cout << " " << v;
            std::cout << "\n";
        }
        auto t_aut1 = std::chrono::steady_clock::now();
        std::cout << "  Aut(F) time: " << std::chrono::duration<double>(t_aut1 - t_aut0).count() << " s\n";

        // Dependency set (exact for n ≤ 20)
        if (params.use_dep_filter) {
            auto dep_sets = compute_dependency_set(tt, n, 500);
            dep_union = 0;
            for (int oi = 0; oi < tt.n_outputs; oi++)
                dep_union |= dep_sets[oi];
            std::cout << "  Dependency union: " << __builtin_popcount(dep_union)
                      << "/" << n << " bits ("
                      << (n - __builtin_popcount(dep_union)) << " don't-care)\n";
        }
    }

    auto t1 = std::chrono::steady_clock::now();

    // Phase 3: Walsh correlations (n ≤ 32)
    std::vector<WalshInfo> walsh;
    if (n <= 32) {
        std::cout << "\nPhase 3: Computing Walsh correlations...\n";
        walsh = compute_walsh_correlations(tt_copy, params.n_threads);
        for (int oi = 0; oi < tt.n_outputs && oi < 8; oi++) {
            std::cout << "  Output " << circ.outputs[output_indices[oi]] << " top-5 bits:";
            // Find top 5 magnitudes
            struct { int idx; int64_t mag; } sorted[32];
            for (int i = 0; i < n; i++) sorted[i] = {i, walsh[oi].walsh_mag[i]};
            std::sort(sorted, sorted + n,
                [](auto& a, auto& b) { return a.mag > b.mag; });
            for (int j = 0; j < std::min(5, n); j++)
                std::cout << " " << sorted[j].idx << "(" << sorted[j].mag << ")";
            std::cout << "\n";
        }
        auto t_walsh = std::chrono::steady_clock::now();
        std::cout << "  Walsh time: " << std::chrono::duration<double>(t_walsh - t1).count() << " s\n";
    } else {
        std::cout << "\nPhase 3: Skipping Walsh correlations (n=" << n << " > 32)\n";
    }

    // Phase 4: Generate and evaluate candidates
    std::cout << "\nPhase 4: Searching for M,b...\n";

#ifndef DIRECTION_TAG
#define DIRECTION_TAG "d1a"
#endif
#ifndef STRATEGY_TAG
#define STRATEGY_TAG "opt1"
#endif

    auto elapsed = [&]() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    };

    // Save best candidate found so far (for time-budget graceful exit)
    auto save_best_on_timeout = [&]() {
        if (params.results_dir.empty() || results.empty()) return;
        // Re-sort: union_T primary, total_T tiebreaker
        std::sort(results.begin(), results.end(),
            [](auto& a, auto& b) {
                if (a.union_T != b.union_T) return a.union_T < b.union_T;
                return a.total_T < b.total_T;
            });
        auto& best = results[0];
        if (best.union_T >= INT64_MAX) return;

        std::cout << "\n  Saving current best (T_union=" << best.union_T << ")...\n";
        namespace fs = std::filesystem;
        fs::create_directories(params.results_dir);

        std::string base = params.results_dir + "/" + params.inst_name
                         + "_" DIRECTION_TAG "_" STRATEGY_TAG;

        // Re-evaluate with save_g_tt=true to get g_tt_raw for .poly
        MbCandidate verified;
        bool can_save = false;
        if (best.m <= 32) {
            verified = evaluate_Mb(tt_copy, best.M_rows.data(), best.b,
                                   best.m, params.n_threads, true);
            can_save = (verified.union_T < INT64_MAX);

            if (can_save && best.m > n) {
                uint64_t cp_masks[16];
                int n_cp = detect_complement_pairs(
                    best.M_rows.data(), best.b, best.m, cp_masks, 16);
                if (n_cp > 0) {
                    for (int oi = 0; oi < tt.n_outputs; oi++) {
                        auto& raw = verified.g_tt_raw[oi];
                        if (raw.empty()) continue;
                        moebius_packed(raw.data(), best.m);
                        filter_pairs(raw.data(), best.m, cp_masks, n_cp);
                        moebius_packed(raw.data(), best.m);
                    }
                }
            }
        }

        if (can_save) {
            save_trans(best, circ, n, base + ".affine");
            save_opt_expr(verified.g_tt_raw, circ, output_indices, best.m, base + ".poly");
            save_opt_T(verified.g_tt_raw, circ, output_indices, best.m, base + "_stats.txt");
            std::cout << "  Saved: " << base << ".*\n";
        } else {
            save_trans(best, circ, n, base + ".affine");
            std::cout << "  Saved: " << base << ".affine (transform only)\n";
        }
    };

    auto time_budget_exhausted = [&]() -> bool {
        if (g_interrupted) return true;
        if (params.time_budget <= 0) return false;
        return elapsed() >= params.time_budget;
    };

    auto has_solution = [&]() -> bool {
        return !results.empty() && results[0].union_T < INT64_MAX;
    };

    // 4a-c: Walsh-guided + random (n ≤ 20 full search; 21-32 permutation only)
    if (n <= 20) {
        CandidateGenerator gen(n, params.max_m, walsh);

        // 4a-pr: Progressive M construction (theory-guided)
        if (params.use_progressive) {
            int prog_max_m = (params.progressive_max_m > 0) ?
                params.progressive_max_m : std::min(n + 8, 32);
            std::cout << "  4a-pr: Progressive M construction (max m=" << prog_max_m << ")...\n";
            auto t_prog0 = std::chrono::steady_clock::now();
            auto prog_cand = progressive_m_search(tt, aut_basis, n, prog_max_m, params.n_threads);
            if (prog_cand.union_T < INT64_MAX) {
                results.push_back(prog_cand);
                std::cout << "    progressive result: m=" << prog_cand.m
                          << " T_union=" << prog_cand.union_T << "\n";
            }
            auto t_prog1 = std::chrono::steady_clock::now();
            std::cout << "    progressive time: "
                      << std::chrono::duration<double>(t_prog1 - t_prog0).count() << " s\n";
        }

        std::cout << "  4a: Single-row Walsh (" << (params.walsh_single_top * 2) << " candidates)\n";
        auto single_rows = gen.gen_walsh_single_rows(params.walsh_single_top);
        for (auto& [row, b] : single_rows) {
            uint32_t M[32] = {0};
            M[0] = row;
            results.push_back(evaluate_Mb(tt_copy, M, b, 1, params.n_threads));
        }

        std::cout << "  4b: Multi-row Walsh (" << params.multi_max_rows
                  << " progressive + pair combinations)\n";
        auto multi_rows = gen.gen_multi_row(params.multi_max_rows);
        for (auto& [rows, b] : multi_rows) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            results.push_back(evaluate_Mb(tt_copy, M, b, m, params.n_threads));
        }

        std::cout << "  4c: Random M,b (" << params.n_random << " candidates)\n";
        auto random_cands = gen.gen_random(params.n_random, params.max_m);
        int rc_idx = 0;
        for (auto& [rows, b] : random_cands) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            results.push_back(evaluate_Mb(tt_copy, M, b, m, params.n_threads));
            if ((++rc_idx) % 1000 == 0 && time_budget_exhausted()) {
                save_best_on_timeout();
                std::cout << "# STATUS: " << (has_solution() ? "has_solution" : "no_solution") << " timeout\n";
                return;
            }
        }

        auto t_search0 = std::chrono::steady_clock::now();
        double walsh_time = std::chrono::duration<double>(t_search0 - t1).count();
        std::cout << "  Search time: " << walsh_time << " s\n";
    } else if (n <= 32 && !walsh.empty()) {
        CandidateGenerator gen(n, params.max_m, walsh);

        // For n>20, the non-bijective evaluate_Mb (m < n) is extremely slow because
        // it iterates the full 2^n truth table per output. Instead, pad Walsh-guided
        // candidates to full-rank n×n matrices and use fast evaluate_Mb_bijective.
        auto eval_as_bijective = [&](uint32_t* M_in, int m_in, uint64_t b_in) -> MbCandidate {
            uint32_t M[32] = {0};
            for (int j = 0; j < m_in; j++) M[j] = M_in[j];

            // Pad to n rows with standard basis vectors
            int r = m_in;
            for (int i = 0; i < n && r < n; i++) {
                M[r] = (1u << i);
                int rank_before = gf2_rank(M, r, n);
                int rank_after = gf2_rank(M, r + 1, n);
                if (rank_after > rank_before) {
                    r++;
                } else {
                    M[r] = 0;
                }
            }

            if (r < n) {
                MbCandidate bad;
                bad.total_T = INT64_MAX;
                bad.union_T = INT64_MAX;
                return bad;
            }
            return evaluate_Mb_bijective(tt_copy, M, b_in, n, n, params.n_threads, false);
        };

        std::cout << "  4a: Single-row Walsh -> m=n permutations ("
                  << (params.walsh_single_top * 2) << " candidates)\n";
        auto single_rows = gen.gen_walsh_single_rows(params.walsh_single_top);
        for (auto& [row, b] : single_rows) {
            uint32_t M[1] = {row};
            results.push_back(eval_as_bijective(M, 1, b));
        }

        std::cout << "  4b: Multi-row Walsh -> m=n matrices"
                  << " (progressive m=1.." << std::min(20, n) << " + XOR pairs)\n";
        auto multi_rows = gen.gen_multi_row(std::min(20, n));
        for (auto& [rows, b] : multi_rows) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            results.push_back(eval_as_bijective(M, m, b));
        }

        std::cout << "  (skipping Phase 4c small-m random for n=" << n << ")\n";
    }

    // 4d: n32 random m=n candidates (for n > 20)
    if (params.n32_random > 0 && n > 20) {
        std::cout << "  4d: n32 random m=n candidates (" << params.n32_random << ")\n";
        std::mt19937_64 rng_n32(42);
        int generated = 0, attempts = 0;
        const int MAX_ATTEMPTS = params.n32_random * 20;

        auto t_n32_0 = std::chrono::steady_clock::now();
        while (generated < params.n32_random && attempts < MAX_ATTEMPTS) {
            attempts++;
            uint32_t M[32] = {0};
            for (int r = 0; r < n; r++)
                for (int i = 0; i < n; i++)
                    if (rng_n32() & 1) M[r] |= (1u << i);

            uint64_t b = 0;
            for (int r = 0; r < n; r++)
                if (rng_n32() & 1) b |= (1ULL << r);

            if (gf2_rank(M, n, n) < n) continue;

            auto t_eval0 = std::chrono::steady_clock::now();
            auto cand = evaluate_Mb(tt_copy, M, b, n, params.n_threads);
            auto t_eval1 = std::chrono::steady_clock::now();
            results.push_back(cand);
            generated++;

            double eval_time = std::chrono::duration<double>(t_eval1 - t_eval0).count();
            double compression = (double)raw.sum_T / std::max(int64_t(1), cand.union_T);
            std::cout << "    n32 #" << generated << "/" << params.n32_random
                      << " (attempt " << attempts << ")"
                      << ": T=" << cand.union_T
                      << " (compression " << compression << "×)"
                      << " time=" << eval_time << "s\n";

            if (time_budget_exhausted()) {
                save_best_on_timeout();
                std::cout << "# STATUS: " << (has_solution() ? "has_solution" : "no_solution") << " timeout\n";
                return;
            }
        }
        auto t_n32_1 = std::chrono::steady_clock::now();
        std::cout << "    n32 random time: " << std::chrono::duration<double>(t_n32_1 - t_n32_0).count() << " s\n";
    }

    // 4f: Complement variable candidates (n ≤ 16)
    if (params.n_complement > 0 && n <= 16) {
        std::cout << "  4f: Complement candidates (" << params.n_complement << " max)\n";
        add_row_pool_candidates(tt_copy, walsh, n, params.n_threads,
                                 params.n_complement, results, dep_union);
    }

    // 4e: Hill climbing from top candidates
    if (params.n_hill_climb > 0) {
        std::sort(results.begin(), results.end(),
            [](auto& a, auto& b) {
                if (a.union_T != b.union_T) return a.union_T < b.union_T;
                return a.total_T < b.total_T;
            });

        int64_t global_best_union = results.empty() ? INT64_MAX : results[0].union_T;
        int64_t global_best_total = results.empty() ? INT64_MAX : results[0].total_T;

        auto is_global_improvement = [&](const MbCandidate& c) -> bool {
            if (c.union_T < global_best_union) return true;
            if (c.union_T == global_best_union && c.total_T < global_best_total) return true;
            return false;
        };

        int n_climb = std::min(params.n_hill_climb, (int)results.size());
        std::cout << "  4e: Hill climbing from top " << n_climb << " candidates (+ identity)\n";

        auto t_hc0 = std::chrono::steady_clock::now();
        for (int ci = 0; ci < n_climb; ci++) {
            if (time_budget_exhausted()) {
                std::cout << "    Hill climb interrupted by time budget\n";
                break;
            }
            auto& base = results[ci];
            if (base.union_T >= INT64_MAX) continue;

            if (base.m == n) {
                // Standard m=n hill climb
                auto improved = hill_climb(tt_copy, base, n, params.n_threads);
                if (improved.union_T < base.union_T) {
                    results.push_back(improved);
                    if (is_global_improvement(improved)) {
                        global_best_union = improved.union_T;
                        global_best_total = improved.total_T;
                        std::cout << "    hill climb #" << ci << ": T_union=" << base.union_T
                                  << " -> " << improved.union_T
                                  << " (sum=" << improved.total_T << ")\n";
                    }
                }
            } else if (base.m > n && base.m <= 32) {
                // Complement m>n hill climb: detect complement pairs dynamically
                uint64_t hc_pair_masks[16];
                int n_hc_pairs = detect_complement_pairs(
                    base.M_rows.data(), base.b, base.m, hc_pair_masks, 16);

                if (n_hc_pairs > 0) {
                    auto improved = hill_climb_complement(tt_copy, base, n, params.n_threads,
                                                           hc_pair_masks, n_hc_pairs);
                    if (improved.union_T < base.union_T) {
                        results.push_back(improved);
                        if (is_global_improvement(improved)) {
                            global_best_union = improved.union_T;
                            global_best_total = improved.total_T;
                            std::cout << "    hill climb (complement) #" << ci
                                      << ": T_union=" << base.union_T
                                      << " -> " << improved.union_T
                                      << " (sum=" << improved.total_T << ")\n";
                        }
                    }
                }
            }

            if (time_budget_exhausted()) {
                save_best_on_timeout();
                std::cout << "# STATUS: " << (has_solution() ? "has_solution" : "no_solution") << " timeout\n";
                return;
            }
        }

        // Hill climb from identity (skip if budget exhausted)
        if (!time_budget_exhausted()) {
            uint32_t ident[32] = {0};
            for (int i = 0; i < n; i++) ident[i] = (1u << i);
            MbCandidate ident_cand = evaluate_Mb(tt_copy, ident, 0, n, params.n_threads);
            auto ident_improved = hill_climb(tt_copy, ident_cand, n, params.n_threads);
            if (ident_improved.union_T < ident_cand.union_T) {
                results.push_back(ident_improved);
                if (is_global_improvement(ident_improved)) {
                    global_best_union = ident_improved.union_T;
                    global_best_total = ident_improved.total_T;
                    std::cout << "    hill climb from identity: T_union=" << ident_cand.union_T
                              << " -> " << ident_improved.union_T
                              << " (sum=" << ident_improved.total_T << ")\n";
                }
            }
        }

        auto t_hc1 = std::chrono::steady_clock::now();
        std::cout << "    Hill climb time: " << std::chrono::duration<double>(t_hc1 - t_hc0).count() << " s\n";
    }

    // Phase 5: Report
    std::sort(results.begin(), results.end(),
        [](auto& a, auto& b) {
            if (a.union_T != b.union_T) return a.union_T < b.union_T;
            return a.total_T < b.total_T;
        });

    std::cout << "\n========================================\n";
    std::cout << "  Results (" << results.size() << " candidates evaluated)\n";
    std::cout << "  Raw ANF Sum T = " << raw.sum_T << "\n";
    std::cout << "========================================\n";

    int n_report = std::min(20, (int)results.size());
    for (int j = 0; j < n_report; j++) {
        auto& r = results[j];
        double compression = (double)raw.sum_T / std::max(int64_t(1), r.union_T);
        std::cout << "  " << (j + 1) << ". m=" << r.m << " T=" << r.total_T << " (union=" << r.union_T
                  << ", compression " << compression << "×) M=[";
        for (int row = 0; row < r.m; row++) {
            if (row > 0) std::cout << "; ";
            for (int i = 0; i < std::min(n, 8); i++)
                std::cout << ((r.M_rows[row] >> i) & 1);
            std::cout << "..";
        }
        std::cout << "]\n";
        if (r.union_T < raw.sum_T && r.m <= 32) {
            std::cout << "    Per-output T:";
            for (int oi = 0; oi < tt.n_outputs; oi++)
                std::cout << " " << circ.outputs[output_indices[oi]] << "=" << r.per_output_T[oi];
            std::cout << "\n";
        }
    }

    // Phase 6: Verify best candidate, save results
    if (!results.empty() && results[0].union_T < INT64_MAX) {
        auto& best = results[0];
        std::cout << "\nBest candidate: m=" << best.m << ", T_sum=" << best.total_T << ", T_union=" << best.union_T << "\n";
    }

    double total_time = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    // Phase 6: verification + file saving
    if (!results.empty() && results[0].union_T < INT64_MAX) {
        auto& best = results[0];

        // Print transform
        std::cout << "// Transform z = Mx + b (GF(2)):\n";
        std::cout << "// m=" << best.m << ", n=" << n << "\n";
        for (int row = 0; row < best.m; row++) {
            std::cout << "// z_" << row << " = ";
            bool first = true;
            for (int i = 0; i < n; i++) {
                if ((best.M_rows[row] >> i) & 1) {
                    if (!first) std::cout << " + ";
                    std::cout << circ.inputs[i];
                    first = false;
                }
            }
            if ((best.b >> row) & 1) {
                if (!first) std::cout << " + ";
                std::cout << "1";
            }
            std::cout << "\n";
        }

        // Detect complement pairs dynamically (works with any row ordering)
        bool is_complement = false;
        uint64_t comp_pair_masks[16] = {0};
        int n_comp_pairs = 0;
        if (best.m > n && best.m <= 32) {
            n_comp_pairs = detect_complement_pairs(
                best.M_rows.data(), best.b, best.m, comp_pair_masks, 16);
            is_complement = (n_comp_pairs > 0);
        }

        // Re-evaluate with save_g_tt=true (with pair filtering for complement)
        MbCandidate verified_cand;
        verified_cand.m = best.m;
        verified_cand.b = best.b;
        verified_cand.M_rows.resize(best.m);
        for (int i = 0; i < best.m; i++) verified_cand.M_rows[i] = best.M_rows[i];
        verified_cand.union_T = best.union_T;
        verified_cand.per_output_T = best.per_output_T;

        if (best.m == n) {
            std::cout << "  Evaluating best candidate for output...\n";
            verified_cand = evaluate_Mb(tt_copy, best.M_rows.data(), best.b, best.m,
                                         params.n_threads, true);
        } else if (is_complement) {
            std::cout << "  Evaluating complement candidate for output...\n";
            verified_cand = evaluate_Mb(tt_copy, best.M_rows.data(), best.b, best.m,
                                         params.n_threads, true);
            // Clean g_tt_raw: Möbius → filter pairs → inverse Möbius
            int64_t n_words_g = (best.m < 6) ? 1 : (int64_t(1) << (best.m - 6));
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                auto& raw = verified_cand.g_tt_raw[oi];
                if (raw.empty()) continue;
                moebius_packed(raw.data(), best.m);
                filter_pairs(raw.data(), best.m, comp_pair_masks, n_comp_pairs);
                moebius_packed(raw.data(), best.m);
            }
        } else {
            std::cout << "  Evaluating best candidate (m=" << best.m << ") for output...\n";
            if (best.m <= 32) {
                verified_cand = evaluate_Mb(tt_copy, best.M_rows.data(), best.b, best.m,
                                             params.n_threads, true);
            }
        }

        // Save results to files
        if (!params.results_dir.empty()) {
            namespace fs = std::filesystem;
            fs::create_directories(params.results_dir);
#ifndef DIRECTION_TAG
#define DIRECTION_TAG "d1a"
#endif
#ifndef STRATEGY_TAG
#define STRATEGY_TAG "opt1"
#endif
            std::string base = params.results_dir + "/" + params.inst_name + "_" + DIRECTION_TAG + "_" + STRATEGY_TAG;

            save_trans(best, circ, n, base + ".affine");

            if (!verified_cand.g_tt_raw.empty()) {
                save_opt_expr(verified_cand.g_tt_raw, circ, output_indices, best.m,
                             base + ".poly");
                save_opt_T(verified_cand.g_tt_raw, circ, output_indices, best.m,
                          base + "_stats.txt");
            }

            // Verification is done separately via verify_anf tool
        }
    }

    bool final_has_solution = !results.empty() && results[0].union_T < INT64_MAX;
    std::cout << "# STATUS: " << (final_has_solution ? "has_solution" : "no_solution") << " no_timeout\n";

    total_time = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::cout << "\nTotal time: " << total_time << " s\n";
}
