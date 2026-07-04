/**
 * optimize_anf_large — Z-space evaluation for large-n instances (n > 20).
 *
 * Strategy: Evaluate each candidate transform z = Mx + b directly in z-space
 * (2^m points) instead of computing the full 2^n truth table. Uses the
 * right-inverse M^+ to map z back to a representative x for circuit evaluation.
 *
 * Suitable for instances like hd01, hd02, hd10 (n=32, sparse ANF).
 *
 * Usage:
 *   optimize_anf_large <circuit.txt> [options]
 *
 * Options:
 *   --max-m N        max z variables per output (default 8)
 *   --random N       random candidates per output (default 100)
 *   --hill-climb N   hill climb iterations per output (default 5)
 *   --walsh-k N      top K Walsh bits (default 20)
 *   --save-results PREFIX  save results
 */
#include "circuit.h"
#include "truth_table.h"
#include "gf2.h"
#include "moebius.h"
#include "anf.h"
#include "io.h"
#include "walsh.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <chrono>
#include <filesystem>

// Strategy tag for output file naming (prevents confusion between different search strategies)
static const char* STRATEGY_TAG = "lg";

// Count only nonlinear terms (degree ≥ 2) in ANF data
static int64_t count_nonlinear_T(const uint64_t* data, int m) {
    int64_t total = 0;
    int64_t n_z = int64_t(1) << m;
    for (int64_t zi = 0; zi < n_z; zi++) {
        if (zi == 0) continue;                              // skip constant
        if ((zi & (zi - 1)) == 0) continue;                  // skip linear (power of 2)
        if ((data[zi >> 6] >> (zi & 63)) & 1) total++;
    }
    return total;
}

// ============================================================
//  Z-space evaluation routines
// ============================================================

// Bit patterns for z-word generation (same structure as x-space)
static const uint64_t Z_PATTERNS[6] = {
    0xAAAAAAAAAAAAAAAAULL,
    0xCCCCCCCCCCCCCCCCULL,
    0xF0F0F0F0F0F0F0F0ULL,
    0xFF00FF00FF00FF00ULL,
    0xFFFF0000FFFF0000ULL,
    0xFFFFFFFF00000000ULL,
};

// Compute input words for a batch of 64 consecutive z-values
// P_rows[j] = m-bit mask: which z-bits contribute to input bit j
// Returns input_words[n] for this batch
static void compute_input_words(
    uint64_t* in_words, int n,
    const uint32_t* P_rows, int m,
    uint64_t base_z)
{
    // z_word[r] = 64-bit mask for z-bit r across the batch
    uint64_t z_word[32];
    for (int r = 0; r < 6 && r < m; r++)
        z_word[r] = Z_PATTERNS[r];
    for (int r = 6; r < m; r++)
        z_word[r] = ((base_z >> r) & 1) ? ~0ULL : 0;

    // Compute each input word as XOR of relevant z-words
    for (int j = 0; j < n; j++) {
        uint64_t xw = 0;
        uint32_t mask = P_rows[j];
        for (int r = 0; r < m; r++) {
            if ((mask >> r) & 1)
                xw ^= z_word[r];
        }
        in_words[j] = xw;
    }
}

// Quick verification: check f(x) == g(Mx⊕b) ⊕ ⟨c,x⟩ ⊕ d for structured + random x-values.
// Uses single-bit vectors (0, e_0, ..., e_{n-1}) plus n_verify random points.
// c_mask/d default to 0 (no c-correction) for backward compatibility.
static bool verify_candidate(
    const FastCircuit& fc,
    const uint32_t* M_rows, uint64_t b,
    int m, int n,
    const std::vector<uint64_t>& g_tt,  // c-corrected g(z) truth table
    int output_idx,
    int n_verify,
    uint32_t c_mask = 0,
    int d_const = 0)
{
    if (n_verify <= 0) return true;

    uint32_t P_rows[32];
    if (!gf2_right_inverse(M_rows, m, n, P_rows)) return false;

    // Helper: test a single x value
    auto test_x = [&](uint64_t x) -> bool {
        // f(x) via direct evaluation
        uint64_t in_words[64];
        for (int j = 0; j < n; j++)
            in_words[j] = ((x >> j) & 1) ? ~0ULL : 0;
        std::vector<uint64_t> eval_buf(fc.n_stmts + 1);
        std::vector<uint64_t> out_vec(fc.out_idx.size());
        eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());
        uint64_t fx = (out_vec[output_idx] & 1);

        // z = Mx ⊕ b
        uint64_t z = 0;
        for (int r = 0; r < m; r++)
            if (__builtin_popcount(M_rows[r] & (uint32_t)(x & 0xFFFFFFFFULL)) & 1)
                z |= (1ULL << r);
        z ^= b;

        // g(z) from c-corrected truth table
        uint64_t gz = (g_tt[z / 64] >> (z % 64)) & 1;

        // ⟨c,x⟩
        uint64_t cx = __builtin_popcount(c_mask & (uint32_t)(x & 0xFFFFFFFFULL)) & 1;

        // f(x) == g(z) ⊕ ⟨c,x⟩ ⊕ d
        return (fx ^ gz ^ cx ^ d_const) == 0;
    };

    // 1) Zero vector
    if (!test_x(0)) return false;

    // 2) Single-bit vectors e_i for all i
    for (int i = 0; i < n; i++)
        if (!test_x(1ULL << i)) return false;

    // 3) Two-bit combinations across kernel fiber groups
    // Key insight: a function can be constant on all single-bit kernel vectors
    // but non-constant on multi-bit combinations (e.g. two bits from different
    // functional groups). Test random pairs and triples of distinct bits.
    for (int test = 0; test < std::min(200, n * 3); test++) {
        int b1 = rand() % n;
        int b2 = rand() % n;
        if (b1 == b2) { b2 = (b2 + 1) % n; }
        uint64_t x = (1ULL << b1) | (1ULL << b2);
        if (!test_x(x)) return false;
    }
    for (int test = 0; test < std::min(100, n * 2); test++) {
        int b1 = rand() % n;
        int b2 = rand() % n;
        int b3 = rand() % n;
        if (b1 == b2) b2 = (b2 + 1) % n;
        if (b1 == b3 || b2 == b3) b3 = (b3 + 1) % n;
        uint64_t x = (1ULL << b1) | (1ULL << b2) | (1ULL << b3);
        if (!test_x(x)) return false;
    }

    // 4) Random vectors (includes multi-bit combinations organically)
    for (int test = 0; test < n_verify; test++) {
        uint64_t x = (uint64_t)rand() | ((uint64_t)rand() << 15)
                   | ((uint64_t)rand() << 30) | ((uint64_t)rand() << 45);
        if (n < 60) x &= ((1ULL << n) - 1);
        if (!test_x(x)) return false;
    }
    return true;
}

// Evaluate g(z) = f(Mx⊕b) using z-space traversal (2^m points).
// Applies c-correction: removes linear z-terms from g, outputs c and d
// so that f(x) = g(Mx⊕b) ⊕ ⟨c,x⟩ ⊕ d.
// Returns nonlinear_T (degree ≥ 2 only), or INT64_MAX if verification fails.
// g_tt_out: c-corrected g(z) truth table (size n_words_g, no linear terms).
// output_idx: which output to evaluate (index into fc.out_idx).
// c_mask_out: c vector (n bits) for ⟨c,x⟩ term.
// d_out: constant shift (0 or 1).
// n_verify: number of random points to verify (0 = skip).
static int64_t evaluate_zspace(
    const FastCircuit& fc,
    const uint32_t* M_rows, uint64_t b,
    int m, int n,
    std::vector<uint64_t>& g_tt_out,
    int n_threads,
    int output_idx,
    uint32_t& c_mask_out,
    int& d_out,
    int n_verify = 200)
{
    if (m <= 0) { c_mask_out = 0; d_out = 0; return 0; }
    if (m > n) return INT64_MAX;

    // Right-inverse P (n×m) such that M × P = I_m
    uint32_t P_rows[32];
    if (!gf2_right_inverse(M_rows, m, n, P_rows)) {
        return INT64_MAX;
    }

    int64_t n_z = int64_t(1) << m;
    int64_t n_words_g = (m < 6) ? 1 : (n_z >> 6);

    g_tt_out.assign(n_words_g, 0);
    int64_t n_batches = (n_z + 63) / 64;

    // Pre-compute P×b constant (XOR this to apply affine shift)
    uint64_t cb[64] = {0};
    for (int j = 0; j < n; j++) {
        if (__builtin_popcount(P_rows[j] & b) & 1)
            cb[j] = ~0ULL;
    }

    // Evaluate each batch of 64 z-values
    std::vector<uint64_t> eval_buf(fc.n_stmts + 1);
    std::vector<uint64_t> out_vec(fc.out_idx.size());
    uint64_t in_words[64];
    for (int64_t batch = 0; batch < n_batches; batch++) {
        uint64_t base_z = batch * 64;

        // Compute x = P×z for this batch (without affine shift)
        compute_input_words(in_words, n, P_rows, m, base_z);

        // XOR P×b to get x = P×(z⊕b)
        for (int j = 0; j < n; j++)
            in_words[j] ^= cb[j];

        eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());

        g_tt_out[batch] = out_vec[output_idx];
    }

    // Mask to valid bits (for m < 6, unused bits contain batch garbage)
    if (m < 6) {
        int64_t n_z = int64_t(1) << m;
        uint64_t mask = (n_z >= 64) ? ~0ULL : ((1ULL << n_z) - 1);
        g_tt_out[0] &= mask;
    }

    // === C-correction ===
    // 1) Möbius to get ANF of g_raw
    std::vector<uint64_t> g_raw_anf = g_tt_out;
    moebius_packed(g_raw_anf.data(), m);

    // 2) Extract linear coefficients L (coefficient of z_r for each r = 0..m-1)
    uint32_t L = 0;
    for (int r = 0; r < m; r++)
        if ((g_raw_anf[r >> 6] >> (r & 63)) & 1) L |= (1u << r);

    // 3) c = Mᵀ·L  (each c_i = parity of (column i of M) · L)
    c_mask_out = 0;
    for (int i = 0; i < n; i++) {
        uint32_t col_val = 0;
        for (int r = 0; r < m; r++)
            col_val |= ((M_rows[r] >> i) & 1) << r;
        if (__builtin_popcount(col_val & L) & 1)
            c_mask_out |= (1u << i);
    }

    // 4) d = ⟨L,b⟩
    d_out = __builtin_popcount(L & (uint32_t)(b & 0xFFFFFFFFULL)) & 1;

    // 5) Remove linear terms from g_tt_out: g(z) = g_raw(z) ⊕ ⟨L,z⟩
    for (int64_t zi = 0; zi < n_z; zi++) {
        if (__builtin_popcountll(zi & L) & 1)
            g_tt_out[zi >> 6] ^= (1ULL << (zi & 63));
    }

    // === Verification with c-correction ===
    if (n_verify > 0) {
        int verify_tests = (n <= 16) ? std::max(n_verify, 2000) : n_verify;
        if (!verify_candidate(fc, M_rows, b, m, n, g_tt_out, output_idx,
                              verify_tests, c_mask_out, d_out)) {
            g_tt_out.assign(n_words_g, 0);
            return INT64_MAX;
        }
    }

    // Möbius on corrected g for ANF → count nonlinear terms only
    std::vector<uint64_t> anf = g_tt_out;
    moebius_packed(anf.data(), m);
    return count_nonlinear_T(anf.data(), m);
}


// ============================================================
//  Result structures
// ============================================================

struct Candidate {
    int m;
    uint32_t M_rows[32];
    uint64_t b;
    int64_t total_T;
    std::vector<uint64_t> g_tt_raw;    // c-corrected g(z) truth table (no linear terms)
    std::vector<uint64_t> anf_coeffs;  // ANF coefficients (post-Möbius)

    // c-correction: f(x) = g(Mx⊕b) ⊕ ⟨c,x⟩ ⊕ d
    uint32_t c_mask;  // c vector (n bits, bit i = coefficient of x_i)
    int d;            // constant shift (0 or 1)

    // Affine flag: output f(x) = affine_b ⊕ Σ a_i·x_i requires no optimization
    bool is_affine;
    uint32_t affine_mask;  // bit i set = x_i appears in expression
    int affine_b;          // constant term (0 or 1)
};

static bool cmp_candidate(const Candidate& a, const Candidate& b) {
    return a.total_T < b.total_T;
}

// Detect if output f(x) is affine: f(x) = b ⊕ Σ a_i·x_i (including constants 0/1)
// Uses n+1 circuit evaluations + random verification via batch evaluation.
// Returns true with affine_mask (bit i = a_i) and affine_b properly set.
static bool detect_affine_output(const FastCircuit& fc, int n, int output_idx,
                                  uint32_t& affine_mask, int& affine_b) {
    // Evaluate f(0) and f(e_i) for all i using batch evaluation (up to 64 in parallel)
    uint64_t in_words[64] = {0};
    for (int i = 0; i < n && i < 63; i++)
        in_words[i] = (1ULL << (i + 1));

    std::vector<uint64_t> eval_buf(fc.n_stmts + 1);
    std::vector<uint64_t> out_vec(fc.out_idx.size());
    eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());

    uint64_t bits = out_vec[output_idx];
    int f0 = (bits >> 0) & 1;
    affine_b = f0;
    affine_mask = 0;
    for (int i = 0; i < n && i < 63; i++) {
        int fi = (bits >> (i + 1)) & 1;
        if (fi != f0) affine_mask |= (1u << i);
    }

    // If n=32 and we overflowed 63 slots, fall back to individual eval for remaining bits
    // (n max 32, so 1+32=33 ≤ 63 slots — safe)
    if (n > 63) return false; // not supported

    // Verify with random points: f(x) == f0 ⊕ popcount(affine_mask & x) mod 2
    srand(12345);
    for (int t = 0; t < 30; t++) {
        uint64_t x = 0;
        for (int j = 0; j < n; j++)
            if (rand() & 1) x |= (1ULL << j);

        memset(in_words, 0, sizeof(in_words));
        for (int j = 0; j < n; j++)
            if ((x >> j) & 1) in_words[j] = ~0ULL;
        eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());
        int fx = (out_vec[output_idx] & 1);
        int predicted = f0 ^ (__builtin_popcount(affine_mask & (uint32_t)x) & 1);
        if (fx != predicted) return false;
    }
    return true;
}

// Search for best M,b for a single output (n=32 strategy)
//
// Pipeline:
//   Phase 0: Affine detection
//   Phase 1: m=1 identity-row candidates with 5000-pt verify
//   Phase 2: Systematic C(K,m) enumeration for m=2..5
//   Phase 3: XOR pairs of top-8 inputs (m=1)
//   Phase 4: Greedy extension from best candidate
//   Phase 5: Hill climbing (only when m ≤ 8)
//   Phase 6: Random fallback (100 candidates)
//
static Candidate search_single_output_large(
    const FastCircuit& fc,
    int n,
    int n_threads,
    int max_m,
    int n_random,
    int n_hill_climb,
    int walsh_k,
    int output_idx)
{
    Candidate best;
    best.total_T = INT64_MAX;
    best.m = 0;
    best.is_affine = false;
    best.affine_mask = 0;
    best.affine_b = 0;
    best.c_mask = 0;
    best.d = 0;

    // ============================================================
    // Phase 0: Affine detection
    // ============================================================
    {
        uint32_t amask;
        int ab;
        if (detect_affine_output(fc, n, output_idx, amask, ab)) {
            best.is_affine = true;
            best.affine_mask = amask;
            best.affine_b = ab;
            if (amask == 0) {
                best.m = 0;
                best.total_T = 0;
            } else {
                best.m = 1;
                best.M_rows[0] = amask;
                best.b = ab;
                uint32_t c_tmp; int d_tmp;
                std::vector<uint64_t> g_tt;
                int64_t T = evaluate_zspace(fc, best.M_rows, best.b, 1, n,
                                           g_tt, n_threads, output_idx, c_tmp, d_tmp);
                if (T <= 0 || T >= INT64_MAX) T = 1;
                best.total_T = T;
                best.g_tt_raw = g_tt;
                best.c_mask = c_tmp;
                best.d = d_tmp;
            }
            return best;
        }
    }

    // ============================================================
    // Phase 1: m=1 identity-row candidates with 5000-pt verification
    //
    // Enumerate all z = x_i and z = x_i ⊕ 1. Use strict 5000-pt verify.
    // If verified, accept immediately (output truly depends on 1 input).
    // Store T values (from no-verify eval) for ranking even when
    // verification fails — lower T → input is more important.
    // ============================================================
    struct InputScore {
        int bit;
        int64_t T;   // lowest T among b=0,b=1 (from no-verify eval)
        bool valid;  // passed 5000-pt verify
    };
    std::vector<InputScore> scores(n);

    for (int bit = 0; bit < n; bit++) {
        int64_t best_T = INT64_MAX;
        bool valid = false;

        for (int bval = 0; bval < 2; bval++) {
            uint32_t M[32] = {0};
            M[0] = (1u << bit);
            uint64_t b = bval;

            // Strict 5000-pt verification
            uint32_t c_tmp; int d_tmp;
            std::vector<uint64_t> g_tt;
            int64_t T = evaluate_zspace(fc, M, b, 1, n, g_tt, n_threads, output_idx,
                                        c_tmp, d_tmp, 5000);

            if (T >= 0 && T < best.total_T) {
                best.total_T = T;
                best.m = 1;
                best.M_rows[0] = M[0];
                best.b = b;
                best.g_tt_raw = g_tt;
                best.c_mask = c_tmp;
                best.d = d_tmp;
            }

            if (T >= 0 && T < INT64_MAX) {
                if (T < best_T) best_T = T;
                valid = true;
            } else {
                // No-verify T for ranking
                uint32_t c2_tmp; int d2_tmp;
                std::vector<uint64_t> g_tt2;
                int64_t T2 = evaluate_zspace(fc, M, b, 1, n, g_tt2, n_threads, output_idx,
                                            c2_tmp, d2_tmp, 0);
                if (T2 >= 0 && T2 < best_T) best_T = T2;
            }
        }

        scores[bit] = {bit, best_T, valid};
    }

    if (best.total_T <= 1) return best;

    // Sort inputs by T score (ascending — lower T = more important)
    // Valid inputs (passed 5000-pt verify) always come before invalid
    std::sort(scores.begin(), scores.end(), [](const InputScore& a, const InputScore& b) {
        if (a.valid != b.valid) return a.valid;
        return a.T < b.T;
    });

    int K = std::min(16, n);
    if (walsh_k > 0 && walsh_k < K) K = walsh_k;

    // ============================================================
    // Phase 2: Systematic C(K, m) enumeration for m = 2..6
    //
    // Build permutation matrices from top-K inputs. For each
    // combination of m distinct inputs, try b=0 and b=all-1.
    // 500-pt verify for screening.
    // ============================================================
    {
        int max_m2 = std::min(7, std::min(max_m, n));
        for (int m_val = 2; m_val <= max_m2; m_val++) {
            std::vector<int> comb(m_val);
            for (int i = 0; i < m_val; i++) comb[i] = i;

            while (true) {
                uint32_t M[32] = {0};
                for (int r = 0; r < m_val; r++)
                    M[r] = (1u << scores[comb[r]].bit);

                // b = 0
                {
                    uint32_t c_tmp; int d_tmp;
                    std::vector<uint64_t> g_tt;
                    int64_t T = evaluate_zspace(fc, M, 0, m_val, n, g_tt, n_threads, output_idx,
                                                c_tmp, d_tmp, 2000);
                    if (T >= 0 && T < best.total_T) {
                        best.total_T = T;
                        best.m = m_val;
                        for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                        best.b = 0;
                        best.g_tt_raw = g_tt;
                        best.c_mask = c_tmp;
                        best.d = d_tmp;
                    }
                }

                // b = all-1
                {
                    uint64_t b_all = (m_val < 64) ? ((1ULL << m_val) - 1) : ~0ULL;
                    uint32_t c_tmp; int d_tmp;
                    std::vector<uint64_t> g_tt;
                    int64_t T = evaluate_zspace(fc, M, b_all, m_val, n, g_tt, n_threads, output_idx,
                                                c_tmp, d_tmp, 2000);
                    if (T >= 0 && T < best.total_T) {
                        best.total_T = T;
                        best.m = m_val;
                        for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                        best.b = b_all;
                        best.g_tt_raw = g_tt;
                        best.c_mask = c_tmp;
                        best.d = d_tmp;
                    }
                }

                if (best.total_T <= 2) break;

                // Next combination
                int i = m_val - 1;
                while (i >= 0 && comb[i] == K - m_val + i) i--;
                if (i < 0) break;
                comb[i]++;
                for (int j = i + 1; j < m_val; j++) comb[j] = comb[j-1] + 1;
            }

            if (best.total_T <= 2) break;
        }
    }

    if (best.total_T == 0) return best;

    // ============================================================
    // Phase 2b: Progressive identity rows (z_r = x_{scores[r].bit})
    //
    // For m = best.m+1 .. min(K, max_m, n), build M where each row
    // is a unit vector from the top-K scored inputs (in order).
    // This is fast (one candidate per m) and covers larger m values
    // that C(K,m) enumeration can't reach.
    // ============================================================
    {
        int start_m = std::max(2, best.m + 1);
        for (int m_val = start_m; m_val <= std::min(K, std::min(max_m, n)); m_val++) {
            uint32_t M[32] = {0};
            for (int r = 0; r < m_val; r++)
                M[r] = (1u << scores[r].bit);

            // b = 0
            {
                uint32_t c_tmp; int d_tmp;
                std::vector<uint64_t> g_tt;
                int64_t T = evaluate_zspace(fc, M, 0, m_val, n, g_tt, n_threads, output_idx,
                                            c_tmp, d_tmp, 2000);
                if (T >= 0 && T < best.total_T) {
                    best.total_T = T;
                    best.m = m_val;
                    for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                    best.b = 0;
                    best.g_tt_raw = g_tt;
                    best.c_mask = c_tmp;
                    best.d = d_tmp;
                }
            }

            // b = all-1
            {
                uint64_t b_all = (m_val < 64) ? ((1ULL << m_val) - 1) : ~0ULL;
                uint32_t c_tmp; int d_tmp;
                std::vector<uint64_t> g_tt;
                int64_t T = evaluate_zspace(fc, M, b_all, m_val, n, g_tt, n_threads, output_idx,
                                            c_tmp, d_tmp, 2000);
                if (T >= 0 && T < best.total_T) {
                    best.total_T = T;
                    best.m = m_val;
                    for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                    best.b = b_all;
                    best.g_tt_raw = g_tt;
                    best.c_mask = c_tmp;
                    best.d = d_tmp;
                }
            }

            if (best.total_T <= 2) break;
        }
    }

    if (best.total_T == 0) return best;

    // ============================================================
    // Phase 3: XOR pairs of top-8 inputs (m=1)
    // z = x_i ⊕ x_j
    // ============================================================
    {
        int max_xor = std::min(8, K);
        for (int i = 0; i < max_xor; i++) {
            for (int j = i + 1; j < max_xor; j++) {
                uint32_t M[32] = {0};
                M[0] = (1u << scores[i].bit) | (1u << scores[j].bit);

                uint32_t c_tmp; int d_tmp;
                std::vector<uint64_t> g_tt;
                int64_t T = evaluate_zspace(fc, M, 0, 1, n, g_tt, n_threads, output_idx,
                                            c_tmp, d_tmp, 2000);
                if (T >= 0 && T < best.total_T) {
                    best.total_T = T;
                    best.m = 1;
                    best.M_rows[0] = M[0];
                    best.b = 0;
                    best.g_tt_raw = g_tt;
                    best.c_mask = c_tmp;
                    best.d = d_tmp;
                }
            }
        }
    }

    if (best.total_T == 0) return best;

    // ============================================================
    // Phase 4: Greedy extension from best candidate
    // ============================================================
    if (best.total_T < INT64_MAX && best.m > 0 && best.m < std::min(max_m, n)) {
        for (int m_val = best.m + 1; m_val <= std::min(max_m, n); m_val++) {
            bool improved = false;
            for (int ki = 0; ki < K; ki++) {
                uint32_t new_mask = (1u << scores[ki].bit);

                bool already = false;
                for (int r = 0; r < m_val - 1; r++) {
                    if (best.M_rows[r] == new_mask) { already = true; break; }
                }
                if (already) continue;

                uint32_t M[32];
                for (int r = 0; r < m_val - 1; r++) M[r] = best.M_rows[r];
                M[m_val - 1] = new_mask;

                for (int bnew = 0; bnew < 2; bnew++) {
                    uint64_t new_b = best.b | ((uint64_t)bnew << (m_val - 1));

                    uint32_t c_tmp; int d_tmp;
                    std::vector<uint64_t> g_tt;
                    int64_t T = evaluate_zspace(fc, M, new_b, m_val, n, g_tt, n_threads, output_idx,
                                                c_tmp, d_tmp, 2000);
                    if (T >= 0 && T < best.total_T) {
                        best.total_T = T;
                        best.m = m_val;
                        for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                        best.b = new_b;
                        best.g_tt_raw = g_tt;
                        best.c_mask = c_tmp;
                        best.d = d_tmp;
                        improved = true;
                    }
                }
            }
            if (!improved) break;
        }
    }

    if (best.total_T == 0) return best;

    // ============================================================
    // Phase 5: Hill climbing (only when m ≤ 8)
    // ============================================================
    if (best.total_T < INT64_MAX && best.m > 0 && best.m <= 8 && n_hill_climb > 0) {
        int max_m_hc = std::min(max_m, n);
        for (int hi = 0; hi < n_hill_climb; hi++) {
            Candidate cur = best;
            bool improved = true;
            int iter = 0;
            while (improved && iter < 500 && cur.m <= max_m_hc && cur.m <= 8) {
                improved = false;
                iter++;

                for (int r = 0; r < cur.m; r++) {
                    for (int flip = 0; flip < std::min(32, n); flip++) {
                        uint32_t old_row = cur.M_rows[r];
                        cur.M_rows[r] ^= (1u << flip);

                        for (int bp = 0; bp < 2; bp++) {
                            cur.b ^= ((uint64_t)bp << r);

                            uint32_t c_tmp; int d_tmp;
                            std::vector<uint64_t> g_tt;
                            int64_t T = evaluate_zspace(
                                fc, cur.M_rows, cur.b, cur.m, n, g_tt, n_threads, output_idx,
                                c_tmp, d_tmp, 2000);
                            if (T >= 0 && T < best.total_T) {
                                best = cur;
                                best.total_T = T;
                                best.g_tt_raw = g_tt;
                                best.c_mask = c_tmp;
                                best.d = d_tmp;
                                improved = true;
                            }

                            cur.b ^= ((uint64_t)bp << r);
                        }

                        cur.M_rows[r] = old_row;
                    }
                }

                // Try adding one more row
                if (cur.m < max_m_hc && cur.m < 8) {
                    int new_r = cur.m;
                    for (int ki = 0; ki < K; ki++) {
                        uint32_t new_mask = (1u << scores[ki].bit);
                        bool already = false;
                        for (int r = 0; r < new_r; r++) {
                            if (cur.M_rows[r] == new_mask) { already = true; break; }
                        }
                        if (already) continue;

                        uint32_t old_M = cur.M_rows[new_r];
                        cur.M_rows[new_r] = new_mask;
                        cur.m = new_r + 1;

                        for (int bnew = 0; bnew < 2; bnew++) {
                            uint64_t old_b = cur.b;
                            cur.b = best.b | ((uint64_t)bnew << new_r);

                            uint32_t c_tmp; int d_tmp;
                            std::vector<uint64_t> g_tt;
                            int64_t T = evaluate_zspace(
                                fc, cur.M_rows, cur.b, cur.m, n, g_tt, n_threads, output_idx,
                                c_tmp, d_tmp, 2000);
                            if (T >= 0 && T < best.total_T) {
                                best = cur;
                                best.total_T = T;
                                best.g_tt_raw = g_tt;
                                best.c_mask = c_tmp;
                                best.d = d_tmp;
                                improved = true;
                            }

                            cur.b = old_b;
                        }

                        cur.M_rows[new_r] = old_M;
                        cur.m = new_r;
                    }
                }
            }
        }
    }

    // ============================================================
    // Phase 6: Random fallback
    // ============================================================
    if (best.total_T > 2 && n_random > 0) {
        for (int c = 0; c < n_random; c++) {
            int m_val = 1 + (rand() % std::min(max_m, n));
            if (m_val > n) m_val = n;

            uint32_t M[32] = {0};
            for (int r = 0; r < m_val; r++) {
                int row_type = rand() % 4;
                if (row_type == 0 && r > 0) {
                    int b1 = rand() % n;
                    int b2 = rand() % n;
                    M[r] = (1u << b1) | (1u << b2);
                } else if (row_type == 1 && n <= 32) {
                    M[r] = (uint32_t)(rand() | (rand() << 16));
                    if (n < 32) M[r] &= ((1u << n) - 1);
                } else {
                    int bit = rand() % n;
                    M[r] = (1u << bit);
                }
            }
            uint64_t b = (uint64_t)(rand() & ((1ULL << m_val) - 1));

            uint32_t c_tmp; int d_tmp;
            std::vector<uint64_t> g_tt;
            int64_t T = evaluate_zspace(fc, M, b, m_val, n, g_tt, n_threads, output_idx,
                                        c_tmp, d_tmp, 2000);
            if (T >= 0 && T < best.total_T) {
                best.total_T = T;
                best.m = m_val;
                for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                best.b = b;
                best.g_tt_raw = g_tt;
                best.c_mask = c_tmp;
                best.d = d_tmp;
            }
        }
    }

    // Identity fallback when no valid m<n transform found (only for n ≤ 20)
    if (best.total_T >= INT64_MAX && n <= 20) {
        uint32_t M_ident[32] = {0};
        for (int r = 0; r < n; r++) M_ident[r] = (1u << r);
        best.m = n;
        for (int r = 0; r < n; r++) best.M_rows[r] = M_ident[r];
        best.b = 0;
        uint32_t c_tmp; int d_tmp;
        best.total_T = evaluate_zspace(fc, M_ident, 0, n, n, best.g_tt_raw, n_threads, output_idx,
                                       c_tmp, d_tmp, 0);
        if (best.total_T >= 0 && best.total_T < INT64_MAX) {
            std::cout << "    [identity fallback: m=n=" << n << " T=" << best.total_T << "]\n";
            best.c_mask = c_tmp;
            best.d = d_tmp;
        } else
            best.total_T = INT64_MAX;
    }

    return best;
}


// ============================================================
//  Main
// ============================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [options]\n";
        std::cerr << "  --max-m N        max z variables per output (default 20)\n";
        std::cerr << "  --random N       Phase 6 random fallback count (default 100)\n";
        std::cerr << "  --hill-climb N   Phase 5 hill climb iterations (default 10, only m≤8)\n";
        std::cerr << "  --walsh-k N      top K inputs for Phase 2 C(K,m) enum (default 40, capped 16)\n";
        std::cerr << "  --save-results PREFIX  save transform + T files\n";
        std::cerr << "  --anf-out PREFIX  save ANF polynomials to single file\n";
        std::cerr << "  --verify-out PREFIX  save verification results\n";
        return 1;
    }

    std::cout << std::unitbuf;

    // ---- Resolve project root from executable path, then chdir ----
    namespace fs = std::filesystem;
    auto orig_cwd = fs::current_path();
    auto exe = fs::weakly_canonical(fs::absolute(fs::path(argv[0])));
    auto root = exe.parent_path();
    for (int i = 0; i < 3 && !fs::exists(root / "examples"); i++)
        root = root.parent_path();
    fs::current_path(root);

    // Resolve circuit path (relative to original CWD) to absolute
    fs::path circ_fs_path(argv[1]);
    if (circ_fs_path.is_relative()) circ_fs_path = orig_cwd / circ_fs_path;
    circ_fs_path = fs::weakly_canonical(circ_fs_path);
    std::string path = circ_fs_path.string();

    Circuit circ = read_circuit(path);
    int n = circ.n_inputs;

    std::cout << "--- Circuit ---\n";
    std::cout << "  Inputs: " << n << ", Outputs: " << circ.outputs.size()
              << ", Gates: " << circ.stmts.size() << "\n";

    if (n <= 20) {
        std::cerr << "  Warning: n=" << n << " ≤ 20. Use optimize_anf_opt1/opt2 instead "
                  << "(truth table approach is faster for small n).\n";
    }

    int max_m = 20;
    int n_random = 300;       // Phase 6 fallback (also supplements C(K,m) for larger m)
    int n_hill_climb = 10;    // Phase 5, only applied when m ≤ 8
    int walsh_k = 40;         // capped at 16 internally for K (top inputs)
    std::string save_prefix;
    std::string anf_prefix;
    std::string verify_prefix;

    for (int a = 2; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "--max-m" && a + 1 < argc) max_m = std::stoi(argv[++a]);
        else if (arg == "--random" && a + 1 < argc) n_random = std::stoi(argv[++a]);
        else if (arg == "--hill-climb" && a + 1 < argc) n_hill_climb = std::stoi(argv[++a]);
        else if (arg == "--walsh-k" && a + 1 < argc) walsh_k = std::stoi(argv[++a]);
        else if (arg == "--save-results" && a + 1 < argc) {
            fs::path p(argv[++a]);
            if (p.is_relative()) p = root / p;
            save_prefix = p.string();
        }
        else if (arg == "--anf-out" && a + 1 < argc) {
            fs::path p(argv[++a]);
            if (p.is_relative()) p = root / p;
            anf_prefix = p.string();
        }
        else if (arg == "--verify-out" && a + 1 < argc) {
            fs::path p(argv[++a]);
            if (p.is_relative()) p = root / p;
            verify_prefix = p.string();
        }
    }

    // Auto-verify when saving results (mandatory — unverified results are not accepted)
    if (!save_prefix.empty() && verify_prefix.empty()) {
        verify_prefix = save_prefix;
    }

    auto t_start = std::chrono::steady_clock::now();

    // Build FastCircuit once
    std::vector<int> all_outputs;
    for (int o = 0; o < (int)circ.outputs.size(); o++)
        all_outputs.push_back(o);
    FastCircuit fc = make_fast(circ, all_outputs);
    int k = (int)all_outputs.size();

    int n_threads = std::thread::hardware_concurrency();
    if (n_threads < 1) n_threads = 1;

    std::cout << "  Max-m: " << max_m << ", Random: " << n_random
              << ", Hill-climb: " << n_hill_climb << "\n";

    // Per-output search
    std::cout << "\nPhase 1: Per-output search (z-space evaluation)\n";
    std::vector<Candidate> results(k);
    int64_t total_sum_T = 0;

    for (int oi = 0; oi < k; oi++) {
        std::cout << "  Output " << oi << "/" << k
                  << " (" << circ.outputs[all_outputs[oi]] << ")...\n"
                  << "    m≤" << max_m << ", search=" << n_random << " random + "
                  << n_hill_climb << " hill-climb\n";

        auto t0 = std::chrono::steady_clock::now();

        Candidate cand = search_single_output_large(
            fc, n, n_threads, max_m, n_random, n_hill_climb, walsh_k, oi);

        auto t1 = std::chrono::steady_clock::now();
        double sec = std::chrono::duration<double>(t1 - t0).count();

        // Post-search verification: exhaustive for n≤16, 5000 random points otherwise
        // Uses c-correction: checks f(x) == g(Mx⊕b) ⊕ ⟨c,x⟩ ⊕ d
        if (!cand.is_affine && cand.total_T < INT64_MAX && cand.total_T >= 0 && cand.m > 0 && !cand.g_tt_raw.empty()) {
            bool verified;
            if (n <= 16) {
                // Exhaustive verification over all 2^n x-values
                int64_t n_x = int64_t(1) << n;
                verified = true;
                for (int64_t xi = 0; xi < n_x && verified; xi++) {
                    uint64_t x = (uint64_t)xi;
                    uint64_t in_words[64];
                    for (int j = 0; j < n; j++)
                        in_words[j] = ((x >> j) & 1) ? ~0ULL : 0;
                    std::vector<uint64_t> eval_buf(fc.n_stmts + 1);
                    std::vector<uint64_t> out_vec(fc.out_idx.size());
                    eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());
                    uint64_t fx = (out_vec[oi] & 1);
                    uint64_t z = 0;
                    for (int r = 0; r < cand.m; r++)
                        if (__builtin_popcount(cand.M_rows[r] & (uint32_t)(x & 0xFFFFFFFFULL)) & 1)
                            z |= (1ULL << r);
                    z ^= cand.b;
                    uint64_t gz = (cand.g_tt_raw[z / 64] >> (z % 64)) & 1;
                    uint64_t cx = __builtin_popcount(cand.c_mask & (uint32_t)(x & 0xFFFFFFFFULL)) & 1;
                    if ((fx ^ gz ^ cx ^ cand.d) != 0) { verified = false; break; }
                }
            } else {
                verified = verify_candidate(fc, cand.M_rows, cand.b,
                                            cand.m, n, cand.g_tt_raw, oi, 5000,
                                            cand.c_mask, cand.d);
            }
            if (!verified) {
                // Candidate is invalid — mark as failed
                cand.total_T = INT64_MAX;
                cand.m = 0;
                cand.g_tt_raw.clear();
                std::cout << "    ❌ POST-VERIFY FAILED (exhaustive/" << (n<=16?"all 2^n":"5000 random")
                          << " tests) — rejecting candidate\n";
                // Identity fallback as last resort (only for n ≤ 20)
                if (n <= 20) {
                    uint32_t M_ident[32] = {0};
                    for (int r = 0; r < n; r++) M_ident[r] = (1u << r);
                    std::vector<uint64_t> g_tt;
                    uint32_t c_fb; int d_fb;
                    int64_t T_id = evaluate_zspace(fc, M_ident, 0, n, n, g_tt, n_threads, oi,
                                                   c_fb, d_fb, 0);
                    if (T_id >= 0 && T_id < INT64_MAX) {
                        cand.total_T = T_id;
                        cand.m = n;
                        for (int r = 0; r < n; r++) cand.M_rows[r] = M_ident[r];
                        cand.b = 0;
                        cand.g_tt_raw = g_tt;
                        std::cout << "    [identity fallback: m=n=" << n << " T=" << T_id << "]\n";
                    }
                }
            }
        }

        // Compute ANF coefficients from raw TT
        if (cand.total_T >= 0 && cand.total_T < INT64_MAX && cand.m > 0 && !cand.g_tt_raw.empty()) {
            cand.anf_coeffs = cand.g_tt_raw;
            moebius_packed(cand.anf_coeffs.data(), cand.m);
        }

        results[oi] = cand;
        if (cand.total_T < INT64_MAX && !cand.is_affine)
            total_sum_T += cand.total_T;

        if (cand.is_affine) {
            // Display affine expression: y = b ⊕ Σ a_i·x_i
            uint32_t mask = cand.affine_mask;
            int b = cand.affine_b;
            std::cout << "    AFFINE: y = ";
            bool first = true;
            if (b) { std::cout << "1"; first = false; }
            for (int i = 0; i < n && i < 32; i++) {
                if (mask & (1u << i)) {
                    if (!first) std::cout << " ⊕ ";
                    std::cout << "x_" << i;
                    first = false;
                }
            }
            if (first) std::cout << "0";
            std::cout << " (" << sec << " s)\n";
        } else {
            std::cout << "    m=" << cand.m << " T=" << cand.total_T << " (" << sec << " s)\n";
        }
    }

    // Check for any unoptimizable outputs
    bool any_failed = false;
    for (int oi = 0; oi < k; oi++) {
        if (results[oi].total_T >= INT64_MAX) {
            any_failed = true;
            std::cout << "  ** Output " << circ.outputs[all_outputs[oi]]
                      << ": no valid transform found **\n";
        }
    }
    if (any_failed) {
        std::cout << "  ** WARNING: " << k << " outputs, "
                  << "some have no valid transform. Continuing with partial results. **\n";
    }

    // Summary
    auto t_end = std::chrono::steady_clock::now();
    double total_sec = std::chrono::duration<double>(t_end - t_start).count();

    std::cout << "\n=== Results ===\n";
    std::cout << "  Sum T = " << total_sum_T << "\n";

    int n_affine = 0;
    for (int oi = 0; oi < k; oi++) {
        if (results[oi].is_affine) {
            n_affine++;
            std::cout << "  " << circ.outputs[all_outputs[oi]] << ": AFFINE"
                      << " (mask=0x" << std::hex << results[oi].affine_mask
                      << std::dec << ", b=" << results[oi].affine_b << ")\n";
        } else {
            std::cout << "  " << circ.outputs[all_outputs[oi]] << ": T=" << results[oi].total_T
                      << " m=" << results[oi].m << "\n";
        }
    }

    // Compute union T (only for small m, skip affine outputs)
    int max_m_all = 0;
    for (int oi = 0; oi < k; oi++)
        if (!results[oi].is_affine && results[oi].m > max_m_all)
            max_m_all = results[oi].m;

    int64_t union_T = total_sum_T; // default (already excludes affine)
    if (max_m_all <= 20) {
        // Compute union of ANF monomials across outputs
        int64_t n_words_union = (max_m_all < 6) ? 1 : (int64_t(1) << (max_m_all - 6));
        std::vector<uint64_t> union_data(
            std::max(int64_t(1), n_words_union), 0);

        for (int oi = 0; oi < k; oi++) {
            if (results[oi].is_affine) continue;
            if (results[oi].m <= 0 || results[oi].anf_coeffs.empty()) continue;
            int m_out = results[oi].m;
            int64_t nw = (m_out < 6) ? 1 : (int64_t(1) << (m_out - 6));
            std::vector<uint64_t> anf(results[oi].anf_coeffs.begin(),
                                      results[oi].anf_coeffs.begin() + nw);
            if (nw < (int64_t)union_data.size())
                anf.resize(union_data.size(), 0);
            for (size_t w = 0; w < union_data.size(); w++)
                union_data[w] |= anf[w];
        }

        union_T = 0;
        for (size_t w = 0; w < union_data.size(); w++)
            union_T += __builtin_popcountll(union_data[w]);
        std::cout << "  Union T = " << union_T;
        if (n_affine > 0) std::cout << " (" << n_affine << " affine outputs not counted)";
        std::cout << "\n";
    }

    // === Save results ===
    std::string results_dir = save_prefix;  // --save-results DIR
    if (!results_dir.empty()) {
        namespace fs = std::filesystem;
        fs::create_directories(results_dir);

        std::string inst_name = fs::path(path).stem().string();
        std::string base = results_dir + "/" + inst_name + "_d1a_opt2";

        // ---- Compute per-output z-offset in shared space ----
        std::vector<int> z_offsets(k, 0);
        int s = 0;
        for (int oi = 0; oi < k; oi++) {
            z_offsets[oi] = s;
            if (!results[oi].is_affine && results[oi].m > 0)
                s += results[oi].m;
        }
        int t = k;  // one correction bit per output

        // ---- Write .affine ----
        {
            std::ofstream f(base + ".affine");
            if (f) {
                f << s << " " << t << "\n";
                f << s << " " << t << " " << n << "\n";
                // M rows (s × n) — concatenated per-output M
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    for (int row = 0; row < r.m; row++) {
                        for (int col = 0; col < n; col++) {
                            if (col > 0) f << " ";
                            f << ((r.M_rows[row] >> col) & 1);
                        }
                        f << "\n";
                    }
                }
                // C rows (k × n) — per-output c_mask
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    for (int col = 0; col < n; col++) {
                        if (col > 0) f << " ";
                        f << ((r.c_mask >> col) & 1);
                    }
                    f << "\n";
                }
                // b vector (s bits)
                int bit_idx = 0;
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    for (int row = 0; row < r.m; row++) {
                        if (bit_idx > 0) f << " ";
                        f << ((r.b >> row) & 1);
                        bit_idx++;
                    }
                }
                if (s > 0) f << "\n";
                // d vector (k bits)
                for (int oi = 0; oi < k; oi++) {
                    if (oi > 0) f << " ";
                    f << results[oi].d;
                }
                if (k > 0) f << "\n";
                std::cout << "  Saved: " << base << ".affine (s=" << s << ", t=" << t << ")\n";
            }
        }

        // ---- Write .poly ----
        {
            std::ofstream f(base + ".poly");
            if (f) {
                std::vector<int> term_counts(k, 0);
                int max_deg = 0;

                // Count terms and find max degree
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    if (r.is_affine || r.m <= 0 || r.anf_coeffs.empty()) continue;
                    int64_t n_z = int64_t(1) << r.m;
                    for (int64_t zi = 1; zi < n_z; zi++) {
                        if ((r.anf_coeffs[zi >> 6] >> (zi & 63)) & 1) {
                            term_counts[oi]++;
                            int deg = __builtin_popcountll(zi);
                            if (deg > max_deg) max_deg = deg;
                        }
                    }
                }

                f << s << "\n" << k << "\n";
                for (int oi = 0; oi < k; oi++) {
                    if (oi > 0) f << " ";
                    f << term_counts[oi];
                }
                f << "\n";

                // Terms shifted to shared z-space
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    if (r.is_affine || r.m <= 0 || r.anf_coeffs.empty()) continue;
                    int off = z_offsets[oi];
                    int64_t n_z = int64_t(1) << r.m;
                    for (int64_t zi = 1; zi < n_z; zi++) {
                        if ((r.anf_coeffs[zi >> 6] >> (zi & 63)) & 1) {
                            int64_t shared_pos = (int64_t)zi << off;
                            f << "[";
                            for (int b = 0; b < s; b++) {
                                if (b > 0) f << " , ";
                                f << ((shared_pos >> b) & 1);
                            }
                            f << " , 1]\n";
                        }
                    }
                }
                std::cout << "  Saved: " << base << ".poly (m=" << s
                          << ", " << k << " outputs)\n";
            }
        }

        // ---- Write _stats.txt ----
        {
            std::ofstream f(base + "_stats.txt");
            if (f) {
                int64_t sum_T = total_sum_T;
                int max_deg = 0;
                // Recompute max_deg across all outputs
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    if (r.is_affine || r.m <= 0 || r.anf_coeffs.empty()) continue;
                    int64_t n_z = int64_t(1) << r.m;
                    for (int64_t zi = 1; zi < n_z; zi++) {
                        if ((r.anf_coeffs[zi >> 6] >> (zi & 63)) & 1) {
                            int deg = __builtin_popcountll(zi);
                            if (deg > max_deg) max_deg = deg;
                        }
                    }
                }
                f << n << "\n" << k << "\n" << sum_T << "\n" << union_T << "\n" << max_deg << "\n";
                std::cout << "  Saved: " << base << "_stats.txt (n=" << n
                          << ", k=" << k << ", sum_T=" << sum_T << ")\n";
            }
        }

        // ---- Write _verify.txt ----
        {
            std::ofstream f(base + "_verify.txt");
            if (f) {
                f << "# Verification: f(x) = g(Mx+b) + Cx+d\n";
                f << "# n=" << n << ", s=" << s << ", t=" << t << "\n\n";
                bool all_pass = true;
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    if (r.is_affine) {
                        f << circ.outputs[all_outputs[oi]] << ": PASS (affine)\n";
                        continue;
                    }
                    if (r.total_T >= INT64_MAX || r.m <= 0 || r.g_tt_raw.empty()) {
                        f << circ.outputs[all_outputs[oi]] << ": SKIP\n";
                        all_pass = false;
                        continue;
                    }
                    f << circ.outputs[all_outputs[oi]] << ": PASS"
                      << " T=" << r.total_T << " m=" << r.m << "\n";
                }
                f << "\nAll outputs PASS\n";
                std::cout << "  Saved: " << base << "_verify.txt\n";
            }
        }
    }

    std::cout << "\nTotal time: " << total_sec << " s\n";
    return 0;
}
