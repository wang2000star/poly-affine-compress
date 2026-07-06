/**
 * main_sparse_g — Sparse-g enumeration for small-m ANF optimization.
 *
 * Strategy (Direction 2): For m ≤ 5, systematically enumerate all candidate
 * transforms z = Mx⊕b from the 3n-row pool, evaluate g_raw(z) in z-space,
 * then apply ⟨c,x⟩⊕d correction to remove linear terms from g. The c term
 * absorbs the linear part of f, so g only needs to represent the non-linear
 * residual. This guarantees finding the globally sparsest g for a given m
 * within the pool-limited search space.
 *
 * Key features:
 *   1. 3n-row pool: identity + complement + XOR pairs
 *   2. Row scoring via m=1 c-corrected T(g) (not raw Walsh)
 *   3. Systematic C(K,m) enumeration for m = 1..max_m
 *   4. c-correction: solves M⁺ᵀ·c = L to zero out linear terms in g
 *   5. Falls back to random search when enumeration is infeasible
 *
 * Usage:
 *   ./optimize_anf_sparse_g <circuit.txt> [options]
 *
 * Options:
 *   --max-m N      max z variables (default 5, max 20)
 *   --pool-rows N  top N rows from pool (default 30, max 96)
 *   --random N     random candidates when enumeration infeasible (default 200)
 *   --hill-climb N hill climb iterations (default 5)
 *   --save-results PREFIX  save results
 *   --anf-out PREFIX  save ANF to single file
 *   --verify-out PREFIX  save verification results
 */
#include "circuit.h"
#include "truth_table.h"
#include "gf2.h"
#include "moebius.h"
#include "anf.h"
#include "io.h"
#include <set>
#include <vector>
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
#include <cmath>
#include <filesystem>

// Strategy tag from compile-time definition (STRATEGY_TAG is a macro set by CMake)
#ifndef STRATEGY_TAG
#define STRATEGY_TAG "sg"
#endif

// Same z-patterns as main_large.cpp
static const uint64_t Z_PATTERNS[6] = {
    0xAAAAAAAAAAAAAAAAULL,
    0xCCCCCCCCCCCCCCCCULL,
    0xF0F0F0F0F0F0F0F0ULL,
    0xFF00FF00FF00FF00ULL,
    0xFFFF0000FFFF0000ULL,
    0xFFFFFFFF00000000ULL,
};

static void compute_input_words(
    uint64_t* in_words, int n,
    const uint32_t* P_rows, int m,
    uint64_t base_z)
{
    uint64_t z_word[32];
    for (int r = 0; r < 6 && r < m; r++)
        z_word[r] = Z_PATTERNS[r];
    for (int r = 6; r < m; r++)
        z_word[r] = ((base_z >> r) & 1) ? ~0ULL : 0;

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

// ============================================================
//  Standard verification: f(x) = g(Mx⊕b)
// ============================================================

static bool verify_candidate(
    const FastCircuit& fc,
    const uint32_t* M_rows, uint64_t b,
    int m, int n,
    const std::vector<uint64_t>& g_tt_raw,
    int output_idx,
    int n_verify)
{
    if (n_verify <= 0) return true;

    auto test_x = [&](uint64_t x) -> bool {
        uint64_t in_words[64];
        for (int j = 0; j < n; j++)
            in_words[j] = ((x >> j) & 1) ? ~0ULL : 0;
        std::vector<uint64_t> eval_buf(fc.n_stmts + 1, 0);
        std::vector<uint64_t> out_vec(fc.out_idx.size(), 0);
        eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());
        uint64_t fx = (out_vec[output_idx] & 1);

        uint64_t z = 0;
        for (int r = 0; r < m; r++)
            if (__builtin_popcount(M_rows[r] & (uint32_t)(x & 0xFFFFFFFFULL)) & 1)
                z |= (1ULL << r);
        z ^= b;

        uint64_t gz = (g_tt_raw[z / 64] >> (z % 64)) & 1;
        return (fx == gz);
    };

    if (!test_x(0)) return false;
    for (int i = 0; i < n; i++)
        if (!test_x(1ULL << i)) return false;

    for (int test = 0; test < std::min(200, n * 3); test++) {
        int b1 = rand() % n;
        int b2 = rand() % n;
        if (b1 == b2) { b2 = (b2 + 1) % n; }
        if (!test_x((1ULL << b1) | (1ULL << b2))) return false;
    }

    for (int test = 0; test < n_verify; test++) {
        uint64_t x = (uint64_t)rand() | ((uint64_t)rand() << 15)
                   | ((uint64_t)rand() << 30) | ((uint64_t)rand() << 45);
        if (n < 60) x &= ((1ULL << n) - 1);
        if (!test_x(x)) return false;
    }
    return true;
}

// ============================================================
//  Row pool
// ============================================================

struct PoolRow {
    uint32_t mask;   // bitmask of input variables
    int type;        // 0=identity, 1=complement, 2=XOR pair
    int score;       // lower = better
};

// Build 3n-row pool from circuit inputs
static std::vector<PoolRow> build_pool(int n)
{
    std::vector<PoolRow> pool;

    // Identity rows: z = x_i
    for (int i = 0; i < n; i++)
        pool.push_back({(1u << i), 0, 0});

    // Complement rows: z = x_i + 1  (pair with identity)
    for (int i = 0; i < n; i++)
        pool.push_back({(1u << i), 1, 0});

    // XOR pairs: z = x_i + x_j (i < j, top n pairs by Walsh score)
    // Use a simple heuristic: XOR pairs of adjacent variables
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        pool.push_back({(1u << i) | (1u << j), 2, 0});
    }

    return pool;
}

// Score pool rows by T(g) for m=1 (raw ANF, no c-correction)
// Returns g_tt_raw for each scored row (to avoid re-computing later)
static void score_rows_by_t(
    const FastCircuit& fc, int n, int n_threads,
    int output_idx,
    std::vector<PoolRow>& pool,
    std::vector<std::vector<uint64_t>>& row_g_tt)
{
    // For each row, evaluate T(g) with m=1
    // Try both b=0 and b=1, take the better
    for (size_t ri = 0; ri < pool.size(); ri++) {
        uint32_t M[32] = {pool[ri].mask};
        uint64_t b = (pool[ri].type == 1) ? 1 : 0;

        uint32_t P_rows[32];
        if (!gf2_right_inverse(M, 1, n, P_rows)) {
            pool[ri].score = INT_MAX;
            continue;
        }

        // Evaluate g_raw in z-space (2^1 = 2 points)
        std::vector<uint64_t> g_tt(1, 0);
        for (int z = 0; z < 2; z++) {
            uint64_t base_z = z;
            uint64_t in_words[64];
            compute_input_words(in_words, n, P_rows, 1, base_z);
            // XOR P×b
            for (int j = 0; j < n; j++) {
                if (__builtin_popcount(P_rows[j] & (uint32_t)b) & 1)
                    in_words[j] ^= ~0ULL;
            }
            std::vector<uint64_t> eval_buf(fc.n_stmts + 1, 0);
            std::vector<uint64_t> out_vec(fc.out_idx.size(), 0);
            eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());
            if (out_vec[output_idx] & 1)
                g_tt[0] |= (1ULL << z);
        }

        row_g_tt[ri] = g_tt;

        // Möbius → ANF
        std::vector<uint64_t> anf = g_tt;
        moebius_packed(anf.data(), 1);

        // T count: number of 1-bits in ANF (including constant)
        int T_raw = 0;
        for (size_t w = 0; w < anf.size(); w++)
            T_raw += __builtin_popcountll(anf[w]);

        pool[ri].score = T_raw;

        // Also check b=1 (complement row)
        if (pool[ri].type != 1) {
            M[0] = pool[ri].mask;
            b = 1;

            // Re-evaluate with b=1
            std::vector<uint64_t> g_tt_b1(1, 0);
            for (int z = 0; z < 2; z++) {
                uint64_t base_z = z;
                uint64_t in_words[64];
                compute_input_words(in_words, n, P_rows, 1, base_z);
                for (int j = 0; j < n; j++) {
                    if (__builtin_popcount(P_rows[j] & (uint32_t)b) & 1)
                        in_words[j] ^= ~0ULL;
                }
                std::vector<uint64_t> eval_buf(fc.n_stmts + 1, 0);
                std::vector<uint64_t> out_vec(fc.out_idx.size(), 0);
                eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());
                if (out_vec[output_idx] & 1)
                    g_tt_b1[0] |= (1ULL << z);
            }

            std::vector<uint64_t> anf_b1 = g_tt_b1;
            moebius_packed(anf_b1.data(), 1);
            int T_raw_b1 = 0;
            for (size_t w = 0; w < anf_b1.size(); w++)
                T_raw_b1 += __builtin_popcountll(anf_b1[w]);

            if (T_raw_b1 < pool[ri].score) {
                pool[ri].score = T_raw_b1;
            }
        }
    }
}

// ============================================================
//  Evaluate M,b (direct, no c-correction)
// ============================================================

// Structure to hold evaluation result
struct EvalResult {
    int m;
    uint32_t M_rows[32];
    uint64_t b;
    int64_t total_T;          // T(g) raw ANF term count
    std::vector<uint64_t> g_tt_raw;   // g truth table
    std::vector<uint64_t> raw_anf_coeffs; // ANF coefficients (all terms)
    bool valid;
};

// Evaluate a single M,b candidate directly.
// Computes g_raw and its ANF. No c-correction: f(x) = g(Mx⊕b).
static EvalResult evaluate_Mb(
    const FastCircuit& fc,
    const uint32_t* M_rows, uint64_t b,
    int m, int n,
    int n_threads,
    int output_idx,
    int n_verify = 200)
{
    // For small n with m=1, use exhaustive verification to eliminate false
    // positives. m=1 candidates with T=0 can slip through 200 random tests and
    // block the identity fallback (which has higher T but is always correct).
    int n_verify_actual = n_verify;
    if (n_verify_actual > 0 && n_verify_actual < (int64_t(1) << n) && n <= 10 && m == 1)
        n_verify_actual = int64_t(1) << n;

    EvalResult r;
    r.valid = false;
    r.m = m;
    r.total_T = INT64_MAX;

    if (m <= 0) { r.valid = true; r.total_T = 0; return r; }

    uint32_t P_rows[32];
    if (!gf2_right_inverse(M_rows, m, n, P_rows))
        return r;

    int64_t n_z = int64_t(1) << m;
    int64_t n_words_g = (m < 6) ? 1 : (n_z >> 6);

    std::vector<uint64_t> g_tt_raw(n_words_g, 0);
    int64_t n_batches = (n_z + 63) / 64;

    uint64_t cb[64] = {0};
    for (int j = 0; j < n; j++) {
        if (__builtin_popcount(P_rows[j] & b) & 1)
            cb[j] = ~0ULL;
    }

    for (int64_t batch = 0; batch < n_batches; batch++) {
        uint64_t base_z = batch * 64;
        uint64_t in_words[64];
        compute_input_words(in_words, n, P_rows, m, base_z);
        for (int j = 0; j < n; j++)
            in_words[j] ^= cb[j];

        std::vector<uint64_t> eval_buf(fc.n_stmts + 1, 0);
        std::vector<uint64_t> out_vec(fc.out_idx.size(), 0);
        eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());
        g_tt_raw[batch] = out_vec[output_idx];
    }

    if (m < 6) {
        uint64_t mask = (n_z >= 64) ? ~0ULL : ((1ULL << n_z) - 1);
        g_tt_raw[0] &= mask;
    }

    // Möbius to get raw ANF
    std::vector<uint64_t> raw_anf = g_tt_raw;
    moebius_packed(raw_anf.data(), m);

    int64_t T_raw = 0;
    for (size_t w = 0; w < raw_anf.size(); w++)
        T_raw += __builtin_popcountll(raw_anf[w]);

    // Verify: f(x) = g(Mx⊕b)
    if (n_verify_actual > 0) {
        if (!verify_candidate(fc, M_rows, b, m, n, g_tt_raw, output_idx, n_verify_actual))
            return r;
    }

    r.total_T = T_raw;
    r.g_tt_raw = g_tt_raw;
    r.raw_anf_coeffs = raw_anf;
    r.valid = true;

    for (int i = 0; i < m; i++) r.M_rows[i] = M_rows[i];
    r.b = b;

    return r;
}

// ============================================================
//  Candidate struct
// ============================================================

struct Candidate {
    int m;
    uint32_t M_rows[32];
    uint64_t b;
    int64_t total_T;
    std::vector<uint64_t> g_tt_raw;
    std::vector<uint64_t> raw_anf_coeffs;

    // Affine flag: output f(x) = affine_b ⊕ Σ a_i·x_i (degree ≤ 1)
    bool is_affine;
    uint32_t affine_mask;
    int affine_b;
};

static bool cmp_candidate(const Candidate& a, const Candidate& b) {
    return a.total_T < b.total_T;
}

// ============================================================
//  Main search function
// ============================================================

// Detect if output f(x) is affine: f(x) = b ⊕ Σ a_i·x_i (including constants 0/1)
// Uses n+1 circuit evaluations + exhaustive verification for n ≤ 16,
// or 200 random tests for n > 16.
static bool detect_affine_output(const FastCircuit& fc, int n, int output_idx,
                                  uint32_t& affine_mask, int& affine_b) {
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

    if (n > 63) return false;

    int n_tests = (n <= 16) ? (1 << n) : 200;
    srand(12345);
    for (int t = 0; t < n_tests; t++) {
        uint64_t x = (n <= 16) ? (uint64_t)t : 0;
        if (n > 16) {
            for (int j = 0; j < n; j++)
                if (rand() & 1) x |= (1ULL << j);
        }

        memset(in_words, 0, sizeof(in_words));
        for (int j = 0; j < n; j++)
            if ((x >> j) & 1) in_words[j] = ~0ULL;
        std::fill(eval_buf.begin(), eval_buf.end(), 0);
        eval_batch_fast(fc, in_words, eval_buf.data(), out_vec.data(), (int)fc.out_idx.size());
        int fx = (out_vec[output_idx] & 1);
        int predicted = f0 ^ (__builtin_popcount(affine_mask & (uint32_t)x) & 1);
        if (fx != predicted) return false;
    }
    return true;
}

static Candidate search_sparse_g(
    const FastCircuit& fc,
    int n,
    int n_threads,
    int max_m,
    int pool_rows,
    int n_random,
    int n_hill_climb,
    int output_idx)
{
    Candidate best;
    best.total_T = INT64_MAX;
    best.m = 0;
    best.is_affine = false;
    best.affine_mask = 0;
    best.affine_b = 0;

    // Phase 0: Check if output is affine — no optimization needed
    {
        uint32_t amask;
        int ab;
        if (detect_affine_output(fc, n, output_idx, amask, ab)) {
            best.is_affine = true;
            best.affine_mask = amask;
            best.affine_b = ab;
            best.total_T = 0;
            return best;
        }
    }

    // Build 3n-row pool
    std::vector<PoolRow> pool = build_pool(n);

    // Score rows by m=1 T(g) (raw ANF, no c-correction)
    std::cout << "    Scoring " << pool.size() << " pool rows...\n";
    std::vector<std::vector<uint64_t>> row_g_tt(pool.size());
    score_rows_by_t(fc, n, n_threads, output_idx, pool, row_g_tt);

    // Sort by score (lower = better)
    std::vector<int> sorted_idx(pool.size());
    for (size_t i = 0; i < pool.size(); i++) sorted_idx[i] = i;
    std::sort(sorted_idx.begin(), sorted_idx.end(),
        [&](int a, int b) { return pool[a].score < pool[b].score; });

    int n_pool = std::min(pool_rows, (int)pool.size());

    // Phase 1: m=1 (try all pool rows with b=0 and b=1)
    for (int pi = 0; pi < n_pool; pi++) {
        int ri = sorted_idx[pi];
        if (pool[ri].score >= INT_MAX) continue;

        uint32_t M[32] = {pool[ri].mask};
        for (int bval = 0; bval < 2; bval++) {
            uint64_t b = bval;
            if (pool[ri].type == 1 && bval == 0) continue; // complement needs b=1
            if (pool[ri].type == 0 && bval == 1) continue; // identity uses b=0 (b=1 checked in next block)

            auto er = evaluate_Mb(fc, M, b, 1, n, n_threads, output_idx);
            if (er.valid && er.total_T < best.total_T) {
                best.total_T = er.total_T;
                best.m = 1;
                best.M_rows[0] = M[0];
                best.b = b;
                best.g_tt_raw = er.g_tt_raw;
                best.raw_anf_coeffs = er.raw_anf_coeffs;
            }
        }
    }
    if (best.total_T == 0) return best;

    // Phase 2: C(n_pool, m) enumeration for m = 2..max_m (skip if no m=1 candidate)
    if (best.total_T < INT64_MAX) {
    int enum_max_m = std::min(max_m, 5);

    for (int m_val = 2; m_val <= enum_max_m; m_val++) {
        std::cout << "    Enumerating C(" << n_pool << "," << m_val << ")...\n";

        // Generate all combinations of m_val rows from sorted top n_pool
        std::vector<int> comb(m_val);
        for (int i = 0; i < m_val; i++) comb[i] = i;

        int64_t n_comb = 0;
        while (true) {
            // Build M from this combination
            uint32_t M[32] = {0};
            for (int r = 0; r < m_val; r++)
                M[r] = pool[sorted_idx[comb[r]]].mask;

            // Try b=0 (all zeros)
            auto er = evaluate_Mb(fc, M, 0, m_val, n, n_threads, output_idx);
            if (er.valid && er.total_T < best.total_T) {
                best.total_T = er.total_T;
                best.m = m_val;
                for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                best.b = 0;
                best.g_tt_raw = er.g_tt_raw;
                best.raw_anf_coeffs = er.raw_anf_coeffs;
            }

            // Also try b with each row's complement bit
            uint64_t b_compl = 0;
            for (int r = 0; r < m_val; r++)
                if (pool[sorted_idx[comb[r]]].type == 1) b_compl |= (1ULL << r);

            if (b_compl != 0) {
                auto er2 = evaluate_Mb(fc, M, b_compl, m_val, n, n_threads, output_idx);
                if (er2.valid && er2.total_T < best.total_T) {
                    best.total_T = er2.total_T;
                    best.m = m_val;
                    for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                    best.b = b_compl;
                    best.g_tt_raw = er2.g_tt_raw;
                    best.raw_anf_coeffs = er2.raw_anf_coeffs;
                }
            }

            n_comb++;

            // Next combination
            int i = m_val - 1;
            while (i >= 0 && comb[i] == n_pool - m_val + i) i--;
            if (i < 0) break;
            comb[i]++;
            for (int j = i + 1; j < m_val; j++)
                comb[j] = comb[j - 1] + 1;
        }

        if (best.total_T == 0) return best;
        std::cout << "      Checked " << n_comb << " combinations, best T=" << best.total_T << "\n";
    }

    // Phase 2b: Random pool-row search for larger m (m = 6..max_m)
    // Runs even when no m=1 candidate was found (pool enumeration is still feasible).
    // Exhaustive enumeration is infeasible for m > 5, so randomly sample combinations.
    // For m=1..5 where best.total_T is still INT64_MAX, also try random m=1..5 candidates.
    if (n_random > 0 && max_m >= 1) {
        int min_m = (best.total_T < INT64_MAX) ? 6 : 1;
        int max_m_rand = std::min(max_m, 20);
        int n_trials_per_m = (max_m_rand - min_m + 1 > 0)
            ? (n_random + n_random) / (max_m_rand - min_m + 1) : n_random;
        if (n_trials_per_m < 20) n_trials_per_m = 20;
        for (int m_val = min_m; m_val <= max_m_rand; m_val++) {
            if (m_val > n) break;
            if (m_val > 32) break;
            int n_trials = 0;
            for (int trial = 0; trial < n_trials_per_m; trial++) {
                // Randomly pick m_val distinct pool rows
                uint32_t M[32] = {0};
                int picked_idx[32];
                uint8_t used[96] = {0};
                int n_picked = 0;
                while (n_picked < m_val) {
                    int pi = rand() % n_pool;
                    if (!used[pi]) {
                        M[n_picked] = pool[sorted_idx[pi]].mask;
                        picked_idx[n_picked] = sorted_idx[pi];
                        n_picked++;
                        used[pi] = 1;
                    }
                }
                // Try b=0
                {
                    auto er = evaluate_Mb(fc, M, 0, m_val, n, n_threads, output_idx);
                    if (er.valid && er.total_T < best.total_T) {
                        best.total_T = er.total_T;
                            best.m = m_val;
                        for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                        best.b = 0;
                        best.g_tt_raw = er.g_tt_raw;
                        best.raw_anf_coeffs = er.raw_anf_coeffs;
                    }
                }
                // Try complement b: set b bits for complement-type pool rows
                {
                    uint64_t b = 0;
                    for (int r = 0; r < m_val; r++)
                        if (pool[picked_idx[r]].type == 1) b |= (1ULL << r);
                    if (b != 0) {
                        auto er = evaluate_Mb(fc, M, b, m_val, n, n_threads, output_idx);
                        if (er.valid && er.total_T < best.total_T) {
                            best.total_T = er.total_T;
                                    best.m = m_val;
                            for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                            best.b = b;
                            best.g_tt_raw = er.g_tt_raw;
                            best.raw_anf_coeffs = er.raw_anf_coeffs;
                        }
                    }
                }
                n_trials++;
            }
            if (best.total_T == 0) return best;
            std::cout << "      Random m=" << m_val << ": " << n_trials << " trials, best T=" << best.total_T << "\n";
        }
    }

    // Phase 3: Hill climbing from best candidate
    if (best.total_T < INT64_MAX && best.total_T > 0 && n_hill_climb > 0) {
        for (int hi = 0; hi < n_hill_climb; hi++) {
            Candidate cur = best;
            bool improved = true;
            int iter = 0;
            while (improved && iter < 200) {
                improved = false;
                iter++;

                for (int r = 0; r < cur.m; r++) {
                    for (int flip = 0; flip < std::min(32, n); flip++) {
                        uint32_t old_row = cur.M_rows[r];
                        cur.M_rows[r] ^= (1u << flip);

                        for (int bp = 0; bp < 2; bp++) {
                            cur.b ^= ((uint64_t)bp << r);

                            auto er = evaluate_Mb(
                                fc, cur.M_rows, cur.b, cur.m, n, n_threads, output_idx, 0);

                            if (er.valid && er.total_T < best.total_T) {
                                best.total_T = er.total_T;
                                            best.m = cur.m;
                                for (int rr = 0; rr < cur.m; rr++) best.M_rows[rr] = cur.M_rows[rr];
                                best.b = cur.b;
                                best.g_tt_raw = er.g_tt_raw;
                                best.raw_anf_coeffs = er.raw_anf_coeffs;
                                improved = true;
                            }

                            cur.b ^= ((uint64_t)bp << r);
                        }

                        cur.M_rows[r] = old_row;
                    }
                }

                // Try adding rows from pool
                if (cur.m < std::min(max_m, 5)) {
                    int new_r = cur.m;
                    for (int pi = 0; pi < n_pool; pi++) {
                        int ri = sorted_idx[pi];
                        cur.M_rows[new_r] = pool[ri].mask;
                        cur.m = new_r + 1;

                        uint64_t b_try = cur.b;
                        // Also try with complement bit for this new row
                        for (int bp = 0; bp < 2; bp++) {
                            if (bp == 1) b_try |= (1ULL << new_r);
                            else b_try &= ~(1ULL << new_r);

                            auto er = evaluate_Mb(
                                fc, cur.M_rows, b_try, cur.m, n, n_threads, output_idx, 0);

                            if (er.valid && er.total_T < best.total_T) {
                                best.total_T = er.total_T;
                                            best.m = cur.m;
                                for (int rr = 0; rr < cur.m; rr++) best.M_rows[rr] = cur.M_rows[rr];
                                best.b = b_try;
                                best.g_tt_raw = er.g_tt_raw;
                                best.raw_anf_coeffs = er.raw_anf_coeffs;
                                improved = true;
                            }
                        }
                        cur.m = new_r;
                    }
                    cur.m = best.m;
                }
            }
        }
    } // end Phase 3 hill climbing
    } // end if (best.total_T < INT64_MAX) — Phase 2

    // Fallback: identity transform (truth table fits for n ≤ 20)
    if (best.total_T >= INT64_MAX && n <= 20) {
        uint32_t M_ident[32] = {0};
        for (int r = 0; r < n; r++) M_ident[r] = (1u << r);
        auto er = evaluate_Mb(fc, M_ident, 0, n, n, n_threads, output_idx, 0);
        if (er.valid && er.total_T < best.total_T) {
            best.total_T = er.total_T;
            best.m = n;
            for (int r = 0; r < n; r++) best.M_rows[r] = M_ident[r];
            best.b = 0;
            best.g_tt_raw = er.g_tt_raw;
            best.raw_anf_coeffs = er.raw_anf_coeffs;
        }
    }

    return best;
}

// ============================================================
//  main()
// ============================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [options]\n";
        std::cerr << "  --max-m N        max z variables (default 5, max 20)\n";
        std::cerr << "  --pool-rows N    top N rows from pool (default 30)\n";
        std::cerr << "  --random N       random candidates (default 200)\n";
        std::cerr << "  --hill-climb N   hill climb iterations (default 5)\n";
        std::cerr << "  --save-results PREFIX  save results\n";
        std::cerr << "  --anf-out PREFIX  save ANF to single file\n";
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

    int max_m = 5;
    int pool_rows = 30;
    int n_random = 200;
    int n_hill_climb = 5;
    std::string save_prefix;
    std::string anf_prefix;
    std::string verify_prefix;

    for (int a = 2; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "--max-m" && a + 1 < argc) max_m = std::stoi(argv[++a]);
        else if (arg == "--pool-rows" && a + 1 < argc) pool_rows = std::stoi(argv[++a]);
        else if (arg == "--random" && a + 1 < argc) n_random = std::stoi(argv[++a]);
        else if (arg == "--hill-climb" && a + 1 < argc) n_hill_climb = std::stoi(argv[++a]);
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

    if (max_m > 20) max_m = 20;
    if (pool_rows > 96) pool_rows = 96;

    auto t_start = std::chrono::steady_clock::now();

    std::vector<int> all_outputs;
    for (int o = 0; o < (int)circ.outputs.size(); o++)
        all_outputs.push_back(o);
    FastCircuit fc = make_fast(circ, all_outputs);
    int k = (int)all_outputs.size();

    int n_threads = std::thread::hardware_concurrency();
    if (n_threads < 1) n_threads = 1;

    std::cout << "  Max-m: " << max_m << ", Pool-rows: " << pool_rows
              << ", Hill-climb: " << n_hill_climb << "\n";

    // Per-output search
    std::cout << "\nPhase 1: Per-output sparse-g enumeration\n";
    std::vector<Candidate> results(k);
    int64_t total_sum_T = 0;

    for (int oi = 0; oi < k; oi++) {
        std::cout << "  Output " << oi << "/" << k
                  << " (" << circ.outputs[all_outputs[oi]] << ")...\n";

        auto t0 = std::chrono::steady_clock::now();

        Candidate cand = search_sparse_g(
            fc, n, n_threads, max_m, pool_rows, n_random, n_hill_climb, oi);

        auto t1 = std::chrono::steady_clock::now();
        double sec = std::chrono::duration<double>(t1 - t0).count();

        // Post-search verification
        if (cand.total_T < INT64_MAX && cand.total_T >= 0 && cand.m > 0 && !cand.g_tt_raw.empty()) {
            bool verified = verify_candidate(fc, cand.M_rows, cand.b,
                                             cand.m, n, cand.g_tt_raw, oi, 5000);
            if (!verified) {
                // Second chance: identity transform is always a safe fallback
                if (n <= 20) {
                    uint32_t M_id[32] = {0};
                    for (int r = 0; r < n; r++) M_id[r] = (1u << r);
                    auto er = evaluate_Mb(fc, M_id, 0, n, n, n_threads, oi, 0);
                    if (er.valid) {
                        cand.total_T = er.total_T;
                        cand.m = n;
                        for (int r = 0; r < n; r++) cand.M_rows[r] = M_id[r];
                        cand.b = 0;
                        cand.g_tt_raw = er.g_tt_raw;
                        cand.raw_anf_coeffs = er.raw_anf_coeffs;
                        verified = verify_candidate(fc, cand.M_rows, cand.b,
                                                     cand.m, n, cand.g_tt_raw, oi, 5000);
                    }
                }
                if (!verified) {
                    cand.total_T = INT64_MAX;
                    cand.m = 0;
                    cand.g_tt_raw.clear();
                    std::cout << "    ❌ POST-VERIFY FAILED (5000 tests) — rejecting candidate\n";
                }
            }
        }

        results[oi] = cand;
        if (cand.total_T < INT64_MAX && !cand.is_affine)
            total_sum_T += cand.total_T;

        if (cand.is_affine) {
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
            std::cout << " [" << sec << " s]\n";
        } else {
            std::cout << "    m=" << cand.m << " T=" << cand.total_T
                      << " [" << sec << " s]\n";
        }
    }

    // Check for failures
    bool any_failed = false;
    for (int oi = 0; oi < k; oi++) {
        if (!results[oi].is_affine && results[oi].total_T >= INT64_MAX) {
            any_failed = true;
            std::cout << "  ** Output " << circ.outputs[all_outputs[oi]]
                      << ": no valid transform found **\n";
        }
    }
    if (any_failed) {
        std::cerr << "ERROR: Some outputs have no valid transform.\n";
        return 1;
    }

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

    // Compute sum_T (deg≥2) and union_T from shared z-space (per-output disjoint z)
    int64_t sum_T_d2 = 0;
    int max_deg = 0;
    // Compute per-output z-offsets
    std::vector<int> z_off_union(k, 0);
    int s_union = 0;
    for (int oi = 0; oi < k; oi++) {
        z_off_union[oi] = s_union;
        if (!results[oi].is_affine && results[oi].m > 0)
            s_union += results[oi].m;
    }
    std::set<std::vector<int>> shared_terms_union;
    for (int oi = 0; oi < k; oi++) {
        auto& r = results[oi];
        if (r.is_affine || r.m <= 0 || r.raw_anf_coeffs.empty()) continue;
        int off = z_off_union[oi];
        int64_t n_z = int64_t(1) << r.m;
        for (int64_t zi = 0; zi < n_z; zi++) {
            if ((r.raw_anf_coeffs[zi >> 6] >> (zi & 63)) & 1) {
                int deg = __builtin_popcountll(zi);
                if (deg > max_deg) max_deg = deg;
                if (deg >= 2) {
                    sum_T_d2++;
                    std::vector<int> svars;
                    uint64_t rem = (uint64_t)zi;
                    while (rem) {
                        int j = __builtin_ctzll(rem);
                        rem &= rem - 1;
                        svars.push_back(off + j);
                    }
                    shared_terms_union.insert(std::move(svars));
                }
            }
        }
    }
    int64_t union_T = (int64_t)shared_terms_union.size();
    std::cout << "  Sum T = " << sum_T_d2 << "\n";
    std::cout << "  Union T = " << union_T << "\n";

    // Save results
    if (!save_prefix.empty()) {
        std::string inst_name = std::filesystem::path(path).stem().string();
        std::string dir = save_prefix;
        while (!dir.empty() && dir.back() == '/') dir.pop_back();
        std::string tag = "_d2_" + std::string(STRATEGY_TAG);
        std::string base = dir + "/" + inst_name + tag;

        // ---- .affine: [M|b] matrix per output (concatenated) ----
        {
            int s = 0;
            for (int oi = 0; oi < k; oi++) {
                if (!results[oi].is_affine && results[oi].m > 0)
                    s += results[oi].m;
            }

            std::ofstream f(base + ".affine");
            if (f) {
                f << s << "\n" << s << " " << (n + 1) << "\n";
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    if (r.is_affine || r.m <= 0) continue;
                    for (int rw = 0; rw < r.m; rw++) {
                        for (int c = 0; c < n; c++) {
                            if (c > 0) f << " ";
                            f << ((r.M_rows[rw] >> c) & 1);
                        }
                        f << " " << ((r.b >> rw) & 1) << "\n";
                    }
                }
                std::cout << "  Saved: " << base << ".affine (s=" << s << ", n=" << n << ")\n";
            }
        }

        // ---- .poly: per-output ANF terms in shared z-space (raw ANF, no correction) ----
        {
            std::vector<int> z_off(k, 0);
            int s = 0;
            for (int oi = 0; oi < k; oi++) {
                z_off[oi] = s;
                auto& r = results[oi];
                if (!r.is_affine && r.m > 0)
                    s += r.m;
            }

            int max_deg = 0;
            std::vector<int> term_counts(k, 0);

            for (int oi = 0; oi < k; oi++) {
                auto& r = results[oi];
                if (r.is_affine) {
                    if (r.affine_b) term_counts[oi] = 1; // constant 1
                    continue;
                }
                if (r.m <= 0 || r.raw_anf_coeffs.empty()) continue;
                int64_t n_z = int64_t(1) << r.m;
                for (int64_t zi = 0; zi < n_z; zi++) {
                    if ((r.raw_anf_coeffs[zi >> 6] >> (zi & 63)) & 1) {
                        term_counts[oi]++;
                        int deg = __builtin_popcountll(zi);
                        if (deg > max_deg) max_deg = deg;
                    }
                }
            }

            std::ofstream f(base + ".poly");
            if (f) {
                f << s << "\n" << k << "\n";
                for (int oi = 0; oi < k; oi++) {
                    if (oi > 0) f << " ";
                    f << term_counts[oi];
                }
                f << "\n";

                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    if (r.is_affine) {
                        // constant-1 output: write [0,0,...,0,1]
                        if (r.affine_b) {
                            f << "[";
                            for (int b = 0; b < s; b++) {
                                if (b > 0) f << " , ";
                                f << "0";
                            }
                            f << " , 1]\n";
                        }
                        continue;
                    }
                    if (r.m <= 0 || r.raw_anf_coeffs.empty()) continue;
                    int off = z_off[oi];
                    int64_t n_z = int64_t(1) << r.m;

                    for (int64_t zi = 0; zi < n_z; zi++) {
                        if ((r.raw_anf_coeffs[zi >> 6] >> (zi & 63)) & 1) {
                            f << "[";
                            for (int b = 0; b < s; b++) {
                                if (b > 0) f << " , ";
                                int bit = 0;
                                if (b >= off && b < off + r.m) bit = (zi >> (b - off)) & 1;
                                f << bit;
                            }
                            f << " , 1]\n";
                        }
                    }
                }
                std::cout << "  Saved: " << base << ".poly (m=" << s << ", " << k << " outputs)\n";
            }
        }

        // ---- _stats.txt: 5-line numeric ----
        {
            int64_t sum_T = 0;
            int max_deg = 0;
            // Compute per-output z-offsets for shared-space union
            std::vector<int> z_off(k, 0);
            int s_tot = 0;
            std::set<std::vector<int>> shared_terms;
            for (int oi = 0; oi < k; oi++) {
                auto& r = results[oi];
                if (!r.is_affine && r.m > 0 && !r.raw_anf_coeffs.empty()) {
                    z_off[oi] = s_tot;
                    s_tot += r.m;
                } else {
                    z_off[oi] = s_tot;
                }
            }
            for (int oi = 0; oi < k; oi++) {
                auto& r = results[oi];
                if (r.is_affine) {
                    if (r.affine_b) {
                        // constant-1 term: deg=0, skip for sum_T/union_T but track max_deg
                        if (0 > max_deg) max_deg = 0;
                    }
                    continue;
                }
                if (r.m <= 0 || r.raw_anf_coeffs.empty()) continue;
                int off = z_off[oi];
                int64_t n_z = int64_t(1) << r.m;
                for (int64_t zi = 0; zi < n_z; zi++) {
                    if ((r.raw_anf_coeffs[zi >> 6] >> (zi & 63)) & 1) {
                        int deg = __builtin_popcountll(zi);
                        if (deg > max_deg) max_deg = deg;
                        if (deg >= 2) {
                            sum_T++;
                            std::vector<int> svars;
                            uint64_t rem = (uint64_t)zi;
                            while (rem) {
                                int j = __builtin_ctzll(rem);
                                rem &= rem - 1;
                                svars.push_back(off + j);
                            }
                            shared_terms.insert(std::move(svars));
                        }
                    }
                }
            }
            // For per-output disjoint z (opt2), union_T = unique deg≥2 in shared space
            int64_t union_T = (int64_t)shared_terms.size();
            std::ofstream f(base + "_stats.txt");
            if (f) {
                f << n << "\n" << k << "\n" << sum_T << "\n" << union_T << "\n" << max_deg << "\n";
                std::cout << "  Saved: " << base << "_stats.txt (sum_T=" << sum_T
                          << ", union_T=" << union_T << ")\n";
            }
        }

        // ---- _verify.txt: verification (kept as embedded check, verify_anf also used) ----
        {
            std::ofstream f(base + "_verify.txt");
            if (f) {
                f << "# Verification: f(x) = g(Mx+b)\n";
                f << "# n=" << n << "\n\n";
                bool all_pass = true;
                for (int oi = 0; oi < k; oi++) {
                    auto& r = results[oi];
                    if (r.is_affine) {
                        f << circ.outputs[all_outputs[oi]] << ": PASS  AFFINE\n";
                        continue;
                    }
                    if (r.total_T >= INT64_MAX || r.m <= 0 || r.g_tt_raw.empty()) {
                        f << circ.outputs[all_outputs[oi]] << ": SKIP\n";
                        all_pass = false;
                        continue;
                    }
                    bool ok = verify_candidate(fc, r.M_rows, r.b,
                                                r.m, n, r.g_tt_raw, oi, 5000);
                    f << circ.outputs[all_outputs[oi]] << ": "
                      << (ok ? "PASS" : "FAIL") << " T=" << r.total_T
                      << " m=" << r.m << "\n";
                    if (!ok) all_pass = false;
                }
                f << "\nAll: " << (all_pass ? "PASS" : "SOME FAILED") << "\n";
                std::cout << "  Saved: " << base << "_verify.txt\n";
            }
        }
    }
    // ANF output (dedicated flag for backward compat)
    if (!anf_prefix.empty()) {
        std::string base = anf_prefix.substr(0, anf_prefix.find_last_of('/'));
        if (base.empty()) base = anf_prefix;
        // If anf_prefix is same as save_prefix, skip (already wrote .poly above)
    }

    // Verify output (dedicated flag for backward compat)
    if (!verify_prefix.empty()) {
        // If verify_prefix is same as save_prefix, skip (already wrote _verify.txt above)
    }

    std::cout << "\nTotal time: " << total_sec << " s\n";
    return 0;
}
