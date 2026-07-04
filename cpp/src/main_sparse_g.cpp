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

// Strategy tag for output file naming (prevents confusion between different search strategies)
static const char* STRATEGY_TAG = "sg";

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
//  Verification with generalized c-correction
//  Checks: f(x) = g(Mx⊕b) ⊕ ⟨c,x⟩ ⊕ d
// ============================================================

static bool verify_candidate_generalized(
    const FastCircuit& fc,
    const uint32_t* M_rows, uint64_t b,
    const std::vector<uint64_t>& c_vec, uint64_t d,
    int m, int n,
    const std::vector<uint64_t>& g_tt_raw,
    int output_idx,
    int n_verify)
{
    if (n_verify <= 0) return true;

    // Pack c_vec into bitmask for fast inner product
    uint64_t c_mask = 0;
    for (int i = 0; i < n && i < 64; i++)
        if (c_vec[i]) c_mask |= (1ULL << i);

    auto test_x = [&](uint64_t x) -> bool {
        // f(x) via circuit evaluation
        uint64_t in_words[64];
        for (int j = 0; j < n; j++)
            in_words[j] = ((x >> j) & 1) ? ~0ULL : 0;
        uint64_t eval_buf[1024], out_vec[32];
        eval_batch_fast(fc, in_words, eval_buf, out_vec, (int)fc.out_idx.size());
        uint64_t fx = (out_vec[output_idx] & 1);

        // z = Mx ⊕ b
        uint64_t z = 0;
        for (int r = 0; r < m; r++)
            if (__builtin_popcount(M_rows[r] & (uint32_t)(x & 0xFFFFFFFFULL)) & 1)
                z |= (1ULL << r);
        z ^= b;

        // g(z) from truth table
        uint64_t gz = (g_tt_raw[z / 64] >> (z % 64)) & 1;

        // ⟨c,x⟩ ⊕ d
        uint64_t cx = (__builtin_popcountll(c_mask & x) & 1);

        // Check: f(x) = g(z) ⊕ ⟨c,x⟩ ⊕ d
        // i.e., fx ⊕ gz ⊕ cx ⊕ d == 0
        return ((fx ^ gz ^ cx ^ (d & 1)) == 0);
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

// Score pool rows by T(g) for m=1 with c-correction
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
            uint64_t eval_buf[1024], out_vec[32];
            eval_batch_fast(fc, in_words, eval_buf, out_vec, (int)fc.out_idx.size());
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

        // Remove linear term (z_0) via c-correction: set c such that
        // the linear term cancels. For m=1, M⁺ is n×1 and M⁺ᵀ·c = L
        // where L = ANF coefficient of z_0 (0 or 1).
        // If L=1, we can cancel it by setting c = M⁺ (any row of M⁺ = P_rows).
        int has_linear = (int)((anf[0] >> 1) & 1);
        int T_corrected = T_raw - has_linear;
        // Also optimize constant via d: T_corrected = min(T_corrected, T_corrected - (has_const?1:0) + 1)
        // Actually: d selects constant = 0 or 1. If T_raw has constant=1, removing it via d
        // reduces T by 1; if constant=0, adding it via d increases T by 1.
        int has_const = (int)(anf[0] & 1);
        int T_opt = std::min(T_corrected, T_corrected - has_const + (1 - has_const));

        pool[ri].score = T_opt;

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
                uint64_t eval_buf[1024], out_vec[32];
                eval_batch_fast(fc, in_words, eval_buf, out_vec, (int)fc.out_idx.size());
                if (out_vec[output_idx] & 1)
                    g_tt_b1[0] |= (1ULL << z);
            }

            std::vector<uint64_t> anf_b1 = g_tt_b1;
            moebius_packed(anf_b1.data(), 1);
            int T_raw_b1 = 0;
            for (size_t w = 0; w < anf_b1.size(); w++)
                T_raw_b1 += __builtin_popcountll(anf_b1[w]);

            int has_linear_b1 = (int)((anf_b1[0] >> 1) & 1);
            int T_corrected_b1 = T_raw_b1 - has_linear_b1;
            int has_const_b1 = (int)(anf_b1[0] & 1);
            int T_opt_b1 = std::min(T_corrected_b1, T_corrected_b1 - has_const_b1 + (1 - has_const_b1));

            if (T_opt_b1 < pool[ri].score) {
                pool[ri].score = T_opt_b1;
                // Don't update g_tt — caller will re-evaluate with chosen b
            }
        }
    }
}

// ============================================================
//  Evaluate M,b with c-correction
// ============================================================

// Structure to hold evaluation result
struct EvalResult {
    int m;
    uint32_t M_rows[32];
    uint64_t b;
    int64_t total_T;          // T(g) after c-correction
    int64_t T_raw;            // T(g) before c-correction
    std::vector<uint64_t> g_tt_raw;   // raw g truth table (pre-Möbius, post c-correction)
    std::vector<uint64_t> anf_coeffs;  // ANF coefficients
    std::vector<uint64_t> c_vec;       // output-side linear coefficients (size n)
    uint64_t d;               // output-side constant
    bool valid;
};

// Evaluate a single M,b candidate with c-correction
// The c-correction removes linear terms from g(z) by solving:
//   M⁺ᵀ · c = L  where L = linear term coefficients of g_raw
// Result: g_corrected(z) = g_raw(z) ⊕ ⟨L, z⊕b⟩
// Verification: f(x) = g_corrected(Mx⊕b) ⊕ ⟨c,x⟩ ⊕ d
static EvalResult evaluate_Mb_ccorrect(
    const FastCircuit& fc,
    const uint32_t* M_rows, uint64_t b,
    int m, int n,
    int n_threads,
    int output_idx,
    int n_verify = 200)
{
    EvalResult r;
    r.valid = false;
    r.m = m;
    r.total_T = INT64_MAX;
    r.T_raw = INT64_MAX;
    r.d = 0;
    r.c_vec.assign(n, 0);

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

        uint64_t eval_buf[1024], out_vec[32];
        eval_batch_fast(fc, in_words, eval_buf, out_vec, (int)fc.out_idx.size());
        g_tt_raw[batch] = out_vec[output_idx];
    }

    if (m < 6) {
        uint64_t mask = (n_z >= 64) ? ~0ULL : ((1ULL << n_z) - 1);
        g_tt_raw[0] &= mask;
    }

    // Step 1: Compute c from linear terms of g_raw's ANF
    // Möbius to get ANF of raw g
    std::vector<uint64_t> raw_anf = g_tt_raw;
    moebius_packed(raw_anf.data(), m);

    // Extract which linear terms exist in g_raw (L[r] = coefficient of z_r)
    uint64_t L_mask = 0;
    for (int r = 0; r < m; r++) {
        int bit_pos = 1 << r;
        if ((raw_anf[bit_pos / 64] >> (bit_pos % 64)) & 1)
            L_mask |= (1ULL << r);
    }

    // Solve M⁺ᵀ · c = L for c
    // M⁺ᵀ[r][j] = bit r of P_rows[j]
    uint32_t selected_cols[32];
    int n_selected = 0;
    for (int j = 0; j < n && n_selected < m; j++) {
        uint32_t test_rows[32];
        for (int k = 0; k < n_selected; k++)
            test_rows[k] = P_rows[selected_cols[k]];
        test_rows[n_selected] = P_rows[j];
        uint32_t work[32];
        for (int k = 0; k <= n_selected; k++) work[k] = test_rows[k];
        int rx = gf2_rank(work, n_selected + 1, m);
        if (rx > n_selected)
            selected_cols[n_selected++] = j;
    }

    if (L_mask != 0 && n_selected == m) {
        uint32_t A[32];
        for (int s = 0; s < m; s++) A[s] = P_rows[selected_cols[s]];
        uint32_t B[32] = {0};
        for (int r = 0; r < m; r++)
            for (int s = 0; s < m; s++)
                if ((A[s] >> r) & 1) B[r] |= (1u << s);
        uint32_t Bw[32];
        uint64_t Lw[32];
        for (int r = 0; r < m; r++) {
            Bw[r] = B[r];
            Lw[r] = (L_mask >> r) & 1;
        }
        for (int col = 0; col < m; col++) {
            int pivot = -1;
            for (int r = col; r < m; r++)
                if ((Bw[r] >> col) & 1) { pivot = r; break; }
            if (pivot < 0) continue;
            std::swap(Bw[col], Bw[pivot]);
            std::swap(Lw[col], Lw[pivot]);
            for (int r = 0; r < m; r++)
                if (r != col && ((Bw[r] >> col) & 1)) {
                    Bw[r] ^= Bw[col];
                    Lw[r] ^= Lw[col];
                }
        }
        for (int s = 0; s < m; s++)
            if (Lw[s]) r.c_vec[selected_cols[s]] = 1;
    }

    // Step 2: Apply c-correction to g_raw:
    // g_corrected(z) = g_raw(z) ⊕ ⟨L, z⊕b⟩
    // This gives g_corrected with NO linear ANF terms.
    std::vector<uint64_t> g_corrected = g_tt_raw;
    if (L_mask != 0) {
        for (int64_t z = 0; z < n_z; z++) {
            uint64_t zb = (uint64_t)z ^ b;
            if (__builtin_popcountll(zb & L_mask) & 1)
                g_corrected[z / 64] ^= (1ULL << (z % 64));
        }
    }

    // Step 3: Verify g_corrected with c-correction
    // Check: f(x) = g_corrected(Mx⊕b) ⊕ ⟨c, x⟩ ⊕ d
    if (n_verify > 0) {
        if (!verify_candidate_generalized(fc, M_rows, b, r.c_vec, 0,
                                           m, n, g_corrected, output_idx, n_verify)) {
            if (!verify_candidate_generalized(fc, M_rows, b, r.c_vec, 1,
                                               m, n, g_corrected, output_idx, n_verify)) {
                return r;
            }
            r.d = 1;
        }
    }

    // Step 4: Compute ANF and T of corrected g
    std::vector<uint64_t> g_anf = g_corrected;
    moebius_packed(g_anf.data(), m);

    int64_t T_corrected = 0;
    for (size_t w = 0; w < g_anf.size(); w++)
        T_corrected += __builtin_popcountll(g_anf[w]);

    int has_const = (int)(g_anf[0] & 1);
    int64_t T_opt = std::min(T_corrected, T_corrected - has_const + (1 - has_const));

    int64_t T_raw_val = 0;
    for (size_t w = 0; w < raw_anf.size(); w++)
        T_raw_val += __builtin_popcountll(raw_anf[w]);

    r.total_T = T_opt;
    r.T_raw = T_raw_val;
    r.g_tt_raw = g_corrected;
    r.anf_coeffs = g_anf;
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
    int64_t T_raw;
    std::vector<uint64_t> g_tt_raw;   // corrected raw TT
    std::vector<uint64_t> anf_coeffs;  // ANF of corrected g
    std::vector<uint64_t> c_vec;       // output-side linear coefficients
    uint64_t d;               // output-side constant

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
// Uses n+1 circuit evaluations + random verification via batch evaluation.
// Returns true with affine_mask (bit i = a_i) and affine_b properly set.
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

    if (n > 63) return false;

    srand(12345);
    for (int t = 0; t < 30; t++) {
        uint64_t x = 0;
        for (int j = 0; j < n; j++)
            if (rand() & 1) x |= (1ULL << j);

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
    best.T_raw = INT64_MAX;
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
            best.T_raw = 0;
            return best;
        }
    }

    // Build 3n-row pool
    std::vector<PoolRow> pool = build_pool(n);

    // Score rows by m=1 T(g) with c-correction
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

            auto er = evaluate_Mb_ccorrect(fc, M, b, 1, n, n_threads, output_idx);
            if (er.valid && er.total_T < best.total_T) {
                best.total_T = er.total_T;
                best.T_raw = er.T_raw;
                best.m = 1;
                best.M_rows[0] = M[0];
                best.b = b;
                best.g_tt_raw = er.g_tt_raw;
                best.c_vec = er.c_vec;
                best.d = er.d;
                best.anf_coeffs = er.anf_coeffs;
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
            auto er = evaluate_Mb_ccorrect(fc, M, 0, m_val, n, n_threads, output_idx);
            if (er.valid && er.total_T < best.total_T) {
                best.total_T = er.total_T;
                best.T_raw = er.T_raw;
                best.m = m_val;
                for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                best.b = 0;
                best.g_tt_raw = er.g_tt_raw;
                best.c_vec = er.c_vec;
                best.d = er.d;
                best.anf_coeffs = er.anf_coeffs;
            }

            // Also try b with each row's complement bit
            uint64_t b_compl = 0;
            for (int r = 0; r < m_val; r++)
                if (pool[sorted_idx[comb[r]]].type == 1) b_compl |= (1ULL << r);

            if (b_compl != 0) {
                auto er2 = evaluate_Mb_ccorrect(fc, M, b_compl, m_val, n, n_threads, output_idx);
                if (er2.valid && er2.total_T < best.total_T) {
                    best.total_T = er2.total_T;
                    best.T_raw = er2.T_raw;
                    best.m = m_val;
                    for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                    best.b = b_compl;
                    best.g_tt_raw = er2.g_tt_raw;
                    best.anf_coeffs = er2.anf_coeffs;
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
                    auto er = evaluate_Mb_ccorrect(fc, M, 0, m_val, n, n_threads, output_idx);
                    if (er.valid && er.total_T < best.total_T) {
                        best.total_T = er.total_T;
                        best.T_raw = er.T_raw;
                        best.m = m_val;
                        for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                        best.b = 0;
                        best.g_tt_raw = er.g_tt_raw;
                        best.c_vec = er.c_vec;
                        best.d = er.d;
                        best.anf_coeffs = er.anf_coeffs;
                    }
                }
                // Try complement b: set b bits for complement-type pool rows
                {
                    uint64_t b = 0;
                    for (int r = 0; r < m_val; r++)
                        if (pool[picked_idx[r]].type == 1) b |= (1ULL << r);
                    if (b != 0) {
                        auto er = evaluate_Mb_ccorrect(fc, M, b, m_val, n, n_threads, output_idx);
                        if (er.valid && er.total_T < best.total_T) {
                            best.total_T = er.total_T;
                            best.T_raw = er.T_raw;
                            best.m = m_val;
                            for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                            best.b = b;
                            best.g_tt_raw = er.g_tt_raw;
                            best.c_vec = er.c_vec;
                            best.d = er.d;
                            best.anf_coeffs = er.anf_coeffs;
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

                            auto er = evaluate_Mb_ccorrect(
                                fc, cur.M_rows, cur.b, cur.m, n, n_threads, output_idx, 0);

                            if (er.valid && er.total_T < best.total_T) {
                                best.total_T = er.total_T;
                                best.T_raw = er.T_raw;
                                best.m = cur.m;
                                for (int rr = 0; rr < cur.m; rr++) best.M_rows[rr] = cur.M_rows[rr];
                                best.b = cur.b;
                                best.g_tt_raw = er.g_tt_raw;
                                best.c_vec = er.c_vec;
                                best.d = er.d;
                                best.anf_coeffs = er.anf_coeffs;
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

                            auto er = evaluate_Mb_ccorrect(
                                fc, cur.M_rows, b_try, cur.m, n, n_threads, output_idx, 0);

                            if (er.valid && er.total_T < best.total_T) {
                                best.total_T = er.total_T;
                                best.T_raw = er.T_raw;
                                best.m = cur.m;
                                for (int rr = 0; rr < cur.m; rr++) best.M_rows[rr] = cur.M_rows[rr];
                                best.b = b_try;
                                best.g_tt_raw = er.g_tt_raw;
                                best.c_vec = er.c_vec;
                                best.d = er.d;
                                best.anf_coeffs = er.anf_coeffs;
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

    // Fallback: identity transform (only for n ≤ 20, where truth table fits)
    if (best.total_T >= INT64_MAX && n <= 20) {
        uint32_t M_ident[32] = {0};
        for (int r = 0; r < n; r++) M_ident[r] = (1u << r);
        auto er = evaluate_Mb_ccorrect(fc, M_ident, 0, n, n, n_threads, output_idx, 0);
        if (er.valid && er.total_T < best.total_T) {
            best.total_T = er.total_T;
            best.T_raw = er.T_raw;
            best.m = n;
            for (int r = 0; r < n; r++) best.M_rows[r] = M_ident[r];
            best.b = 0;
            best.g_tt_raw = er.g_tt_raw;
            best.c_vec = er.c_vec;
            best.d = er.d;
            best.anf_coeffs = er.anf_coeffs;
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
    std::string path = argv[1];
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
        else if (arg == "--save-results" && a + 1 < argc) save_prefix = argv[++a];
        else if (arg == "--anf-out" && a + 1 < argc) anf_prefix = argv[++a];
        else if (arg == "--verify-out" && a + 1 < argc) verify_prefix = argv[++a];
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
            bool verified = verify_candidate_generalized(fc, cand.M_rows, cand.b,
                                             cand.c_vec, cand.d, cand.m, n, cand.g_tt_raw, oi, 5000);
            if (!verified) {
                cand.total_T = INT64_MAX;
                cand.m = 0;
                cand.g_tt_raw.clear();
                std::cout << "    ❌ POST-VERIFY FAILED (5000 tests) — rejecting candidate\n";
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
                      << " (raw T=" << cand.T_raw << ") [" << sec << " s]\n";
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
                      << " (raw=" << results[oi].T_raw << ") m=" << results[oi].m << "\n";
        }
    }

    // Union T (skip affine outputs)
    int max_m_all = 0;
    for (int oi = 0; oi < k; oi++)
        if (!results[oi].is_affine && results[oi].m > max_m_all) max_m_all = results[oi].m;

    int64_t union_T = total_sum_T;
    if (max_m_all <= 20) {
        int64_t n_words_union = (max_m_all < 6) ? 1 : (int64_t(1) << (max_m_all - 6));
        std::vector<uint64_t> union_data(std::max(int64_t(1), n_words_union), 0);

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
        std::cout << "  Union T = " << union_T << "\n";
    }

    // Save results
    if (!save_prefix.empty()) {
        {
            std::ofstream f(save_prefix + "_" + STRATEGY_TAG + "_trans.poly");
            f << "# sparse-g enumeration: z_i = M_i x + b_i (per-output)\n";
            for (int oi = 0; oi < k; oi++) {
                auto& r = results[oi];
                if (r.is_affine) {
                    f << "# " << circ.outputs[all_outputs[oi]] << ": AFFINE (constant)\n";
                    continue;
                }
                if (r.m <= 0) continue;
                f << "# " << circ.outputs[all_outputs[oi]] << " (m=" << r.m
                  << ", T=" << r.total_T << ", raw_T=" << r.T_raw << ")\n";
                for (int row = 0; row < r.m; row++) {
                    f << "z_" << oi << "_" << row << " = ";
                    bool first = true;
                    for (int i = 0; i < n; i++) {
                        if ((r.M_rows[row] >> i) & 1) {
                            if (!first) f << " + ";
                            f << circ.inputs[i];
                            first = false;
                        }
                    }
                    if ((r.b >> row) & 1) {
                        if (!first) f << " + ";
                        f << "1";
                    } else if (first) {
                        f << "0";
                    }
                    f << "\n";
                }
            }
        }

        {
            std::ofstream f(save_prefix + "_" + STRATEGY_TAG + "_T.txt");
            f << "# T(g) results for " << path << "\n";
            f << "# direction=1, strategy=opt2\n";
            int a1 = 0, a2 = 0, a3 = 0;
            for (int oi = 0; oi < k; oi++) {
                if (results[oi].is_affine && results[oi].affine_mask == 0 && results[oi].affine_b == 0)
                    a1++;
                else if (results[oi].is_affine)
                    a2++;
                else
                    a3++;
            }
            f << "# a1=" << a1 << ", a2=" << a2 << ", a3=" << a3 << "\n";
            f << "sum_T = " << total_sum_T << "\n";
            f << "union_T = " << union_T << "\n";
            f << "---\n";
            for (int oi = 0; oi < k; oi++) {
                auto& r = results[oi];
                if (r.is_affine) {
                    if (r.affine_mask == 0 && r.affine_b == 0)
                        f << circ.outputs[all_outputs[oi]] << ": AFFINE  (constant 0)\n";
                    else if (r.affine_mask == 0 && r.affine_b == 1)
                        f << circ.outputs[all_outputs[oi]] << ": AFFINE  (constant 1)\n";
                    else
                        f << circ.outputs[all_outputs[oi]] << ": AFFINE  m=1  nonlinear_terms=0\n";
                } else {
                    // Count nonlinear terms (degree ≥ 2)
                    int nl_terms = 0;
                    if (!r.anf_coeffs.empty() && r.m > 0) {
                        int64_t n_z = int64_t(1) << r.m;
                        for (int64_t zi = 3; zi < n_z; zi++) {
                            // degree ≥ 2 means popcount ≥ 2
                            if (__builtin_popcountll(zi) >= 2 && ((r.anf_coeffs[zi >> 6] >> (zi & 63)) & 1))
                                nl_terms++;
                        }
                    }
                    f << circ.outputs[all_outputs[oi]] << ":  T=" << r.total_T
                      << "  m=" << r.m << "  nonlinear_terms=" << nl_terms << "\n";
                }
            }
        }

        std::cout << "  Saved: " << save_prefix << "_" << STRATEGY_TAG << "_trans.poly\n";
        std::cout << "  Saved: " << save_prefix << "_" << STRATEGY_TAG << "_T.txt\n";
    }

    // ANF output
    if (!anf_prefix.empty()) {
        std::ofstream f(anf_prefix + "_" + STRATEGY_TAG + "_anf.poly");
        f << "# ANF for " << path << " (n=" << n << ", k=" << k << ")\n";
        for (int oi = 0; oi < k; oi++) {
            auto& r = results[oi];
            if (r.is_affine) {
                f << "# " << circ.outputs[all_outputs[oi]] << "  (AFFINE)\n";
                continue;
            }
            if (r.m <= 0 || r.anf_coeffs.empty()) continue;
            f << "# " << circ.outputs[all_outputs[oi]] << "\n";
            int64_t n_z = int64_t(1) << r.m;
            for (int64_t zi = 1; zi < n_z; zi++) {
                if ((r.anf_coeffs[zi / 64] >> (zi % 64)) & 1) {
                    f << "  ";
                    bool first = true;
                    for (int b = 0; b < r.m; b++) {
                        if ((zi >> b) & 1) {
                            if (!first) f << " * ";
                            f << "z_" << oi << "_" << b;
                            first = false;
                        }
                    }
                    f << "\n";
                }
            }
        }
        std::cout << "  Saved: " << anf_prefix << "_" << STRATEGY_TAG << "_anf.poly\n";
    }

    // Verify output
    if (!verify_prefix.empty()) {
        std::ofstream f(verify_prefix + "_" + STRATEGY_TAG + "_verify.txt");
        f << "# Verification results for " << path << "\n\n";
        bool all_pass = true;
        for (int oi = 0; oi < k; oi++) {
            auto& r = results[oi];
            if (r.is_affine) {
                f << circ.outputs[all_outputs[oi]] << ": PASS  AFFINE"
                  << "  y=";
                bool first = true;
                if (r.affine_b) { f << "1"; first = false; }
                for (int i = 0; i < n && i < 32; i++) {
                    if (r.affine_mask & (1u << i)) {
                        if (!first) f << " ⊕ ";
                        f << "x_" << i;
                        first = false;
                    }
                }
                if (first) f << "0";
                f << "\n";
                continue;
            }
            if (r.total_T >= INT64_MAX || r.m <= 0 || r.g_tt_raw.empty()) {
                f << circ.outputs[all_outputs[oi]] << ": SKIP\n";
                all_pass = false;
                continue;
            }
            bool ok = verify_candidate_generalized(fc, r.M_rows, r.b, r.c_vec, r.d,
                                        r.m, n, r.g_tt_raw, oi, 5000);
            f << circ.outputs[all_outputs[oi]] << ": "
              << (ok ? "PASS" : "FAIL") << " T=" << r.total_T
              << " m=" << r.m << "\n";
            if (!ok) all_pass = false;
        }
        f << "\nAll: " << (all_pass ? "PASS" : "SOME FAILED") << "\n";
        std::cout << "  Saved: " << verify_prefix << "_" << STRATEGY_TAG << "_verify.txt\n";
    }

    std::cout << "\nTotal time: " << total_sec << " s\n";
    return 0;
}
