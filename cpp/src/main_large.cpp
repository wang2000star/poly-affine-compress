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

// Quick verification: check f(x) == g(Mx⊕b) for structured + random x-values.
// Uses single-bit vectors (0, e_0, ..., e_{n-1}) plus n_verify random points.
// For sparse functions, also tests random pairs/triples of kernel vectors.
static bool verify_candidate(
    const FastCircuit& fc,
    const uint32_t* M_rows, uint64_t b,
    int m, int n,
    const std::vector<uint64_t>& g_tt_raw,
    int output_idx,
    int n_verify)
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
        uint64_t eval_buf[1024], out_vec[32];
        eval_batch_fast(fc, in_words, eval_buf, out_vec, (int)fc.out_idx.size());
        uint64_t fx = (out_vec[output_idx] & 1);

        // z = Mx ⊕ b
        uint64_t z = 0;
        for (int r = 0; r < m; r++)
            if (__builtin_popcount(M_rows[r] & (uint32_t)(x & 0xFFFFFFFFULL)) & 1)
                z |= (1ULL << r);
        z ^= b;

        // g(z) from raw truth table
        uint64_t gz = (g_tt_raw[z / 64] >> (z % 64)) & 1;
        return (fx == gz);
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
// Returns total_T, or INT64_MAX if rank(M) < m (invalid) or verification fails.
// g_tt_raw_out: pre-Möbius g(z) truth table (size n_words_g).
// output_idx: which output to evaluate (index into fc.out_idx).
// n_verify: number of random points to verify consistency (0 = skip).
static int64_t evaluate_zspace(
    const FastCircuit& fc,
    const uint32_t* M_rows, uint64_t b,
    int m, int n,
    std::vector<uint64_t>& g_tt_raw_out,
    int n_threads,
    int output_idx,
    int n_verify = 200)
{
    if (m <= 0) return 0;
    if (m > n) return INT64_MAX;

    // Right-inverse P (n×m) such that M × P = I_m
    uint32_t P_rows[32];
    if (!gf2_right_inverse(M_rows, m, n, P_rows)) {
        return INT64_MAX;
    }

    int64_t n_z = int64_t(1) << m;
    int64_t n_words_g = (m < 6) ? 1 : (n_z >> 6);

    g_tt_raw_out.assign(n_words_g, 0);
    int64_t n_batches = (n_z + 63) / 64;

    // Pre-compute P×b constant (XOR this to apply affine shift)
    uint64_t cb[64] = {0};
    for (int j = 0; j < n; j++) {
        if (__builtin_popcount(P_rows[j] & b) & 1)
            cb[j] = ~0ULL;
    }

    // Evaluate each batch of 64 z-values
    for (int64_t batch = 0; batch < n_batches; batch++) {
        uint64_t base_z = batch * 64;

        // Compute x = P×z for this batch (without affine shift)
        uint64_t in_words[64];
        compute_input_words(in_words, n, P_rows, m, base_z);

        // XOR P×b to get x = P×(z⊕b)
        for (int j = 0; j < n; j++)
            in_words[j] ^= cb[j];

        uint64_t eval_buf[1024];
        uint64_t out_vec[32];
        eval_batch_fast(fc, in_words, eval_buf, out_vec, (int)fc.out_idx.size());

        g_tt_raw_out[batch] = out_vec[output_idx];
    }

    // Mask to valid bits (for m < 6, unused bits contain batch garbage)
    if (m < 6) {
        int64_t n_z = int64_t(1) << m;
        uint64_t mask = (n_z >= 64) ? ~0ULL : ((1ULL << n_z) - 1);
        g_tt_raw_out[0] &= mask;
    }
    // Verify consistency before Möbius (raw TT is still intact)
    if (n_verify > 0) {
        if (!verify_candidate(fc, M_rows, b, m, n, g_tt_raw_out, output_idx, n_verify)) {
            g_tt_raw_out.assign(n_words_g, 0);
            return INT64_MAX;
        }
    }

    // Möbius → T count
    moebius_packed(g_tt_raw_out.data(), m);
    return count_T(g_tt_raw_out.data(), m);
}


// ============================================================
//  Result structures
// ============================================================

struct Candidate {
    int m;
    uint32_t M_rows[32];
    uint64_t b;
    int64_t total_T;
    std::vector<uint64_t> g_tt_raw;
};

static bool cmp_candidate(const Candidate& a, const Candidate& b) {
    return a.total_T < b.total_T;
}

// Search for best M,b for a single output
// For n=32 large instances: uses z-space evaluation
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

    // Phase 1: Identity-row candidates (z_i = x_j or z_i = x_j + 1)
    // Each row selects a single input variable
    for (int bit = 0; bit < n; bit++) {
        for (int bval = 0; bval < 2; bval++) {
            uint32_t M[32] = {0};
            M[0] = (1u << bit);
            uint64_t b = bval;

            std::vector<uint64_t> g_tt;
            int64_t T = evaluate_zspace(fc, M, b, 1, n, g_tt, n_threads, output_idx);
            if (T < best.total_T) {
                best.total_T = T;
                best.m = 1;
                best.M_rows[0] = M[0];
                best.b = b;
                best.g_tt_raw = g_tt;
            }
        }
    }

    if (best.total_T == 0) return best;  // T=0 is optimal, skip remaining phases

    // Phase 2: Walsh-guided candidates (for n ≤ 20 for Walsh)
    // For n=32, Walsh is expensive. Use a simple approach instead.
    // Try pairs of identity rows (XOR combinations)
    if (n <= 20) {
        // Compute Walsh for this output's circuit (not available with FastCircuit alone)
        // Skip Walsh for n > 20
    } else {
        // Simple multi-row: combine top identity rows
        // Each new row selects one more input variable
        for (int m_val = 2; m_val <= std::min(max_m, n); m_val++) {
            uint32_t M[32] = {0};
            for (int r = 0; r < m_val && r < n; r++)
                M[r] = (1u << r);
            uint64_t b = 0;

            std::vector<uint64_t> g_tt;
            int64_t T = evaluate_zspace(fc, M, b, m_val, n, g_tt, n_threads, output_idx);
            if (T < best.total_T) {
                best.total_T = T;
                best.m = m_val;
                for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                best.b = b;
                best.g_tt_raw = g_tt;
            }
        }

        // Try XOR pairs of top 8 bits as single-row candidates
        for (int i = 0; i < std::min(8, n); i++) {
            for (int j = i + 1; j < std::min(8, n); j++) {
                uint32_t M[32] = {0};
                M[0] = (1u << i) | (1u << j);
                for (int bval = 0; bval < 2; bval++) {
                    uint64_t b = bval;
                    std::vector<uint64_t> g_tt;
                    int64_t T = evaluate_zspace(fc, M, b, 1, n, g_tt, n_threads, output_idx);
                    if (T < best.total_T) {
                        best.total_T = T;
                        best.m = 1;
                        best.M_rows[0] = M[0];
                        best.b = b;
                        best.g_tt_raw = g_tt;
                    }
                }
            }
        }
    }

    // Short-circuit: if already optimal, skip random search
    if (best.total_T == 0) return best;

    // Phase 3: Random permutation-like candidates
    // Each row selects a random input bit (or small XOR combination)
    // Also try dense rows where each bit is independently set
    for (int c = 0; c < n_random; c++) {
        int m_val = 1 + (rand() % std::min(max_m, n));
        if (m_val > n) m_val = n;

        // Generate random rows
        uint32_t M[32] = {0};
        uint32_t used = 0;
        for (int r = 0; r < m_val; r++) {
            int row_type = rand() % 4;
            if (row_type == 0 && r > 0) {
                // XOR of two random bits
                int b1 = rand() % n;
                int b2 = rand() % n;
                M[r] = (1u << b1) | (1u << b2);
            } else if (row_type == 1 && n <= 32) {
                // Dense row: each bit independently set with 50% probability
                M[r] = (uint32_t)(rand() | (rand() << 16));
                if (n < 32) M[r] &= ((1u << n) - 1);
            } else {
                // Single random bit (prefer unused bits)
                int bit;
                int attempts = 0;
                do {
                    bit = rand() % n;
                    attempts++;
                } while (((used >> bit) & 1) && attempts < 20);
                M[r] = (1u << bit);
                used |= (1u << bit);
            }
        }
        uint64_t b = (uint64_t)(rand() & ((1ULL << m_val) - 1));

        std::vector<uint64_t> g_tt;
        int64_t T = evaluate_zspace(fc, M, b, m_val, n, g_tt, n_threads, output_idx);
        if (T < best.total_T) {
            best.total_T = T;
            best.m = m_val;
            for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
            best.b = b;
            best.g_tt_raw = g_tt;
        }
    }

    // Phase 3b: Special structured candidates with dense random rows
    if (n <= 32) {
        for (int trial = 0; trial < std::min(100, n_random / 2); trial++) {
            int m_val = 1 + (rand() % std::min(max_m, n));
            if (m_val > n) m_val = n;
            uint32_t M[32] = {0};
            for (int r = 0; r < m_val; r++) {
                uint32_t dense = (uint32_t)(rand() | (rand() << 16));
                if (n < 32) dense &= ((1u << n) - 1);
                M[r] = dense;
            }
            uint64_t b = (uint64_t)(rand() & ((1ULL << m_val) - 1));
            std::vector<uint64_t> g_tt;
            int64_t T = evaluate_zspace(fc, M, b, m_val, n, g_tt, n_threads, output_idx);
            if (T < best.total_T) {
                best.total_T = T;
                best.m = m_val;
                for (int r = 0; r < m_val; r++) best.M_rows[r] = M[r];
                best.b = b;
                best.g_tt_raw = g_tt;
            }
        }
    }

    if (best.total_T == 0) return best;

    // Phase 4: Hill climbing from best candidate
    if (best.total_T < INT64_MAX && n_hill_climb > 0) {
        for (int hi = 0; hi < n_hill_climb; hi++) {
            Candidate cur = best;
            bool improved = true;
            int iter = 0;
            while (improved && iter < 500) {
                improved = false;
                iter++;

                // Try changing one row's bits
                for (int r = 0; r < cur.m; r++) {
                    for (int flip = 0; flip < std::min(32, n); flip++) {
                        uint32_t old_row = cur.M_rows[r];
                        cur.M_rows[r] ^= (1u << flip);

                        // Also try toggling b bits
                        for (int bp = 0; bp < 2; bp++) {
                            cur.b ^= ((uint64_t)bp << r);

                            std::vector<uint64_t> g_tt;
                            int64_t T = evaluate_zspace(
                                fc, cur.M_rows, cur.b, cur.m, n, g_tt, n_threads, output_idx);

                            if (T < best.total_T) {
                                best = cur;
                                best.total_T = T;
                                best.g_tt_raw = g_tt;
                                improved = true;
                            }

                            cur.b ^= ((uint64_t)bp << r);
                        }

                        cur.M_rows[r] = old_row;
                    }
                }

                // Try adding one more row
                if (cur.m < std::min(max_m, n)) {
                    int new_r = cur.m;
                    for (int bit = 0; bit < n; bit++) {
                        uint32_t old_M = cur.M_rows[new_r];
                        cur.M_rows[new_r] = (1u << bit);
                        cur.m = new_r + 1;

                        std::vector<uint64_t> g_tt;
                        int64_t T = evaluate_zspace(
                            fc, cur.M_rows, cur.b, cur.m, n, g_tt, n_threads, output_idx);

                        if (T < best.total_T) {
                            best = cur;
                            best.total_T = T;
                            best.g_tt_raw = g_tt;
                            improved = true;
                        }

                        cur.M_rows[new_r] = old_M;
                        cur.m = cur.m;
                    }
                    cur.m = best.m;
                }
            }
        }
    }

    // Fallback: identity transform when no valid m<n transform found
    if (best.total_T >= INT64_MAX && n <= 24) {
        uint32_t M_ident[32] = {0};
        for (int r = 0; r < n; r++) M_ident[r] = (1u << r);
        best.m = n;
        for (int r = 0; r < n; r++) best.M_rows[r] = M_ident[r];
        best.b = 0;
        best.total_T = evaluate_zspace(fc, M_ident, 0, n, n, best.g_tt_raw, n_threads, output_idx, 0);
        if (best.total_T >= 0 && best.total_T < INT64_MAX)
            std::cout << "    [identity fallback: m=n=" << n << " T=" << best.total_T << "]\n";
        else
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
        std::cerr << "  --max-m N      max z variables per output (default 8)\n";
        std::cerr << "  --random N     random candidates per output (default 100)\n";
        std::cerr << "  --hill-climb N hill climb iter per output (default 5)\n";
        std::cerr << "  --walsh-k N    top K Walsh bits (default 20)\n";
        std::cerr << "  --save-results PREFIX  save combined result\n";
        return 1;
    }

    std::cout << std::unitbuf;
    std::string path = argv[1];
    Circuit circ = read_circuit(path);
    int n = circ.n_inputs;

    std::cout << "--- Circuit ---\n";
    std::cout << "  Inputs: " << n << ", Outputs: " << circ.outputs.size()
              << ", Gates: " << circ.stmts.size() << "\n";

    if (n <= 20) {
        std::cerr << "  Warning: n=" << n << " ≤ 20. Use optimize_anf_opt1/opt2 instead "
                  << "(truth table approach is faster for small n).\n";
    }

    int max_m = 8;
    int n_random = 100;
    int n_hill_climb = 5;
    int walsh_k = 20;
    std::string save_prefix;

    for (int a = 2; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "--max-m" && a + 1 < argc) max_m = std::stoi(argv[++a]);
        else if (arg == "--random" && a + 1 < argc) n_random = std::stoi(argv[++a]);
        else if (arg == "--hill-climb" && a + 1 < argc) n_hill_climb = std::stoi(argv[++a]);
        else if (arg == "--walsh-k" && a + 1 < argc) walsh_k = std::stoi(argv[++a]);
        else if (arg == "--save-results" && a + 1 < argc) save_prefix = argv[++a];
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

        // Post-search verification: re-verify best candidate with many random points
        if (cand.total_T < INT64_MAX && cand.total_T >= 0 && cand.m > 0 && !cand.g_tt_raw.empty()) {
            bool verified = verify_candidate(fc, cand.M_rows, cand.b,
                                             cand.m, n, cand.g_tt_raw, oi, 5000);
            if (!verified) {
                // Candidate is invalid — mark as failed
                cand.total_T = INT64_MAX;
                cand.m = 0;
                cand.g_tt_raw.clear();
                std::cout << "    ❌ POST-VERIFY FAILED (5000 tests) — rejecting candidate\n";
            }
        }

        results[oi] = cand;
        if (cand.total_T < INT64_MAX)
            total_sum_T += cand.total_T;

        std::cout << "    m=" << cand.m << " T=" << cand.total_T << " (" << sec << " s)\n";
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
        std::cerr << "ERROR: " << k << " outputs, some have no valid transform. Exiting.\n";
        return 1;
    }

    // Summary
    auto t_end = std::chrono::steady_clock::now();
    double total_sec = std::chrono::duration<double>(t_end - t_start).count();

    std::cout << "\n=== Results ===\n";
    std::cout << "  Sum T = " << total_sum_T << "\n";

    for (int oi = 0; oi < k; oi++) {
        std::cout << "  " << circ.outputs[all_outputs[oi]] << ": T=" << results[oi].total_T
                  << " m=" << results[oi].m << "\n";
    }

    // Compute union T (only for small m)
    int max_m_all = 0;
    for (int oi = 0; oi < k; oi++)
        if (results[oi].m > max_m_all) max_m_all = results[oi].m;

    int64_t union_T = total_sum_T; // default
    if (max_m_all <= 20) {
        // Compute union of ANF monomials across outputs
        int64_t n_words_union = (max_m_all < 6) ? 1 : (int64_t(1) << (max_m_all - 6));
        std::vector<uint64_t> union_data(
            std::max(int64_t(1), n_words_union), 0);

        for (int oi = 0; oi < k; oi++) {
            if (results[oi].m <= 0 || results[oi].g_tt_raw.empty()) continue;
            int m_out = results[oi].m;
            int64_t nw = (m_out < 6) ? 1 : (int64_t(1) << (m_out - 6));
            std::vector<uint64_t> anf(results[oi].g_tt_raw.begin(),
                                      results[oi].g_tt_raw.begin() + nw);
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
        // Trans file
        {
            std::ofstream f(save_prefix + "_trans.poly");
            f << "# large-n transform z_i = M_i x + b_i (per-output)\n";
            for (int oi = 0; oi < k; oi++) {
                auto& r = results[oi];
                if (r.m <= 0) continue;
                f << "# " << circ.outputs[all_outputs[oi]] << " (m=" << r.m
                  << ", T=" << r.total_T << ")\n";
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

        // T file
        {
            std::ofstream f(save_prefix + "_T.poly");
            f << "# large-n per-output T(g) results\n";
            f << "# sum T = " << total_sum_T << "\n";
            for (int oi = 0; oi < k; oi++)
                f << circ.outputs[all_outputs[oi]] << ": T=" << results[oi].total_T
                  << " m=" << results[oi].m << "\n";
        }

        std::cout << "  Saved: " << save_prefix << "_trans.poly\n";
        std::cout << "  Saved: " << save_prefix << "_T.poly\n";
    }

    std::cout << "\nTotal time: " << total_sec << " s\n";
    return 0;
}
