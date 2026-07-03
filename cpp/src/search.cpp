#include "search.h"
#include "affine.h"
#include "anf.h"
#include "gf2.h"
#include "io.h"
#include "verify.h"
#include "walsh.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <thread>

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
            uint32_t b = 0;
            for (int j = 0; j < m; j++) {
                if (rng() & 1) b |= (1u << j);
            }
            candidates.emplace_back(rows, b);
        }
        return candidates;
    }
};

// ============================================================
//  Hill climbing
// ============================================================

MbCandidate hill_climb(const TruthTable& tt, const MbCandidate& start,
                       int n, int n_threads)
{
    MbCandidate best = start;
    int m = start.m;

    auto eval = [&](uint32_t M[32], uint32_t b) -> MbCandidate {
        return evaluate_Mb(tt, M, b, m, n_threads);
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
                MbCandidate cand = eval(cur.M_rows, cur.b);
                if (cand.total_T < best.total_T && cand.total_T != INT64_MAX)
                    best = cand;
                if (cand.total_T < cur.total_T && cand.total_T != INT64_MAX) {
                    cur = cand;
                    improved = true;
                } else {
                    cur.M_rows[r] ^= (1u << c);
                }
            }
        }

        for (int r = 0; r < m && !improved; r++) {
            cur.b ^= (1u << r);
            MbCandidate cand = eval(cur.M_rows, cur.b);
            if (cand.total_T < best.total_T && cand.total_T != INT64_MAX)
                best = cand;
            if (cand.total_T < cur.total_T && cand.total_T != INT64_MAX) {
                cur = cand;
                improved = true;
            } else {
                cur.b ^= (1u << r);
            }
        }
    } while (improved && pass < MAX_PASSES);

    return best;
}

// ============================================================
//  run_search — full pipeline (Phase 1..6)
// ============================================================

void run_search(const TruthTable& tt, const Circuit& circ,
                const std::vector<int>& output_indices,
                const SearchParams& params)
{
    int n = tt.n;
    auto t0 = std::chrono::steady_clock::now();

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
    if (!params.anf_out_prefix.empty() && n <= 16) {
        save_raw_anf(tt_copy, circ, output_indices, params.anf_out_prefix + "_expr.poly");
        save_raw_T(tt_copy, circ, output_indices, params.anf_out_prefix + "_T.poly");
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
    std::vector<MbCandidate> results;

    // 4a-c: Walsh-guided + random (n ≤ 20 full search; 21-32 permutation only)
    if (n <= 20) {
        CandidateGenerator gen(n, params.max_m, walsh);

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
        for (auto& [rows, b] : random_cands) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            results.push_back(evaluate_Mb(tt_copy, M, b, m, params.n_threads));
        }

        auto t_search0 = std::chrono::steady_clock::now();
        double walsh_time = std::chrono::duration<double>(t_search0 - t1).count();
        std::cout << "  Search time: " << walsh_time << " s\n";
    } else if (n <= 32 && !walsh.empty()) {
        CandidateGenerator gen(n, params.max_m, walsh);

        std::cout << "  4a: Single-row Walsh (" << (params.walsh_single_top * 2) << " candidates)\n";
        auto single_rows = gen.gen_walsh_single_rows(params.walsh_single_top);
        for (auto& [row, b] : single_rows) {
            uint32_t M[32] = {0};
            M[0] = row;
            results.push_back(evaluate_Mb(tt_copy, M, b, 1, params.n_threads));
        }

        std::cout << "  4b: Multi-row Walsh (progressive m=1.." << std::min(20, n)
                  << " + XOR pairs)\n";
        auto multi_rows = gen.gen_multi_row(std::min(20, n));
        for (auto& [rows, b] : multi_rows) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            results.push_back(evaluate_Mb(tt_copy, M, b, m, params.n_threads));
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

            uint32_t b = 0;
            for (int r = 0; r < n; r++)
                if (rng_n32() & 1) b |= (1u << r);

            if (gf2_rank(M, n, n) < n) continue;

            auto t_eval0 = std::chrono::steady_clock::now();
            auto cand = evaluate_Mb(tt_copy, M, b, n, params.n_threads);
            auto t_eval1 = std::chrono::steady_clock::now();
            results.push_back(cand);
            generated++;

            double eval_time = std::chrono::duration<double>(t_eval1 - t_eval0).count();
            double compression = (double)raw.sum_T / std::max(int64_t(1), cand.total_T);
            std::cout << "    n32 #" << generated << "/" << params.n32_random
                      << " (attempt " << attempts << ")"
                      << ": T=" << cand.total_T
                      << " (compression " << compression << "×)"
                      << " time=" << eval_time << "s\n";
        }
        auto t_n32_1 = std::chrono::steady_clock::now();
        std::cout << "    n32 random time: " << std::chrono::duration<double>(t_n32_1 - t_n32_0).count() << " s\n";
    }

    // 4e: Hill climbing from top candidates (m=n only)
    if (params.n_hill_climb > 0) {
        std::sort(results.begin(), results.end(),
            [](auto& a, auto& b) { return a.total_T < b.total_T; });

        int n_climb = std::min(params.n_hill_climb, (int)results.size());
        std::cout << "  4e: Hill climbing from top " << n_climb << " candidates (+ identity)\n";

        auto t_hc0 = std::chrono::steady_clock::now();
        for (int ci = 0; ci < n_climb; ci++) {
            auto& base = results[ci];
            if (base.total_T >= INT64_MAX) continue;
            if (base.m != n) continue;

            auto improved = hill_climb(tt_copy, base, n, params.n_threads);
            if (improved.total_T < base.total_T) {
                results.push_back(improved);
                std::cout << "    hill climb #" << ci << ": T=" << base.total_T << " -> T=" << improved.total_T << "\n";
            }
        }

        // Hill climb from identity
        uint32_t ident[32] = {0};
        for (int i = 0; i < n; i++) ident[i] = (1u << i);
        MbCandidate ident_cand = evaluate_Mb(tt_copy, ident, 0, n, params.n_threads);
        auto ident_improved = hill_climb(tt_copy, ident_cand, n, params.n_threads);
        if (ident_improved.total_T < ident_cand.total_T) {
            results.push_back(ident_improved);
            std::cout << "    hill climb from identity: T=" << ident_cand.total_T << " -> T=" << ident_improved.total_T << "\n";
        }

        auto t_hc1 = std::chrono::steady_clock::now();
        std::cout << "    Hill climb time: " << std::chrono::duration<double>(t_hc1 - t_hc0).count() << " s\n";
    }

    // Phase 5: Report
    std::sort(results.begin(), results.end(),
        [](auto& a, auto& b) { return a.total_T < b.total_T; });

    std::cout << "\n========================================\n";
    std::cout << "  Results (" << results.size() << " candidates evaluated)\n";
    std::cout << "  Raw ANF Sum T = " << raw.sum_T << "\n";
    std::cout << "========================================\n";

    int n_report = std::min(20, (int)results.size());
    for (int j = 0; j < n_report; j++) {
        auto& r = results[j];
        double compression = (double)raw.sum_T / std::max(int64_t(1), r.total_T);
        std::cout << "  " << (j + 1) << ". m=" << r.m << " T=" << r.total_T
                  << " (compression " << compression << "×) M=[";
        for (int row = 0; row < r.m; row++) {
            if (row > 0) std::cout << "; ";
            for (int i = 0; i < std::min(n, 8); i++)
                std::cout << ((r.M_rows[row] >> i) & 1);
            std::cout << "..";
        }
        std::cout << "]\n";
        if (r.total_T < raw.sum_T && r.m <= 32) {
            std::cout << "    Per-output T:";
            for (int oi = 0; oi < tt.n_outputs; oi++)
                std::cout << " " << circ.outputs[output_indices[oi]] << "=" << r.per_output_T[oi];
            std::cout << "\n";
        }
    }

    // Phase 6: Verify best candidate, save results
    if (!results.empty() && results[0].total_T < INT64_MAX) {
        auto& best = results[0];
        std::cout << "\nBest candidate: m=" << best.m << ", T=" << best.total_T << "\n";
    }

    double total_time = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    // Phase 6: verification + file saving
    if (!results.empty() && results[0].total_T < INT64_MAX) {
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

        // Re-evaluate with save_g_tt=true
        bool all_verified = true;
        MbCandidate verified_cand;
        verified_cand.m = best.m;
        verified_cand.b = best.b;
        for (int i = 0; i < best.m; i++) verified_cand.M_rows[i] = best.M_rows[i];
        verified_cand.total_T = best.total_T;
        verified_cand.per_output_T = best.per_output_T;

        if (best.m == n) {
            std::cout << "Verifying best candidate (5000 random tests per output)...\n";
            verified_cand = evaluate_Mb(tt_copy, best.M_rows, best.b, best.m,
                                         params.n_threads, true);
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                bool ok = verify_transform(tt_copy, verified_cand.g_tt_raw[oi], oi,
                    best.M_rows, best.b, best.m, n, 5000);
                std::cout << "  " << circ.outputs[output_indices[oi]]
                          << (ok ? " ✅ Verified (5000 tests)" : " ❌ FAILED") << "\n";
                if (!ok) all_verified = false;
            }
        } else {
            std::cout << "Skipping verification (non-bijective m=" << best.m << " < n=" << n << ")\n";
        }

        if (all_verified)
            std::cout << "✅ All outputs verified!\n";

        // Save results to files
        if (!params.save_results_prefix.empty()) {
            std::string prefix = params.save_results_prefix;

            save_trans(best, circ, n, prefix + "_trans.poly");

            if (!verified_cand.g_tt_raw.empty()) {
                save_opt_expr(verified_cand.g_tt_raw, circ, output_indices, best.m,
                             prefix + "_expr.poly");
                save_opt_T(verified_cand.g_tt_raw, circ, output_indices, best.m,
                          prefix + "_T.poly");
            }

            save_verify(tt_copy, verified_cand.g_tt_raw, best, n, output_indices, circ,
                       prefix + "_verify.txt");
        }
    }

    total_time = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::cout << "\nTotal time: " << total_time << " s\n";
}
