/**
 * optimize_anf_opt2 — Per-output affine transform ANF optimization.
 *
 * Strategy: For each output independently, find the best affine transform
 * z_i = M_i x + b_i that minimizes T(g_i) where f_i(x) = g_i(z_i).
 * Then merge all per-output transforms into a block-diagonal M and
 * combined g(z) expression.
 *
 * Usage:
 *   optimize_anf_opt2 <circuit.txt> [options]
 *
 * Options:
 *   --max-m N     max z variables per output (default 12)
 *   --walsh-k N   top K Walsh bits (default 30)
 *   --random N    random candidates per output (default 40)
 *   --hill-climb N hill climb iterations per output (default 10)
 *   --save-results PREFIX  save combined result to PREFIX_opt2_*
 *
 * Output files (with --save-results PREFIX):
 *   PREFIX_opt2_trans.poly   — block-diagonal transform
 *   PREFIX_opt2_expr.poly    — combined ANF expression
 *   PREFIX_opt2_T.poly       — per-output and total term counts
 */
#include "circuit.h"
#include "truth_table.h"
#include "affine.h"
#include "anf.h"
#include "moebius.h"
#include "io.h"
#include "search.h"
#include "walsh.h"
#include "gf2.h"
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <thread>
#include <algorithm>
#include <map>
#include <set>
#include <filesystem>

struct Opt2Result {
    int m;                              // total z variables
    int n;                              // input variables
    std::vector<uint32_t> M_rows;       // block-diagonal M (m rows)
    uint64_t b;                         // combined b
    int64_t sum_T;                      // sum of per-output T
    int64_t union_T;                    // unique monomials
    std::vector<std::vector<uint64_t>> g_tt_raw; // per-output g truth tables (post-merge index)
};

// Run search for a single output
static MbCandidate search_single_output(
    const TruthTable& tt, int oi, int n, const SearchParams& params, int& best_m)
{
    // Create single-output truth table
    TruthTable tt1;
    tt1.n = tt.n;
    tt1.n_outputs = 1;
    tt1.n_words = tt.n_words;
    tt1.output_indices = {tt.output_indices[oi]};
    tt1.tt = {tt.tt[oi]};  // share data (no copy)

    // Need to copy for Möbius in run_search
    TruthTable tt1_copy = tt1;

    // We'll do a manual search instead of using run_search
    // since run_search prints results and expects multi-output

    int n_threads = params.n_threads;
    int64_t n_words = tt1.n_words;

    std::vector<MbCandidate> candidates;
    int64_t raw_T = 0;
    {
        TruthTable tc = tt1;
        moebius_packed(tc.tt[0].data(), n);
        for (int64_t w = 0; w < n_words; w++)
            raw_T += __builtin_popcountll(tc.tt[0][w]);
    }

    int max_m_local = std::min(params.max_m, n);
    if (max_m_local < 1) max_m_local = 1;

    // Phase 1: Identity-row candidates (always available)
    for (int bit = 0; bit < n; bit++) {
        uint32_t M[32] = {0};
        M[0] = (1u << bit);
        uint64_t b_cand = 0;
        auto cand = evaluate_Mb(tt1, M, b_cand, 1, n_threads);
        cand.m = 1;
        candidates.push_back(cand);

        // Also try with b=1
        b_cand = 1;
        cand = evaluate_Mb(tt1, M, b_cand, 1, n_threads);
        cand.m = 1;
        candidates.push_back(cand);
    }

    // Phase 1b: d1c product-based exhaustive b-search (n ≤ 32, support t ≤ 20)
    if (n <= 32) {
        int64_t n_words = tt1.n_words;
        // Compute ANF from truth table copy
        TruthTable tc_d1c = tt1;
        moebius_packed(tc_d1c.tt[0].data(), n);
        const uint64_t* anf = tc_d1c.tt[0].data();

        // Compute support mask (same logic as compute_output_info in search.cpp)
        uint64_t support = 0, lower = 0;
        for (int64_t w = 0; w < n_words; w++) {
            uint64_t val = anf[w];
            if (val == 0) continue;
            support |= ((uint64_t)w << 6);
            if (val == UINT64_MAX) { lower |= 0x3F; }
            else { for (uint64_t v = val; v; v &= v - 1)
                       lower |= (uint64_t)__builtin_ctzll(v); }
        }
        uint64_t support_mask = support | lower;
        if (n < 64) support_mask &= (1ULL << n) - 1;
        int t = __builtin_popcountll(support_mask);

        if (t > 0 && t <= 20) {
            uint64_t best_b = 0;
            int64_t best_T = exhaustive_search_best_b(anf, n, support_mask,
                                                       20, best_b);
            if (best_T >= 0) {
                int64_t raw_T = 0;
                for (int64_t w = 0; w < n_words; w++)
                    raw_T += __builtin_popcountll(anf[w]);
                std::cerr << "    d1c: support=" << support_mask << " t=" << t
                          << " raw_T=" << raw_T << " best_T=" << best_T
                          << " best_b=" << best_b << "\n";
                {
                    // best_b from exhaustive_search_best_b is already a full n-bit b vector
                    uint32_t M_d1c[32] = {0};
                    for (int i = 0; i < n && i < 32; i++)
                        M_d1c[i] = (1u << i);
                    auto cand = evaluate_Mb(tt1, M_d1c, best_b, n, n_threads);
                    cand.m = n;
                    if (cand.total_T >= 0 && cand.total_T <= raw_T)
                        candidates.push_back(cand);
                }
            }
        }
    }

    // Phase 2: Walsh correlations (n ≤ 32, skip if walsh-k=0)
    if (n <= 32 && params.walsh_single_top > 0) {
        auto walsh_vec = compute_walsh_correlations(tt1, n_threads);
        WalshInfo walsh = walsh_vec[0];

        std::vector<int> order(n);
        for (int i = 0; i < n; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return std::abs(walsh.walsh_mag[a]) > std::abs(walsh.walsh_mag[b]);
        });

        int top_k = std::min(n, params.walsh_single_top);

        if (n <= 20) {
            // Small n: direct evaluate_Mb with m < n (fast consistency check)
            for (int m_rows = 2; m_rows <= std::min(8, n); m_rows++) {
                uint32_t M[32] = {0};
                for (int r = 0; r < m_rows && r < top_k; r++)
                    M[r] = (1u << order[r]);
                uint64_t b_cand = 0;
                auto cand = evaluate_Mb(tt1, M, b_cand, m_rows, n_threads);
                cand.m = m_rows;
                if (cand.total_T <= raw_T)
                    candidates.push_back(cand);
            }

            // XOR pairs of top bits (as single combined row)
            for (int i = 0; i < std::min(16, top_k); i++) {
                for (int j = i + 1; j < std::min(16, top_k); j++) {
                    uint32_t M[32] = {0};
                    M[0] = (1u << order[i]) | (1u << order[j]);
                    uint64_t b_cand = 0;
                    auto cand = evaluate_Mb(tt1, M, b_cand, 1, n_threads);
                    cand.m = 1;
                    if (cand.total_T <= raw_T)
                        candidates.push_back(cand);
                }
            }
        } else {
            // n=21..32: bijective evaluation path (pad permutation rows to full rank)
            auto eval_padded = [&](const uint32_t* M_in, int m_in, uint64_t b_in) -> MbCandidate {
                uint32_t M[32] = {0};
                for (int j = 0; j < m_in; j++) M[j] = M_in[j];
                int r = m_in;
                for (int i = 0; i < n && r < n; i++) {
                    M[r] = (1u << i);
                    int ra_before = gf2_rank(M, r, n);
                    int ra_after  = gf2_rank(M, r + 1, n);
                    if (ra_after > ra_before) r++;
                    else M[r] = 0;
                }
                if (r < n) { MbCandidate bad; bad.total_T = INT64_MAX; return bad; }
                return evaluate_Mb_bijective(tt1, M, b_in, n, n, n_threads, false);
            };

            // Single-row Walsh: each top bit as a permutation row, padded to full rank
            for (int i = 0; i < top_k; i++) {
                uint32_t M_row[1] = {1u << order[i]};
                auto c0 = eval_padded(M_row, 1, 0);
                if (c0.total_T < INT64_MAX && c0.total_T <= raw_T)
                    candidates.push_back(c0);
                auto c1 = eval_padded(M_row, 1, 1);
                if (c1.total_T < INT64_MAX && c1.total_T <= raw_T)
                    candidates.push_back(c1);
            }

            // Multi-row progressive: top Walsh bits as separate rows
            for (int m_rows = 2; m_rows <= std::min(8, n); m_rows++) {
                uint32_t M[32] = {0};
                for (int r = 0; r < m_rows && r < top_k; r++)
                    M[r] = (1u << order[r]);
                auto c = eval_padded(M, m_rows, 0);
                if (c.total_T < INT64_MAX && c.total_T <= raw_T)
                    candidates.push_back(c);
            }

            // XOR pairs of top bits (as single combined row, padded)
            for (int i = 0; i < std::min(16, top_k); i++) {
                for (int j = i + 1; j < std::min(16, top_k); j++) {
                    uint32_t M[32] = {0};
                    M[0] = (1u << order[i]) | (1u << order[j]);
                    auto c = eval_padded(M, 1, 0);
                    if (c.total_T < INT64_MAX && c.total_T <= raw_T)
                        candidates.push_back(c);
                }
            }
        }
    }

    // Phase 3: Random candidates
    if (params.n_random > 0) {
        if (n <= 20) {
            for (int i = 0; i < params.n_random; i++) {
                int m_cand = 1 + (rand() % max_m_local);
                uint32_t M[32] = {0};
                for (int r = 0; r < m_cand; r++)
                    M[r] = (uint32_t)(rand() & ((1u << n) - 1));
                uint64_t b_cand = (uint64_t)(rand() & ((1ULL << m_cand) - 1));
                auto cand = evaluate_Mb(tt1, M, b_cand, m_cand, n_threads);
                cand.m = m_cand;
                if (cand.total_T <= raw_T)
                    candidates.push_back(cand);
            }
        } else {
            // n=21..32: full-rank random M only (bijective path is much faster)
            std::mt19937_64 rng_n32(42 + oi);
            int generated = 0, attempts = 0;
            const int MAX_ATTEMPTS = params.n_random * 20;
            while (generated < params.n_random && attempts < MAX_ATTEMPTS) {
                attempts++;
                uint32_t M[32] = {0};
                for (int r = 0; r < n; r++)
                    for (int i = 0; i < n; i++)
                        if (rng_n32() & 1) M[r] |= (1u << i);
                uint64_t b = 0;
                for (int r = 0; r < n; r++)
                    if (rng_n32() & 1) b |= (1ULL << r);
                if (gf2_rank(M, n, n) < n) continue;
                auto cand = evaluate_Mb(tt1, M, b, n, n_threads);
                if (cand.total_T < INT64_MAX && cand.total_T <= raw_T) {
                    candidates.push_back(cand);
                    generated++;
                } else {
                    generated++; // count even if not better, to limit runtime
                }
            }
        }
    }

    // Sort by total_T ascending
    std::sort(candidates.begin(), candidates.end(),
              [](const MbCandidate& a, const MbCandidate& b) {
                  return a.total_T < b.total_T;
              });

    if (candidates.empty()) {
        // Fallback: identity transform
        MbCandidate fallback;
        fallback.m = n;
        fallback.total_T = raw_T;
        fallback.b = 0;
        fallback.M_rows.resize(n);
        for (int r = 0; r < n; r++) {
            fallback.M_rows[r] = (r < 32) ? (1u << r) : 0;
        }
        return fallback;
    }

    MbCandidate best_candidate;
    bool found = false;

    // Try hill climbing from top candidates
    int n_hc = std::min(params.n_hill_climb, (int)candidates.size());
    for (int hi = 0; hi < n_hc; hi++) {
        auto hc = hill_climb(tt1, candidates[hi], n, n_threads);
        if (!found || hc.total_T < best_candidate.total_T) {
            best_candidate = hc;
            found = true;
        }
    }

    // Also try hill climbing from identity (only if we have hill climb iterations)
    if (params.n_hill_climb > 0) {
        TruthTable tc = tt1;
        moebius_packed(tc.tt[0].data(), n);
        int64_t ident_T = 0;
        for (int64_t w = 0; w < n_words; w++)
            ident_T += __builtin_popcountll(tc.tt[0][w]);
        MbCandidate ident_start;
        ident_start.m = n;
        ident_start.total_T = ident_T;
        ident_start.b = 0;
        ident_start.M_rows.resize(n);
        for (int r = 0; r < n; r++)
            ident_start.M_rows[r] = (1u << r);
        auto hc = hill_climb(tt1, ident_start, n, n_threads);
        if (!found || hc.total_T < best_candidate.total_T) {
            best_candidate = hc;
            found = true;
        }
    }

    if (!found) best_candidate = candidates[0];

    // Evaluate the best candidate
    best_candidate = evaluate_Mb(tt1, best_candidate.M_rows.data(), best_candidate.b,
                                  best_candidate.m, n_threads);
    return best_candidate;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [options]\n";
        std::cerr << "  --max-m N      max z variables per output (default 12)\n";
        std::cerr << "  --walsh-k N    top K Walsh bits (default 30)\n";
        std::cerr << "  --random N     random candidates per output (default 40)\n";
        std::cerr << "  --hill-climb N hill climb per output (default 10)\n";
        std::cerr << "  --save-results PREFIX  save combined result\n";
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

    // Resolve circuit path (relative to original CWD), then make relative to project root
    fs::path circ_fs_path(argv[1]);
    if (circ_fs_path.is_relative()) circ_fs_path = orig_cwd / circ_fs_path;
    circ_fs_path = fs::weakly_canonical(circ_fs_path);
    std::string path = circ_fs_path.string();

    Circuit circ = read_circuit(path);
    int n = circ.n_inputs;
    std::vector<int> output_indices;
    for (int o = 0; o < (int)circ.outputs.size(); o++)
        output_indices.push_back(o);
    if (output_indices.empty()) { std::cerr << "No outputs\n"; return 1; }

    SearchParams params;
    params.n_threads = std::thread::hardware_concurrency();
    if (params.n_threads < 1) params.n_threads = 1;
    std::string save_prefix;

    for (int a = 2; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "--max-m" && a + 1 < argc) params.max_m = std::stoi(argv[++a]);
        else if (arg == "--walsh-k" && a + 1 < argc) params.walsh_single_top = std::stoi(argv[++a]);
        else if (arg == "--random" && a + 1 < argc) params.n_random = std::stoi(argv[++a]);
        else if (arg == "--hill-climb" && a + 1 < argc) params.n_hill_climb = std::stoi(argv[++a]);
        else if (arg == "--save-results" && a + 1 < argc) {
            fs::path p(argv[++a]);
            if (p.is_relative()) p = root / p;
            save_prefix = p.string();
        }
    }

    if (n > 32) {
        std::cerr << "  n=" << n << " > 32 not yet supported for opt2\n";
        return 1;
    }

    std::cout << "--- Circuit ---\n";
    std::cout << "  Inputs: " << n << ", Outputs: " << output_indices.size()
              << ", Strategy: opt2 (per-output transforms)\n";

    int k = (int)output_indices.size();
    std::vector<MbCandidate> per_output(k);
    Opt2Result result;
    result.n = n;
    result.m = 0;
    result.sum_T = 0;

    // Phase 1: Compute truth table (all outputs at once)
    std::cout << "\nPhase 1: Computing truth table...\n";
    TruthTable tt = compute_truth_table(circ, output_indices, params.n_threads);
    std::cout << "  Truth table: 2^" << n << " = " << (int64_t(1) << n)
              << " inputs, " << tt.n_words << " batches, "
              << k << " outputs\n";

    // Phase 2: Search per output
    std::cout << "\nPhase 2: Per-output search\n";
    for (int oi = 0; oi < k; oi++) {
        std::cout << "  Output " << oi << "/" << k << " (" << circ.outputs[output_indices[oi]] << ")...\n";
        auto cand = search_single_output(tt, oi, n, params, result.m);
        per_output[oi] = cand;
        result.m += cand.m;
        result.sum_T += cand.total_T;
        std::cout << "    m=" << cand.m << " T=" << cand.total_T << "\n";
    }

    // Phase 3: Build merged transform, deduplicate rows, remap ANF
    std::cout << "\nPhase 3: Merging and deduplicating transforms\n";

    // Step 1: Build sig_to_idx mapping from all per-output (M_row, b_bit)
    std::map<std::pair<uint32_t, int>, int> sig_to_idx;
    std::vector<std::vector<int>> out_var_map(k);
    int original_m = result.m;

    for (int oi = 0; oi < k; oi++) {
        int mi = per_output[oi].m;
        out_var_map[oi].resize(mi);
        for (int r = 0; r < mi; r++) {
            uint32_t M_row = per_output[oi].M_rows[r];
            int b_bit = (per_output[oi].b >> r) & 1;
            auto key = std::make_pair(M_row, b_bit);
            auto it = sig_to_idx.find(key);
            if (it == sig_to_idx.end())
                it = sig_to_idx.insert({key, (int)sig_to_idx.size()}).first;
            out_var_map[oi][r] = it->second;
        }
    }
    int shared_m = (int)sig_to_idx.size();

    // Step 2: Build deduplicated M and b
    result.M_rows.resize(shared_m);
    result.b = 0;
    for (auto& [key, idx] : sig_to_idx) {
        result.M_rows[idx] = key.first;
        if (key.second)
            result.b |= (1ULL << idx);
    }
    result.m = shared_m;

    // Step 3: For each output, evaluate and remap ANF into shared z-space
    // Store ALL remapped terms for .poly, but only count deg≥2 for sum_T/union_T
    std::vector<std::vector<std::vector<int>>> out_all_terms(k);
    std::set<std::vector<int>> shared_terms;    // unique deg≥2 across all outputs
    std::vector<int> out_deg2_count(k, 0);      // per-output deg≥2 counts
    int max_deg = 0;

    for (int oi = 0; oi < k; oi++) {
        const auto& cand = per_output[oi];
        if (cand.m == 0) continue;

        TruthTable tt_oi;
        tt_oi.n = n;
        tt_oi.n_outputs = 1;
        tt_oi.n_words = tt.n_words;
        tt_oi.tt = {tt.tt[oi]};

        auto ev = evaluate_Mb(tt_oi, cand.M_rows.data(), cand.b,
                              cand.m, params.n_threads, true);
        if (ev.g_tt_raw.empty()) continue;

        std::vector<uint64_t> anf(ev.g_tt_raw[0]);
        moebius_packed(anf.data(), cand.m);

        int64_t nw = (cand.m < 6) ? 1 : (int64_t(1) << (cand.m - 6));
        for (int64_t w = 0; w < nw; w++) {
            uint64_t word = anf[w];
            while (word) {
                int bit_idx = __builtin_ctzll(word);
                word &= word - 1;
                int mono_mask = (int)((w << 6) | bit_idx);

                std::vector<int> svars;
                int rem = mono_mask;
                while (rem) {
                    int j = __builtin_ctzll(rem);
                    rem &= rem - 1;
                    if (j < (int)out_var_map[oi].size())
                        svars.push_back(out_var_map[oi][j]);
                }
                std::sort(svars.begin(), svars.end());
                int deg = (int)svars.size();
                if (deg > max_deg) max_deg = deg;

                // Store ALL terms for .poly
                out_all_terms[oi].push_back(svars);

                // Only deg≥2 counts for sum_T / union_T
                if (deg >= 2) {
                    out_deg2_count[oi]++;
                    shared_terms.insert(std::move(svars));
                }
            }
        }
    }

    // sum_T = sum of per-output unique deg≥2 terms
    // union_T = unique deg≥2 terms across all outputs
    result.sum_T = 0;
    for (int oi = 0; oi < k; oi++)
        result.sum_T += out_deg2_count[oi];
    result.union_T = (int64_t)shared_terms.size();

    std::cout << "  Dedup: " << original_m << " → " << shared_m << " z variables\n";
    std::cout << "  Sum T = " << result.sum_T << "\n";
    std::cout << "  Union T = " << result.union_T << "\n";

    // ---- Save results ----
    if (!save_prefix.empty()) {
        std::string inst_name = std::filesystem::path(path).stem().string();
        std::string dir = save_prefix;
        while (!dir.empty() && dir.back() == '/') dir.pop_back();
#ifndef STRATEGY_TAG
#define STRATEGY_TAG "opt2"
#endif
#ifndef DIRECTION_TAG
#define DIRECTION_TAG "d3"
#endif
        std::string tag = std::string("_") + DIRECTION_TAG + "_" + std::string(STRATEGY_TAG);
        std::string prefix = dir + "/" + inst_name + tag;

        // ---- .affine: s × (n+1) matrix [M | b] (deduplicated) ----
        {
            std::ofstream f(prefix + ".affine");
            if (f) {
                f << shared_m << "\n" << shared_m << " " << (n + 1) << "\n";
                for (int r = 0; r < shared_m; r++) {
                    for (int c = 0; c < n; c++) {
                        if (c > 0) f << " ";
                        f << ((result.M_rows[r] >> c) & 1);
                    }
                    f << " " << ((result.b >> r) & 1) << "\n";
                }
                std::cout << "  Saved: " << prefix << ".affine (s=" << shared_m << ")\n";
            }
        }

        // ---- .poly: per-output ANF in shared z-space (ALL terms including deg 0/1) ----
        {
            std::vector<int> term_counts(k, 0);
            for (int oi = 0; oi < k; oi++)
                term_counts[oi] = (int)out_all_terms[oi].size();

            std::ofstream f(prefix + ".poly");
            if (f) {
                f << shared_m << "\n" << k << "\n";
                for (int oi = 0; oi < k; oi++) {
                    if (oi > 0) f << " ";
                    f << term_counts[oi];
                }
                f << "\n";
                for (int oi = 0; oi < k; oi++) {
                    for (const auto& svars : out_all_terms[oi]) {
                        f << "[";
                        int svi = 0;
                        for (int b = 0; b < shared_m; b++) {
                            if (b > 0) f << " , ";
                            int bit = (svi < (int)svars.size() && svars[svi] == b) ? 1 : 0;
                            if (bit) svi++;
                            f << bit;
                        }
                        f << " , 1]\n";
                    }
                }
                std::cout << "  Saved: " << prefix << ".poly (m=" << shared_m
                          << ", " << k << " outputs)\n";
            }
        }

        // ---- _stats.txt: 5-line numeric ----
        {
            std::ofstream f(prefix + "_stats.txt");
            if (f) {
                f << n << "\n" << k << "\n" << result.sum_T << "\n"
                  << result.union_T << "\n" << max_deg << "\n";
                std::cout << "  Saved: " << prefix << "_stats.txt (sum_T="
                          << result.sum_T << ", union_T=" << result.union_T << ")\n";
            }
        }
    }

    std::cout << "\nTotal time: skipped\n";
    return 0;
}
