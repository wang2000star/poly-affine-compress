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
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <algorithm>

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

    // Phase 2: Walsh multi-row (for n ≤ 20)
    if (n <= 20) {
        auto walsh_vec = compute_walsh_correlations(tt1, n_threads);
        WalshInfo walsh = walsh_vec[0];

        std::vector<int> order(n);
        for (int i = 0; i < n; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return std::abs(walsh.walsh_mag[a]) > std::abs(walsh.walsh_mag[b]);
        });

        int top_k = std::min(n, params.walsh_single_top);
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

        // XOR pairs of top bits
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
    }

    // Phase 3: Random candidates
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

    // Also try hill climbing from identity
    {
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

    // Evaluate the best candidate to get g_tt_raw
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
    std::string path = argv[1];

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
        else if (arg == "--save-results" && a + 1 < argc) save_prefix = argv[++a];
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

    // Phase 3: Build merged transform and compute union T
    std::cout << "\nPhase 3: Merging transforms (total m=" << result.m << ")\n";
    result.M_rows.resize(result.m, 0);
    result.b = 0;
    int offset = 0;
    for (int oi = 0; oi < k; oi++) {
        int mi = per_output[oi].m;
        for (int r = 0; r < mi; r++) {
            result.M_rows[offset + r] = per_output[oi].M_rows[r];
        }
        uint64_t bi = per_output[oi].b;
        for (int r = 0; r < mi; r++) {
            if ((bi >> r) & 1)
                result.b |= (1ULL << (offset + r));
        }
        offset += mi;
    }

    result.sum_T = 0;
    for (int oi = 0; oi < k; oi++)
        result.sum_T += per_output[oi].total_T;

    if (result.m <= 20) {
        // Small m: full merged evaluation
        TruthTable tt_combined = tt;
        auto merged = evaluate_Mb(tt_combined, result.M_rows.data(), result.b,
                                   result.m, params.n_threads, true);
        result.g_tt_raw = merged.g_tt_raw;

        int64_t n_words_g = (result.m < 6) ? 1 : (int64_t(1) << (result.m - 6));
        std::vector<uint64_t> union_data(n_words_g, 0);
        for (int oi = 0; oi < k; oi++) {
            std::vector<uint64_t> anf(result.g_tt_raw[oi]);
            moebius_packed(anf.data(), result.m);
            for (int64_t w = 0; w < n_words_g; w++)
                union_data[w] |= anf[w];
        }
        result.union_T = 0;
        for (int64_t w = 0; w < n_words_g; w++)
            result.union_T += __builtin_popcountll(union_data[w]);
    } else {
        // Large m: per-output z-variables are disjoint → union_T = sum_T
        result.union_T = result.sum_T;
    }

    std::cout << "  Sum T = " << result.sum_T << "\n";
    std::cout << "  Union T = " << result.union_T << "\n";

    // Save results
    if (!save_prefix.empty()) {
        std::string p = save_prefix;

        // Trans file
        {
            std::ofstream f(p + "_opt2_trans.poly");
            f << "# opt2 transform z_i = M_i x + b_i (per-output, block-diagonal M)\n";
            f << "# m=" << result.m << ", n=" << n << "\n";
            int off = 0;
            for (int oi = 0; oi < k; oi++) {
                int mi = per_output[oi].m;
                for (int r = 0; r < mi; r++) {
                    f << "z_" << (off + r) << " = ";
                    bool first = true;
                    for (int i = 0; i < n; i++) {
                        if ((per_output[oi].M_rows[r] >> i) & 1) {
                            if (!first) f << " + ";
                            f << circ.inputs[i];
                            first = false;
                        }
                    }
                    if ((per_output[oi].b >> r) & 1) {
                        if (!first) f << " + ";
                        f << "1";
                    } else if (first) {
                        f << "0";
                    }
                    f << "\n";
                }
                off += mi;
            }
        }

        // Expr file (only when merged g_tt_raw is available)
        if (!result.g_tt_raw.empty()) {
            save_opt_expr(result.g_tt_raw, circ, output_indices, result.m, p + "_opt2_expr.poly");
        }

        // T file
        {
            std::ofstream f(p + "_opt2_T.poly");
            f << "# opt2 T(g) (per-output, m=" << result.m << ")\n";
            f << "# sum T = " << result.sum_T << "\n";
            f << "# union T = " << result.union_T << "\n";
            for (int oi = 0; oi < k; oi++)
                f << circ.outputs[output_indices[oi]] << ": T=" << per_output[oi].total_T << "\n";
        }

        std::cout << "\n  Saved: " << p << "_opt2_trans.poly\n";
        if (!result.g_tt_raw.empty())
            std::cout << "  Saved ANF: " << p << "_opt2_expr.poly (" << result.sum_T << " terms)\n";
        std::cout << "  Saved T: " << p << "_opt2_T.poly (sum=" << result.sum_T
                  << ", union=" << result.union_T << ")\n";
    }

    std::cout << "\nTotal time: skipped\n";
    return 0;
}
