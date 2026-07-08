/**
 * optimize_anf_gatebuilder — Direction 1b (gate-level linear mapping construction).
 *
 * Builds (g, M, b) incrementally from circuit gates, then applies
 * c-correction and outputs the result.
 *
 * Usage:
 *   optimize_anf_gatebuilder <circuit.txt> [--threshold N] [--out-dir DIR]
 *
 * Output files (in --out-dir):
 *   {inst}_gatebuilder_trans.poly — M, b
 *   {inst}_gatebuilder_anf.poly  — corrected g (degree ≥ 2 only)
 *   {inst}_gatebuilder_T.txt     — per-output T values
 */
#include "circuit.h"
#include "gate_builder.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <filesystem>
#include <climits>

namespace fs = std::filesystem;

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}

static std::string strip_ext(const std::string& path) {
    auto pos = path.rfind('.');
    return (pos == std::string::npos) ? path : path.substr(0, pos);
}

// ====================================================================
//  Per-output compression: complement selection only.
//  Greedy merge is NOT used (see compress_output for explanation).
// ====================================================================

// Run complement selection on one output's ANF, updating b in place.
// This is the only safe per-output compression: it XORs z_j XOR 1 into
// selected variables to reduce T, without changing the function.
static void compress_output(SparseANF& g, std::vector<uint32_t>& M_rows,
                             uint64_t& b, int& m_used) {
    if (g.T() == 0) return;
    if (M_rows.size() <= 1) return;

    // Complement selection: for each z_j, try XOR with 1 (z_j -> not z_j).
    // This reduces T without changing m or introducing verification error.
    std::vector<uint32_t> comp_M;  // unit-vector rows, irrelevant for M
    uint64_t comp_b = 0;
    g.complement_search_greedy(comp_M, comp_b);
    b ^= comp_b;

    // IMPORTANT: Greedy merge (z_i = z_i XOR z_j, then remove z_i) is DISABLED
    // because it only produces correct results when M_i = M_j and b_i = b_j.
    // In the gate builder's per-output z-space, different z-variables nearly
    // always have different (M,b), so the merge silently changes f(x).
    // Use truth-table directions d1a, d2, or d3 for safe variable reduction.

    m_used = (int)M_rows.size();  // must match affine row count
}

// ====================================================================
//  Per-output processing
// ====================================================================
struct OutputResult {
    std::string name;
    SparseANF g_all;          // full ANF in z-space (all terms, no correction)
    std::vector<uint32_t> M_rows;  // per-output M (raw, no c/d rows)
    uint64_t b;
    int m_used;                // active z variables used
};

static std::vector<OutputResult> process_outputs(
    const Circuit& circ,
    const GateBuilder& gb,
    bool verbose = false)
{
    std::vector<OutputResult> results;
    for (const auto& out_name : circ.outputs) {
        auto* si = gb.get_signal(out_name);
        if (!si) {
            printf("  WARNING: output %s not found\n", out_name.c_str());
            continue;
        }

        if (verbose)
            printf("  Output %s: m=%d T=%d T_nl=%d\n",
                   out_name.c_str(), si->m, si->g.T(), si->g.T_nonlinear());

        // Determine how many z variables are actually used
        int m_used = 0;
        for (auto& [mask, v] : si->g.terms()) {
            if (v) {
                int max_bit = 63 - __builtin_clzll(mask);
                if (max_bit + 1 > m_used) m_used = max_bit + 1;
            }
        }

        results.push_back({out_name, si->g, si->M_rows, si->b, m_used});

        if (verbose) {
            printf("    m_used=%d T=%d\n", m_used, si->g.T());
        }
    }
    return results;
}

// ====================================================================
//  File writers — matrix format (.affine, .poly, _stats.txt)
// ====================================================================

static void write_affine(const std::string& path,
                          const std::vector<OutputResult>& results,
                          const Circuit& circ)
{
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }

    int n = circ.n_inputs;

    // Total z vars = sum of per-output m_used
    int s = 0;
    for (const auto& r : results) s += r.m_used;

    f << s << "\n" << s << " " << (n + 1) << "\n";

    // s x (n+1) matrix: [M | b] -- write only m_used rows per output
    for (const auto& r : results) {
        for (int ri = 0; ri < r.m_used; ri++) {
            for (int c = 0; c < n; c++) {
                if (c > 0) f << " ";
                f << ((r.M_rows[ri] >> c) & 1);
            }
            f << " " << ((r.b >> ri) & 1) << "\n";
        }
    }

    std::cout << "  Saved: " << path << " (s=" << s << ", n=" << n << ")\n";
}

static void write_poly(const std::string& path,
                        const std::vector<OutputResult>& results)
{
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }

    int k = (int)results.size();
    if (k == 0) return;

    // Compute z offsets and total s
    std::vector<int> z_offsets(k, 0);
    int s = 0;
    for (int oi = 0; oi < k; oi++) {
        z_offsets[oi] = s;
        s += results[oi].m_used;
    }

    // Count terms per output
    std::vector<int> term_counts(k, 0);
    for (int oi = 0; oi < k; oi++) {
        if (results[oi].g_all.T() > 0) {
            for (auto& [mask, v] : results[oi].g_all.terms()) {
                if (v) term_counts[oi]++;
            }
        }
    }

    f << s << "\n" << k << "\n";
    for (int oi = 0; oi < k; oi++) {
        if (oi > 0) f << " ";
        f << term_counts[oi];
    }
    f << "\n";

    for (int oi = 0; oi < k; oi++) {
        int off = z_offsets[oi];

        // g_all terms (shifted by off)
        if (results[oi].g_all.T() > 0) {
            for (auto& [mask, v] : results[oi].g_all.terms()) {
                if (!v) continue;
                f << "[";
                for (int b = 0; b < s; b++) {
                    if (b > 0) f << " , ";
                    int bit = 0;
                    int mi = results[oi].m_used;
                    if (b >= off && b < off + mi) bit = (mask >> (b - off)) & 1;
                    f << bit;
                }
                f << " , 1]\n";
            }
        }
    }
    std::cout << "  Saved: " << path << " (m=" << s << ", " << k << " outputs)\n";
}

static void write_stats(const std::string& path,
                         const std::vector<OutputResult>& results,
                         int n)
{
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }

    int k = (int)results.size();
    int64_t sum_T = 0;
    int max_deg = 0;
    for (int oi = 0; oi < k; oi++) {
        sum_T += results[oi].g_all.T_nonlinear();
        if (results[oi].g_all.T() > 0) {
            for (auto& [mask, v] : results[oi].g_all.terms()) {
                int deg = __builtin_popcountll(mask);
                if (deg > max_deg) max_deg = deg;
            }
        }
    }

    // Each output has disjoint z → no cross-output overlap, union_T == sum_T
    int64_t union_T = sum_T;
    f << n << "\n" << k << "\n" << sum_T << "\n" << union_T << "\n" << max_deg << "\n";
    std::cout << "  Saved: " << path << " (sum_T=" << sum_T << ", union_T=" << union_T << ")\n";
}

// ====================================================================
//  Main
// ====================================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [--threshold N] [--out-dir DIR]\n";
        return 1;
    }

    std::string circ_path;
    int threshold = INT32_MAX;  // no inline compression; done per-output below
    std::string out_dir = "";

    // ---- Resolve project root from executable path, then chdir ----
    {
        auto orig_cwd = fs::current_path();
        auto exe = fs::weakly_canonical(fs::absolute(fs::path(argv[0])));
        auto root = exe.parent_path();
        for (int i = 0; i < 3 && !fs::exists(root / "examples"); i++)
            root = root.parent_path();
        fs::current_path(root);

        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "--threshold" && i + 1 < argc)
                threshold = std::stoi(argv[++i]);
            else if (arg == "--out-dir" && i + 1 < argc) {
                fs::path p(argv[++i]);
                if (p.is_relative()) p = root / p;
                out_dir = p.string();
            }
            else if (arg[0] != '-')
                circ_path = arg;
        }

        if (circ_path.empty()) {
            std::cerr << "  ERROR: no circuit file specified\n";
            return 1;
        }

        // Resolve circuit path
        fs::path cp(circ_path);
        if (cp.is_relative()) cp = orig_cwd / cp;
        circ_path = fs::weakly_canonical(cp).string();
    }

    // Derive instance name from circuit path
    std::string inst = strip_ext(fs::path(circ_path).filename().string());

    // Default: results/{inst}/
    if (out_dir.empty())
        out_dir = "results/" + inst;

    printf("GateBuilder: %s\n", inst.c_str());
    printf("  threshold=%d out_dir=%s\n", threshold, out_dir.c_str());

    // Read circuit
    Circuit circ = read_circuit(circ_path);
    printf("  Circuit: %d inputs, %d outputs, %d gates\n",
           circ.n_inputs, (int)circ.outputs.size(), (int)circ.stmts.size());

    // Build (g, M, b)
    GateBuilder gb;
    gb.set_threshold(threshold);
    gb.build(circ);
    printf("  Build complete\n");

    // Extract outputs
    auto results = process_outputs(circ, gb, true);

    // Per-output compression (each output in its own z-space)
    for (auto& r : results) {
        if (r.g_all.T() > 0) {
            compress_output(r.g_all, r.M_rows, r.b, r.m_used);
        }
    }
    printf("  Per-output compression complete\n");
    for (const auto& r : results)
        printf("    %s: m=%d T=%d\n", r.name.c_str(), r.m_used, r.g_all.T());

    // ---- Deduplicate M rows across outputs, remap ANF, recalculate stats ----
    // (M_row, b_bit) → shared z index
    std::map<std::pair<uint32_t, int>, int> sig_to_idx;
    for (auto& r : results) {
        for (int ri = 0; ri < r.m_used; ri++) {
            auto key = std::make_pair(r.M_rows[ri], (int)((r.b >> ri) & 1));
            if (sig_to_idx.find(key) == sig_to_idx.end())
                sig_to_idx[key] = (int)sig_to_idx.size();
        }
    }
    int shared_m = (int)sig_to_idx.size();

    // Build deduplicated M, b
    std::vector<uint32_t> dedup_M(shared_m);
    uint64_t dedup_b = 0;
    for (auto& [key, idx] : sig_to_idx) {
        dedup_M[idx] = key.first;
        if (key.second) dedup_b |= (1ULL << idx);
    }

    // Per-output local z → shared z mapping
    std::vector<std::vector<int>> var_map(results.size());
    for (int oi = 0; oi < (int)results.size(); oi++) {
        auto& r = results[oi];
        var_map[oi].resize(r.m_used);
        for (int ri = 0; ri < r.m_used; ri++) {
            auto key = std::make_pair(r.M_rows[ri], (int)((r.b >> ri) & 1));
            var_map[oi][ri] = sig_to_idx[key];
        }
    }

    // Remap ANF terms and compute stats
    struct OutTerm { std::vector<int> svars; };
    std::vector<std::vector<OutTerm>> out_terms(results.size());
    std::set<std::vector<int>> shared_deg2_terms;
    int64_t sum_T = 0;
    int max_deg = 0;

    for (int oi = 0; oi < (int)results.size(); oi++) {
        auto& r = results[oi];
        for (auto& [mask, v] : r.g_all.terms()) {
            if (!v) continue;
            std::vector<int> svars;
            uint64_t m = mask;
            while (m) {
                int j = __builtin_ctzll(m);
                m &= m - 1;
                if (j < (int)var_map[oi].size())
                    svars.push_back(var_map[oi][j]);
            }
            std::sort(svars.begin(), svars.end());
            int deg = (int)svars.size();
            if (deg > max_deg) max_deg = deg;
            out_terms[oi].push_back({std::move(svars)});
            if (deg >= 2) {
                sum_T++;
                shared_deg2_terms.insert(out_terms[oi].back().svars);
            }
        }
    }
    int64_t union_T = (int64_t)shared_deg2_terms.size();

    printf("  Dedup: %d outputs, %d unique z (from %d raw rows)\n",
           (int)results.size(), shared_m, (int)sig_to_idx.size());
    printf("  sum_T=%ld union_T=%ld max_deg=%d\n", (long)sum_T, (long)union_T, max_deg);

    // ---- Write output files ----
    fs::create_directories(out_dir);
#ifndef STRATEGY_TAG
#define STRATEGY_TAG "opt2"
#endif
    std::string tag = "_d1b_" + std::string(STRATEGY_TAG);
    std::string base = out_dir + "/" + inst + tag;

    // .affine: s × (n+1) matrix [M | b] (deduplicated)
    {
        std::ofstream f(base + ".affine");
        if (f) {
            f << shared_m << "\n" << shared_m << " " << (circ.n_inputs + 1) << "\n";
            for (int r = 0; r < shared_m; r++) {
                for (int c = 0; c < circ.n_inputs; c++) {
                    if (c > 0) f << " ";
                    f << ((dedup_M[r] >> c) & 1);
                }
                f << " " << ((dedup_b >> r) & 1) << "\n";
            }
            std::cout << "  Saved: " << base << ".affine (s=" << shared_m << ", n=" << circ.n_inputs << ")\n";
        }
    }

    // .poly: per-output ANF in shared z-space
    {
        int k = (int)results.size();
        std::vector<int> term_counts(k, 0);
        for (int oi = 0; oi < k; oi++)
            term_counts[oi] = (int)out_terms[oi].size();

        std::ofstream f(base + ".poly");
        if (f) {
            f << shared_m << "\n" << k << "\n";
            for (int oi = 0; oi < k; oi++) {
                if (oi > 0) f << " ";
                f << term_counts[oi];
            }
            f << "\n";
            for (int oi = 0; oi < k; oi++) {
                for (const auto& term : out_terms[oi]) {
                    f << "[";
                    int svi = 0;
                    for (int b = 0; b < shared_m; b++) {
                        if (b > 0) f << " , ";
                        int bit = (svi < (int)term.svars.size() && term.svars[svi] == b) ? 1 : 0;
                        if (bit) svi++;
                        f << bit;
                    }
                    f << " , 1]\n";
                }
            }
            std::cout << "  Saved: " << base << ".poly (m=" << shared_m
                      << ", " << k << " outputs)\n";
        }
    }

    // _stats.txt: 5 lines
    {
        std::ofstream f(base + "_stats.txt");
        if (f) {
            f << circ.n_inputs << "\n" << results.size() << "\n"
              << sum_T << "\n" << union_T << "\n" << max_deg << "\n";
            std::cout << "  Saved: " << base << "_stats.txt (sum_T="
                      << sum_T << ", union_T=" << union_T << ")\n";
        }
    }

    bool has_solution = (shared_m > 0);
    std::cout << "# STATUS: " << (has_solution ? "has_solution" : "no_solution") << " no_timeout\n";

    std::cout << "  Done. Files written to " << base << ".*\n";
    return 0;
}
