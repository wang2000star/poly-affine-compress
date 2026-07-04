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
#include <filesystem>

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

    // s × (n+1) matrix: [M | b]
    for (const auto& r : results) {
        for (int ri = 0; ri < (int)r.M_rows.size(); ri++) {
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
                uint64_t shared_pos = (uint64_t)mask << off;
                f << "[";
                for (int b = 0; b < s; b++) {
                    if (b > 0) f << " , ";
                    f << ((shared_pos >> b) & 1);
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
        int t = results[oi].g_all.T();
        sum_T += t;
        if (t > 0) {
            for (auto& [mask, v] : results[oi].g_all.terms()) {
                if (v) {
                    int deg = __builtin_popcountll(mask);
                    if (deg > max_deg) max_deg = deg;
                }
            }
        }
    }

    // Each output has disjoint z, so sum_T == union_T
    f << n << "\n" << k << "\n" << sum_T << "\n" << sum_T << "\n" << max_deg << "\n";
    std::cout << "  Saved: " << path << " (sum_T=" << sum_T << ")\n";
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
    int threshold = 4096;
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

    // Prepare z variable names
    std::vector<std::string> z_names;
    for (int j = 0; j < 64; j++)
        z_names.push_back("z_" + std::to_string(j));

    // Write output files
    fs::create_directories(out_dir);
    int k_tag = (int)results.size();
#ifndef STRATEGY_TAG
#define STRATEGY_TAG "opt2"
#endif
    std::string tag = "_d1b_" + std::string(STRATEGY_TAG);
    std::string base = out_dir + "/" + inst + tag;
    write_affine(base + ".affine", results, circ);
    write_poly(base + ".poly", results);
    write_stats(base + "_stats.txt", results, circ.n_inputs);

    std::cout << "  Done. Files written to " << (out_dir + "/" + inst + tag + ".*") << "\n";
    return 0;
}
