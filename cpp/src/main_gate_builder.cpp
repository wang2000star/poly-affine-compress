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
//  C-correction: remove linear terms from g, compute c and d
// ====================================================================
struct CCorrection {
    std::vector<int> c;      // c vector (length n)
    int d;                    // d scalar
    SparseANF g_corrected;    // g with linear terms removed
};

static CCorrection apply_correction(
    const SparseANF& g,
    const std::vector<uint32_t>& M_rows,
    uint64_t b,
    int n_inputs)
{
    SparseANF g_copy = g;
    auto L = g_copy.extract_linear_terms();

    // c = M^T · L: XOR of M_rows[j] for each j where L[j]=1
    std::vector<int> c_vec(n_inputs, 0);
    for (int j = 0; j < (int)L.size(); j++) {
        if (L[j] && j < (int)M_rows.size()) {
            uint32_t row = M_rows[j];
            for (int i = 0; i < n_inputs; i++) {
                if ((row >> i) & 1) c_vec[i] ^= 1;
            }
        }
    }

    // d = <L, b> = XOR of L[j]*b_bit_j
    int d_val = 0;
    for (int j = 0; j < (int)L.size(); j++) {
        if (L[j] && ((b >> j) & 1)) d_val ^= 1;
    }

    return {c_vec, d_val, g_copy};
}

// ====================================================================
//  Per-output processing
// ====================================================================
struct OutputResult {
    std::string name;
    SparseANF g_corrected;    // nonlinear ANF in z-space
    std::vector<uint32_t> M_rows;  // per-output M
    uint64_t b;
    std::vector<int> c;
    int d;
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

        auto cc = apply_correction(si->g, si->M_rows, si->b, circ.n_inputs);

        // Determine how many z variables are actually used
        int m_used = 0;
        for (auto& [mask, v] : cc.g_corrected.terms()) {
            if (v) {
                int max_bit = 63 - __builtin_clzll(mask);
                if (max_bit + 1 > m_used) m_used = max_bit + 1;
            }
        }

        results.push_back({out_name, cc.g_corrected, si->M_rows, si->b,
                           cc.c, cc.d, m_used});

        if (verbose) {
            printf("    after c-corr: T_nl=%d m_used=%d\n",
                   cc.g_corrected.T_nonlinear(), m_used);
        }
    }
    return results;
}

// ====================================================================
//  File writers
// ====================================================================

static void write_trans(const std::string& path, const std::string& inst,
                         const std::vector<std::string>& inputs,
                         const std::vector<OutputResult>& results)
{
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    f << "# Transform z = Mx + b for " << inst << "\n";
    f << "# direction=1, strategy=opt2, variant=1b (gate-level)\n";

    int z_idx = 0;
    for (auto& r : results) {
        for (int j = 0; j < (int)r.M_rows.size(); j++) {
            uint32_t row = r.M_rows[j];
            if (row == 0 && !((r.b >> j) & 1)) continue;
            f << "z_" << z_idx << " = ";
            bool first = true;
            for (int i = 0; i < (int)inputs.size(); i++) {
                if ((row >> i) & 1) {
                    if (!first) f << " + ";
                    f << inputs[i];
                    first = false;
                }
            }
            if ((r.b >> j) & 1) {
                if (!first) f << " + ";
                f << "1";
                first = false;
            }
            f << "  # " << r.name << "\n";
            z_idx++;
        }
        // c rows: ⟨c, x⟩
        for (int i = 0; i < (int)r.c.size(); i++) {
            if (r.c[i]) {
                f << "z_" << z_idx << " = " << inputs[i] << "  # c for " << r.name << "\n";
                z_idx++;
            }
        }
        // d row
        if (r.d) {
            f << "z_" << z_idx << " = 1  # d for " << r.name << "\n";
            z_idx++;
        }
    }
}

static void write_anf(const std::string& path, const std::string& inst,
                       const std::vector<OutputResult>& results,
                       const std::vector<std::string>& z_names)
{
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    f << "# ANF for " << inst << " (g corrected, degree >= 2 only)\n";
    f << "# direction=1, strategy=opt2, variant=1b\n";

    for (auto& r : results) {
        if (r.g_corrected.T_nonlinear() == 0) continue;
        f << "# " << r.name << "\n";
        auto& terms = r.g_corrected.terms();
        std::vector<std::pair<uint64_t, int>> sorted(terms.begin(), terms.end());
        std::sort(sorted.begin(), sorted.end());
        for (auto& [mask, v] : sorted) {
            if (!v || __builtin_popcountll(mask) < 2) continue;
            f << "  ";
            bool first = true;
            for (int j = 0; j < 64; j++) {
                if ((mask >> j) & 1) {
                    if (!first) f << " * ";
                    f << (j < (int)z_names.size() ? z_names[j] : ("z_" + std::to_string(j)));
                    first = false;
                }
            }
            f << "\n";
        }
    }
}

static void write_T(const std::string& path, const std::string& inst,
                     const std::vector<OutputResult>& results)
{
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    int64_t sum_T = 0;
    for (auto& r : results) sum_T += r.g_corrected.T_nonlinear();
    f << "# T(g) results for " << inst << "\n";
    f << "# direction=1, strategy=opt2, variant=1b\n";
    f << "sum_T = " << sum_T << "\n";
    f << "---\n";
    for (auto& r : results) {
        int T_nl = r.g_corrected.T_nonlinear();
        f << r.name << ":  T=" << T_nl << "  m=" << r.m_used
          << "  nonlinear_terms=" << T_nl << "\n";
    }
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

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--threshold" && i + 1 < argc)
            threshold = std::stoi(argv[++i]);
        else if (arg == "--out-dir" && i + 1 < argc)
            out_dir = argv[++i];
        else if (arg[0] != '-')
            circ_path = arg;
    }

    if (circ_path.empty()) {
        std::cerr << "  ERROR: no circuit file specified\n";
        return 1;
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
    write_trans(out_dir + "/" + inst + "_gatebuilder_trans.poly", inst,
                circ.inputs, results);
    write_anf(out_dir + "/" + inst + "_gatebuilder_anf.poly", inst,
              results, z_names);
    write_T(out_dir + "/" + inst + "_gatebuilder_T.txt", inst, results);

    printf("  Done. Files written to %s/%s_gatebuilder_*\n",
           out_dir.c_str(), inst.c_str());
    return 0;
}
