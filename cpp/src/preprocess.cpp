/**
 * preprocess.cpp — One-time circuit preprocessing
 *
 * Reads raw circuit, remaps variables, classifies outputs, saves fixed files.
 *
 * Usage:
 *   ./preprocess <circuit.txt> [--out-dir DIR]
 *
 * Remapping:
 *   inputs  (INORDER)  →  x_0, x_1, ..., x_{n-1}
 *   outputs (OUTORDER) →  y_0, y_1, ..., y_{k-1}
 *   internal signals    →  t_0, t_1, ... (stmt order)
 *
 * Output:
 *   {out_dir}/{inst}.txt        — remapped circuit
 *   {out_dir}/{inst}_stats.txt  — statistics and output classification
 */
#include "circuit.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <random>

namespace fs = std::filesystem;

// ====================================================================
//  Helpers
// ====================================================================

static std::string strip_ext(const std::string& path) {
    auto pos = path.rfind('.');
    return (pos == std::string::npos) ? path : path.substr(0, pos);
}

// Möbius transform in-place (n ≤ 16)
static void mobius_inplace(std::vector<uint64_t>& tt, int n) {
    int size = 1 << n;
    for (int i = 0; i < n; i++) {
        int step = 1 << i;
        for (int j = 0; j < size; j += 2 * step)
            for (int k = j; k < j + step; k++)
                tt[k + step] ^= tt[k];
    }
}

// Count ANF terms with degree ≥ 2
static int count_nonlinear(const std::vector<uint64_t>& anf, int n) {
    int cnt = 0;
    for (int mask = 0; mask < (1 << n); mask++)
        if (anf[mask] && __builtin_popcount(mask) >= 2) cnt++;
    return cnt;
}

// Maximum degree in ANF
static int max_degree(const std::vector<uint64_t>& anf, int n) {
    int d = 0;
    for (int mask = 0; mask < (1 << n); mask++)
        if (anf[mask]) d = std::max(d, __builtin_popcount(mask));
    return d;
}

// Affine detection via BLR test (any n ≤ 32)
// Returns coefficients (a, b) such that f(x) = ⟨a, x⟩ ⊕ b
// Returns false if nonlinear
static bool detect_affine(const Circuit& circ, int output_idx,
                           uint32_t& a, int& b) {
    int n = circ.n_inputs;
    auto& eval = eval_circuit_point;

    // b = f(0)
    auto r0 = eval(circ, 0);
    b = r0[output_idx];

    // a_i = f(e_i) ⊕ f(0)
    a = 0;
    for (int i = 0; i < n; i++) {
        auto ri = eval(circ, 1U << i);
        if (ri[output_idx] ^ r0[output_idx])
            a |= (1U << i);
    }

    // BLR test: 100 random triples
    std::mt19937 rng(12345);
    uint32_t mask = (n == 32) ? ~0U : ((1U << n) - 1);
    for (int t = 0; t < 100; t++) {
        uint32_t x = (uint32_t)rng() & mask;
        uint32_t u = (uint32_t)rng() & mask;
        uint32_t v = (uint32_t)rng() & mask;

        auto fx    = eval(circ, x);
        auto fx_u  = eval(circ, x ^ u);
        auto fx_v  = eval(circ, x ^ v);
        auto fx_uv = eval(circ, x ^ u ^ v);

        if (fx[output_idx] ^ fx_u[output_idx] ^
            fx_v[output_idx] ^ fx_uv[output_idx])
            return false;
    }

    // Verify coefficients: 200 random points
    for (int t = 0; t < 200; t++) {
        uint32_t x = (uint32_t)rng() & mask;
        auto fx = eval(circ, x)[output_idx];
        int dot = __builtin_popcount(x & a) & 1;
        if (fx != (dot ^ b)) return false;
    }
    return true;
}

// ====================================================================
//  Variable remapping
// ====================================================================

struct NameMapper {
    std::unordered_map<std::string, std::string> orig_to_new;
    std::vector<std::string> x_names;  // x_0 ...
    std::vector<std::string> y_names;  // y_0 ...
    std::vector<std::string> t_names;  // t_0 ...
    int t_count{0};

    void build(const Circuit& circ) {
        // Map inputs → x_i
        for (int i = 0; i < (int)circ.inputs.size(); i++) {
            std::string xn = "x_" + std::to_string(i);
            x_names.push_back(xn);
            orig_to_new[circ.inputs[i]] = xn;
        }
        // Map outputs → y_j
        for (int j = 0; j < (int)circ.outputs.size(); j++) {
            std::string yn = "y_" + std::to_string(j);
            y_names.push_back(yn);
            orig_to_new[circ.outputs[j]] = yn;
        }
        // Map internal signals → t_k
        for (auto& st : circ.stmts) {
            if (orig_to_new.find(st.name) != orig_to_new.end()) continue;
            std::string tn = "t_" + std::to_string(t_count);
            t_names.push_back(tn);
            orig_to_new[st.name] = tn;
            t_count++;
        }
    }

    std::string map(const std::string& name) const {
        auto it = orig_to_new.find(name);
        return (it != orig_to_new.end()) ? it->second : name;
    }
};

// ====================================================================
//  Write remapped circuit
// ====================================================================

static void write_remapped_circuit(const std::string& path,
                                    const Circuit& circ,
                                    const NameMapper& mapper) {
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }

    // INORDER
    f << "INORDER =";
    for (auto& xn : mapper.x_names) f << " " << xn;
    f << ";\n";

    // OUTORDER
    f << "OUTORDER =";
    for (auto& yn : mapper.y_names) f << " " << yn;
    f << ";\n";

    // Statements
    for (auto& st : circ.stmts) {
        std::string new_name = mapper.map(st.name);
        switch (st.op) {
            case Op::CONST0: f << new_name << " = 0;\n"; break;
            case Op::CONST1: f << new_name << " = 1;\n"; break;
            case Op::INPUT:  f << new_name << " = " << mapper.map(st.arg1) << ";\n"; break;
            case Op::NOT:    f << new_name << " = !" << mapper.map(st.arg1) << ";\n"; break;
            case Op::AND:    f << new_name << " = " << mapper.map(st.arg1) << " * " << mapper.map(st.arg2) << ";\n"; break;
            case Op::XOR:    f << new_name << " = " << mapper.map(st.arg1) << " + " << mapper.map(st.arg2) << ";\n"; break;
        }
    }
}

// ====================================================================
//  Write stats
// ====================================================================

static void write_stats(const std::string& path,
                         const std::string& inst,
                         const Circuit& circ,
                         const NameMapper& mapper) {
    int n = circ.n_inputs;
    int n_out = (int)circ.outputs.size();
    int n_internal = mapper.t_count;

    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }

    f << "# Instance: " << inst << "\n";
    f << "n_inputs = " << n << "\n";
    f << "n_outputs = " << n_out << "\n";
    f << "n_internal = " << n_internal << "\n";
    f << "n_and = " << count_and_gates(circ) << "\n";
    f << "n_xor = " << count_xor_gates(circ) << "\n";
    f << "n_not = " << count_not_gates(circ) << "\n";
    f << "\n";

    int n_const = 0, n_affine = 0, n_nonlinear = 0;

    for (int j = 0; j < n_out; j++) {
        auto& yname = mapper.y_names[j];
        auto& orig = circ.outputs[j];

        f << "[" << yname << "]\n";
        f << "original = " << orig << "\n";

        // Check if it's a direct constant
        bool is_const = false;
        int const_val = 0;
        // Try to trace through aliases to find CONST0/CONST1
        {
            // First check if the output name is assigned to a constant directly
            for (auto& st : circ.stmts) {
                if (st.name == orig && (st.op == Op::CONST0 || st.op == Op::CONST1)) {
                    is_const = true;
                    const_val = (st.op == Op::CONST1) ? 1 : 0;
                    break;
                }
            }
        }

        if (is_const) {
            f << "type = constant\n";
            f << "value = " << const_val << "\n";
            n_const++;
        } else if (n <= 16) {
            // Exact ANF via truth table + Möbius
            int size = 1 << n;
            std::vector<uint64_t> tt(size, 0);
            for (int x = 0; x < size; x++) {
                auto r = eval_circuit_point(circ, x);
                tt[x] = r[j];
            }
            mobius_inplace(tt, n);
            int T = count_nonlinear(tt, n);
            int deg = max_degree(tt, n);

            if (T == 0) {
                // Affine (or constant, but we already checked constant)
                // Compute (a, b)
                uint32_t a = 0;
                int b = (int)tt[0];  // constant term
                for (int i = 0; i < n; i++)
                    if (tt[1 << i]) a |= (1U << i);
                f << "type = affine\n";
                f << "mask = 0b";
                for (int i = n - 1; i >= 0; i--)
                    f << ((a >> i) & 1);
                f << "\n";
                f << "b = " << b << "\n";
                n_affine++;
            } else {
                f << "type = nonlinear\n";
                f << "T_raw = " << T << "\n";
                f << "deg = " << deg << "\n";
                n_nonlinear++;
            }
        } else {
            // n = 32: BLR test for affine detection
            uint32_t a = 0;
            int b = 0;
            if (detect_affine(circ, j, a, b)) {
                f << "type = affine\n";
                f << "mask = 0b";
                for (int i = n - 1; i >= 0; i--)
                    f << ((a >> i) & 1);
                f << "\n";
                f << "b = " << b << "\n";
                n_affine++;
            } else {
                f << "type = nonlinear\n";
                f << "T_raw = N/A (n=32)\n";
                n_nonlinear++;
            }
        }
        f << "\n";
    }

    f << "# Summary\n";
    f << "constant = " << n_const << "\n";
    f << "affine = " << n_affine << "\n";
    f << "nonlinear = " << n_nonlinear << "\n";
}

// ====================================================================
//  Main
// ====================================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [--out-dir DIR]\n";
        return 1;
    }

    std::string circ_path;
    std::string out_dir = "";

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--out-dir" && i + 1 < argc)
            out_dir = argv[++i];
        else if (arg[0] != '-')
            circ_path = arg;
    }

    if (circ_path.empty()) {
        std::cerr << "  ERROR: no circuit file specified\n";
        return 1;
    }

    std::string inst = strip_ext(fs::path(circ_path).filename().string());
    if (out_dir.empty())
        out_dir = "preprocessed/" + inst;

    printf("Preprocess: %s\n", inst.c_str());
    printf("  input:  %s\n", circ_path.c_str());
    printf("  output: %s/\n", out_dir.c_str());

    Circuit circ = read_circuit(circ_path);
    printf("  Circuit: n=%d, outputs=%d, stmts=%zu\n",
           circ.n_inputs, (int)circ.outputs.size(), circ.stmts.size());

    NameMapper mapper;
    mapper.build(circ);
    printf("  Internal signals: %d\n", mapper.t_count);

    fs::create_directories(out_dir);

    // Write remapped circuit
    std::string circ_out = out_dir + "/" + inst + ".txt";
    write_remapped_circuit(circ_out, circ, mapper);
    printf("  Remapped circuit: %s\n", circ_out.c_str());

    // Write stats
    std::string stats_out = out_dir + "/" + inst + "_stats.txt";
    write_stats(stats_out, inst, circ, mapper);
    printf("  Stats: %s\n", stats_out.c_str());

    printf("  Done.\n");
    return 0;
}
