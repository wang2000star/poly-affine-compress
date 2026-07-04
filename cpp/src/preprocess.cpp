/**
 * preprocess.cpp — One-time circuit preprocessing
 *
 * Reads raw circuit, remaps variables, classifies outputs, saves fixed files.
 *
 * Usage:
 *   ./preprocess <circuit.txt> [--out-dir DIR]
 *
 * Output (in --out-dir, default preprocessed/{inst}/):
 *   {inst}.eqn           — remapped circuit
 *   {inst}_stats.txt     — 5-line statistics
 *   {inst}_const.txt     — constant outputs
 *   {inst}_affine.txt    — affine outputs
 *   {inst}_nonlinear.txt — nonlinear outputs
 *   {inst}.poly          — raw polynomial (n≤16 only)
 */
#include "circuit.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
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

// Build truth table for a single output (n ≤ 16)
static std::vector<uint64_t> build_tt(const Circuit& circ, int output_idx, int n) {
    int size = 1 << n;
    std::vector<uint64_t> tt(size, 0);
    for (int x = 0; x < size; x++) {
        auto r = eval_circuit_point(circ, x);
        tt[x] = r[output_idx];
    }
    return tt;
}

// Affine detection via BLR test (any n ≤ 32)
// Returns coefficients (a, b) such that f(x) = ⟨a, x⟩ ⊕ b
// Returns false if nonlinear
static bool detect_affine(const Circuit& circ, int output_idx,
                           uint32_t& a, int& b) {
    int n = circ.n_inputs;
    // b = f(0)
    auto r0 = eval_circuit_point(circ, 0);
    b = r0[output_idx];
    // a_i = f(e_i) ⊕ f(0)
    a = 0;
    for (int i = 0; i < n; i++) {
        auto ri = eval_circuit_point(circ, 1U << i);
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
        auto fx    = eval_circuit_point(circ, x);
        auto fx_u  = eval_circuit_point(circ, x ^ u);
        auto fx_v  = eval_circuit_point(circ, x ^ v);
        auto fx_uv = eval_circuit_point(circ, x ^ u ^ v);
        if (fx[output_idx] ^ fx_u[output_idx] ^
            fx_v[output_idx] ^ fx_uv[output_idx])
            return false;
    }
    // Verify coefficients: 200 random points
    for (int t = 0; t < 200; t++) {
        uint32_t x = (uint32_t)rng() & mask;
        auto fx = eval_circuit_point(circ, x)[output_idx];
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
    std::vector<std::string> x_names;
    std::vector<std::string> y_names;
    std::vector<std::string> t_names;      // new names t_0, t_1, ...
    std::vector<std::string> t_orig_names; // original internal variable names
    int t_count{0};

    void build(const Circuit& circ) {
        for (int i = 0; i < (int)circ.inputs.size(); i++) {
            std::string xn = "x_" + std::to_string(i);
            x_names.push_back(xn);
            orig_to_new[circ.inputs[i]] = xn;
        }
        for (int j = 0; j < (int)circ.outputs.size(); j++) {
            std::string yn = "y_" + std::to_string(j);
            y_names.push_back(yn);
            orig_to_new[circ.outputs[j]] = yn;
        }
        for (auto& st : circ.stmts) {
            if (orig_to_new.find(st.name) != orig_to_new.end()) continue;
            std::string tn = "t_" + std::to_string(t_count);
            t_names.push_back(tn);
            t_orig_names.push_back(st.name);
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
//  Write remapped circuit (.eqn)
// ====================================================================

static void write_circuit_eqn(const std::string& path,
                               const Circuit& circ,
                               const NameMapper& mapper) {
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    f << "INORDER =";
    for (auto& xn : mapper.x_names) f << " " << xn;
    f << ";\n";
    f << "OUTORDER =";
    for (auto& yn : mapper.y_names) f << " " << yn;
    f << ";\n";
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
//  Write variable mapping (.bij)
//  Three sections: inputs (orig → x_i), outputs (orig → y_j), internal (orig → t_k)
// ====================================================================

static void write_bij(const std::string& path,
                       const Circuit& circ,
                       const NameMapper& mapper) {
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    // Inputs
    f << (int)circ.inputs.size() << "\n";
    for (auto& name : circ.inputs)
        f << name << "\n";
    // Outputs
    f << (int)circ.outputs.size() << "\n";
    for (auto& name : circ.outputs)
        f << name << "\n";
    // Internal (t_0, t_1, ...)
    f << mapper.t_count << "\n";
    for (auto& name : mapper.t_orig_names)
        f << name << "\n";
}

// ====================================================================
//  Write stats (5 lines)
// ====================================================================

static void write_stats(const std::string& path,
                         const Circuit& circ,
                         const NameMapper& mapper,
                         int sum_T, int union_T, int max_deg) {
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    f << circ.n_inputs << "\n";
    f << (int)circ.outputs.size() << "\n";
    f << sum_T << "\n";
    f << union_T << "\n";
    f << max_deg << "\n";
}

// ====================================================================
//  Write const file
// ====================================================================

static void write_const(const std::string& path,
                         const std::vector<int>& const_indices,
                         const std::vector<int>& const_values) {
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    int k = (int)const_indices.size();
    f << k << "\n";
    if (k == 0) return;
    for (int i = 0; i < k; i++) {
        if (i > 0) f << " ";
        f << const_indices[i];
    }
    f << "\n";
    for (int i = 0; i < k; i++) {
        if (i > 0) f << " ";
        f << const_values[i];
    }
    f << "\n";
}

// ====================================================================
//  Write affine file
// ====================================================================

static void write_affine(const std::string& path,
                          const std::vector<int>& affine_indices,
                          const std::vector<uint32_t>& affine_masks,
                          const std::vector<int>& affine_b,
                          int n) {
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    int k = (int)affine_indices.size();
    f << k << "\n";
    if (k == 0) return;
    // Line 2: indices
    for (int i = 0; i < k; i++) {
        if (i > 0) f << " ";
        f << affine_indices[i];
    }
    f << "\n";
    // Line 3: rows cols
    f << k << " " << (n + 1) << "\n";
    // Matrix rows: [x_0, x_1, ..., x_{n-1}, b]
    for (int i = 0; i < k; i++) {
        for (int b = 0; b < n; b++)
            f << ((affine_masks[i] >> b) & 1) << " ";
        f << affine_b[i] << "\n";
    }
}

// ====================================================================
//  Write nonlinear file
// ====================================================================

static void write_nonlinear(const std::string& path,
                             const std::vector<int>& nonlinear_indices,
                             const std::vector<int>& nonlinear_T,
                             const std::vector<int>& nonlinear_deg) {
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    int k = (int)nonlinear_indices.size();
    f << k << "\n";
    if (k == 0) return;
    for (int i = 0; i < k; i++) {
        if (i > 0) f << " ";
        f << nonlinear_indices[i];
    }
    f << "\n";
    // n=32: no T/deg info available, stop after indices
    if (nonlinear_T.empty()) return;
    for (int i = 0; i < k; i++) {
        if (i > 0) f << " ";
        f << nonlinear_T[i];
    }
    f << "\n";
    for (int i = 0; i < k; i++) {
        if (i > 0) f << " ";
        f << nonlinear_deg[i];
    }
    f << "\n";
}

// ====================================================================
//  Write .poly (n ≤ 16 only)
// ====================================================================

static void write_poly(const std::string& path,
                        const Circuit& circ, int n) {
    std::ofstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot write " << path << "\n"; return; }
    int n_out = (int)circ.outputs.size();
    f << n << "\n";
    f << n_out << "\n";

    // Compute ANF for each output, collect term counts and terms
    std::vector<std::vector<uint64_t>> terms_list;
    for (int j = 0; j < n_out; j++) {
        auto tt = build_tt(circ, j, n);
        mobius_inplace(tt, n);
        // Collect terms with coeff=1
        int size = 1 << n;
        int cnt = 0;
        std::vector<uint64_t> terms;
        for (int x = 0; x < size; x++) {
            if (tt[x]) {
                cnt++;
                terms.push_back((uint64_t)x | (1ULL << n)); // append coeff=1 bit
            }
        }
        f << cnt;
        if (j < n_out - 1) f << " ";
        terms_list.push_back(terms);
    }
    f << "\n";

    // Write all monomial rows
    for (auto& terms : terms_list) {
        for (uint64_t entry : terms) {
            uint64_t mask = entry & ((1ULL << n) - 1);
            int coeff = (entry >> n) & 1;
            if (!coeff) continue;
            f << "[";
            for (int b = 0; b < n; b++) {
                f << ((mask >> b) & 1);
                f << ",";
            }
            f << coeff << "]\n";
        }
    }
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
    int n = circ.n_inputs;
    int n_out = (int)circ.outputs.size();

    NameMapper mapper;
    mapper.build(circ);

    printf("  Circuit: n=%d, outputs=%d, internal=%d\n",
           n, n_out, mapper.t_count);

    fs::create_directories(out_dir);

    // ── Write .eqn (remapped circuit) ──────────────────────────────
    std::string eqn_path = out_dir + "/" + inst + ".eqn";
    write_circuit_eqn(eqn_path, circ, mapper);
    printf("  .eqn: %s\n", eqn_path.c_str());

    // ── Write .bij (variable mapping) ─────────────────────────────
    std::string bij_path = out_dir + "/" + inst + ".bij";
    write_bij(bij_path, circ, mapper);
    printf("  .bij: %s\n", bij_path.c_str());

    // ── Classify outputs ──────────────────────────────────────────
    // Determine constant/affine/nonlinear for each output

    std::vector<int> const_indices, const_values;
    std::vector<int> affine_indices;
    std::vector<uint32_t> affine_masks;
    std::vector<int> affine_b;
    std::vector<int> nonlinear_indices;

    // For n≤16: full stats via Möbius
    std::vector<int> nonlinear_T, nonlinear_deg;
    int sum_T = 0, union_T = 0, max_deg_out = 0;

    if (n <= 16) {
        // Compute union_T: XOR all output truth tables, then Möbius
        int union_size = 1 << n;
        std::vector<uint64_t> union_tt(union_size, 0);

        for (int j = 0; j < n_out; j++) {
            auto tt = build_tt(circ, j, n);
            mobius_inplace(tt, n);
            int T = count_nonlinear(tt, n);
            int deg = max_degree(tt, n);
            sum_T += T;
            if (deg > max_deg_out) max_deg_out = deg;

            // Check if constant
            bool is_const = true;
            for (int x = 0; x < union_size; x++) {
                if (tt[x] && __builtin_popcount(x) >= 1) {
                    is_const = false;
                    break;
                }
            }
            if (is_const) {
                const_indices.push_back(j);
                const_values.push_back((int)tt[0]); // constant value
                continue;
            }

            // Check if affine (T==0 means only constant + linear terms)
            if (T == 0) {
                uint32_t a = 0;
                for (int i = 0; i < n; i++)
                    if (tt[1U << i]) a |= (1U << i);
                affine_indices.push_back(j);
                affine_masks.push_back(a);
                affine_b.push_back((int)tt[0]);
                continue;
            }

            nonlinear_indices.push_back(j);
            nonlinear_T.push_back(T);
            nonlinear_deg.push_back(deg);

            // Union TT: OR of all output ANF masks (nonlinear only)
            for (int x = 0; x < union_size; x++)
                if (tt[x] && __builtin_popcount(x) >= 2)
                    union_tt[x] = 1;
        }
        for (int x = 0; x < union_size; x++)
            if (union_tt[x]) union_T++;
    } else {
        // n = 32: use BLR for affine detection
        for (int j = 0; j < n_out; j++) {
            // Check constant by checking stmts
            bool is_const = false;
            int const_val = 0;
            for (auto& st : circ.stmts) {
                if (st.name == circ.outputs[j] &&
                    (st.op == Op::CONST0 || st.op == Op::CONST1)) {
                    is_const = true;
                    const_val = (st.op == Op::CONST1) ? 1 : 0;
                    break;
                }
            }
            if (is_const) {
                const_indices.push_back(j);
                const_values.push_back(const_val);
                continue;
            }

            uint32_t a = 0;
            int b = 0;
            if (detect_affine(circ, j, a, b)) {
                affine_indices.push_back(j);
                affine_masks.push_back(a);
                affine_b.push_back(b);
                continue;
            }

            nonlinear_indices.push_back(j);
        }
        sum_T = -1;
        union_T = -1;
        max_deg_out = -1;
    }

    // ── Write _stats.txt ──────────────────────────────────────────
    std::string stats_path = out_dir + "/" + inst + "_stats.txt";
    write_stats(stats_path, circ, mapper, sum_T, union_T, max_deg_out);
    printf("  stats: %s\n", stats_path.c_str());

    // ── Write const.txt ───────────────────────────────────────────
    std::string const_path = out_dir + "/" + inst + "_const.txt";
    write_const(const_path, const_indices, const_values);
    printf("  const: %s (%d)\n", const_path.c_str(), (int)const_indices.size());

    // ── Write affine.txt ──────────────────────────────────────────
    std::string affine_path = out_dir + "/" + inst + "_affine.txt";
    write_affine(affine_path, affine_indices, affine_masks, affine_b, n);
    printf("  affine: %s (%d)\n", affine_path.c_str(), (int)affine_indices.size());

    // ── Write nonlinear.txt ───────────────────────────────────────
    std::string nonlinear_path = out_dir + "/" + inst + "_nonlinear.txt";
    write_nonlinear(nonlinear_path, nonlinear_indices, nonlinear_T, nonlinear_deg);
    printf("  nonlinear: %s (%d)\n", nonlinear_path.c_str(), (int)nonlinear_indices.size());

    // ── Write .poly (n ≤ 16 only) ─────────────────────────────────
    if (n <= 16) {
        std::string poly_path = out_dir + "/" + inst + ".poly";
        write_poly(poly_path, circ, n);
        printf("  .poly: %s\n", poly_path.c_str());
    }

    printf("  Done.\n");
    return 0;
}
