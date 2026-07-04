#include "circuit.h"
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

// ============================================================
//  ANF term: a product (AND) of a subset of z-variables
// ============================================================

struct AnfTerm {
    std::vector<int> vars;  // indices of z_i; empty = constant 1
};

struct AnfOutput {
    std::string name;
    std::vector<AnfTerm> terms;
};

// ============================================================
//  Affine transform z = Mx + b
// ============================================================

struct Transform {
    int m = 0, n = 0;
    uint32_t M_rows[32] = {0};
    uint64_t b = 0;
};

// ============================================================
//  Helpers
// ============================================================

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}

// ============================================================
//  Parse expr.poly — extract output ANF expressions
// ============================================================

static std::vector<AnfOutput> parse_expr_poly(const std::string& path, int& m) {
    std::vector<AnfOutput> outputs;
    std::ifstream f(path);
    if (!f) { std::cerr << "ERROR: cannot open " << path << "\n"; return {}; }

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (m == 0) {
                auto p = line.find("m=");
                if (p != std::string::npos) m = std::stoi(line.substr(p + 2));
            }
            continue;
        }
        auto eq = line.find(" = ");
        if (eq == std::string::npos) continue;

        AnfOutput out;
        out.name = line.substr(0, eq);
        std::string expr = line.substr(eq + 3);

        std::string trimmed_expr = trim(expr);
        if (trimmed_expr == "0") {
            // "= 0" means constant zero → no terms, eval returns false
        } else {
            std::stringstream ss(expr);
            std::string term_str;
            while (std::getline(ss, term_str, '+')) {
                term_str = trim(term_str);
                if (term_str.empty()) continue;

                AnfTerm term;
                if (term_str != "1") {
                    size_t pos = 0;
                    while (pos < term_str.size()) {
                        auto zpos = term_str.find("z_", pos);
                        if (zpos == std::string::npos) break;
                        pos = zpos + 2;
                        size_t end = pos;
                        while (end < term_str.size() && std::isdigit((unsigned char)term_str[end]))
                            end++;
                        term.vars.push_back(std::stoi(term_str.substr(pos, end - pos)));
                        pos = end;
                    }
                }
                out.terms.push_back(term);
            }
        }
        outputs.push_back(out);
    }
    return outputs;
}

// ============================================================
//  Parse trans.poly — extract M, b, m, n
// ============================================================

static Transform parse_trans_poly(const std::string& path,
                                  const std::vector<std::string>& input_names)
{
    Transform t;
    std::ifstream f(path);
    if (!f) { std::cerr << "ERROR: cannot open " << path << "\n"; return t; }

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.find("m=") != std::string::npos && t.m == 0) {
                auto pm = line.find("m=");
                auto pn = line.find("n=");
                t.m = std::stoi(line.substr(pm + 2));
                if (pn != std::string::npos) t.n = std::stoi(line.substr(pn + 2));
            }
            continue;
        }
        auto eq = line.find(" = ");
        if (eq == std::string::npos) continue;

        // lhs = "z_i"
        int row = std::stoi(line.substr(2, eq - 2));

        std::string rhs = line.substr(eq + 3);
        std::stringstream ss(rhs);
        std::string token;
        while (ss >> token) {
            if (token == "+") continue;
            if (token == "1") {
                t.b |= (1ULL << row);
            } else {
                for (int i = 0; i < (int)input_names.size(); i++) {
                    if (token == input_names[i]) {
                        t.M_rows[row] |= (1u << i);
                        break;
                    }
                }
            }
        }
    }
    return t;
}

// ============================================================
//  Evaluate ANF at a given z
// ============================================================

static bool eval_anf(const AnfOutput& out, uint32_t z) {
    bool result = false;
    for (const auto& term : out.terms) {
        bool val = true;
        for (int v : term.vars)
            if (!((z >> v) & 1)) { val = false; break; }
        result ^= val;
    }
    return result;
}

// ============================================================
//  Main: standalone ANF verification
//
//  Usage:
//    verify_anf <circuit.txt> <trans.poly> <expr.poly> [n_tests] [--output verify.txt]
//
//  Verifies f(x) == g(Mx+b) for n_tests random inputs,
//  using the original circuit and the optimized ANF expression.
//
//  Output columns:
//    x (n bits)  f(x) (k bits)  z (m bits)  g(z) (k bits)
// ============================================================

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <circuit.txt> <trans.poly> <expr.poly> [n_tests] [--output file]\n";
        return 1;
    }

    std::string circ_path = argv[1];
    std::string trans_path = argv[2];
    std::string expr_path = argv[3];

    int n_tests = 5000;
    std::string output_path;
    for (int a = 4; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "--output" && a + 1 < argc)
            output_path = argv[++a];
        else
            n_tests = std::stoi(arg);
    }

    // ---- Read circuit ----
    Circuit circ = read_circuit(circ_path);
    int n = circ.n_inputs;

    // ---- Read expr ----
    int m = 0;
    auto outputs = parse_expr_poly(expr_path, m);
    if (outputs.empty()) {
        std::cerr << "ERROR: no ANF outputs found in " << expr_path << "\n";
        return 1;
    }

    // ---- Read transform ----
    Transform t = parse_trans_poly(trans_path, circ.inputs);
    if (t.m == 0) { std::cerr << "ERROR: failed to parse " << trans_path << "\n"; return 1; }
    m = (m > 0) ? m : t.m;

    // ---- Map output names to circuit indices ----
    std::vector<int> out_to_circ(outputs.size(), -1);
    for (size_t oi = 0; oi < outputs.size(); oi++) {
        for (int ci = 0; ci < (int)circ.outputs.size(); ci++) {
            if (circ.outputs[ci] == outputs[oi].name) {
                out_to_circ[oi] = ci;
                break;
            }
        }
        if (out_to_circ[oi] < 0)
            std::cerr << "WARNING: output " << outputs[oi].name << " not found in circuit\n";
    }

    // ---- Run verification ----
    uint32_t mask = (n == 32) ? 0xFFFFFFFF : ((1u << n) - 1);
    int k = (int)outputs.size();
    int n_mismatch = 0;

    std::ofstream f_out;
    std::ostream* out = &std::cout;
    if (!output_path.empty()) {
        f_out.open(output_path);
        if (!f_out) { std::cerr << "ERROR: cannot write " << output_path << "\n"; return 1; }
        out = &f_out;
    }

    *out << "# Verification for (n=" << n << ", k=" << k << " outputs)\n";
    *out << "# Transform: z = Mx + b (m=" << t.m << ")\n";
    *out << "# Tests: " << n_tests << "\n\n";
    *out << "# x  f(x)  z  g(z)\n";

    auto bits_to_str = [](const std::vector<bool>& bits) -> std::string {
        std::string s(bits.size(), '0');
        for (size_t i = 0; i < bits.size(); i++)
            if (bits[i]) s[bits.size() - 1 - i] = '1';
        return s;
    };
    auto int_to_bits = [bits_to_str](uint32_t val, int nbits) -> std::string {
        std::vector<bool> bits(nbits, false);
        for (int i = 0; i < nbits; i++) if ((val >> i) & 1) bits[i] = true;
        return bits_to_str(bits);
    };

    for (int test = 0; test < n_tests; test++) {
        uint32_t x = (uint32_t)(((uint64_t)test * 0x9e3779b97f4a7c15ULL) & mask);
        uint64_t z = t.b;
        for (int row = 0; row < t.m; row++) {
            if (__builtin_popcount(t.M_rows[row] & x) & 1)
                z ^= (1ULL << row);
        }

        auto circ_vals = eval_circuit_point(circ, x);
        std::vector<bool> f_vec(k, false), g_vec(k, false);
        for (int oi = 0; oi < k; oi++) {
            if (out_to_circ[oi] >= 0 && (size_t)out_to_circ[oi] < circ_vals.size())
                f_vec[oi] = circ_vals[out_to_circ[oi]];
            if (eval_anf(outputs[oi], z))
                g_vec[oi] = true;
        }

        bool match = (f_vec == g_vec);
        *out << int_to_bits(x, n) << "  "
             << bits_to_str(f_vec) << "  "
             << int_to_bits(z, t.m) << "  "
             << bits_to_str(g_vec) << "\n";

        if (!match) n_mismatch++;
    }

    *out << "\n";
    if (n_mismatch == 0) {
        *out << "✅ All outputs verified! (" << n_tests << " tests, 0 mismatches)\n";
    } else {
        *out << "❌ " << n_mismatch << " mismatches out of " << n_tests << " tests\n";
    }

    if (!output_path.empty())
        std::cout << "  Saved: " << output_path << " (" << n_mismatch << "/" << n_tests << " mismatches)\n";

    return (n_mismatch == 0) ? 0 : 1;
}
