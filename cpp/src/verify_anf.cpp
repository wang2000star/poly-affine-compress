/**
 * verify_anf.cpp — Verification: f(x) = g(Mx+b) + Cx+d
 *
 * Reads circuit (.eqn), transform (.affine), and ANF (.poly),
 * then verifies the identity for random test points.
 *
 * Usage:
 *   ./verify_anf <circuit.eqn> <transform.affine> <anf.poly> [n_tests] [--output file]
 */
#include "circuit.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

// ====================================================================
//  Data structures
// ====================================================================

struct AffineTransform {
    int s{0}, t{0}, n{0};
    std::vector<uint32_t> M_rows; // s rows, each an n-bit mask (A = Mx+b)
    std::vector<uint32_t> C_rows; // t rows, each an n-bit mask (B = Cx+d)
    std::vector<uint8_t> b;       // s bias bits
    std::vector<uint8_t> d;       // t bias bits
};

struct PolyTerm {
    uint32_t mask;   // which z-variables appear (bits 0..m-1)
    int coeff;       // 0 or 1
};

struct PolyOutput {
    std::vector<PolyTerm> terms;
};

// ====================================================================
//  Helpers
// ====================================================================

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}

// ====================================================================
//  Parse .affine
// ====================================================================

static AffineTransform parse_affine(const std::string& path) {
    AffineTransform tfm;
    std::ifstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot open " << path << "\n"; return tfm; }

    f >> tfm.s >> tfm.t;
    int s_check = 0, t_check = 0;
    f >> s_check >> t_check >> tfm.n;
    if (s_check != tfm.s || t_check != tfm.t) {
        std::cerr << "  ERROR: affine header mismatch (s=" << tfm.s
                  << " t=" << tfm.t << ", but matrix claims "
                  << s_check << " " << t_check << ")\n";
        return tfm;
    }

    tfm.M_rows.resize(tfm.s);
    for (int i = 0; i < tfm.s; i++) {
        uint32_t row = 0;
        for (int b = 0; b < tfm.n; b++) {
            int bit; f >> bit;
            if (bit) row |= (1u << b);
        }
        tfm.M_rows[i] = row;
    }

    tfm.C_rows.resize(tfm.t);
    for (int i = 0; i < tfm.t; i++) {
        uint32_t row = 0;
        for (int b = 0; b < tfm.n; b++) {
            int bit; f >> bit;
            if (bit) row |= (1u << b);
        }
        tfm.C_rows[i] = row;
    }

    tfm.b.resize(tfm.s);
    for (int i = 0; i < tfm.s; i++) {
        int bit; f >> bit;
        tfm.b[i] = (uint8_t)bit;
    }

    tfm.d.resize(tfm.t);
    for (int i = 0; i < tfm.t; i++) {
        int bit; f >> bit;
        tfm.d[i] = (uint8_t)bit;
    }

    return tfm;
}

// ====================================================================
//  Parse .poly (matrix format)
// ====================================================================

static std::vector<PolyOutput> parse_poly(const std::string& path, int& m) {
    std::vector<PolyOutput> outputs;
    std::ifstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot open " << path << "\n"; return outputs; }

    f >> m;
    int k; f >> k;
    std::vector<int> term_counts(k);
    for (int i = 0; i < k; i++) f >> term_counts[i];

    outputs.resize(k);
    for (int oi = 0; oi < k; oi++) {
        outputs[oi].terms.resize(term_counts[oi]);
        for (int ti = 0; ti < term_counts[oi]; ti++) {
            // Parse [e_0, e_1, ..., e_{m-1}, coeff]
            char c; f >> c; // skip '['
            uint32_t mask = 0;
            for (int b = 0; b < m; b++) {
                int bit; f >> bit; f >> c; // read bit and comma/']
                if (bit) mask |= (1u << b);
            }
            // Last value is coefficient
            int coeff; f >> coeff;
            f >> c; // skip ']'
            outputs[oi].terms[ti] = {mask, coeff};
        }
    }
    return outputs;
}

// ====================================================================
//  Evaluate ANF at given z (z = A, s-bit value)
// ====================================================================

static bool eval_poly(const PolyOutput& out, uint32_t z) {
    bool result = false;
    for (const auto& term : out.terms) {
        if (term.coeff == 0) continue;
        if ((z & term.mask) == term.mask)
            result ^= true;
    }
    return result;
}

// ====================================================================
//  Compute A = Mx+b, B = Cx+d
// ====================================================================

static uint32_t compute_A(const AffineTransform& tfm, uint32_t x) {
    uint32_t A = 0;
    for (int j = 0; j < tfm.s; j++) {
        if ((__builtin_popcount(tfm.M_rows[j] & x) & 1) ^ tfm.b[j])
            A |= (1u << j);
    }
    return A;
}

static void compute_B(const AffineTransform& tfm, uint32_t x,
                       std::vector<uint8_t>& B) {
    B.resize(tfm.t);
    for (int j = 0; j < tfm.t; j++) {
        B[j] = (__builtin_popcount(tfm.C_rows[j] & x) & 1) ^ tfm.d[j];
    }
}

// ====================================================================
//  Main
// ====================================================================

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <circuit.eqn> <transform.affine> <anf.poly>"
                  << " [n_tests] [--output file]\n";
        return 1;
    }

    std::string circ_path = argv[1];
    std::string affine_path = argv[2];
    std::string poly_path = argv[3];

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
    int k = (int)circ.outputs.size();
    printf("  Circuit: %s (n=%d, outputs=%d)\n",
           circ_path.c_str(), n, k);

    // ---- Read affine transform ----
    AffineTransform tfm = parse_affine(affine_path);
    if (tfm.s == 0 && tfm.t == 0) {
        std::cerr << "  ERROR: failed to parse " << affine_path << "\n";
        return 1;
    }
    if (tfm.n != n) {
        std::cerr << "  ERROR: transform n=" << tfm.n
                  << " != circuit n=" << n << "\n";
        return 1;
    }
    if (tfm.t != k) {
        std::cerr << "  ERROR: transform t=" << tfm.t
                  << " != circuit outputs=" << k << "\n";
        return 1;
    }
    printf("  Transform: s=%d, t=%d, n=%d\n", tfm.s, tfm.t, tfm.n);

    // ---- Read ANF poly ----
    int m = 0;
    auto poly_outputs = parse_poly(poly_path, m);
    if (poly_outputs.empty()) {
        std::cerr << "  ERROR: failed to parse " << poly_path << "\n";
        return 1;
    }
    if (m != tfm.s) {
        std::cerr << "  ERROR: poly m=" << m
                  << " != transform s=" << tfm.s << "\n";
        return 1;
    }
    if ((int)poly_outputs.size() != k) {
        std::cerr << "  ERROR: poly outputs=" << poly_outputs.size()
                  << " != circuit outputs=" << k << "\n";
        return 1;
    }
    printf("  ANF: m=%d, %zu outputs\n", m, poly_outputs.size());

    // ---- Run verification ----
    uint32_t mask = (n == 32) ? 0xFFFFFFFF : ((1u << n) - 1);
    std::vector<int> mismatch_per_output(k, 0);

    std::ofstream f_out;
    std::ostream* out = &std::cout;
    if (!output_path.empty()) {
        f_out.open(output_path);
        if (!f_out) { std::cerr << "  ERROR: cannot write " << output_path << "\n"; return 1; }
        out = &f_out;
    }

    *out << "# Verification: f(x) = g(Mx+b) + Cx+d\n";
    *out << "# n=" << n << ", s=" << tfm.s << ", t=" << tfm.t << "\n";
    *out << "# Tests: " << n_tests << "\n\n";

    std::vector<uint8_t> B;
    for (int test = 0; test < n_tests; test++) {
        uint32_t x = (uint32_t)(((uint64_t)test * 0x9e3779b97f4a7c15ULL) & mask);
        uint32_t A = compute_A(tfm, x);
        compute_B(tfm, x, B);

        auto f_vals = eval_circuit_point(circ, x);

        for (int j = 0; j < k; j++) {
            uint8_t f_val = f_vals[j];
            uint8_t g_val = eval_poly(poly_outputs[j], A) ? 1 : 0;
            uint8_t rhs = g_val ^ B[j];
            if (f_val != rhs)
                mismatch_per_output[j]++;
        }
    }

    // ---- Report results ----
    int n_fail = 0;
    for (int j = 0; j < k; j++) {
        if (mismatch_per_output[j] == 0) {
            *out << "output " << j << ": PASS (" << n_tests << " tests)\n";
        } else {
            *out << "output " << j << ": FAIL (" << mismatch_per_output[j]
                 << "/" << n_tests << ")\n";
            n_fail++;
        }
    }
    *out << "\n";

    if (n_fail == 0) {
        *out << "All outputs PASS (" << n_tests << " tests)\n";
        printf("  Verification: ✅ PASS (%d outputs, %d tests)\n", k, n_tests);
    } else {
        *out << n_fail << "/" << k << " outputs FAIL\n";
        printf("  Verification: ❌ FAIL (%d/%d outputs failed, %d tests)\n",
               n_fail, k, n_tests);
        for (int j = 0; j < k; j++) {
            if (mismatch_per_output[j] > 0)
                printf("    output %d: %d/%d FAIL\n",
                       j, mismatch_per_output[j], n_tests);
        }
    }

    if (!output_path.empty())
        printf("  Saved: %s\n", output_path.c_str());

    return (n_fail == 0) ? 0 : 1;
}
