/**
 * verify_anf.cpp — Verification: f(x) = g(Mx+b)
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
    int s{0}, n{0};
    std::vector<uint32_t> M_rows; // s rows, each an n-bit mask
    std::vector<uint8_t> b;       // s bias bits
};

struct PolyTerm {
    uint32_t mask;   // which z-variables appear (bits 0..m-1)
    int coeff;       // 0 or 1
};

struct PolyOutput {
    std::vector<PolyTerm> terms;
};

// ====================================================================
//  Parse .affine
// ====================================================================

static AffineTransform parse_affine(const std::string& path) {
    AffineTransform tfm;
    std::ifstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot open " << path << "\n"; return tfm; }

    f >> tfm.s;
    int rows = 0, cols = 0;
    f >> rows >> cols;
    tfm.n = cols - 1;
    if (rows != tfm.s || cols != tfm.n + 1) {
        std::cerr << "  ERROR: affine header mismatch (s=" << tfm.s
                  << ", matrix claims " << rows << "×" << cols
                  << ", expected " << tfm.s << "×" << (tfm.n + 1) << ")\n";
        return tfm;
    }

    // s × (n+1) matrix: [M | b]
    tfm.M_rows.resize(tfm.s);
    tfm.b.resize(tfm.s);
    for (int i = 0; i < tfm.s; i++) {
        uint32_t row = 0;
        for (int col = 0; col < tfm.n; col++) {
            int bit; f >> bit;
            if (bit) row |= (1u << col);
        }
        tfm.M_rows[i] = row;
        int b_bit; f >> b_bit;
        tfm.b[i] = (uint8_t)b_bit;
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
            int coeff; f >> coeff;
            f >> c; // skip ']'
            outputs[oi].terms[ti] = {mask, coeff};
        }
    }
    return outputs;
}

// ====================================================================
//  Evaluate ANF at given z-value
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
    if (tfm.s == 0) {
        std::cerr << "  ERROR: failed to parse " << affine_path << "\n";
        return 1;
    }
    if (tfm.n != n) {
        std::cerr << "  ERROR: transform n=" << tfm.n
                  << " != circuit n=" << n << "\n";
        return 1;
    }
    printf("  Transform: s=%d, n=%d\n", tfm.s, tfm.n);

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

    // ---- Verification: f(x) == g(Mx+b) ----
    // For n ≤ 16, exhaustive enumeration is cheap and guarantees correctness.
    int n_total = (n <= 16) ? (1 << n) : n_tests;
    uint32_t mask = (n == 32) ? 0xFFFFFFFF : ((1u << n) - 1);
    std::vector<int> mismatch_per_output(k, 0);

    std::ofstream f_out;
    std::ostream* out = &std::cout;
    if (!output_path.empty()) {
        f_out.open(output_path);
        if (!f_out) { std::cerr << "  ERROR: cannot write " << output_path << "\n"; return 1; }
        out = &f_out;
    }

    *out << "# Verification: f(x) = g(Mx+b)\n";
    *out << "# n=" << n << ", s=" << tfm.s << "\n";
    if (n <= 16)
        *out << "# Exhaustive: 2^" << n << " = " << n_total << " inputs\n\n";
    else
        *out << "# Tests: " << n_total << "\n\n";

    for (int test = 0; test < n_total; test++) {
        uint32_t x = (n <= 16) ? (uint32_t)test
                               : (uint32_t)(((uint64_t)test * 0x9e3779b97f4a7c15ULL) & mask);
        // z = Mx+b
        uint32_t z = 0;
        for (int j = 0; j < tfm.s; j++) {
            if ((__builtin_popcount(tfm.M_rows[j] & x) & 1) ^ tfm.b[j])
                z |= (1u << j);
        }

        auto f_vals = eval_circuit_point(circ, x);

        for (int j = 0; j < k; j++) {
            uint8_t f_val = f_vals[j];
            uint8_t g_val = eval_poly(poly_outputs[j], z) ? 1 : 0;
            if (f_val != g_val)
                mismatch_per_output[j]++;
        }
    }

    // ---- Report results ----
    int n_fail = 0;
    int n_display = n_total;
    for (int j = 0; j < k; j++) {
        if (mismatch_per_output[j] == 0) {
            *out << "output " << j << ": PASS (" << n_display << " tests)\n";
        } else {
            *out << "output " << j << ": FAIL (" << mismatch_per_output[j]
                 << "/" << n_display << ")\n";
            n_fail++;
        }
    }
    *out << "\n";

    if (n_fail == 0) {
        *out << "All outputs PASS (" << n_display << " tests)\n";
        printf("  Verification: ✅ PASS (%d outputs, %d tests)\n", k, n_display);
    } else {
        *out << n_fail << "/" << k << " outputs FAIL\n";
        printf("  Verification: ❌ FAIL (%d/%d outputs failed, %d tests)\n",
               n_fail, k, n_display);
        for (int j = 0; j < k; j++) {
            if (mismatch_per_output[j] > 0)
                printf("    output %d: %d/%d FAIL\n",
                       j, mismatch_per_output[j], n_display);
        }
    }

    if (!output_path.empty())
        printf("  Saved: %s\n", output_path.c_str());

    return (n_fail == 0) ? 0 : 1;
}
