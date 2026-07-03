#include "io.h"
#include "anf.h"
#include "moebius.h"
#include "verify.h"
#include <fstream>
#include <iostream>
#include <iomanip>

// ============================================================
//  Save raw ANF — _expr.poly
// ============================================================

void save_raw_anf(const TruthTable& tt, const Circuit& circ,
                  const std::vector<int>& output_indices,
                  const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    f << "# Raw ANF for circuit (n=" << tt.n << ", k=" << tt.n_outputs << " outputs)\n";
    if (tt.n <= 16) {
        f << "# Variables:";
        for (int i = tt.n - 1; i >= 0; i--)
            f << " " << circ.inputs[i];
        f << "\n";
    }

    int total_terms = 0;
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        std::string out_name = circ.outputs[output_indices[oi]];
        f << out_name << " = ";

        int64_t n_words_g = (tt.n < 6) ? 1 : (int64_t(1) << (tt.n - 6));
        int count = 0;
        for (int64_t w = 0; w < n_words_g; w++) {
            uint64_t word = tt.tt[oi][w];
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                int64_t pos = (w << 6) | bit;
                if (count > 0) f << " + ";
                if (pos == 0) { f << "1"; count++; continue; }
                bool first = true;
                for (int j = 0; j < tt.n; j++) {
                    if ((pos >> j) & 1) {
                        if (!first) f << " * ";
                        f << circ.inputs[j];
                        first = false;
                    }
                }
                count++;
            }
        }
        if (count == 0) f << "0";
        f << "\n";
        total_terms += count;
    }
    std::cout << "  Saved ANF: " << fname << " (" << total_terms << " terms)\n";
}

// ============================================================
//  Save raw term counts — _T.poly
// ============================================================

static int64_t count_union_T(const TruthTable& tt) {
    if (tt.n_outputs == 0) return 0;
    int64_t n_words = tt.n_words;
    int64_t union_T = 0;
    for (int64_t w = 0; w < n_words; w++) {
        uint64_t or_val = 0;
        for (int oi = 0; oi < tt.n_outputs; oi++)
            or_val |= tt.tt[oi][w];
        union_T += __builtin_popcountll(or_val);
    }
    return union_T;
}

void save_raw_T(const TruthTable& tt, const Circuit& circ,
                const std::vector<int>& output_indices,
                const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    int64_t sum_T = 0;
    for (int oi = 0; oi < tt.n_outputs; oi++)
        sum_T += count_T(tt.tt[oi].data(), tt.n);
    int64_t union_T = count_union_T(tt);

    f << "# Raw ANF term counts for circuit (n=" << tt.n << ", k=" << tt.n_outputs << " outputs)\n";
    f << "# sum T = " << sum_T << "\n";
    f << "# union T = " << union_T << "\n";
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        int64_t T = count_T(tt.tt[oi].data(), tt.n);
        f << circ.outputs[output_indices[oi]] << ": T=" << T << "\n";
    }
    std::cout << "  Saved T: " << fname << " (sum=" << sum_T << ", union=" << union_T << ")\n";
}

// ============================================================
//  Save optimized ANF — _expr.poly (z-space)
// ============================================================

void save_opt_expr(const std::vector<std::vector<uint64_t>>& g_tt_raw,
                   const Circuit& circ,
                   const std::vector<int>& output_indices,
                   int m, const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    f << "# opt1 ANF (z-space, m=" << m << ")\n";
    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    int total_terms = 0;

    for (int oi = 0; oi < (int)g_tt_raw.size(); oi++) {
        std::string out_name = circ.outputs[output_indices[oi]];
        f << out_name << " = ";

        std::vector<uint64_t> anf(g_tt_raw[oi]);
        moebius_packed(anf.data(), m);

        int count = 0;
        for (int64_t w = 0; w < n_words_g; w++) {
            uint64_t word = anf[w];
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                int64_t pos = (w << 6) | bit;
                if (count > 0) f << " + ";
                if (pos == 0) { f << "1"; count++; continue; }
                bool first = true;
                for (int j = 0; j < m; j++) {
                    if ((pos >> j) & 1) {
                        if (!first) f << " * ";
                        f << "z_" << j;
                        first = false;
                    }
                }
                count++;
            }
        }
        if (count == 0) f << "0";
        f << "\n";
        total_terms += count;
    }
    std::cout << "  Saved ANF: " << fname << " (" << total_terms << " terms)\n";
}

// ============================================================
//  Save optimized term counts — _T.poly
// ============================================================

void save_opt_T(const std::vector<std::vector<uint64_t>>& g_tt_raw,
                const Circuit& circ,
                const std::vector<int>& output_indices,
                int m, const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    int64_t sum_T = 0;
    std::vector<int64_t> per_T;

    std::vector<uint64_t> union_data(n_words_g, 0);
    for (int oi = 0; oi < (int)g_tt_raw.size(); oi++) {
        std::vector<uint64_t> anf(g_tt_raw[oi]);
        moebius_packed(anf.data(), m);
        int64_t T = count_T(anf.data(), m);
        per_T.push_back(T);
        sum_T += T;
        for (int64_t w = 0; w < n_words_g; w++)
            union_data[w] |= anf[w];
    }

    int64_t union_T = 0;
    for (int64_t w = 0; w < n_words_g; w++)
        union_T += __builtin_popcountll(union_data[w]);

    f << "# opt1 T(g) (m=" << m << ")\n";
    f << "# sum T = " << sum_T << "\n";
    f << "# union T = " << union_T << "\n";
    for (int oi = 0; oi < (int)g_tt_raw.size(); oi++) {
        f << circ.outputs[output_indices[oi]] << ": T=" << per_T[oi] << "\n";
    }
    std::cout << "  Saved T: " << fname << " (sum=" << sum_T << ", union=" << union_T << ")\n";
}

// ============================================================
//  Save transform — _trans.poly
// ============================================================

void save_trans(const MbCandidate& best, const Circuit& circ,
                int n, const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    f << "# opt1 transform z = Mx + b (GF(2))\n";
    f << "# m=" << best.m << ", n=" << n << "\n";
    for (int row = 0; row < best.m; row++) {
        f << "z_" << row << " = ";
        bool first = true;
        for (int i = 0; i < n; i++) {
            if ((best.M_rows[row] >> i) & 1) {
                if (!first) f << " + ";
                f << circ.inputs[i];
                first = false;
            }
        }
        if ((best.b >> row) & 1) {
            if (!first) f << " + ";
            f << "1";
        } else if (first) {
            f << "0";
        }
        f << "\n";
    }
    std::cout << "  Saved: " << fname << "\n";
}

// ============================================================
//  Save verification — _verify.txt
// ============================================================

void save_verify(const TruthTable& tt,
                 const std::vector<std::vector<uint64_t>>& g_tt_raw,
                 const MbCandidate& best, int n,
                 const std::vector<int>& output_indices,
                 const Circuit& circ,
                 const std::string& fname)
{
    std::ofstream f_ver(fname);
    if (!f_ver) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    f_ver << "# Verification for (n=" << n << ", k=" << tt.n_outputs << " outputs)\n";
    f_ver << "# Strategy: opt1\n";
    f_ver << "# Transform: z = Mx + b (m=" << best.m << ")\n";
    f_ver << "# Tests: 5000\n\n";

    int n_tests = 5000;
    uint32_t mask = (n == 32) ? 0xFFFFFFFF : ((1u << n) - 1);
    int n_mismatch = 0;
    int first_mismatch_oi = -1;

    f_ver << "# x  f(x)  z  g(z)\n";

    bool has_g_data = !g_tt_raw.empty() && g_tt_raw.size() == (size_t)tt.n_outputs;

    for (int t = 0; t < n_tests; t++) {
        uint32_t x = (uint32_t)(((uint64_t)t * 0x9e3779b97f4a7c15ULL) & mask);
        uint32_t z = best.b;
        for (int row = 0; row < best.m; row++) {
            if (__builtin_popcount(best.M_rows[row] & x) & 1)
                z ^= (1u << row);
        }

        uint64_t f_vec = pack_output_vector(tt, x);
        uint64_t g_vec = 0;
        if (has_g_data) {
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                const auto& g_tt = g_tt_raw[oi];
                uint64_t bit = (g_tt[z >> 6] >> (z & 63)) & 1;
                g_vec |= (bit << oi);
            }
        }

        auto to_bits = [](uint32_t val, int bits) -> std::string {
            std::string s(bits, '0');
            for (int b = bits - 1; b >= 0; b--, val >>= 1)
                if (val & 1) s[b] = '1';
            return s;
        };
        f_ver << to_bits(x, n) << "  "
              << to_bits(f_vec, tt.n_outputs) << "  "
              << to_bits(z, best.m) << "  "
              << to_bits(g_vec, tt.n_outputs) << "\n";

        if (f_vec != g_vec) {
            n_mismatch++;
            if (first_mismatch_oi < 0) {
                uint64_t diff = f_vec ^ g_vec;
                first_mismatch_oi = __builtin_ctzll(diff);
            }
        }
    }

    f_ver << "\n";
    if (n_mismatch == 0) {
        f_ver << "✅ All outputs verified! (5000 tests, 0 mismatches)\n";
    } else {
        f_ver << "❌ " << n_mismatch << " mismatches out of " << n_tests << " tests\n";
        f_ver << "   First mismatch at output: "
              << circ.outputs[output_indices[first_mismatch_oi]] << "\n";
    }
    std::cout << "  Saved: " << fname << " (" << n_mismatch << "/" << n_tests << " mismatches)\n";
}
