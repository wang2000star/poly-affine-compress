#include "io.h"
#include "anf.h"
#include "moebius.h"
#include "verify.h"
#include <fstream>
#include <iostream>

// ============================================================
//  Save raw ANF — {inst}_raw.poly (matrix format)
//  m = n (original x-variables), per-output term rows
// ============================================================

void save_raw_anf(const TruthTable& tt, const Circuit& circ,
                  const std::vector<int>& output_indices,
                  const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    int m = tt.n;
    int k = tt.n_outputs;
    int64_t n_words = tt.n_words;

    // Count terms per output and overall max degree
    std::vector<int> term_counts(k, 0);
    int max_deg = 0;
    for (int oi = 0; oi < k; oi++) {
        std::vector<uint64_t> anf(tt.tt[oi]);
        moebius_packed(anf.data(), m);
        for (int64_t w = 0; w < n_words; w++) {
            uint64_t word = anf[w];
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                term_counts[oi]++;
                int deg = __builtin_popcountll((w << 6) | bit);
                if (deg > max_deg) max_deg = deg;
            }
        }
    }

    // Header
    f << m << "\n" << k << "\n";
    for (int oi = 0; oi < k; oi++) {
        if (oi > 0) f << " ";
        f << term_counts[oi];
    }
    f << "\n";

    // Terms: [e_0 , e_1 , ... , e_{m-1} , coeff]
    for (int oi = 0; oi < k; oi++) {
        std::vector<uint64_t> anf(tt.tt[oi]);
        moebius_packed(anf.data(), m);
        for (int64_t w = 0; w < n_words; w++) {
            uint64_t word = anf[w];
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                int64_t pos = (w << 6) | bit;
                f << "[";
                for (int b = 0; b < m; b++) {
                    if (b > 0) f << " , ";
                    f << ((pos >> b) & 1);
                }
                f << " , 1]\n";
            }
        }
    }

    int total = 0;
    for (int c : term_counts) total += c;
    std::cout << "  Saved raw poly: " << fname << " (m=" << m << ", "
              << k << " outputs, " << total << " terms)\n";
}

// ============================================================
//  Save raw stats — {inst}_raw_stats.txt (5-line numeric)
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

    int m = tt.n;
    int k = tt.n_outputs;
    int64_t n_words = tt.n_words;

    int64_t sum_T = 0;
    int max_deg = 0;
    for (int oi = 0; oi < k; oi++) {
        std::vector<uint64_t> anf(tt.tt[oi]);
        moebius_packed(anf.data(), m);
        sum_T += count_T(anf.data(), m);
        for (int64_t w = 0; w < n_words; w++) {
            uint64_t word = anf[w];
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                int64_t pos = (w << 6) | bit;
                int deg = __builtin_popcountll(pos);
                if (deg > max_deg) max_deg = deg;
            }
        }
    }

    // Raw ANF union = union of original truth tables (pre-Möbius)
    int64_t union_T = count_union_T(tt);

    // 5-line pure numeric
    f << tt.n << "\n";
    f << k << "\n";
    f << sum_T << "\n";
    f << union_T << "\n";
    f << max_deg << "\n";

    std::cout << "  Saved raw stats: " << fname << " (n=" << tt.n
              << ", k=" << k << ", sum_T=" << sum_T
              << ", union_T=" << union_T << ")\n";
}

// ============================================================
//  Save optimized ANF — {inst}_{combo}.poly (matrix format)
//  m = number of z-variables, per-output term rows
// ============================================================

void save_opt_expr(const std::vector<std::vector<uint64_t>>& g_tt_raw,
                   const Circuit& circ,
                   const std::vector<int>& output_indices,
                   int m, const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    int k = (int)g_tt_raw.size();
    if (k == 0) return;
    int64_t n_words = (m < 6) ? 1 : (int64_t(1) << (m - 6));

    // Count terms per output
    std::vector<int> term_counts(k, 0);
    for (int oi = 0; oi < k; oi++) {
        std::vector<uint64_t> anf(g_tt_raw[oi]);
        moebius_packed(anf.data(), m);
        for (int64_t w = 0; w < n_words; w++)
            term_counts[oi] += __builtin_popcountll(anf[w]);
    }

    // Header
    f << m << "\n" << k << "\n";
    for (int oi = 0; oi < k; oi++) {
        if (oi > 0) f << " ";
        f << term_counts[oi];
    }
    f << "\n";

    // Terms
    for (int oi = 0; oi < k; oi++) {
        std::vector<uint64_t> anf(g_tt_raw[oi]);
        moebius_packed(anf.data(), m);
        for (int64_t w = 0; w < n_words; w++) {
            uint64_t word = anf[w];
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                int64_t pos = (w << 6) | bit;
                f << "[";
                for (int b = 0; b < m; b++) {
                    if (b > 0) f << " , ";
                    f << ((pos >> b) & 1);
                }
                f << " , 1]\n";
            }
        }
    }

    int total = 0;
    for (int c : term_counts) total += c;
    std::cout << "  Saved poly: " << fname << " (m=" << m << ", "
              << k << " outputs, " << total << " terms)\n";
}

// ============================================================
//  Save optimized stats — {inst}_{combo}_stats.txt (5-line numeric)
// ============================================================

void save_opt_T(const std::vector<std::vector<uint64_t>>& g_tt_raw,
                const Circuit& circ,
                const std::vector<int>& output_indices,
                int m, const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    int k = (int)g_tt_raw.size();
    if (k == 0) return;
    int64_t n_words = (m < 6) ? 1 : (int64_t(1) << (m - 6));

    int64_t sum_T = 0;
    int max_deg = 0;

    std::vector<uint64_t> union_data(n_words, 0);
    for (int oi = 0; oi < k; oi++) {
        std::vector<uint64_t> anf(g_tt_raw[oi]);
        moebius_packed(anf.data(), m);
        int64_t T = count_T(anf.data(), m);
        sum_T += T;
        for (int64_t w = 0; w < n_words; w++) {
            union_data[w] |= anf[w];
            uint64_t word = anf[w];
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                int64_t pos = (w << 6) | bit;
                int deg = __builtin_popcountll(pos);
                if (deg > max_deg) max_deg = deg;
            }
        }
    }

    int64_t union_T = 0;
    for (int64_t w = 0; w < n_words; w++)
        union_T += __builtin_popcountll(union_data[w]);

    // 5-line pure numeric
    f << circ.n_inputs << "\n";
    f << k << "\n";
    f << sum_T << "\n";
    f << union_T << "\n";
    f << max_deg << "\n";

    std::cout << "  Saved stats: " << fname << " (n=" << circ.n_inputs
              << ", k=" << k << ", sum_T=" << sum_T
              << ", union_T=" << union_T << ")\n";
}

// ============================================================
//  Save affine transform — {inst}_{combo}.affine (matrix format)
//  Header: s t
//          s t n
//  Then: M (s×n), C (t×n), b (s bits), d (t bits)
//  When t=0: no C rows, no d line
// ============================================================

void save_trans(const MbCandidate& best, const Circuit& circ,
                int n, const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) { std::cerr << "  ERROR: cannot write " << fname << "\n"; return; }

    int s = best.m;

    f << s << "\n";
    f << s << " " << (n + 1) << "\n";

    // s × (n+1) matrix: [M | b]
    for (int row = 0; row < s; row++) {
        for (int col = 0; col < n; col++) {
            if (col > 0) f << " ";
            f << ((best.M_rows[row] >> col) & 1);
        }
        f << " " << ((best.b >> row) & 1) << "\n";
    }

    std::cout << "  Saved affine: " << fname << " (s=" << s << ", n=" << n << ")\n";
}

// ============================================================
//  Save verification — {inst}_{combo}_verify.txt
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

    f_ver << "# Verification: f(x) = g(Mx+b) + Cx+d\n";
    f_ver << "# n=" << n << ", s=" << best.m << ", t=0\n";
    f_ver << "# Tests: 5000\n\n";

    int n_tests = 5000;
    uint32_t mask = (n == 32) ? 0xFFFFFFFF : ((1u << n) - 1);
    int n_mismatch = 0;
    int first_mismatch_oi = -1;

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
        f_ver << "All outputs PASS (" << n_tests << " tests)\n";
    } else {
        f_ver << n_mismatch << "/" << n_tests << " mismatches\n";
        f_ver << "First mismatch at output: "
              << circ.outputs[output_indices[first_mismatch_oi]] << "\n";
    }
    std::cout << "  Saved verify: " << fname << " (" << n_mismatch << "/" << n_tests << " mismatches)\n";
}
