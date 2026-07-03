/**
 * optimize_anf — Truth-table-based ANF optimization for n ≤ 32.
 *
 * Core idea: given f(x) computed via truth table, search for affine
 * transform z = Mx + b (over GF(2)) such that g(z) = f(M^{-1}(z⊕b))
 * has fewer ANF terms (T(g) < T(f)).
 *
 * This works with the BIT-PACKED TRUTH TABLE directly, avoiding the need
 * for a sparse ANF representation.  This makes it suitable for dense ANF
 * functions (e.g. hd09, hd11, hd12) where Python's SparseANF explodes.
 */

#include <cstdint>
#include <chrono>
#include <optional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <cmath>
#include <cassert>
#include <cstring>
#include <random>
#include <iomanip>

// ============================================================
//  Circuit representation (same as compute_raw_anf)
// ============================================================

enum class Op : uint8_t { INPUT, CONST0, CONST1, NOT, AND, XOR };

struct Stmt {
    std::string name;
    Op op;
    std::string arg1, arg2;
};

struct Circuit {
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    std::vector<Stmt> stmts;
    std::unordered_map<std::string, int> name_to_idx;
    int n_inputs;
};

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}

static std::optional<Circuit> parse_circuit(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot open " << path << "\n"; return std::nullopt; }
    std::stringstream buf;
    buf << f.rdbuf();
    std::string text = buf.str();
    for (auto& c : text) if (c == '\n' || c == '\r') c = ' ';
    while (text.find("  ") != std::string::npos)
        text.replace(text.find("  "), 2, " ");

    Circuit circ;
    auto p = text.find("INORDER");
    if (p == std::string::npos) return std::nullopt;
    auto eq = text.find('=', p);
    auto semi = text.find(';', p);
    if (eq == std::string::npos || semi == std::string::npos) return std::nullopt;
    {
        std::string list = trim(text.substr(eq + 1, semi - eq - 1));
        std::stringstream ss(list); std::string name;
        while (ss >> name) circ.inputs.push_back(name);
    }
    circ.n_inputs = (int)circ.inputs.size();

    p = text.find("OUTORDER");
    if (p == std::string::npos) return std::nullopt;
    eq = text.find('=', p);
    semi = text.find(';', p);
    if (eq == std::string::npos || semi == std::string::npos) return std::nullopt;
    {
        std::string list = trim(text.substr(eq + 1, semi - eq - 1));
        std::stringstream ss(list); std::string name;
        while (ss >> name) circ.outputs.push_back(name);
    }

    for (auto& inp : circ.inputs)
        circ.name_to_idx[inp] = -(int)circ.name_to_idx.size() - 1;

    auto body_start = text.find(';', std::max(
        text.find("INORDER"), text.find("OUTORDER"))) + 1;
    std::string body = text.substr(body_start);
    std::stringstream ss(body);
    std::string token;
    while (std::getline(ss, token, ';')) {
        token = trim(token);
        if (token.empty()) continue;
        auto m = token.find('=');
        if (m == std::string::npos) continue;
        std::string lhs = trim(token.substr(0, m));
        std::string rhs = trim(token.substr(m + 1));
        if (rhs.empty()) continue;

        Stmt s; s.name = lhs;
        if (rhs == "0" || rhs == "false") { s.op = Op::CONST0; }
        else if (rhs == "1" || rhs == "true") { s.op = Op::CONST1; }
        else if (rhs[0] == '!') { s.op = Op::NOT; s.arg1 = trim(rhs.substr(1)); }
        else if (rhs.find('*') != std::string::npos) {
            s.op = Op::AND;
            auto star = rhs.find('*');
            s.arg1 = trim(rhs.substr(0, star));
            s.arg2 = trim(rhs.substr(star + 1));
        } else if (rhs.find('+') != std::string::npos) {
            s.op = Op::XOR;
            auto plus = rhs.find('+');
            s.arg1 = trim(rhs.substr(0, plus));
            s.arg2 = trim(rhs.substr(plus + 1));
        } else {
            s.op = Op::INPUT;
            s.arg1 = rhs;
        }
        circ.stmts.push_back(s);
        circ.name_to_idx[lhs] = (int)circ.stmts.size();
    }
    return circ;
}

// ============================================================
//  Input word generation (same as compute_raw_anf)
// ============================================================

static const uint64_t LOW_PATTERNS[6] = {
    0xAAAAAAAAAAAAAAAAULL,  // bit 0
    0xCCCCCCCCCCCCCCCCULL,  // bit 1
    0xF0F0F0F0F0F0F0F0ULL,  // bit 2
    0xFF00FF00FF00FF00ULL,  // bit 3
    0xFFFF0000FFFF0000ULL,  // bit 4
    0xFFFFFFFF00000000ULL,  // bit 5
};

static uint64_t input_word_for_bit(int i, int64_t base) {
    if (i < 6) return LOW_PATTERNS[i];
    int64_t cycle = int64_t(1) << i;
    int64_t offset = base % cycle;
    int64_t rem = cycle - offset;
    int bit_val = (base >> i) & 1;
    if (rem >= 64) return bit_val ? ~0ULL : 0;
    uint64_t ones = (1ULL << rem) - 1;
    return bit_val ? ones : ~ones;
}

// ============================================================
//  Bit-parallel circuit evaluation
// ============================================================

static void eval_batch(
    const Circuit& circ,
    const std::vector<uint64_t>& input_words,
    std::vector<uint64_t>& eval_buf,
    uint64_t* out_results, int n_outputs)
{
    auto get_val = [&](const std::string& name) -> uint64_t {
        auto it = circ.name_to_idx.find(name);
        if (it == circ.name_to_idx.end()) return 0;
        int idx = it->second;
        if (idx < 0) {
            int inp_idx = -idx - 1;
            return (inp_idx < (int)input_words.size()) ? input_words[inp_idx] : 0;
        }
        return eval_buf[idx];
    };

    for (int i = 0; i < (int)circ.stmts.size(); i++) {
        const auto& st = circ.stmts[i];
        uint64_t r;
        switch (st.op) {
            case Op::CONST0: r = 0; break;
            case Op::CONST1: r = ~0ULL; break;
            case Op::INPUT:  r = get_val(st.arg1); break;
            case Op::NOT:    r = ~get_val(st.arg1); break;
            case Op::AND:    r = get_val(st.arg1) & get_val(st.arg2); break;
            case Op::XOR:    r = get_val(st.arg1) ^ get_val(st.arg2); break;
        }
        eval_buf[i + 1] = r;
    }

    for (int o = 0; o < n_outputs; o++) {
        auto it = circ.name_to_idx.find(circ.outputs[o]);
        if (it != circ.name_to_idx.end()) {
            int idx = it->second;
            if (idx < 0) {
                int inp_idx = -idx - 1;
                out_results[o] = (inp_idx < (int)input_words.size()) ? input_words[inp_idx] : 0;
            } else {
                out_results[o] = eval_buf[idx];
            }
        }
    }
}

// ============================================================
//  Truth table computation (bit-packed, multi-threaded)
// ============================================================

struct TruthTable {
    int n;                          // number of input variables
    int n_outputs;                  // number of outputs
    int64_t n_words;                // 2^n / 64
    std::vector<std::vector<uint64_t>> tt;  // tt[output][word]
    std::vector<int> output_indices; // which original output indices
};

static TruthTable compute_truth_table(
    const Circuit& circ,
    const std::vector<int>& output_indices,
    int n_threads)
{
    int n = circ.n_inputs;
    int64_t N = int64_t(1) << n;
    int64_t n_words = N / 64;
    int n_out = (int)output_indices.size();

    TruthTable result;
    result.n = n;
    result.n_outputs = n_out;
    result.n_words = n_words;
    result.output_indices = output_indices;
    result.tt.resize(n_out, std::vector<uint64_t>(n_words, 0));

    auto worker = [&](int64_t b_start, int64_t b_end) {
        std::vector<uint64_t> eval_buf(circ.stmts.size() + 1, 0);
        std::vector<uint64_t> in_words(circ.n_inputs, 0);
        std::vector<uint64_t> out_vec(n_out, 0);
        for (int64_t b = b_start; b < b_end; b++) {
            int64_t base = b * 64;
            for (int i = 0; i < circ.n_inputs; i++)
                in_words[i] = input_word_for_bit(i, base);
            eval_batch(circ, in_words, eval_buf, out_vec.data(), n_out);
            for (int oi = 0; oi < n_out; oi++)
                result.tt[oi][b] = out_vec[oi];
        }
    };

    if (n_threads <= 1) {
        worker(0, n_words);
    } else {
        int64_t chunk = std::max(int64_t(1), (n_words + n_threads - 1) / n_threads);
        std::vector<std::thread> threads;
        for (int t = 0; t < n_threads; t++) {
            int64_t start = t * chunk;
            int64_t end = std::min(start + chunk, n_words);
            if (start < end)
                threads.emplace_back(worker, start, end);
        }
        for (auto& th : threads) th.join();
    }

    return result;
}

// ============================================================
//  Möbius transform (same as compute_raw_anf)
// ============================================================

static void moebius_packed(uint64_t* data, int n) {
    if (n <= 0) return;
    int64_t words = (n < 6) ? 1 : (int64_t(1) << (n - 6));
    for (int i = 0; i < n; i++) {
        int64_t step = int64_t(1) << i;
        if (step >= 64) {
            int ws = (int)(step / 64);
            for (int64_t w = ws; w < words; w++) {
                if ((w / ws) & 1)
                    data[w] ^= data[w - ws];
            }
        } else {
            int b = (int)step;
            uint64_t mask_lo = (1ULL << b) - 1;
            for (int64_t w = 0; w < words; w++) {
                uint64_t val = data[w];
                uint64_t res = 0;
                for (int g = 0; g < 64; g += 2 * b) {
                    uint64_t lo = (val >> g) & mask_lo;
                    uint64_t hi = (val >> (g + b)) & mask_lo;
                    hi ^= lo;
                    res |= lo << g;
                    res |= hi << (g + b);
                }
                data[w] = res;
            }
        }
    }
}

// Parallel Möbius: each word/group is independent, partition among threads.
// Bit-level passes (i<6): partition words.  Word-level passes (i≥6): partition groups.
// For n=32 with n_threads=8, this reduces 21s → ~3s per output.
static void moebius_packed_mt(uint64_t* data, int n, int n_threads) {
    if (n <= 0) return;
    int64_t words = (n < 6) ? 1 : (int64_t(1) << (n - 6));
    int use_threads = std::min(n_threads, 8);  // cap at 8 (memory bandwidth bound)

    for (int i = 0; i < n; i++) {
        int64_t step = int64_t(1) << i;
        if (step >= 64) {
            int ws = (int)(step / 64);
            int64_t group_words = int64_t(ws) * 2;
            int64_t n_groups = words / group_words;
            int n_workers = std::min(use_threads, (int)n_groups);

            auto worker = [&](int tid) {
                int64_t g_per = (n_groups + n_workers - 1) / n_workers;
                int64_t g0 = tid * g_per;
                int64_t g1 = std::min(g0 + g_per, n_groups);
                for (int64_t g = g0; g < g1; g++) {
                    int64_t w0 = g * group_words + ws;
                    int64_t w1 = g * group_words + group_words;
                    for (int64_t w = w0; w < w1; w++)
                        data[w] ^= data[w - ws];
                }
            };

            std::vector<std::thread> threads;
            for (int t = 0; t < n_workers; t++)
                threads.emplace_back(worker, t);
            for (auto& th : threads) th.join();

            // Remaining partial group (single-threaded)
            int64_t rem_start = n_groups * group_words + ws;
            for (int64_t w = rem_start; w < words; w++) {
                if ((w / ws) & 1) data[w] ^= data[w - ws];
            }
        } else {
            // Bit-level: each word is independent
            int n_workers = std::min(use_threads, (int)words);
            int b = (int)step;
            uint64_t mask_lo = (1ULL << b) - 1;

            auto worker = [&](int tid) {
                int64_t w_per = (words + n_workers - 1) / n_workers;
                int64_t w0 = tid * w_per;
                int64_t w1 = std::min(w0 + w_per, words);
                for (int64_t w = w0; w < w1; w++) {
                    uint64_t val = data[w];
                    uint64_t res = 0;
                    for (int gb = 0; gb < 64; gb += 2 * b) {
                        uint64_t lo = (val >> gb) & mask_lo;
                        uint64_t hi = (val >> (gb + b)) & mask_lo;
                        hi ^= lo;
                        res |= lo << gb;
                        res |= hi << (gb + b);
                    }
                    data[w] = res;
                }
            };

            std::vector<std::thread> threads;
            for (int t = 0; t < n_workers; t++)
                threads.emplace_back(worker, t);
            for (auto& th : threads) th.join();
        }
    }
}

// ============================================================
//  Count T from ANF (after Möbius)
// ============================================================

static int64_t count_T(const uint64_t* data, int n) {
    if (n <= 0) return 0;
    int64_t words = (n < 6) ? 1 : (int64_t(1) << (n - 6));
    int64_t T = 0;
    for (int64_t w = 0; w < words; w++)
        T += __builtin_popcountll(data[w]);
    return T;
}

// ============================================================
//  Walsh single-bit correlations from truth table
// ============================================================

struct WalshInfo {
    int output_idx;
    int64_t walsh_mag[32];  // |W_i| for each input bit (n ≤ 32)
    int64_t walsh_raw[32];  // raw W_i (positive = correlated with x_i=1)
};

static std::vector<WalshInfo> compute_walsh_correlations(
    const TruthTable& tt, int n_threads)
{
    int n = tt.n;
    int n_out = tt.n_outputs;
    int64_t n_words = tt.n_words;

    std::vector<WalshInfo> results(n_out);
    for (int oi = 0; oi < n_out; oi++) {
        results[oi].output_idx = oi;
        for (int i = 0; i < 32; i++) {
            results[oi].walsh_mag[i] = 0;
            results[oi].walsh_raw[i] = 0;
        }
    }

    // For each output, iterate truth table and accumulate counts.
    // W_i = count(f=1 ∧ x_i=1) - count(f=1 ∧ x_i=0)
    // We compute count_f1_x1[i] and count_f1_x0[i] in one pass per output.
    auto worker = [&](int oi, int64_t w_start, int64_t w_end) {
        const uint64_t* data = tt.tt[oi].data();
        int64_t cnt_f1_x0[32] = {0};
        int64_t cnt_f1_x1[32] = {0};

        for (int64_t w = w_start; w < w_end; w++) {
            uint64_t word = data[w];
            if (word == 0) continue;
            int64_t base = w * 64;
            int pop = __builtin_popcountll(word);

            // Bits 0-5: x_i toggles within word
            for (int i = 0; i < 6 && i < n; i++) {
                uint64_t ones_mask = LOW_PATTERNS[i];
                int f1_x1 = __builtin_popcountll(word & ones_mask);
                int f1_x0 = __builtin_popcountll(word & ~ones_mask);
                cnt_f1_x1[i] += f1_x1;
                cnt_f1_x0[i] += f1_x0;
            }

            // Bits 6..n-1: x_i is constant across word
            for (int i = 6; i < n; i++) {
                if ((base >> i) & 1) {
                    cnt_f1_x1[i] += pop;  // x_i = 1
                } else {
                    cnt_f1_x0[i] += pop;  // x_i = 0
                }
            }
        }

        // Accumulate into results (thread-safe since each output is independent)
        for (int i = 0; i < 32 && i < n; i++) {
            int64_t raw = cnt_f1_x1[i] - cnt_f1_x0[i];
            results[oi].walsh_raw[i] = raw;
            results[oi].walsh_mag[i] = raw >= 0 ? raw : -raw;
        }
    };

    if (n_threads <= 1) {
        for (int oi = 0; oi < n_out; oi++)
            worker(oi, 0, n_words);
    } else {
        // Process each output in parallel
        for (int oi = 0; oi < n_out; oi++) {
            int64_t chunk = std::max(int64_t(1), (n_words + n_threads - 1) / n_threads);
            std::vector<std::thread> threads;
            for (int t = 0; t < n_threads; t++) {
                int64_t start = t * chunk;
                int64_t end = std::min(start + chunk, n_words);
                if (start < end)
                    threads.emplace_back(worker, oi, start, end);
            }
            for (auto& th : threads) th.join();
        }
    }

    return results;
}

// ============================================================
//  Evaluate a single M,b candidate for ONE output
// ============================================================
//  GF(2) matrix rank (Gaussian elimination on m×n binary matrix)
// ============================================================

static int gf2_rank(const uint32_t* rows, int m, int n) {
    if (m <= 0 || n <= 0) return 0;
    std::vector<uint32_t> tmp(rows, rows + m);
    for (int i = 0; i < m; i++) {
        if (n < 32) tmp[i] &= (1u << n) - 1;
    }
    int rank = 0;
    for (int col = n - 1; col >= 0 && rank < m; col--) {
        int pivot = -1;
        for (int i = rank; i < m; i++) {
            if ((tmp[i] >> col) & 1) { pivot = i; break; }
        }
        if (pivot < 0) continue;
        std::swap(tmp[rank], tmp[pivot]);
        uint32_t pivot_val = tmp[rank];
        for (int i = 0; i < m; i++) {
            if (i != rank && ((tmp[i] >> col) & 1))
                tmp[i] ^= pivot_val;
        }
        rank++;
    }
    return rank;
}

// ============================================================
//  Verify M,b transform: check f(x) = g(Mx+b) for random x
// ============================================================

static bool verify_transform(
    const TruthTable& tt,
    const std::vector<uint64_t>& g_tt,   // truth table of g (BEFORE Möbius)
    int output_idx,
    const uint32_t* M_rows,
    uint32_t b,
    int m,
    int n,
    int n_tests)
{
    uint32_t M_offsets[64];
    for (int d = 0; d < 64; d++) {
        uint32_t z = 0;
        for (int row = 0; row < m; row++) {
            if (__builtin_popcount(M_rows[row] & d) & 1)
                z |= (1u << row);
        }
        M_offsets[d] = z;
    }

    uint32_t mask = (n == 32) ? 0xFFFFFFFF : ((1u << n) - 1);
    const uint64_t* f_tt = tt.tt[output_idx].data();

    for (int t = 0; t < n_tests; t++) {
        uint32_t x = (uint32_t)(((uint64_t)t * 0x9e3779b97f4a7c15ULL) & mask);
        // compute z = Mx + b
        uint32_t z = b;
        for (int row = 0; row < m; row++) {
            if (__builtin_popcount(M_rows[row] & x) & 1)
                z ^= (1u << row);
        }
        uint64_t f_val = (f_tt[x >> 6] >> (x & 63)) & 1;
        uint64_t g_val = (g_tt[z >> 6] >> (z & 63)) & 1;
        if (f_val != g_val) return false;
    }
    return true;
}


// ============================================================
//  Updated evaluate_Mb_single_output with coset consistency check.
//  For m < n (non-invertible M), must verify f is constant on each
//  coset of ker(M).  Old code XOR'd over cosets, giving PARITY not
//  COMMON VALUE — leading to spurious "T=0" results.
// ============================================================

struct MbResult {
    int64_t T;       // ANF term count (INT64_MAX if inconsistent)
    bool consistent; // whether f is constant on each coset of ker(M)
};

//  State per coset: 0=UNSEEN, 1=SAW_ZERO, 2=SAW_ONE, 3=CONFLICT
static uint8_t merge_states(uint8_t a, uint8_t b) {
    if (a == 0) return b;
    if (b == 0) return a;
    if (a == 3 || b == 3 || a != b) return 3;
    return a;  // a == b ∈ {1,2}
}

static MbResult evaluate_Mb_single_output(
    const uint64_t* f_tt,
    int n,
    int m,
    const uint32_t* M_rows,
    uint32_t b,
    int64_t n_words_f,
    int n_threads)
{
    if (m <= 0) return {0, true};

    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));

    uint32_t M_offsets[64];
    for (int d = 0; d < 64; d++) {
        uint32_t z = 0;
        for (int row = 0; row < m; row++) {
            if (__builtin_popcount(M_rows[row] & d) & 1)
                z |= (1u << row);
        }
        M_offsets[d] = z;
    }

    // Worker: for each x with f(x)=1, compute z=Mx+b and update state.
    // Since we iterate over SET bits of f_tt, f(x) is always 1.
    // SAW_ONE → state=2.
    // If we also encounter a 0 → would be CONFLICT.  But since all
    // entries here have f(x)=1, we only see SAW_ONE or UNSEEN.
    // CONFLICT can only arise if some x in the SAME coset has f(x)=0
    // while another has f(x)=1.  We only iterate over f(x)=1 entries,
    // so CONFLICT means we saw two different 1s for the same coset
    // (which is fine, still consistent — all 1s).  A TRUE conflict
    // (mix of 0 and 1) can only be detected by ALSO iterating over
    // f(x)=0 entries (where f_tt has 0 bits).
    //
    // Therefore the current approach is INCOMPLETE: we only see
    // f(x)=1 entries, so we can't detect conflicts with f(x)=0.
    //
    // We need to either:
    //   (a) iterate over ALL f_tt entries (0 and 1), or
    //   (b) use a different consistency check
    //
    // For now, (a): iterate over f(x)=0 entries too.

    // Extended worker that processes BOTH 0 and 1 entries.
    // Since n_words_f can be huge, we make 2 passes per coset:
    // pass 0: set state from any f(x)=1 entry
    // pass 1: verify no f(x)=0 in same coset (if state==2)

    // Actually, more efficient: process all f_tt entries.
    // For words with mixed bits, some bits are 0 and some 1.
    // For each z, we track whether we've seen both a 0 and a 1.

    // Unified compute: for each f_tt word, process both 1-bits (f=1)
    // and 0-bits (f=0), tracking per-coset state (UNSEEN/SAW_ZERO/SAW_ONE/CONFLICT).
    // g_data[z] = 1 iff state[z] == SAW_ONE after processing.
    auto compute = [&](int64_t w_start, int64_t w_end,
                       uint8_t* state, uint64_t* g_data) {
        for (int64_t w = w_start; w < w_end; w++) {
            uint64_t word = f_tt[w];
            // NOTE: can't skip all-zero words — must process zero-bits for consistency
            int64_t base = w * 64;

            uint32_t z_base = b;
            for (int row = 0; row < m; row++) {
                if (__builtin_popcount(M_rows[row] & (uint32_t)base) & 1)
                    z_base ^= (1u << row);
            }

            // Helper: process a word's bits (1s = f=1, mapped to value 'fx_val')
            // After processing fx_val=1 bits, process fx_val=0 bits.
            uint64_t ones = word;
            uint64_t zeros = ~word;
            // Mask: only first 2^n bits are valid; for n<64, mask off the rest.
            // Actually, n_words_f * 64 may have extra bits in the last word.
            // For simplicity, we process all 64 bits — the state vector has
            // enough entries for all z values, and bits beyond the truth table
            // are always 0 (initialized as such).

            // 1-bits (f=1)
            uint64_t bits = ones;
            while (bits) {
                int bpos = __builtin_ctzll(bits);
                bits &= bits - 1;
                uint32_t z = z_base ^ M_offsets[bpos];
                uint8_t& st = state[z];
                if (st == 0)      { st = 2; g_data[z >> 6] |= (1ULL << (z & 63)); }
                else if (st == 1) { st = 3; }
                // st == 2: already SAW_ONE, consistent
                // st == 3: already CONFLICT
            }

            // 0-bits (f=0)
            bits = zeros;
            while (bits) {
                int bpos = __builtin_ctzll(bits);
                bits &= bits - 1;
                uint32_t z = z_base ^ M_offsets[bpos];
                uint8_t& st = state[z];
                if (st == 2)      { st = 3; }  // mix of 1 and 0 → conflict
                else if (st == 0) { st = 1; }
                // st == 1: already SAW_ZERO, consistent
                // st == 3: already CONFLICT
            }
        }
    };

    if (m <= 20) {
        int n_z = 1 << m;
        int64_t n_thr = std::min(n_threads, 64);
        int64_t chunk = std::max(int64_t(1), (n_words_f + n_thr - 1) / n_thr);

        std::vector<std::vector<uint8_t>> per_thread_state(n_thr, std::vector<uint8_t>(n_z, 0));
        std::vector<std::vector<uint64_t>> per_thread_g(n_thr,
            std::vector<uint64_t>(std::max(int64_t(1), n_words_g), 0));

        std::vector<std::thread> threads;
        for (int t = 0; t < n_thr; t++) {
            int64_t start = t * chunk;
            int64_t end = std::min(start + chunk, n_words_f);
            if (start < end)
                threads.emplace_back(compute, start, end,
                    per_thread_state[t].data(), per_thread_g[t].data());
        }
        for (auto& th : threads) th.join();

        // Merge per-thread state + g_tt into thread 0
        bool consistent = true;
        for (int t = 1; t < n_thr && consistent; t++) {
            for (int zi = 0; zi < n_z; zi++) {
                uint8_t s0 = per_thread_state[0][zi];
                uint8_t st = per_thread_state[t][zi];
                if (st == 0) continue;
                if (s0 == 0) {
                    per_thread_state[0][zi] = st;
                    if (st == 2) per_thread_g[0][zi >> 6] |= (1ULL << (zi & 63));
                } else if (s0 != st) {
                    consistent = false;
                    break;
                }
            }
            // OR g_tt for threads that had values
            for (int64_t w = 0; w < std::max(int64_t(1), n_words_g); w++)
                per_thread_g[0][w] |= per_thread_g[t][w];
        }

        // Check for any CONFLICT in thread 0's state
        for (int zi = 0; zi < n_z; zi++) {
            if (per_thread_state[0][zi] == 3) { consistent = false; break; }
        }

        if (!consistent) return {INT64_MAX, false};

        moebius_packed(per_thread_g[0].data(), m);
        return {count_T(per_thread_g[0].data(), m), true};
    } else {
        // m > 20: single-threaded
        if (m > 25) {
            // m > 25: state vector too large, fall back to XOR (may be wrong)
            // This path should rarely be reached for practical m
            int64_t n_words_g_val = (m < 6) ? 1 : (int64_t(1) << (m - 6));
            std::vector<uint64_t> g_tt(std::max(int64_t(1), n_words_g_val), 0);
            for (int64_t w = 0; w < n_words_f; w++) {
                uint64_t word = f_tt[w];
                if (word == 0) continue;
                int64_t base = w * 64;
                uint32_t z_base = b;
                for (int row = 0; row < m; row++) {
                    if (__builtin_popcount(M_rows[row] & (uint32_t)base) & 1)
                        z_base ^= (1u << row);
                }
                while (word) {
                    int bit = __builtin_ctzll(word);
                    word &= word - 1;
                    uint32_t z = z_base ^ M_offsets[bit];
                    g_tt[z >> 6] ^= (1ULL << (z & 63));
                }
            }
            moebius_packed(g_tt.data(), m);
            return {count_T(g_tt.data(), m), true};
        }

        int n_z = 1 << m;
        std::vector<uint8_t> state(n_z, 0);
        int64_t n_words_g_val = (m < 6) ? 1 : (int64_t(1) << (m - 6));
        std::vector<uint64_t> g_tt(std::max(int64_t(1), n_words_g_val), 0);

        compute(0, n_words_f, state.data(), g_tt.data());

        for (int zi = 0; zi < n_z; zi++) {
            if (state[zi] == 3) return {INT64_MAX, false};
        }

        moebius_packed(g_tt.data(), m);
        return {count_T(g_tt.data(), m), true};
    }
}

// ============================================================
//  Evaluate M,b for ALL outputs
// ============================================================

struct MbCandidate {
    int m;
    uint32_t M_rows[32];  // each row = 32 bits for n ≤ 32
    uint32_t b;           // packed: bit i = b_i
    int64_t total_T;       // sum of T(g) across all outputs
    std::vector<int64_t> per_output_T;
    std::vector<std::vector<uint64_t>> g_tt_raw;  // pre-Möbius g_tt per output (empty if not saved)
};

// ============================================================
//  Evaluate M,b for bijective transform (m=n, full rank).
//
//  Strategy: iterate over f_tt set bits (sequential scan, cache-friendly)
//  and map x → z = Mx+b.  Uses per-thread g_tt arrays (no atomic needed
//  since threads write to disjoint z ranges for bijective M).
//  Merged after all threads finish.
// ============================================================

static MbCandidate evaluate_Mb_bijective(
    const TruthTable& tt,
    const uint32_t* M_rows,
    uint32_t b,
    int m,
    int n,
    int n_threads,
    bool save_g_tt = false)
{
    MbCandidate cand;
    cand.m = m;
    cand.b = b;
    for (int i = 0; i < m; i++) cand.M_rows[i] = M_rows[i];
    cand.total_T = 0;
    cand.per_output_T.resize(tt.n_outputs, 0);

    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    int64_t n_words_f = tt.n_words;

    // M_offsets[d] = M * d for d in [0, 64)
    uint32_t M_offsets[64];
    for (int d = 0; d < 64; d++) {
        uint32_t z = 0;
        for (int row = 0; row < m; row++) {
            if (__builtin_popcount(M_rows[row] & d) & 1)
                z |= (1u << row);
        }
        M_offsets[d] = z;
    }

    // Use limited threads for bijective path (each thread needs its own g_tt).
    // All threads share the f_tt scan (read-only) but each writes to a private g_tt.
    int use_threads = std::min(n_threads, 8);
    int64_t chunk = std::max(int64_t(1), (n_words_f + use_threads - 1) / use_threads);

    // Pre-alloc scratch buffers for the entire run (reused across outputs)
    std::vector<uint64_t*> per_thread_g(use_threads, nullptr);
    for (int t = 0; t < use_threads; t++)
        per_thread_g[t] = (uint64_t*)calloc(n_words_g, sizeof(uint64_t));

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        const uint64_t* f_tt = tt.tt[oi].data();

        // Exact constant-zero check (scan all words — cheap for n≤32)
        bool is_constant_zero = true;
        for (int64_t w = 0; w < n_words_f && is_constant_zero; w++) {
            if (f_tt[w] != 0) { is_constant_zero = false; }
        }
        if (is_constant_zero) {
            cand.per_output_T[oi] = 0;
            if (save_g_tt) {
                std::vector<uint64_t> g_zero(n_words_g, 0);
                cand.g_tt_raw.push_back(std::move(g_zero));
            }
            continue;
        }

        // Clear per-thread g_tt buffers: each thread clears its ENTIRE own buffer.
        // (NOT partitioning the address space — each thread owns a full-size buffer.)
        auto clear_fn = [&](int tid) {
            memset(per_thread_g[tid], 0, n_words_g * sizeof(uint64_t));
        };
        {
            std::vector<std::thread> clear_thr;
            for (int t = 0; t < use_threads; t++)
                clear_thr.emplace_back(clear_fn, t);
            for (auto& th : clear_thr) th.join();
        }

        // Worker: process f_tt words in chunks
        auto worker = [&](int tid) {
            int64_t w_start = tid * chunk;
            int64_t w_end = std::min(w_start + chunk, n_words_f);
            uint64_t* g_local = per_thread_g[tid];

            for (int64_t w = w_start; w < w_end; w++) {
                uint64_t word = f_tt[w];
                if (word == 0) continue;
                int64_t base = w * 64;

                uint32_t z_base = b;
                for (int row = 0; row < m; row++) {
                    if (__builtin_popcount(M_rows[row] & (uint32_t)base) & 1)
                        z_base ^= (1u << row);
                }

                while (word) {
                    int bit = __builtin_ctzll(word);
                    word &= word - 1;
                    uint32_t z = z_base ^ M_offsets[bit];
                    g_local[z >> 6] |= (1ULL << (z & 63));
                }
            }
        };

        {
            std::vector<std::thread> threads;
            for (int t = 0; t < use_threads; t++)
                threads.emplace_back(worker, t);
            for (auto& th : threads) th.join();
        }

        // Merge: OR all per-thread g_tt into thread 0's array (parallel)
        if (use_threads > 1) {
            auto merge_fn = [&](int tid) {
                int64_t ms = tid * chunk;
                int64_t me = std::min(ms + chunk, n_words_g);
                uint64_t* dst = per_thread_g[0];
                for (int t = 1; t < use_threads; t++) {
                    uint64_t* src = per_thread_g[t];
                    for (int64_t w = ms; w < me; w++)
                        dst[w] |= src[w];
                }
            };
            std::vector<std::thread> merge_thr;
            for (int t = 0; t < use_threads; t++)
                merge_thr.emplace_back(merge_fn, t);
            for (auto& th : merge_thr) th.join();
        }

        // Save pre-Möbius g_tt if requested (for verification)
        if (save_g_tt) {
            std::vector<uint64_t> g_copy(n_words_g, 0);
            memcpy(g_copy.data(), per_thread_g[0], n_words_g * sizeof(uint64_t));
            cand.g_tt_raw.push_back(std::move(g_copy));
        }

        moebius_packed_mt(per_thread_g[0], m, n_threads);
        int64_t T = count_T(per_thread_g[0], m);
        cand.per_output_T[oi] = T;
        cand.total_T += T;
    }

    // Free scratch buffers
    for (int t = 0; t < use_threads; t++)
        free(per_thread_g[t]);

    return cand;
}

static MbCandidate evaluate_Mb(
    const TruthTable& tt,
    const uint32_t* M_rows,
    uint32_t b,
    int m,
    int n_threads,
    bool save_g_tt = false)
{
    MbCandidate cand;
    cand.m = m;
    cand.b = b;
    for (int i = 0; i < m; i++) cand.M_rows[i] = M_rows[i];
    cand.total_T = 0;
    cand.per_output_T.resize(tt.n_outputs, 0);

    // For m=n: use bijective evaluation (fast XOR path, needed for save_g_tt)
    // When n > 20, always use bijective path.
    // For n ≤ 20, use bijective path when save_g_tt is requested.
    if (m == tt.n && (tt.n > 20 || save_g_tt)) {
        int rank = gf2_rank(M_rows, m, tt.n);
        if (rank < m) {
            cand.total_T = INT64_MAX;  // singular M — skip
            return cand;
        }
        return evaluate_Mb_bijective(tt, M_rows, b, m, tt.n, n_threads, save_g_tt);
    }

    bool any_inconsistent = false;
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        MbResult res = evaluate_Mb_single_output(
            tt.tt[oi].data(), tt.n, m, M_rows, b,
            tt.n_words, n_threads);
        if (!res.consistent) { any_inconsistent = true; break; }
        cand.per_output_T[oi] = res.T;
        cand.total_T += res.T;
    }

    if (any_inconsistent) {
        cand.total_T = INT64_MAX;  // signal invalid to caller
    }
    return cand;
}

// ============================================================
//  Generate candidate M,b
// ============================================================

struct CandidateGenerator {
    int n;                       // input variables
    int max_m_search;            // max m to search
    const std::vector<WalshInfo>& walsh;
    std::mt19937_64 rng;

    CandidateGenerator(int n, int max_m, const std::vector<WalshInfo>& w)
        : n(n), max_m_search(max_m), walsh(w), rng(42) {}

    // Generate single-row candidates from Walsh (each bit with high |W_i|)
    std::vector<std::pair<uint32_t, uint32_t>> gen_walsh_single_rows(int top_k = 40) {
        // Score each input bit by sum of |W_i| across outputs (weighted)
        struct BitScore {
            int bit;
            double score;
        };
        std::vector<BitScore> scores(n);
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (auto& w : walsh)
                s += (double)w.walsh_mag[i] / (1 << (n - 1));  // normalize
            scores[i] = {i, s};
        }
        std::sort(scores.begin(), scores.end(),
            [](auto& a, auto& b) { return a.score > b.score; });

        std::vector<std::pair<uint32_t, uint32_t>> candidates;
        int k = std::min(top_k, n);
        for (int j = 0; j < k; j++) {
            uint32_t row = 1u << scores[j].bit;
            // Try both b=0 and b=1
            candidates.emplace_back(row, 0u);
            candidates.emplace_back(row, 1u);
        }
        return candidates;
    }

    // Generate multi-row candidates by combining top Walsh bits
    std::vector<std::pair<std::vector<uint32_t>, uint32_t>> gen_multi_row(int max_rows = 12) {
        // Score bits as above
        struct BitScore { int bit; double score; };
        std::vector<BitScore> scores(n);
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (auto& w : walsh)
                s += (double)w.walsh_mag[i] / (1 << (n - 1));
            scores[i] = {i, s};
        }
        std::sort(scores.begin(), scores.end(),
            [](auto& a, auto& b) { return a.score > b.score; });

        std::vector<std::pair<std::vector<uint32_t>, uint32_t>> candidates;

        // Build progressive M: add top bits one by one
        for (int m = 1; m <= std::min(max_rows, n); m++) {
            std::vector<uint32_t> rows;
            for (int j = 0; j < m; j++)
                rows.push_back(1u << scores[j].bit);
            candidates.emplace_back(rows, 0u);
        }

        // Try combining each pair of top 8 bits (XOR pairs)
        int top_n = std::min(8, n);
        for (int i = 0; i < top_n; i++) {
            for (int j = i + 1; j < top_n; j++) {
                uint32_t row = (1u << scores[i].bit) | (1u << scores[j].bit);
                candidates.emplace_back(std::vector<uint32_t>{row}, 0u);
                candidates.emplace_back(std::vector<uint32_t>{row}, 1u);
            }
        }

        return candidates;
    }

    // Generate random M,b candidates
    std::vector<std::pair<std::vector<uint32_t>, uint32_t>> gen_random(
        int n_candidates, int max_m)
    {
        std::vector<std::pair<std::vector<uint32_t>, uint32_t>> candidates;
        for (int c = 0; c < n_candidates; c++) {
            int m = (rng() % max_m) + 1;
            std::vector<uint32_t> rows(m, 0);
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    if (rng() & 1)
                        rows[j] |= (1u << i);
                }
            }
            uint32_t b = 0;
            for (int j = 0; j < m; j++) {
                if (rng() & 1)
                    b |= (1u << j);
            }
            candidates.emplace_back(rows, b);
        }
        return candidates;
    }
};

// ============================================================
//  Compute raw ANF T and degree (baseline)
// ============================================================

struct RawANFInfo {
    int64_t sum_T;
    std::vector<int64_t> per_output_T;
    std::vector<int> per_output_deg;
    int overall_max_deg;
};

static RawANFInfo compute_raw_anf_info(TruthTable& tt) {
    RawANFInfo info;
    info.per_output_T.resize(tt.n_outputs, 0);
    info.per_output_deg.resize(tt.n_outputs, 0);
    info.sum_T = 0;
    info.overall_max_deg = 0;

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        moebius_packed(tt.tt[oi].data(), tt.n);
        info.per_output_T[oi] = count_T(tt.tt[oi].data(), tt.n);

        // Max degree
        int n = tt.n;
        int64_t words = (n < 6) ? 1 : (int64_t(1) << (n - 6));
        int max_deg = 0;
        for (int64_t w = words - 1; w >= 0 && max_deg < n; w--) {
            uint64_t val = tt.tt[oi][w];
            if (val == 0) continue;
            int base_pop = __builtin_popcountll(w);
            if (base_pop + 6 <= max_deg) break;
            while (val) {
                int bit = __builtin_ctzll(val);
                val &= val - 1;
                int d = base_pop + __builtin_popcountll(bit);
                if (d > max_deg) { max_deg = d; if (max_deg == n) break; }
            }
        }
        info.per_output_deg[oi] = max_deg;
        info.sum_T += info.per_output_T[oi];
        if (max_deg > info.overall_max_deg) info.overall_max_deg = max_deg;
    }
    // Restore truth table from ANF (Möbius is self-inverse)
    for (int oi = 0; oi < tt.n_outputs; oi++)
        moebius_packed(tt.tt[oi].data(), tt.n);
    return info;
}

// ============================================================
//  Hill climbing: iteratively flip M bits to reduce T
// ============================================================
//  For m=n, tries flipping each bit of M and b, accepting
//  changes that reduce total_T.  Repeats until convergence.
//  This is O(m*(n+1)) evaluations per full pass.
// ============================================================

static MbCandidate hill_climb(
    const TruthTable& tt,
    const MbCandidate& start,
    int n,
    int n_threads)
{
    MbCandidate best = start;
    int m = start.m;

    auto eval = [&](uint32_t M[32], uint32_t b) -> MbCandidate {
        return evaluate_Mb(tt, M, b, m, n_threads);
    };

    // Make a working copy
    MbCandidate cur = best;
    bool improved;
    int pass = 0;
    const int MAX_PASSES = 50;

    do {
        improved = false;
        pass++;

        // Try flipping each bit of M
        for (int r = 0; r < m && !improved; r++) {
            for (int c = 0; c < n && !improved; c++) {
                cur.M_rows[r] ^= (1u << c);
                MbCandidate cand = eval(cur.M_rows, cur.b);
                if (cand.total_T < best.total_T && cand.total_T != INT64_MAX) {
                    best = cand;
                }
                if (cand.total_T < cur.total_T && cand.total_T != INT64_MAX) {
                    cur = cand;
                    improved = true;
                } else {
                    cur.M_rows[r] ^= (1u << c);  // revert
                }
            }
        }

        // Try flipping each bit of b
        for (int r = 0; r < m && !improved; r++) {
            cur.b ^= (1u << r);
            MbCandidate cand = eval(cur.M_rows, cur.b);
            if (cand.total_T < best.total_T && cand.total_T != INT64_MAX) {
                best = cand;
            }
            if (cand.total_T < cur.total_T && cand.total_T != INT64_MAX) {
                cur = cand;
                improved = true;
            } else {
                cur.b ^= (1u << r);  // revert
            }
        }
    } while (improved && pass < MAX_PASSES);

    return best;
}


// ============================================================
//  Save raw ANF to _expr.poly file
// ============================================================
//  Format:  output_name = term1 + term2 + ...
//  Terms are monomials: var1 * var2 * ...
//  '0' means constant false, '1' means constant true.
// ============================================================

static void save_raw_anf(
    const TruthTable& tt,
    const Circuit& circ,
    const std::vector<int>& output_indices,
    const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) {
        std::cerr << "  ERROR: cannot write " << fname << "\n";
        return;
    }

    std::string circ_name = circ.inputs.empty() ? "circuit" : circ.inputs[0].substr(0, circ.inputs[0].find('_'));
    f << "# Raw ANF for circuit (n=" << tt.n << ", k=" << tt.n_outputs << " outputs)\n";
    f << "# Variables:";
    for (auto& inp : circ.inputs) f << " " << inp;
    f << "\n";

    const auto& inputs = circ.inputs;
    int64_t n_words = tt.n_words;
    int total_terms = 0;

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        std::string out_name = circ.outputs[output_indices[oi]];
        f << out_name << " = ";

        const uint64_t* data = tt.tt[oi].data();
        int count = 0;

        // Check for constant 0 (all zeros)
        bool all_zero = true;
        for (int64_t w = 0; w < n_words && all_zero; w++) {
            if (data[w] != 0) all_zero = false;
        }
        if (all_zero) {
            f << "0\n";
            continue;
        }

        for (int64_t w = 0; w < n_words; w++) {
            uint64_t word = data[w];
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                int64_t pos = (w << 6) | bit;

                if (count > 0) f << " + ";

                if (pos == 0) {
                    f << "1";
                    count++;
                    continue;
                }

                bool first = true;
                for (int i = 0; i < tt.n && i < 32; i++) {
                    if ((pos >> i) & 1) {
                        if (!first) f << " * ";
                        f << inputs[i];
                        first = false;
                    }
                }
                count++;
            }
        }
        f << "\n";
        total_terms += count;
    }
    std::cout << "  Saved ANF: " << fname << " (" << total_terms << " terms)\n";
}

// ============================================================
//  Save term counts to _T.poly
// ============================================================

static int64_t count_union_T(const TruthTable& tt) {
    int64_t n_words = tt.n_words;
    std::vector<uint64_t> union_data(n_words, 0);
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        for (int64_t w = 0; w < n_words; w++)
            union_data[w] |= tt.tt[oi][w];
    }
    return count_T(union_data.data(), tt.n);
}

static void save_raw_T(
    const TruthTable& tt,
    const Circuit& circ,
    const std::vector<int>& output_indices,
    const std::string& fname,
    int64_t sum_T)
{
    std::ofstream f(fname);
    if (!f) {
        std::cerr << "  ERROR: cannot write " << fname << "\n";
        return;
    }

    int64_t union_T = count_union_T(tt);

    f << "# Raw ANF term counts for circuit (n=" << tt.n
      << ", k=" << tt.n_outputs << " outputs)\n";
    f << "# sum T = " << sum_T << "\n";
    f << "# union T = " << union_T << "\n";
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        int64_t T = count_T(tt.tt[oi].data(), tt.n);
        f << circ.outputs[output_indices[oi]] << ": T=" << T << "\n";
    }
    std::cout << "  Saved T: " << fname << " (sum=" << sum_T << ", union=" << union_T << ")\n";
}

// ============================================================
//  Save optimized ANF (z-space) to _expr.poly
// ============================================================

static void save_opt_anf(
    const std::vector<std::vector<uint64_t>>& g_tt_raw,  // pre-Möbius per output
    const Circuit& circ,
    const std::vector<int>& output_indices,
    int m,
    const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) {
        std::cerr << "  ERROR: cannot write " << fname << "\n";
        return;
    }

    f << "# opt1 ANF (z-space, m=" << m << ")\n";
    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    int total_terms = 0;

    for (int oi = 0; oi < (int)g_tt_raw.size(); oi++) {
        std::string out_name = circ.outputs[output_indices[oi]];
        f << out_name << " = ";

        // Copy and apply Möbius to get ANF
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
//  Save optimized term counts + union T to _T.poly
// ============================================================

static void save_opt_T(
    const std::vector<std::vector<uint64_t>>& g_tt_raw,
    const Circuit& circ,
    const std::vector<int>& output_indices,
    int m,
    const std::string& fname)
{
    std::ofstream f(fname);
    if (!f) {
        std::cerr << "  ERROR: cannot write " << fname << "\n";
        return;
    }

    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    int64_t sum_T = 0;
    std::vector<int64_t> per_T;

    // Apply Möbius to each output and compute T + union
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
    int64_t union_T = count_T(union_data.data(), m);

    f << "# opt1 T(g) (m=" << m << ")\n";
    f << "# sum T = " << sum_T << "\n";
    f << "# union T = " << union_T << "\n";
    for (int oi = 0; oi < (int)g_tt_raw.size(); oi++)
        f << circ.outputs[output_indices[oi]] << ": T=" << per_T[oi] << "\n";
    std::cout << "  Saved T: " << fname << " (sum=" << sum_T << ", union=" << union_T << ")\n";
}

// ============================================================
//  Pack output vector from truth table for a given x
// ============================================================

static uint64_t pack_output_vector(
    const TruthTable& tt,
    uint32_t x)
{
    uint64_t result = 0;
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        uint64_t bit = (tt.tt[oi][x >> 6] >> (x & 63)) & 1;
        result |= (bit << oi);
    }
    return result;
}


// ============================================================
//  Main search function
// ============================================================

struct SearchParams {
    int max_m = 12;              // maximum number of z variables to try
    int walsh_single_top = 30;   // top K Walsh bits for single-row candidates
    int multi_max_rows = 10;     // max rows for multi-row Walsh candidates
    int n_random = 40;           // number of random M,b candidates
    int n32_random = 0;          // n=32 random m=n candidates (full rank only)
    int n_hill_climb = 10;       // number of top candidates to hill climb from
    int n_threads = 104;
    bool verbose = true;
    std::string anf_out_prefix;  // if non-empty, save raw ANF to PREFIX.poly (n ≤ 16)
    std::string save_results_prefix;  // if non-empty, save best candidate results to PREFIX_*
};

static void run_search(
    Circuit& circ,
    const std::vector<int>& output_indices,
    const SearchParams& params)
{
    int n = circ.n_inputs;

    std::cout << std::unitbuf;
    std::cout << "--- " << "Circuit" << " ---\n";
    std::cout << "  Inputs: " << n << ", Outputs to optimize: " << output_indices.size() << "\n";

    auto t0 = std::chrono::steady_clock::now();

    // Phase 1: Compute truth table
    std::cout << "Phase 1: Computing truth table...\n";
    auto tt = compute_truth_table(circ, output_indices, params.n_threads);
    auto t1 = std::chrono::steady_clock::now();
    double dt = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "  Truth table: 2^" << n << " = " << (int64_t(1) << n)
              << " inputs, " << tt.n_words << " batches, "
              << tt.n_outputs << " output(s), " << params.n_threads << " thread(s)\n";
    std::cout << "  Time: " << dt << " s\n";

    // Phase 2: Compute raw ANF baseline (Möbius transform on truth table)
    // Skipped for n > 20 (Möbius on 2^32 = 16 GB data is too slow single-threaded)
    RawANFInfo raw;
    if (n > 20) {
        std::cout << "\nPhase 2: Skipping raw ANF baseline (n=" << n << " > 20, too expensive)\n";
        raw.sum_T = int64_t(1) << n;  // worst-case estimate
        raw.overall_max_deg = n;
        raw.per_output_T.assign(tt.n_outputs, int64_t(1) << n);
        raw.per_output_deg.assign(tt.n_outputs, n);
    } else {
        std::cout << "\nPhase 2: Raw ANF baseline...\n";
        TruthTable tt_copy_baseline = tt;
        raw = compute_raw_anf_info(tt_copy_baseline);
        std::cout << "  Sum T = " << raw.sum_T << "\n";
        std::cout << "  Max deg = " << raw.overall_max_deg << "\n";
        for (int oi = 0; oi < tt.n_outputs; oi++) {
            std::cout << "    " << circ.outputs[output_indices[oi]]
                      << ": T=" << raw.per_output_T[oi]
                      << ", m=" << raw.per_output_deg[oi] << "\n";
        }
    }

    // Save raw ANF to _expr.poly + _T.poly if requested (n ≤ 16)
    if (!params.anf_out_prefix.empty() && n <= 16) {
        std::string prefix = params.anf_out_prefix;
        // Copy tt, apply Möbius to get ANF data
        TruthTable tt_copy(tt);
        for (int oi = 0; oi < tt_copy.n_outputs; oi++)
            moebius_packed(tt_copy.tt[oi].data(), tt_copy.n);
        save_raw_anf(tt_copy, circ, output_indices, prefix + "_expr.poly");
        save_raw_T(tt_copy, circ, output_indices, prefix + "_T.poly", raw.sum_T);
    }

    // Phase 3: Compute Walsh correlations (for n ≤ 32: used for permutation search)
    std::vector<WalshInfo> walsh;
    if (n <= 32) {
        std::cout << "\nPhase 3: Computing Walsh correlations...\n";
        walsh = compute_walsh_correlations(tt, params.n_threads);
        for (int oi = 0; oi < tt.n_outputs; oi++) {
            std::cout << "  Output " << circ.outputs[output_indices[oi]] << " top-5 bits: ";
            std::vector<std::pair<int, int64_t>> bits;
            for (int i = 0; i < n; i++)
                bits.emplace_back(i, walsh[oi].walsh_mag[i]);
            std::sort(bits.begin(), bits.end(),
                [](auto& a, auto& b) { return a.second > b.second; });
            for (int j = 0; j < std::min(5, (int)bits.size()); j++)
                std::cout << bits[j].first << "(" << bits[j].second << ") ";
            std::cout << "\n";
        }

        auto t_walsh = std::chrono::steady_clock::now();
        std::cout << "  Walsh time: " << std::chrono::duration<double>(t_walsh - t1).count() << " s\n";
    } else {
        std::cout << "\nPhase 3: Skipping Walsh correlations (n=" << n << " > 32)\n";
    }

    // Phase 4: Generate and evaluate candidates
    std::cout << "\nPhase 4: Searching for M,b...\n";

    std::vector<MbCandidate> results;

    // 4a-c: Walsh-guided + random candidates (for n ≤ 20, full search; for n=21-32, permutation search only)
    if (n <= 20) {
        CandidateGenerator gen(n, params.max_m, walsh);

        // 4a: Single-row Walsh candidates
        std::cout << "  4a: Single-row Walsh (" << params.walsh_single_top * 2 << " candidates)\n";
        auto single_rows = gen.gen_walsh_single_rows(params.walsh_single_top);
        for (auto& [row, b] : single_rows) {
            uint32_t M[32] = {0};
            M[0] = row;
            auto cand = evaluate_Mb(tt, M, b, 1, params.n_threads);
            results.push_back(cand);
        }

        // 4b: Multi-row Walsh candidates (linear combinations)
        std::cout << "  4b: Multi-row Walsh (" << params.multi_max_rows
                  << " progressive + pair combinations)\n";
        auto multi_rows = gen.gen_multi_row(params.multi_max_rows);
        for (auto& [rows, b] : multi_rows) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            auto cand = evaluate_Mb(tt, M, b, m, params.n_threads);
            results.push_back(cand);
        }

        // 4c: Random M,b
        std::cout << "  4c: Random M,b (" << params.n_random << " candidates)\n";
        auto random_cands = gen.gen_random(params.n_random, params.max_m);
        for (auto& [rows, b] : random_cands) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            auto cand = evaluate_Mb(tt, M, b, m, params.n_threads);
            results.push_back(cand);
        }

        auto t_search0 = std::chrono::steady_clock::now();
        double walsh_time = std::chrono::duration<double>(t_search0 - t1).count();
        std::cout << "  Search time: " << walsh_time << " s\n";
    } else if (n <= 32 && !walsh.empty()) {
        // n=21-32: permutation search via Walsh (unit-vector rows only, no small-m random)
        CandidateGenerator gen(n, params.max_m, walsh);

        // 4a: Single-row Walsh (unit vectors)
        std::cout << "  4a: Single-row Walsh (" << params.walsh_single_top * 2 << " candidates)\n";
        auto single_rows = gen.gen_walsh_single_rows(params.walsh_single_top);
        for (auto& [row, b] : single_rows) {
            uint32_t M[32] = {0};
            M[0] = row;
            auto cand = evaluate_Mb(tt, M, b, 1, params.n_threads);
            results.push_back(cand);
        }

        // 4b: Multi-row Walsh (progressive permutation matrices)
        std::cout << "  4b: Multi-row Walsh (progressive m=1.." << std::min(20, n)
                  << " + XOR pairs)\n";
        auto multi_rows = gen.gen_multi_row(std::min(20, n));
        for (auto& [rows, b] : multi_rows) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            auto cand = evaluate_Mb(tt, M, b, m, params.n_threads);
            results.push_back(cand);
        }

        std::cout << "  (skipping Phase 4c small-m random for n=" << n << ")\n";
    }

    // 4d: n32 random m=n candidates (for n > 20)
    if (params.n32_random > 0 && n > 20) {
        std::cout << "  4d: n32 random m=n candidates (" << params.n32_random << ")\n";
        std::mt19937_64 rng_n32(42);
        int generated = 0;
        int attempts = 0;
        const int MAX_ATTEMPTS = params.n32_random * 20;

        auto t_n32_0 = std::chrono::steady_clock::now();

        while (generated < params.n32_random && attempts < MAX_ATTEMPTS) {
            attempts++;

            // Generate random m=n M matrix
            uint32_t M[32] = {0};
            for (int r = 0; r < n; r++) {
                for (int i = 0; i < n; i++) {
                    if (rng_n32() & 1) M[r] |= (1u << i);
                }
            }

            // Generate random b vector
            uint32_t b = 0;
            for (int r = 0; r < n; r++) {
                if (rng_n32() & 1) b |= (1u << r);
            }

            // Skip singular M (fast rank check)
            if (gf2_rank(M, n, n) < n) continue;

            // Use evaluate_Mb (which dispatches to evaluate_Mb_bijective for m=n, n>20)
            auto t_eval0 = std::chrono::steady_clock::now();
            auto cand = evaluate_Mb(tt, M, b, n, params.n_threads);
            auto t_eval1 = std::chrono::steady_clock::now();

            results.push_back(cand);
            generated++;

            double eval_time = std::chrono::duration<double>(t_eval1 - t_eval0).count();
            double compression = (double)raw.sum_T / std::max(int64_t(1), cand.total_T);
            std::cout << "    n32 #" << generated << "/" << params.n32_random
                      << " (attempt " << attempts << ")"
                      << ": T=" << cand.total_T
                      << " (compression " << compression << "×)"
                      << " time=" << eval_time << "s\n";
        }

        auto t_n32_1 = std::chrono::steady_clock::now();
        std::cout << "    n32 random time: " << std::chrono::duration<double>(t_n32_1 - t_n32_0).count() << " s\n";
    }

    // 4e: Hill climbing from top candidates (only m=n, skips candidates where m != n)
    if (params.n_hill_climb > 0) {
        // Sort current results to find top candidates
        std::sort(results.begin(), results.end(),
            [](auto& a, auto& b) { return a.total_T < b.total_T; });

        int n_climb = std::min(params.n_hill_climb, (int)results.size());
        std::cout << "  4e: Hill climbing from top " << n_climb << " candidates (+ identity)\n";

        auto t_hc0 = std::chrono::steady_clock::now();

        for (int ci = 0; ci < n_climb; ci++) {
            auto& base = results[ci];
            if (base.total_T >= INT64_MAX) continue;
            if (base.m != n) continue;  // only m=n (invertible) for now

            auto improved = hill_climb(tt, base, n, params.n_threads);
            if (improved.total_T < base.total_T) {
                results.push_back(improved);
                std::cout << "    hill climb #" << ci << ": T=" << base.total_T << " -> T=" << improved.total_T << "\n";
            }
        }

        // Also hill climb from identity (m=n)
        uint32_t ident[32] = {0};
        for (int i = 0; i < n; i++) ident[i] = (1u << i);
        MbCandidate ident_cand = evaluate_Mb(tt, ident, 0, n, params.n_threads);
        auto ident_improved = hill_climb(tt, ident_cand, n, params.n_threads);
        if (ident_improved.total_T < ident_cand.total_T) {
            results.push_back(ident_improved);
            std::cout << "    hill climb from identity: T=" << ident_cand.total_T << " -> T=" << ident_improved.total_T << "\n";
        }

        auto t_hc1 = std::chrono::steady_clock::now();
        std::cout << "    Hill climb time: " << std::chrono::duration<double>(t_hc1 - t_hc0).count() << " s\n";
    }

    // Phase 5: Report
    std::cout << "\n========================================\n";
    std::cout << "  Results (" << results.size() << " candidates evaluated)\n";
    std::cout << "  Raw ANF Sum T = " << raw.sum_T << "\n";
    std::cout << "========================================\n";

    // Sort by total T
    std::sort(results.begin(), results.end(),
        [](auto& a, auto& b) { return a.total_T < b.total_T; });

    int n_report = std::min(20, (int)results.size());
    for (int j = 0; j < n_report; j++) {
        auto& r = results[j];
        double compression = (double)raw.sum_T / std::max(int64_t(1), r.total_T);
        std::cout << "  " << (j + 1) << ". m=" << r.m
                  << " T=" << r.total_T
                  << " (compression " << compression << "×)";

        // Show M matrix
        std::cout << " M=[";
        for (int row = 0; row < r.m; row++) {
            if (row > 0) std::cout << "; ";
            for (int i = 0; i < std::min(n, 8); i++)
                std::cout << ((r.M_rows[row] >> i) & 1);
            std::cout << "..";
        }
        std::cout << "]\n";

        if (r.total_T < raw.sum_T && r.m <= 32) {
            // Print per-output T
            std::cout << "    Per-output T:";
            for (int oi = 0; oi < tt.n_outputs; oi++)
                std::cout << " " << circ.outputs[output_indices[oi]] << "=" << r.per_output_T[oi];
            std::cout << "\n";
        }
    }

    // Phase 6: For the best candidate, compute full ANF of g and verify
    if (!results.empty() && results[0].total_T < INT64_MAX) {
        auto& best = results[0];
        std::cout << "\nBest candidate: m=" << best.m << ", T=" << best.total_T << "\n";

        // Save the transformation for inspection
        std::cout << "// Transform z = Mx + b (GF(2)):\n";
        std::cout << "// m=" << best.m << ", n=" << n << "\n";
        for (int row = 0; row < best.m; row++) {
            std::cout << "// z_" << row << " = ";
            bool first = true;
            for (int i = 0; i < n; i++) {
                if ((best.M_rows[row] >> i) & 1) {
                    if (!first) std::cout << " + ";
                    std::cout << circ.inputs[i];
                    first = false;
                }
            }
            if ((best.b >> row) & 1) {
                if (!first) std::cout << " + ";
                std::cout << "1";
            }
            std::cout << "\n";
        }

        // Verify: re-evaluate best candidate saving pre-Möbius g_tt, then call verify_transform
        bool all_verified = true;
        MbCandidate verified_cand;
        verified_cand.m = best.m;
        verified_cand.b = best.b;
        for (int i = 0; i < best.m; i++) verified_cand.M_rows[i] = best.M_rows[i];
        verified_cand.total_T = best.total_T;
        verified_cand.per_output_T = best.per_output_T;

        if (best.m == n) {
            std::cout << "Verifying best candidate (5000 random tests per output)...\n";
            verified_cand = evaluate_Mb(tt, best.M_rows, best.b, best.m, params.n_threads, true);
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                bool ok = verify_transform(tt, verified_cand.g_tt_raw[oi], oi,
                    best.M_rows, best.b, best.m, n, 5000);
                std::cout << "  " << circ.outputs[output_indices[oi]]
                          << (ok ? " ✅ Verified (5000 tests)" : " ❌ FAILED") << "\n";
                if (!ok) all_verified = false;
            }
        } else {
            std::cout << "Skipping verification (non-bijective m=" << best.m << " < n=" << n << ")\n";
        }
        if (all_verified) {
            std::cout << "✅ All outputs verified!\n";
        }

        // Save results to files if requested
        if (!params.save_results_prefix.empty()) {
            std::string prefix = params.save_results_prefix;

            // (a) Save transform: _trans.poly
            std::string fname_trans = prefix + "_trans.poly";
            std::ofstream f_trans(fname_trans);
            if (f_trans) {
                f_trans << "# opt1 transform z = Mx + b (GF(2))\n";
                f_trans << "# m=" << best.m << ", n=" << n << "\n";
                for (int row = 0; row < best.m; row++) {
                    f_trans << "z_" << row << " = ";
                    bool first = true;
                    for (int i = 0; i < n; i++) {
                        if ((best.M_rows[row] >> i) & 1) {
                            if (!first) f_trans << " + ";
                            f_trans << circ.inputs[i];
                            first = false;
                        }
                    }
                    if ((best.b >> row) & 1) {
                        if (!first) f_trans << " + ";
                        f_trans << "1";
                    } else if (first) {
                        f_trans << "0";  // z_i = 0 (constant false)
                    }
                    f_trans << "\n";
                }
                std::cout << "  Saved: " << fname_trans << "\n";
            }

            // (b) Save optimized ANF + T if we have g_tt_raw
            if (!verified_cand.g_tt_raw.empty()) {
                save_opt_anf(verified_cand.g_tt_raw, circ, output_indices, best.m,
                             prefix + "_expr.poly");
                save_opt_T(verified_cand.g_tt_raw, circ, output_indices, best.m,
                           prefix + "_T.poly");
            }

            // (c) Save verification data: _verify.txt with hex test vectors
            std::string fname_ver = prefix + "_verify.txt";
            std::ofstream f_ver(fname_ver);
            if (f_ver) {
                f_ver << "# Verification for (n=" << n << ", k=" << tt.n_outputs << " outputs)\n";
                f_ver << "# Strategy: opt1\n";
                f_ver << "# Transform: z = Mx + b (m=" << best.m << ")\n";
                f_ver << "# Tests: 5000\n\n";

                int n_tests = 5000;
                uint32_t mask = (n == 32) ? 0xFFFFFFFF : ((1u << n) - 1);
                int hex_bytes = (tt.n_outputs + 7) / 8;
                int n_mismatch = 0;
                int first_mismatch_oi = -1;

                f_ver << "# x  f(x)  g(z)";
                if (hex_bytes > 1) f_ver << "  (hex, " << hex_bytes << " bytes)";
                f_ver << "\n";

                for (int t = 0; t < n_tests; t++) {
                    uint32_t x = (uint32_t)(((uint64_t)t * 0x9e3779b97f4a7c15ULL) & mask);
                    // Compute z = Mx + b
                    uint32_t z = best.b;
                    for (int row = 0; row < best.m; row++) {
                        if (__builtin_popcount(best.M_rows[row] & x) & 1)
                            z ^= (1u << row);
                    }

                    uint64_t f_vec = pack_output_vector(tt, x);
                    uint64_t g_vec = 0;
                    for (int oi = 0; oi < tt.n_outputs; oi++) {
                        const auto& g_tt = verified_cand.g_tt_raw[oi];
                        uint64_t bit = (g_tt[z >> 6] >> (z & 63)) & 1;
                        g_vec |= (bit << oi);
                    }

                    if (hex_bytes <= 1) {
                        f_ver << std::hex << std::setw(2) << std::setfill('0') << (int)x << "  "
                              << std::setw(2) << (int)f_vec << "  "
                              << std::setw(2) << (int)g_vec << std::dec << "\n";
                    } else {
                        f_ver << std::hex << std::setw(2 * hex_bytes) << std::setfill('0') << x << "  "
                              << std::setw(2 * hex_bytes) << f_vec << "  "
                              << std::setw(2 * hex_bytes) << g_vec << std::dec << "\n";
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
                    f_ver << "✅ All outputs verified! (5000 tests, 0 mismatches)\n";
                } else {
                    f_ver << "❌ " << n_mismatch << " mismatches out of " << n_tests << " tests\n";
                    f_ver << "   First mismatch at output: "
                          << circ.outputs[output_indices[first_mismatch_oi]] << "\n";
                }
                std::cout << "  Saved: " << fname_ver
                          << " (" << n_mismatch << "/" << n_tests << " mismatches)\n";
            }
        }
    }

    double total_time = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::cout << "\nTotal time: " << total_time << " s\n";
}

// ============================================================
//  Main
// ============================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [options]\n";
        std::cerr << "  --max-m N      max z variables to search (default 12)\n";
        std::cerr << "  --walsh-k N    top K Walsh bits (default 30)\n";
        std::cerr << "  --random N     random candidates (default 40)\n";
        std::cerr << "  --n32-random N n>20 full-rank random m=n candidates (default 0)\n";
        std::cerr << "  --hill-climb N hill climb from top N candidates (default 10)\n";
        std::cerr << "  --anf-out PREF save raw ANF to PREFIX_expr.poly + PREFIX_T.poly (n≤16 only)\n";
        std::cerr << "  --save-results PREFIX  save best candidate: PREFIX_best_trans.txt,\n";
        std::cerr << "                          PREFIX_best_anf.poly, PREFIX_verify.txt\n";
        return 1;
    }

    std::cout << std::unitbuf;
    std::string path = argv[1];
    auto circ_opt = parse_circuit(path);
    if (!circ_opt) { std::cerr << "Failed to parse\n"; return 1; }
    auto& circ = *circ_opt;

    // All non-constant outputs
    std::vector<int> output_indices;
    for (int o = 0; o < (int)circ.outputs.size(); o++)
        output_indices.push_back(o);

    if (output_indices.empty()) {
        std::cerr << "No outputs to process\n";
        return 1;
    }

    SearchParams params;
    params.n_threads = std::thread::hardware_concurrency();
    if (params.n_threads < 1) params.n_threads = 1;

    // Parse options
    for (int a = 2; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "--max-m" && a + 1 < argc)
            params.max_m = std::stoi(argv[++a]);
        else if (arg == "--walsh-k" && a + 1 < argc)
            params.walsh_single_top = std::stoi(argv[++a]);
        else if (arg == "--random" && a + 1 < argc)
            params.n_random = std::stoi(argv[++a]);
        else if (arg == "--n32-random" && a + 1 < argc)
            params.n32_random = std::stoi(argv[++a]);
        else if (arg == "--hill-climb" && a + 1 < argc)
            params.n_hill_climb = std::stoi(argv[++a]);
        else if (arg == "--anf-out" && a + 1 < argc)
            params.anf_out_prefix = argv[++a];
        else if (arg == "--save-results" && a + 1 < argc)
            params.save_results_prefix = argv[++a];
    }

    run_search(circ, output_indices, params);

    return 0;
}
