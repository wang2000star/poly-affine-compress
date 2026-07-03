/**
 * compute_raw_anf — Bit-parallel exhaustive ANF computation for n ≤ 32.
 *
 * Reads a .txt circuit, computes truth table via bit-parallel simulation,
 * performs in-place Möbius transform on bit-packed arrays, and outputs
 * T(g) per output plus union T.
 *
 * Usage: compute_raw_anf <circuit.txt> [output_name ...]
 *   If no output names given, all non-constant outputs are processed.
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
#include <optional>
#include <cmath>
#include <cassert>

// ============================================================
//  Circuit representation
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
    std::unordered_map<std::string, int> name_to_idx; // >0 → stmt+1, <0 → input
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

    // INORDER
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

    // OUTORDER
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

    // Input indices (negative, -1-based)
    for (auto& inp : circ.inputs)
        circ.name_to_idx[inp] = -(int)circ.name_to_idx.size() - 1;

    // Body: after last header semicolon
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
            s.op = Op::XOR;  // '+' = XOR in this codebase (xor_semantics)
            auto plus = rhs.find('+');
            s.arg1 = trim(rhs.substr(0, plus));
            s.arg2 = trim(rhs.substr(plus + 1));
        } else {
            s.op = Op::INPUT;  // alias
            s.arg1 = rhs;
        }
        circ.stmts.push_back(s);
        circ.name_to_idx[lhs] = (int)circ.stmts.size();
    }
    return circ;
}

// ============================================================
//  Input word generation
// ============================================================

// Precomputed 64-bit patterns for input bits i=0..5
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
    int64_t rem = cycle - offset;     // inputs until next toggle
    int bit_val = (base >> i) & 1;
    if (rem >= 64) return bit_val ? ~0ULL : 0;
    uint64_t ones = (1ULL << rem) - 1;
    return bit_val ? ones : ~ones;
}

// ============================================================
//  Bit-parallel circuit evaluation
// ============================================================

static std::vector<uint64_t> eval_batch(
    const Circuit& circ,
    const std::vector<uint64_t>& input_words,
    std::vector<uint64_t>& eval_buf)  // reusable buffer, size = stmts.size()+1
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

    std::vector<uint64_t> result(circ.outputs.size(), 0);
    for (int o = 0; o < (int)circ.outputs.size(); o++) {
        auto it = circ.name_to_idx.find(circ.outputs[o]);
        if (it != circ.name_to_idx.end()) {
            int idx = it->second;
            if (idx < 0) {
                int inp_idx = -idx - 1;
                result[o] = (inp_idx < (int)input_words.size()) ? input_words[inp_idx] : 0;
            } else {
                result[o] = eval_buf[idx];
            }
        }
    }
    return result;
}

// ============================================================
//  Bit-packed Möbius transform (in-place)
// ============================================================

/** In-place Möbius transform on bit-packed array (N = 2^n bits). */
static void moebius_packed(uint64_t* data, int n) {
    int64_t words = int64_t(1) << (n - 6);  // N / 64

    for (int i = 0; i < n; i++) {
        int64_t step = int64_t(1) << i;

        if (step >= 64) {
            // Word-aligned: step/64 words per block
            int ws = (int)(step / 64);
            for (int64_t w = ws; w < words; w++) {
                if ((w / ws) & 1)
                    data[w] ^= data[w - ws];
            }
        } else {
            // Within-word: step < 64, process each word independently
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

// ============================================================
//  Main computation
// ============================================================

struct RawANFResult {
    std::vector<int64_t> output_T;   // T per output
    std::vector<int> output_max_deg; // max monomial degree per output
    int64_t sum_T;
    int64_t union_T;                 // -1 = not computed
    int overall_max_deg;             // max degree across all outputs
};

static std::optional<RawANFResult> compute_raw_anf(
    const Circuit& circ,
    const std::vector<int>& output_indices,
    int n_threads,
    const std::string& out_prefix = "")
{
    int n = circ.n_inputs;
    if (n > 32) {
        std::cerr << "  SKIP: n=" << n << " > 32\n";
        return std::nullopt;
    }

    int64_t N = int64_t(1) << n;
    int64_t n_batches = N / 64;
    int n_out = (int)output_indices.size();
    int words_per_output = (int)(N / 64);

    // Allocate truth tables: each output = bit-packed array
    std::vector<std::vector<uint64_t>> tts(n_out,
        std::vector<uint64_t>(words_per_output, 0));

    std::cout << "  Truth table: 2^" << n << " = " << N
              << " inputs, " << n_batches << " batches"
              << ", " << n_out << " output(s), " << n_threads << " thread(s)\n";

    // Parallel batch processing
    auto worker = [&](int64_t b_start, int64_t b_end) {
        std::vector<uint64_t> eval_buf(circ.stmts.size() + 1, 0);
        std::vector<uint64_t> in_words(circ.n_inputs, 0);
        for (int64_t b = b_start; b < b_end; b++) {
            int64_t base = b * 64;
            for (int i = 0; i < circ.n_inputs; i++)
                in_words[i] = input_word_for_bit(i, base);
            auto out = eval_batch(circ, in_words, eval_buf);
            for (int oi = 0; oi < n_out; oi++)
                tts[oi][b] = out[output_indices[oi]];
        }
    };

    if (n_threads <= 1) {
        worker(0, n_batches);
    } else {
        int64_t chunk = std::max(int64_t(1), (n_batches + n_threads - 1) / n_threads);
        std::vector<std::thread> threads;
        for (int t = 0; t < n_threads; t++) {
            int64_t start = t * chunk;
            int64_t end = std::min(start + chunk, n_batches);
            if (start < end)
                threads.emplace_back(worker, start, end);
        }
        for (auto& th : threads) th.join();
    }

    // Möbius transform + counting for each output
    RawANFResult result;
    result.output_T.resize(n_out, 0);
    result.output_max_deg.resize(n_out, 0);
    result.sum_T = 0;
    result.overall_max_deg = 0;

    for (int oi = 0; oi < n_out; oi++) {
        moebius_packed(tts[oi].data(), n);
        int64_t T = 0;
        int max_deg = 0;
        int64_t words = int64_t(1) << (n - 6);
        // Count T in forward pass (fast)
        for (int64_t w = 0; w < words; w++)
            T += __builtin_popcountll(tts[oi][w]);
        // Compute max degree: scan from high to low for early exit.
        // deg(packed_idx) = popcount(packed_idx).
        // Each word w covers indices [w*64, w*64+63]. For bit b in [0,63],
        // deg = popcount(w) + popcount(b)  (since w*64+b = (w<<6)|b).
        // Max possible degree in word w is popcount(w) + 6.
        for (int64_t w = words - 1; w >= 0 && max_deg < n; w--) {
            uint64_t val = tts[oi][w];
            if (val == 0) continue;
            // Skip if no monomial in this word can beat current max_deg
            int base_pop = __builtin_popcountll(w);
            if (base_pop + 6 <= max_deg) break;
            while (val) {
                int bit = __builtin_ctzll(val);
                val &= val - 1;
                int d = base_pop + __builtin_popcountll(bit);
                if (d > max_deg) max_deg = d;
                if (max_deg == n) break;
            }
        }
        result.output_T[oi] = T;
        result.output_max_deg[oi] = max_deg;
        result.sum_T += T;
        if (max_deg > result.overall_max_deg) result.overall_max_deg = max_deg;
        std::cout << "    " << circ.outputs[output_indices[oi]]
                  << ": T=" << T << ", m=" << max_deg << "\n";
    }
    std::cout << "  Sum T = " << result.sum_T << "\n";

    // Save raw ANF files if out_prefix specified
    if (!out_prefix.empty() && n <= 16) {
        // expr.poly
        std::string expr_path = out_prefix + "_expr.poly";
        std::ofstream fe(expr_path);
        if (fe) {
            fe << "# Raw ANF for circuit (n=" << n << ", k=" << n_out << " outputs)\n";
            if (n <= 16) {
                fe << "# Variables:";
                for (int i = n - 1; i >= 0; i--)
                    fe << " " << circ.inputs[i];
                fe << "\n";
            }
            for (int oi = 0; oi < n_out; oi++) {
                fe << circ.outputs[output_indices[oi]] << " = ";
                int count = 0;
                int64_t words = int64_t(1) << (n - 6);
                for (int64_t w = 0; w < words; w++) {
                    uint64_t word = tts[oi][w];
                    while (word) {
                        int bit = __builtin_ctzll(word);
                        word &= word - 1;
                        int64_t pos = (w << 6) | bit;
                        if (count > 0) fe << " + ";
                        if (pos == 0) { fe << "1"; count++; continue; }
                        bool first = true;
                        for (int j = 0; j < n; j++) {
                            if ((pos >> j) & 1) {
                                if (!first) fe << " * ";
                                fe << circ.inputs[j];
                                first = false;
                            }
                        }
                        count++;
                    }
                }
                if (count == 0) fe << "0";
                fe << "\n";
            }
        }
        std::cout << "  Saved: " << expr_path << "\n";

        // T.poly
        std::string T_path = out_prefix + "_T.poly";
        std::ofstream fT(T_path);
        if (fT) {
            fT << "# Raw ANF term counts for circuit (n=" << n << ", k=" << n_out << " outputs)\n";
            fT << "# sum T = " << result.sum_T << "\n";
            for (int oi = 0; oi < n_out; oi++)
                fT << circ.outputs[output_indices[oi]] << ": T=" << result.output_T[oi] << "\n";
        }
        std::cout << "  Saved: " << T_path << "\n";
    }

    // Union T: OR all arrays then popcount
    // For n ≤ 24: exact (needs ≤ 256 MB bitmap)
    // For n = 25..32: exact (OR in-place on one array, needs 512 MB peak)
    if (n_out == 1) {
        result.union_T = result.sum_T;
    } else {
        // OR all but first into first, then popcount
        for (int oi = 1; oi < n_out; oi++)
            for (int w = 0; w < words_per_output; w++)
                tts[0][w] |= tts[oi][w];
        result.union_T = 0;
        for (auto w : tts[0])
            result.union_T += __builtin_popcountll(w);
        std::cout << "  Union T = " << result.union_T << "\n";
    }

    return result;
}

// ============================================================
//  Main
// ============================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [output_name ...] [--out-prefix PREFIX]\n";
        std::cerr << "  --out-prefix PREFIX  save raw ANF to PREFIX_expr.poly and PREFIX_T.poly\n";
        return 1;
    }

    std::string path = argv[1];
    auto circ_opt = parse_circuit(path);
    if (!circ_opt) { std::cerr << "Failed to parse\n"; return 1; }
    const auto& circ = *circ_opt;

    // Determine output indices and output prefix
    std::vector<int> output_indices;
    std::string out_prefix;
    for (int a = 2; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "--out-prefix" && a + 1 < argc) {
            out_prefix = argv[++a];
        } else {
            for (int o = 0; o < (int)circ.outputs.size(); o++) {
                if (circ.outputs[o] == arg) { output_indices.push_back(o); break; }
            }
        }
    }
    if (output_indices.empty()) {
        for (int o = 0; o < (int)circ.outputs.size(); o++)
            output_indices.push_back(o);
    }

    if (output_indices.empty()) {
        std::cerr << "No outputs to process\n";
        return 1;
    }

    int n_threads = std::thread::hardware_concurrency();
    if (n_threads < 1) n_threads = 1;

    std::cout << "--- " << path << " ---\n";
    std::cout << "  Inputs: " << circ.n_inputs << ", Outputs: "
              << output_indices.size() << "\n";

    auto t0 = std::chrono::steady_clock::now();
    auto result = compute_raw_anf(circ, output_indices, n_threads, out_prefix);
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    if (result)
        std::cout << "  Time: " << elapsed << " s\n";

    return 0;
}
