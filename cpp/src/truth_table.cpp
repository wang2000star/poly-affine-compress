#include "truth_table.h"
#include <algorithm>
#include <thread>
#include <vector>

// Bit patterns for input word generation
static const uint64_t LOW_PATTERNS[6] = {
    0xAAAAAAAAAAAAAAAAULL,
    0xCCCCCCCCCCCCCCCCULL,
    0xF0F0F0F0F0F0F0F0ULL,
    0xFF00FF00FF00FF00ULL,
    0xFFFF0000FFFF0000ULL,
    0xFFFFFFFF00000000ULL,
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

TruthTable compute_truth_table(
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
