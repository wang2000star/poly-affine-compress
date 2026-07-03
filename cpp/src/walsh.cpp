#include "walsh.h"
#include <algorithm>
#include <thread>
#include <vector>

// Shared bit patterns (same as truth_table.cpp)
static const uint64_t LOW_PATTERNS[6] = {
    0xAAAAAAAAAAAAAAAAULL,
    0xCCCCCCCCCCCCCCCCULL,
    0xF0F0F0F0F0F0F0F0ULL,
    0xFF00FF00FF00FF00ULL,
    0xFFFF0000FFFF0000ULL,
    0xFFFFFFFF00000000ULL,
};

std::vector<WalshInfo> compute_walsh_correlations(const TruthTable& tt, int n_threads) {
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

    auto worker = [&](int oi, int64_t w_start, int64_t w_end) {
        const uint64_t* data = tt.tt[oi].data();
        int64_t cnt_f1_x0[32] = {0};
        int64_t cnt_f1_x1[32] = {0};

        for (int64_t w = w_start; w < w_end; w++) {
            uint64_t word = data[w];
            if (word == 0) continue;
            int64_t base = w * 64;
            int pop = __builtin_popcountll(word);

            for (int i = 0; i < 6 && i < n; i++) {
                uint64_t ones_mask = LOW_PATTERNS[i];
                int f1_x1 = __builtin_popcountll(word & ones_mask);
                int f1_x0 = __builtin_popcountll(word & ~ones_mask);
                cnt_f1_x1[i] += f1_x1;
                cnt_f1_x0[i] += f1_x0;
            }

            for (int i = 6; i < n; i++) {
                if ((base >> i) & 1)
                    cnt_f1_x1[i] += pop;
                else
                    cnt_f1_x0[i] += pop;
            }
        }

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
