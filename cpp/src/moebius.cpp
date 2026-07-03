#include "moebius.h"
#include <algorithm>
#include <thread>
#include <vector>

// In-word Möbius transform (bits 0..5) on a single 64-bit word.
// For each i: XOR bit at position (p - 2^i) into position p, for all p where bit i of p is 1.
// Parallel bitwise implementation:
//   i=0: v ^= ((v & 0x5555555555555555) << 1)   -- XOR adjacent bits upward
//   i=1: v ^= ((v & 0x3333333333333333) << 2)
//   i=2: v ^= ((v & 0x0F0F0F0F0F0F0F0F) << 4)
//   i=3: v ^= ((v & 0x00FF00FF00FF00FF) << 8)
//   i=4: v ^= ((v & 0x0000FFFF0000FFFF) << 16)
//   i=5: v ^= ((v & 0x00000000FFFFFFFF) << 32)
static uint64_t moebius_word(uint64_t v) {
    v ^= ((v & 0x5555555555555555ULL) << 1);
    v ^= ((v & 0x3333333333333333ULL) << 2);
    v ^= ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
    v ^= ((v & 0x00FF00FF00FF00FFULL) << 8);
    v ^= ((v & 0x0000FFFF0000FFFFULL) << 16);
    v ^= ((v & 0x00000000FFFFFFFFULL) << 32);
    return v;
}

void moebius_packed(uint64_t* data, int n) {
    int64_t n_words = int64_t(1) << (n - 6);

    // In-word pass: bits 0..5 (independent per word)
    for (int64_t w = 0; w < n_words; w++)
        data[w] = moebius_word(data[w]);

    // Cross-word pass: bits 6..n-1
    for (int i = 6; i < n; i++) {
        int64_t step_words = int64_t(1) << (i - 6);  // word-level distance
        for (int64_t j = 0; j < n_words; j += step_words * 2) {
            for (int64_t k = j; k < j + step_words; k++)
                data[k + step_words] ^= data[k];
        }
    }
}

void moebius_packed_mt(uint64_t* data, int n, int n_threads) {
    if (n <= 16 || n_threads <= 1) {
        moebius_packed(data, n);
        return;
    }
    int64_t n_words = int64_t(1) << (n - 6);
    int use_threads = std::min(n_threads, 8);
    int64_t chunk = std::max(int64_t(1), (n_words + use_threads - 1) / use_threads);

    // In-word pass: parallel by word
    auto word_worker = [&](int64_t start, int64_t end) {
        for (int64_t w = start; w < end; w++)
            data[w] = moebius_word(data[w]);
    };
    {
        std::vector<std::thread> threads;
        for (int t = 0; t < use_threads; t++) {
            int64_t start = t * chunk;
            int64_t end = std::min(start + chunk, n_words);
            if (start < end)
                threads.emplace_back(word_worker, start, end);
        }
        for (auto& th : threads) th.join();
    }

    // Cross-word pass: bits 6..n-1 (sequential, cross-word dependencies)
    for (int i = 6; i < n; i++) {
        int64_t step_words = int64_t(1) << (i - 6);
        for (int64_t j = 0; j < n_words; j += step_words * 2) {
            for (int64_t k = j; k < j + step_words; k++)
                data[k + step_words] ^= data[k];
        }
    }
}
