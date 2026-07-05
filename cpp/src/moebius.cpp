#include "moebius.h"
#include <algorithm>
#include <thread>
#include <vector>

// In-word Moebius transform for i-th bit (one pass).
// Separated so moebius_packed can apply only the first n passes (n < 6).
static uint64_t moebius_pass(uint64_t v, int i) {
    static const uint64_t MASKS[6] = {
        0x5555555555555555ULL,
        0x3333333333333333ULL,
        0x0F0F0F0F0F0F0F0FULL,
        0x00FF00FF00FF00FFULL,
        0x0000FFFF0000FFFFULL,
        0x00000000FFFFFFFFULL,
    };
    return v ^ ((v & MASKS[i]) << (1 << i));
}

// Full 6-bit in-word Moebius (kept for backward compatibility)
static uint64_t moebius_word(uint64_t v) {
    for (int i = 0; i < 6; i++) v = moebius_pass(v, i);
    return v;
}

void moebius_packed(uint64_t* data, int n) {
    int64_t n_words = (n < 6) ? 1 : (int64_t(1) << (n - 6));

    // In-word pass: only apply the first min(n, 6) passes.
    // When n < 6, higher passes would XOR data across 2^i boundaries
    // into the valid region, corrupting it if unused bits are non-zero.
    int in_word_bits = (n < 6) ? n : 6;
    for (int64_t w = 0; w < n_words; w++) {
        uint64_t v = data[w];
        for (int i = 0; i < in_word_bits; i++)
            v = moebius_pass(v, i);
        data[w] = v;
    }

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
    int in_word_bits = (n < 6) ? n : 6;
    auto word_worker = [&](int64_t start, int64_t end) {
        for (int64_t w = start; w < end; w++) {
            uint64_t v = data[w];
            for (int i = 0; i < in_word_bits; i++)
                v = moebius_pass(v, i);
            data[w] = v;
        }
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

void moebius_upward_packed(uint64_t* data, int n) {
    int64_t n_words = (n < 6) ? 1 : (int64_t(1) << (n - 6));

    // In-word pass: reverse direction (lo ^= hi instead of hi ^= lo)
    int in_word_bits = (n < 6) ? n : 6;
    for (int i = 0; i < in_word_bits; i++) {
        int b = 1 << i;
        uint64_t mask = (1ULL << b) - 1;
        for (int64_t w = 0; w < n_words; w++) {
            uint64_t v = data[w];
            uint64_t res = 0;
            for (int g = 0; g < 64; g += 2 * b) {
                uint64_t lo = (v >> g) & mask;
                uint64_t hi = (v >> (g + b)) & mask;
                lo ^= hi;          // upward: low gets XOR of low and high
                res |= lo << g;
                res |= hi << (g + b);
            }
            data[w] = res;
        }
    }

    // Cross-word pass: upward direction
    for (int i = 6; i < n; i++) {
        int64_t step_words = int64_t(1) << (i - 6);
        for (int64_t w = 0; w + step_words < n_words; w++) {
            if (!((w / step_words) & 1))
                data[w] ^= data[w + step_words];  // upward: low ^= high
        }
    }
}
