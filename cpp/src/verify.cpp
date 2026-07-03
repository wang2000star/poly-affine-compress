#include "verify.h"

bool verify_transform(const TruthTable& tt,
                      const std::vector<uint64_t>& g_tt,
                      int output_idx,
                      const uint32_t* M_rows,
                      uint32_t b, int m, int n,
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

uint64_t pack_output_vector(const TruthTable& tt, uint32_t x) {
    uint64_t result = 0;
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        uint64_t bit = (tt.tt[oi][x >> 6] >> (x & 63)) & 1;
        result |= (bit << oi);
    }
    return result;
}
