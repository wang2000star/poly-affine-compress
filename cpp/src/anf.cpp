#include "anf.h"
#include "moebius.h"
#include <algorithm>

int64_t count_T(const uint64_t* data, int n) {
    int64_t total = 0;
    int64_t words = (n < 6) ? 1 : (int64_t(1) << (n - 6));
    for (int64_t w = 0; w < words; w++) {
        uint64_t word = data[w];
        // Exclude constant (index 0) and degree-1 terms (indices 2^i)
        if (w == 0) {
            // Constant at bit 0, degree-1 at bits 1,2,4,8,16,32
            word &= ~0x100010117ULL;
        } else if ((w & (w - 1)) == 0) {
            // Power-of-2 word index → degree-1 term at bit 0 (x_{k} where k=6+log2(w))
            word &= ~1ULL;
        }
        total += __builtin_popcountll(word);
    }
    return total;
}

int64_t count_T_paired(const uint64_t* data, int m,
                       const uint64_t* pair_masks, int n_pairs)
{
    int64_t words = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    int64_t total = 0;
    for (int64_t w = 0; w < words; w++) {
        uint64_t val = data[w];
        while (val) {
            int bit = __builtin_ctzll(val);
            val &= val - 1;
            int64_t pos = (w << 6) | bit;
            bool skip = false;
            for (int p = 0; p < n_pairs; p++) {
                if ((pos & pair_masks[p]) == pair_masks[p]) {
                    skip = true;
                    break;
                }
            }
            if (!skip) total++;
        }
    }
    return total;
}

void filter_pairs(uint64_t* data, int m,
                  const uint64_t* pair_masks, int n_pairs)
{
    int64_t words = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    for (int64_t w = 0; w < words; w++) {
        uint64_t val = data[w];
        while (val) {
            int bit = __builtin_ctzll(val);
            val &= val - 1;
            int64_t pos = (w << 6) | bit;
            for (int p = 0; p < n_pairs; p++) {
                if ((pos & pair_masks[p]) == pair_masks[p]) {
                    data[w] &= ~(1ULL << bit);
                    break;
                }
            }
        }
    }
}

RawANFInfo compute_raw_anf_info(TruthTable& tt) {
    RawANFInfo info;
    info.per_output_T.resize(tt.n_outputs, 0);
    info.per_output_deg.resize(tt.n_outputs, 0);
    info.sum_T = 0;
    info.overall_max_deg = 0;

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        moebius_packed(tt.tt[oi].data(), tt.n);
        info.per_output_T[oi] = count_T(tt.tt[oi].data(), tt.n);

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
