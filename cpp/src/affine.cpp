#include "affine.h"
#include "anf.h"
#include "gf2.h"
#include "moebius.h"
#include <algorithm>
#include <cstring>
#include <thread>
#include <vector>

// ============================================================
//  evaluate_Mb_single_output — non-bijective path (m < n or m=n with save_g_tt=false)
// ============================================================

MbResult evaluate_Mb_single_output(
    const uint64_t* f_tt, int n, int m,
    const uint32_t* M_rows, uint64_t b,
    int64_t n_words_f, int n_threads,
    std::vector<uint64_t>* g_tt_raw_out)
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

    auto compute = [&](int64_t w_start, int64_t w_end,
                       uint8_t* state, uint64_t* g_data) {
        for (int64_t w = w_start; w < w_end; w++) {
            uint64_t word = f_tt[w];
            int64_t base = w * 64;

            uint32_t z_base = b;
            for (int row = 0; row < m; row++) {
                if (__builtin_popcount(M_rows[row] & (uint32_t)base) & 1)
                    z_base ^= (1u << row);
            }

            uint64_t ones = word;
            uint64_t zeros = ~word;

            // 1-bits (f=1)
            uint64_t bits = ones;
            while (bits) {
                int bpos = __builtin_ctzll(bits);
                bits &= bits - 1;
                uint32_t z = z_base ^ M_offsets[bpos];
                uint8_t& st = state[z];
                if (st == 0)      { st = 2; g_data[z >> 6] |= (1ULL << (z & 63)); }
                else if (st == 1) { st = 3; }
            }

            // 0-bits (f=0)
            bits = zeros;
            while (bits) {
                int bpos = __builtin_ctzll(bits);
                bits &= bits - 1;
                uint32_t z = z_base ^ M_offsets[bpos];
                uint8_t& st = state[z];
                if (st == 2)      { st = 3; }
                else if (st == 0) { st = 1; }
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
            for (int64_t w = 0; w < std::max(int64_t(1), n_words_g); w++)
                per_thread_g[0][w] |= per_thread_g[t][w];
        }

        for (int zi = 0; zi < n_z; zi++) {
            if (per_thread_state[0][zi] == 3) { consistent = false; break; }
        }

        if (!consistent) return {INT64_MAX, false};

        if (g_tt_raw_out) {
            g_tt_raw_out->assign(per_thread_g[0].begin(), per_thread_g[0].end());
        }
        moebius_packed(per_thread_g[0].data(), m);
        return {count_T(per_thread_g[0].data(), m), true};
    } else {
        // m > 20: single-threaded
        if (m > 25) {
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
            if (g_tt_raw_out) *g_tt_raw_out = g_tt;
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

        if (g_tt_raw_out) *g_tt_raw_out = g_tt;
        moebius_packed(g_tt.data(), m);
        return {count_T(g_tt.data(), m), true};
    }
}

// ============================================================
//  evaluate_Mb_bijective — with optional save_g_tt for verification
// ============================================================

MbCandidate evaluate_Mb_bijective(
    const TruthTable& tt,
    const uint32_t* M_rows, uint64_t b,
    int m, int n, int n_threads,
    bool save_g_tt)
{
    MbCandidate cand;
    cand.m = m;
    cand.b = b;
    cand.M_rows.resize(m);
    for (int i = 0; i < m; i++) cand.M_rows[i] = M_rows[i];
    cand.total_T = 0;
    cand.union_T = 0;
    cand.per_output_T.resize(tt.n_outputs, 0);

    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    int64_t n_words_f = tt.n_words;

    uint32_t M_offsets[64];
    for (int d = 0; d < 64; d++) {
        uint32_t z = 0;
        for (int row = 0; row < m; row++) {
            if (__builtin_popcount(M_rows[row] & d) & 1)
                z |= (1u << row);
        }
        M_offsets[d] = z;
    }

    int use_threads = std::min(n_threads, 8);
    int64_t chunk = std::max(int64_t(1), (n_words_f + use_threads - 1) / use_threads);

    std::vector<uint64_t*> per_thread_g(use_threads, nullptr);
    for (int t = 0; t < use_threads; t++)
        per_thread_g[t] = (uint64_t*)calloc(n_words_g, sizeof(uint64_t));

    std::vector<uint64_t> union_anf(n_words_g, 0);  // union accumulator

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        const uint64_t* f_tt = tt.tt[oi].data();

        // Exact constant-zero check
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

        auto clear_fn = [&](int tid) {
            memset(per_thread_g[tid], 0, n_words_g * sizeof(uint64_t));
        };
        {
            std::vector<std::thread> clear_thr;
            for (int t = 0; t < use_threads; t++)
                clear_thr.emplace_back(clear_fn, t);
            for (auto& th : clear_thr) th.join();
        }

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

        // Merge: OR all per-thread g_tt into thread 0
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

        if (save_g_tt) {
            std::vector<uint64_t> g_copy(n_words_g, 0);
            memcpy(g_copy.data(), per_thread_g[0], n_words_g * sizeof(uint64_t));
            cand.g_tt_raw.push_back(std::move(g_copy));
        }

        moebius_packed_mt(per_thread_g[0], m, n_threads);
        int64_t T = count_T(per_thread_g[0], m);
        cand.per_output_T[oi] = T;
        cand.total_T += T;
        for (int64_t w = 0; w < n_words_g; w++)
            union_anf[w] |= per_thread_g[0][w];
    }

    cand.union_T = count_T(union_anf.data(), m);

    for (int t = 0; t < use_threads; t++)
        free(per_thread_g[t]);

    return cand;
}

// ============================================================
//  evaluate_Mb — top-level dispatcher
// ============================================================

MbCandidate evaluate_Mb(
    const TruthTable& tt,
    const uint32_t* M_rows, uint64_t b,
    int m, int n_threads,
    bool save_g_tt)
{
    MbCandidate cand;
    cand.m = m;
    cand.b = b;
    cand.M_rows.resize(m);
    for (int i = 0; i < m; i++) cand.M_rows[i] = M_rows[i];
    cand.total_T = 0;
    cand.union_T = 0;
    cand.per_output_T.resize(tt.n_outputs, 0);

    if (m == tt.n) {
        int rank = gf2_rank(M_rows, m, tt.n);
        if (rank < m) {
            cand.total_T = INT64_MAX;
            cand.union_T = INT64_MAX;
            return cand;
        }
        return evaluate_Mb_bijective(tt, M_rows, b, m, tt.n, n_threads, save_g_tt);
    }

    bool any_inconsistent = false;
    cand.g_tt_raw.clear();
    if (save_g_tt) cand.g_tt_raw.resize(tt.n_outputs);
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        std::vector<uint64_t>* g_out = save_g_tt ? &cand.g_tt_raw[oi] : nullptr;
        MbResult res = evaluate_Mb_single_output(
            tt.tt[oi].data(), tt.n, m, M_rows, b,
            tt.n_words, n_threads, g_out);
        if (!res.consistent) { any_inconsistent = true; break; }
        cand.per_output_T[oi] = res.T;
        cand.total_T += res.T;
    }

    if (any_inconsistent) {
        cand.total_T = INT64_MAX;
        cand.union_T = INT64_MAX;
    } else {
        int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
        if (save_g_tt) {
            std::vector<uint64_t> union_anf(n_words_g, 0);
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                std::vector<uint64_t> anf = cand.g_tt_raw[oi];
                moebius_packed(anf.data(), m);
                for (int64_t w = 0; w < n_words_g; w++)
                    union_anf[w] |= anf[w];
            }
            cand.union_T = count_T(union_anf.data(), m);
        } else {
            cand.union_T = cand.total_T;
        }
    }
    return cand;
}
