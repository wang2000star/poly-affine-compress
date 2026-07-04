#include "search.h"
#include "affine.h"
#include "anf.h"
#include "gf2.h"
#include "io.h"
#include "moebius.h"
#include "verify.h"
#include "walsh.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <thread>

// ============================================================
//  CandidateGenerator
// ============================================================

struct CandidateGenerator {
    int n;
    int max_m_search;
    const std::vector<WalshInfo>& walsh;
    std::mt19937_64 rng;

    CandidateGenerator(int n, int max_m, const std::vector<WalshInfo>& w)
        : n(n), max_m_search(max_m), walsh(w), rng(42) {}

    // Single-row candidates from Walsh (each bit with high |W_i|)
    std::vector<std::pair<uint32_t, uint32_t>> gen_walsh_single_rows(int top_k = 40) {
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

        std::vector<std::pair<uint32_t, uint32_t>> candidates;
        int k = std::min(top_k, n);
        for (int j = 0; j < k; j++) {
            uint32_t row = 1u << scores[j].bit;
            candidates.emplace_back(row, 0u);
            candidates.emplace_back(row, 1u);
        }
        return candidates;
    }

    // Multi-row candidates by combining top Walsh bits
    std::vector<std::pair<std::vector<uint32_t>, uint32_t>> gen_multi_row(int max_rows = 12) {
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

        for (int m = 1; m <= std::min(max_rows, n); m++) {
            std::vector<uint32_t> rows;
            for (int j = 0; j < m; j++)
                rows.push_back(1u << scores[j].bit);
            candidates.emplace_back(rows, 0u);
        }

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

    // Random M,b candidates
    std::vector<std::pair<std::vector<uint32_t>, uint32_t>> gen_random(
        int n_candidates, int max_m)
    {
        std::vector<std::pair<std::vector<uint32_t>, uint32_t>> candidates;
        for (int c = 0; c < n_candidates; c++) {
            int m = (rng() % max_m) + 1;
            std::vector<uint32_t> rows(m, 0);
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    if (rng() & 1) rows[j] |= (1u << i);
                }
            }
            uint64_t b = 0;
            for (int j = 0; j < m; j++) {
                if (rng() & 1) b |= (1ULL << j);
            }
            candidates.emplace_back(rows, b);
        }
        return candidates;
    }
};

// ============================================================
//  Complement pair detection and row-pool evaluation
//
//  Build a 3n-row pool (W,E) as suggested by the user:
//    Rows 0..n-1:   z_i = x_i          (identity, b=0)
//    Rows n..2n-1:  z_{n+i} = x_i ⊕ 1  (complement, b=1)
//    Rows 2n..3n-1: XOR of ≥2 input bits (b arbitrary)
//
//  From this pool, select any s rows as M,b.
//  Complement pairs: if row i (single-bit e_c, b=0) AND row j (single-bit e_c, b=1)
//  are both selected, their monomials are filtered post-Möbius.
// ============================================================

// Detect complement pairs in M,b: pairs of single-bit rows with same bit but opposite b.
// Returns number of pairs found, fills pair_masks with (1<<row_a)|(1<<row_b) for each.
static int detect_complement_pairs(const uint32_t* M_rows, uint64_t b, int m,
                                    uint64_t* pair_masks_out, int max_pairs)
{
    int np = 0;
    for (int i = 0; i < m && np < max_pairs; i++) {
        uint32_t ri = M_rows[i];
        if (ri == 0) continue;
        int bit_i = __builtin_ctz(ri);
        if (ri != (1u << bit_i)) continue;
        for (int j = i + 1; j < m && np < max_pairs; j++) {
            uint32_t rj = M_rows[j];
            if (rj == 0) continue;
            int bit_j = __builtin_ctz(rj);
            if (rj != (1u << bit_j)) continue;
            if (bit_i == bit_j && ((b >> i) & 1) != ((b >> j) & 1))
                pair_masks_out[np++] = (1ULL << i) | (1ULL << j);
        }
    }
    return np;
}

// Evaluate M,b with complement-pair filtering (if any pairs detected).
static MbCandidate evaluate_pool_rowset(
    const TruthTable& tt,
    const uint32_t* M_rows, uint64_t b, int m,
    int n, int n_threads)
{
    // Detect complement pairs
    uint64_t pair_masks[16];
    int np = detect_complement_pairs(M_rows, b, m, pair_masks, 16);

    MbCandidate cand = evaluate_Mb(tt, M_rows, b, m, n_threads, np > 0);
    if (cand.union_T >= INT64_MAX) return cand;

    if (np == 0) return cand;  // no filtering needed

    // Filter: Möbius → remove pair monomials → inverse Möbius
    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    std::vector<uint64_t> union_anf(n_words_g, 0);
    int64_t total_T = 0;

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        auto& raw = cand.g_tt_raw[oi];
        if (raw.empty()) { cand.per_output_T[oi] = 0; continue; }
        moebius_packed(raw.data(), m);
        filter_pairs(raw.data(), m, pair_masks, np);
        int64_t T_filt = count_T(raw.data(), m);
        cand.per_output_T[oi] = T_filt;
        total_T += T_filt;
        for (int64_t w = 0; w < n_words_g; w++)
            union_anf[w] |= raw[w];
        moebius_packed(raw.data(), m);
    }
    cand.total_T = total_T;
    cand.union_T = 0;
    for (int64_t w = 0; w < n_words_g; w++)
        cand.union_T += __builtin_popcountll(union_anf[w]);
    return cand;
}

// Build the 3n-row pool and score each row by Walsh magnitude.
// Row types:
//   0..n-1:   identity, b=0
//   n..2n-1:  complement, b=1
//   2n..3n-1: XOR of two top-Walsh bits (each row has ≥2 ones)
// Returns a vector of {M_row, b_row, score}.
static void build_row_pool(int n,
    const std::vector<WalshInfo>& walsh,
    std::vector<std::pair<uint32_t, uint32_t>>& pool_rows,
    std::vector<double>& scores)
{
    pool_rows.clear();
    scores.clear();

    // Score each input bit by Walsh magnitude
    struct BitScore { int bit; double score; };
    std::vector<BitScore> bit_scores(n);
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (auto& w : walsh)
            s += (double)w.walsh_mag[i] / (1 << (n - 1));
        bit_scores[i] = {i, s};
    }
    std::sort(bit_scores.begin(), bit_scores.end(),
        [](auto& a, auto& b) { return a.score > b.score; });

    // Group 1: identity rows (b=0) and complement rows (b=1)
    for (int i = 0; i < n; i++) {
        uint32_t row = (1u << i);
        pool_rows.emplace_back(row, 0u);
        scores.push_back(bit_scores[i].score);
        pool_rows.emplace_back(row, 1u);
        scores.push_back(bit_scores[i].score);
    }

    // Group 2: XOR pairs of top Walsh bits (≥2 ones per row)
    int top_n = std::min(n, 8);
    for (int a = 0; a < top_n; a++) {
        for (int b_idx = a + 1; b_idx < top_n; b_idx++) {
            uint32_t row = (1u << bit_scores[a].bit) | (1u << bit_scores[b_idx].bit);
            double score = bit_scores[a].score + bit_scores[b_idx].score;
            pool_rows.emplace_back(row, 0u);
            scores.push_back(score);
            pool_rows.emplace_back(row, 1u);
            scores.push_back(score);
        }
    }
}

// Generate candidates from the 3n-row pool.
// Strategy: start from identity transform, try swapping/replacing rows.
static void add_row_pool_candidates(
    const TruthTable& tt,
    const std::vector<WalshInfo>& walsh,
    int n, int n_threads, int max_candidates,
    std::vector<MbCandidate>& results)
{
    if (walsh.empty() || max_candidates <= 0) return;
    int pool_size = 2 * n + 2 * std::min(n, 8) * (std::min(n, 8) - 1) / 2;
    if (pool_size <= 0) return;

    std::vector<std::pair<uint32_t, uint32_t>> pool_rows;
    std::vector<double> scores;
    build_row_pool(n, walsh, pool_rows, scores);

    // Index each bit's identity row and complement row in the pool.
    // Pool layout: first n = identity(b=0), next n = identity(b=1), then XOR-pair rows.
    // So for bit i: id_row = i, comp_row = i + n.
    std::vector<int> id_idx(n), comp_idx(n);
    for (int i = 0; i < n; i++) {
        id_idx[i] = i;
        comp_idx[i] = i + n;
    }

    int generated = 0;
    auto try_Mb = [&](const std::vector<int>& selected, const char* label) {
        if (generated >= max_candidates) return;
        int m = (int)selected.size();
        uint32_t M[32] = {0};
        uint64_t bv = 0;
        for (int j = 0; j < m; j++) {
            M[j] = pool_rows[selected[j]].first;
            if (pool_rows[selected[j]].second) bv |= (1ULL << j);
        }
        auto cand = evaluate_pool_rowset(tt, M, bv, m, n, n_threads);
        if (cand.union_T < INT64_MAX) {
            results.push_back(cand);
            generated++;
            std::cout << "      pool " << label << " m=" << m
                      << " T_union=" << cand.union_T << "\n";
        }
    };

    // 1) Identity transform (baseline)
    {
        std::vector<int> sel(n);
        for (int i = 0; i < n; i++) sel[i] = id_idx[i];
        try_Mb(sel, "identity");
    }

    // 2) Full complement set: identity + complement rows (m=2n)
    if (generated < max_candidates) {
        std::vector<int> sel;
        for (int i = 0; i < n; i++) {
            sel.push_back(id_idx[i]);
            sel.push_back(comp_idx[i]);
        }
        try_Mb(sel, "full-complement");
    }

    // 3) Try swapping each identity row with its complement (one at a time)
    //    M has identity rows, but bit i is swapped to complement (b=1).
    if (generated < max_candidates) {
        for (int i = 0; i < n && generated < max_candidates; i++) {
            std::vector<int> sel;
            for (int j = 0; j < n; j++)
                sel.push_back(j == i ? comp_idx[j] : id_idx[j]);
            try_Mb(sel, "swap-comp");
        }
    }

    // 4) Try swapping pairs of identity rows with complements
    if (generated < max_candidates) {
        for (int i = 0; i < n && generated < max_candidates; i++) {
            for (int j = i + 1; j < n && generated < max_candidates; j++) {
                std::vector<int> sel;
                for (int k = 0; k < n; k++)
                    sel.push_back((k == i || k == j) ? comp_idx[k] : id_idx[k]);
                try_Mb(sel, "swap2-comp");
            }
        }
    }

    // 5) Greedy: start from best complement mix, try flipping more bits
    if (generated < max_candidates) {
        // Find best combination so far by scanning results
        int64_t best_union = INT64_MAX;
        std::vector<int> best_config(n, 0); // 0 = identity, 1 = complement
        for (auto& r : results) {
            if (r.m == n && r.union_T < best_union) {
                best_union = r.union_T;
            }
        }
        // Greedy: start from identity, flip each bit to complement if it helps
        std::vector<int> config(n, 0);
        bool improved = true;
        while (improved && generated < max_candidates) {
            improved = false;
            for (int i = 0; i < n && generated < max_candidates; i++) {
                if (config[i]) continue; // already flipped
                config[i] = 1; // try flip
                std::vector<int> sel;
                for (int j = 0; j < n; j++)
                    sel.push_back(config[j] ? comp_idx[j] : id_idx[j]);
                uint32_t M[32] = {0};
                uint64_t bv = 0;
                for (int j = 0; j < n; j++) {
                    M[j] = pool_rows[sel[j]].first;
                    if (pool_rows[sel[j]].second) bv |= (1ULL << j);
                }
                auto cand = evaluate_pool_rowset(tt, M, bv, n, n, n_threads);
                if (cand.union_T < best_union) {
                    best_union = cand.union_T;
                    best_config = config;
                    improved = true;
                    results.push_back(cand);
                    generated++;
                    std::cout << "      pool greedy-flip bit " << i
                              << " T_union=" << cand.union_T << "\n";
                } else {
                    config[i] = 0; // revert
                }
            }
        }
    }

    // 6) Try adding XOR-pair rows on top of identity (m > n)
    if (generated < max_candidates) {
        int n_xor_pairs = std::min(12, (int)pool_rows.size() - 2 * n);
        for (int xi = 0; xi < n_xor_pairs && generated < max_candidates; xi++) {
            int pool_idx = 2 * n + xi; // XOR-pair rows start at 2n
            std::vector<int> sel;
            for (int j = 0; j < n; j++) sel.push_back(id_idx[j]);
            sel.push_back(pool_idx); // add XOR-pair row
            try_Mb(sel, "add-xor");
        }
    }

    std::cout << "    Pool candidates: " << generated << "\n";
}

// ============================================================
//  Hill climbing
// ============================================================

// ============================================================
//  Hill climbing
// ============================================================

static MbCandidate evaluate_with_pair_filter(
    const TruthTable& tt,
    const uint32_t* M_rows, uint64_t b, int m, int n_threads,
    const uint64_t* pair_masks, int n_pairs)
{
    MbCandidate cand = evaluate_Mb(tt, M_rows, b, m, n_threads, true);
    if (cand.union_T >= INT64_MAX) return cand;

    int64_t n_words_g = (m < 6) ? 1 : (int64_t(1) << (m - 6));
    std::vector<uint64_t> union_anf(n_words_g, 0);
    int64_t total_T = 0;

    for (int oi = 0; oi < tt.n_outputs; oi++) {
        auto& raw = cand.g_tt_raw[oi];
        if (raw.empty()) { cand.per_output_T[oi] = 0; continue; }
        moebius_packed(raw.data(), m);
        filter_pairs(raw.data(), m, pair_masks, n_pairs);
        int64_t T_filt = count_T(raw.data(), m);
        cand.per_output_T[oi] = T_filt;
        total_T += T_filt;
        for (int64_t w = 0; w < n_words_g; w++)
            union_anf[w] |= raw[w];
        moebius_packed(raw.data(), m);
    }
    cand.total_T = total_T;
    cand.union_T = 0;
    for (int64_t w = 0; w < n_words_g; w++)
        cand.union_T += __builtin_popcountll(union_anf[w]);
    return cand;
}

// Hill climb from a complement candidate: modify non-paired identity rows
// but keep complement rows and their paired identity rows fixed.
static MbCandidate hill_climb_complement(
    const TruthTable& tt, const MbCandidate& start,
    int n, int n_threads,
    const uint64_t* pair_masks, int n_pairs)
{
    MbCandidate best = start;
    int m = start.m;
    int k = m - n;

    // Identify which identity rows are paired with complement rows (must stay fixed)
    bool row_locked[32] = {false};
    for (int j = 0; j < k; j++) {
        // Extract the complement bit from pair_masks[j]:
        // pair_masks[j] = (1<<comp_bit) | (1<<(n+j))
        int comp_bit = __builtin_ctzll(pair_masks[j]);
        row_locked[comp_bit] = true;
    }

    MbCandidate cur = best;
    bool improved;
    int pass = 0;
    const int MAX_PASSES = 50;

    do {
        improved = false;
        pass++;

        // Flip bits in non-locked identity rows only
        for (int r = 0; r < n && !improved; r++) {
            if (row_locked[r]) continue;
            for (int c = 0; c < n && !improved; c++) {
                cur.M_rows[r] ^= (1u << c);
                MbCandidate cand = evaluate_with_pair_filter(
                    tt, cur.M_rows.data(), cur.b, m, n_threads, pair_masks, n_pairs);
                if (cand.union_T < best.union_T && cand.union_T != INT64_MAX)
                    best = cand;
                if (cand.union_T < cur.union_T && cand.union_T != INT64_MAX) {
                    cur = cand;
                    improved = true;
                } else {
                    cur.M_rows[r] ^= (1u << c);
                }
            }
        }

        // Flip b-bit in non-locked rows only
        for (int r = 0; r < n && !improved; r++) {
            if (row_locked[r]) continue;
            cur.b ^= (1ULL << r);
            MbCandidate cand = evaluate_with_pair_filter(
                tt, cur.M_rows.data(), cur.b, m, n_threads, pair_masks, n_pairs);
            if (cand.union_T < best.union_T && cand.union_T != INT64_MAX)
                best = cand;
            if (cand.union_T < cur.union_T && cand.union_T != INT64_MAX) {
                cur = cand;
                improved = true;
            } else {
                cur.b ^= (1ULL << r);
            }
        }
    } while (improved && pass < MAX_PASSES);

    return best;
}

// Old hill_climb_complement kept for reference, redirect to new version
// (The new version above replaces the old one.)

MbCandidate hill_climb(const TruthTable& tt, const MbCandidate& start,
                       int n, int n_threads)
{
    MbCandidate best = start;
    int m = start.m;

    auto eval = [&](uint32_t* M, uint64_t b) -> MbCandidate {
        return evaluate_Mb(tt, M, b, m, n_threads);
    };

    MbCandidate cur = best;
    bool improved;
    int pass = 0;
    const int MAX_PASSES = 50;

    do {
        improved = false;
        pass++;

        for (int r = 0; r < m && !improved; r++) {
            for (int c = 0; c < n && !improved; c++) {
                cur.M_rows[r] ^= (1u << c);
                MbCandidate cand = eval(cur.M_rows.data(), cur.b);
                if (cand.union_T < best.union_T && cand.union_T != INT64_MAX)
                    best = cand;
                if (cand.union_T < cur.union_T && cand.union_T != INT64_MAX) {
                    cur = cand;
                    improved = true;
                } else {
                    cur.M_rows[r] ^= (1u << c);
                }
            }
        }

        for (int r = 0; r < m && !improved; r++) {
            cur.b ^= (1ULL << r);
            MbCandidate cand = eval(cur.M_rows.data(), cur.b);
            if (cand.union_T < best.union_T && cand.union_T != INT64_MAX)
                best = cand;
            if (cand.union_T < cur.union_T && cand.union_T != INT64_MAX) {
                cur = cand;
                improved = true;
            } else {
                cur.b ^= (1ULL << r);
            }
        }
    } while (improved && pass < MAX_PASSES);

    return best;
}

// ============================================================
//  run_search — full pipeline (Phase 1..6)
// ============================================================

void run_search(const TruthTable& tt, const Circuit& circ,
                const std::vector<int>& output_indices,
                const SearchParams& params)
{
    int n = tt.n;
    auto t0 = std::chrono::steady_clock::now();

    // Phase 2: Raw ANF baseline
    std::cout << "\nPhase 2: Raw ANF baseline...\n";
    TruthTable tt_copy = tt;  // mutable copy for ANF transform
    RawANFInfo raw = compute_raw_anf_info(tt_copy);
    std::cout << "  Sum T = " << raw.sum_T << "\n";
    std::cout << "  Max deg = " << raw.overall_max_deg << "\n";
    for (int oi = 0; oi < tt.n_outputs; oi++) {
        std::cout << "    " << circ.outputs[output_indices[oi]]
                  << ": T=" << raw.per_output_T[oi]
                  << ", m=" << raw.per_output_deg[oi] << "\n";
    }

    // Save raw ANF if requested
    if (!params.anf_out_prefix.empty() && n <= 16) {
        save_raw_anf(tt_copy, circ, output_indices, params.anf_out_prefix + "_expr.poly");
        save_raw_T(tt_copy, circ, output_indices, params.anf_out_prefix + "_T.poly");
    }

    auto t1 = std::chrono::steady_clock::now();

    // Phase 3: Walsh correlations (n ≤ 32)
    std::vector<WalshInfo> walsh;
    if (n <= 32) {
        std::cout << "\nPhase 3: Computing Walsh correlations...\n";
        walsh = compute_walsh_correlations(tt_copy, params.n_threads);
        for (int oi = 0; oi < tt.n_outputs && oi < 8; oi++) {
            std::cout << "  Output " << circ.outputs[output_indices[oi]] << " top-5 bits:";
            // Find top 5 magnitudes
            struct { int idx; int64_t mag; } sorted[32];
            for (int i = 0; i < n; i++) sorted[i] = {i, walsh[oi].walsh_mag[i]};
            std::sort(sorted, sorted + n,
                [](auto& a, auto& b) { return a.mag > b.mag; });
            for (int j = 0; j < std::min(5, n); j++)
                std::cout << " " << sorted[j].idx << "(" << sorted[j].mag << ")";
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

    // 4a-c: Walsh-guided + random (n ≤ 20 full search; 21-32 permutation only)
    if (n <= 20) {
        CandidateGenerator gen(n, params.max_m, walsh);

        std::cout << "  4a: Single-row Walsh (" << (params.walsh_single_top * 2) << " candidates)\n";
        auto single_rows = gen.gen_walsh_single_rows(params.walsh_single_top);
        for (auto& [row, b] : single_rows) {
            uint32_t M[32] = {0};
            M[0] = row;
            results.push_back(evaluate_Mb(tt_copy, M, b, 1, params.n_threads));
        }

        std::cout << "  4b: Multi-row Walsh (" << params.multi_max_rows
                  << " progressive + pair combinations)\n";
        auto multi_rows = gen.gen_multi_row(params.multi_max_rows);
        for (auto& [rows, b] : multi_rows) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            results.push_back(evaluate_Mb(tt_copy, M, b, m, params.n_threads));
        }

        std::cout << "  4c: Random M,b (" << params.n_random << " candidates)\n";
        auto random_cands = gen.gen_random(params.n_random, params.max_m);
        for (auto& [rows, b] : random_cands) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            results.push_back(evaluate_Mb(tt_copy, M, b, m, params.n_threads));
        }

        auto t_search0 = std::chrono::steady_clock::now();
        double walsh_time = std::chrono::duration<double>(t_search0 - t1).count();
        std::cout << "  Search time: " << walsh_time << " s\n";
    } else if (n <= 32 && !walsh.empty()) {
        CandidateGenerator gen(n, params.max_m, walsh);

        std::cout << "  4a: Single-row Walsh (" << (params.walsh_single_top * 2) << " candidates)\n";
        auto single_rows = gen.gen_walsh_single_rows(params.walsh_single_top);
        for (auto& [row, b] : single_rows) {
            uint32_t M[32] = {0};
            M[0] = row;
            results.push_back(evaluate_Mb(tt_copy, M, b, 1, params.n_threads));
        }

        std::cout << "  4b: Multi-row Walsh (progressive m=1.." << std::min(20, n)
                  << " + XOR pairs)\n";
        auto multi_rows = gen.gen_multi_row(std::min(20, n));
        for (auto& [rows, b] : multi_rows) {
            int m = (int)rows.size();
            uint32_t M[32] = {0};
            for (int j = 0; j < m; j++) M[j] = rows[j];
            results.push_back(evaluate_Mb(tt_copy, M, b, m, params.n_threads));
        }

        std::cout << "  (skipping Phase 4c small-m random for n=" << n << ")\n";
    }

    // 4d: n32 random m=n candidates (for n > 20)
    if (params.n32_random > 0 && n > 20) {
        std::cout << "  4d: n32 random m=n candidates (" << params.n32_random << ")\n";
        std::mt19937_64 rng_n32(42);
        int generated = 0, attempts = 0;
        const int MAX_ATTEMPTS = params.n32_random * 20;

        auto t_n32_0 = std::chrono::steady_clock::now();
        while (generated < params.n32_random && attempts < MAX_ATTEMPTS) {
            attempts++;
            uint32_t M[32] = {0};
            for (int r = 0; r < n; r++)
                for (int i = 0; i < n; i++)
                    if (rng_n32() & 1) M[r] |= (1u << i);

            uint64_t b = 0;
            for (int r = 0; r < n; r++)
                if (rng_n32() & 1) b |= (1ULL << r);

            if (gf2_rank(M, n, n) < n) continue;

            auto t_eval0 = std::chrono::steady_clock::now();
            auto cand = evaluate_Mb(tt_copy, M, b, n, params.n_threads);
            auto t_eval1 = std::chrono::steady_clock::now();
            results.push_back(cand);
            generated++;

            double eval_time = std::chrono::duration<double>(t_eval1 - t_eval0).count();
            double compression = (double)raw.sum_T / std::max(int64_t(1), cand.union_T);
            std::cout << "    n32 #" << generated << "/" << params.n32_random
                      << " (attempt " << attempts << ")"
                      << ": T=" << cand.union_T
                      << " (compression " << compression << "×)"
                      << " time=" << eval_time << "s\n";
        }
        auto t_n32_1 = std::chrono::steady_clock::now();
        std::cout << "    n32 random time: " << std::chrono::duration<double>(t_n32_1 - t_n32_0).count() << " s\n";
    }

    // 4f: Complement variable candidates (n ≤ 16)
    if (params.n_complement > 0 && n <= 16) {
        std::cout << "  4f: Complement candidates (" << params.n_complement << " max)\n";
        add_row_pool_candidates(tt_copy, walsh, n, params.n_threads,
                                 params.n_complement, results);
    }

    // 4e: Hill climbing from top candidates
    if (params.n_hill_climb > 0) {
        std::sort(results.begin(), results.end(),
            [](auto& a, auto& b) { return a.union_T < b.union_T; });

        int n_climb = std::min(params.n_hill_climb, (int)results.size());
        std::cout << "  4e: Hill climbing from top " << n_climb << " candidates (+ identity)\n";

        auto t_hc0 = std::chrono::steady_clock::now();
        for (int ci = 0; ci < n_climb; ci++) {
            auto& base = results[ci];
            if (base.union_T >= INT64_MAX) continue;

            if (base.m == n) {
                // Standard m=n hill climb
                auto improved = hill_climb(tt_copy, base, n, params.n_threads);
                if (improved.union_T < base.union_T) {
                    results.push_back(improved);
                    std::cout << "    hill climb #" << ci << ": T=" << base.union_T << " -> T=" << improved.union_T << "\n";
                }
            } else if (base.m > n && base.m <= 32) {
                // Complement m>n hill climb: detect complement pairs dynamically
                uint64_t hc_pair_masks[16];
                int n_hc_pairs = detect_complement_pairs(
                    base.M_rows.data(), base.b, base.m, hc_pair_masks, 16);

                if (n_hc_pairs > 0) {
                    auto improved = hill_climb_complement(tt_copy, base, n, params.n_threads,
                                                           hc_pair_masks, n_hc_pairs);
                    if (improved.union_T < base.union_T) {
                        results.push_back(improved);
                        std::cout << "    hill climb (complement) #" << ci
                                  << ": T=" << base.union_T << " -> T=" << improved.union_T << "\n";
                    }
                }
            }
        }

        // Hill climb from identity
        uint32_t ident[32] = {0};
        for (int i = 0; i < n; i++) ident[i] = (1u << i);
        MbCandidate ident_cand = evaluate_Mb(tt_copy, ident, 0, n, params.n_threads);
        auto ident_improved = hill_climb(tt_copy, ident_cand, n, params.n_threads);
        if (ident_improved.union_T < ident_cand.union_T) {
            results.push_back(ident_improved);
            std::cout << "    hill climb from identity: T=" << ident_cand.union_T << " -> T=" << ident_improved.union_T << "\n";
        }

        auto t_hc1 = std::chrono::steady_clock::now();
        std::cout << "    Hill climb time: " << std::chrono::duration<double>(t_hc1 - t_hc0).count() << " s\n";
    }

    // Phase 5: Report
    std::sort(results.begin(), results.end(),
        [](auto& a, auto& b) { return a.union_T < b.union_T; });

    std::cout << "\n========================================\n";
    std::cout << "  Results (" << results.size() << " candidates evaluated)\n";
    std::cout << "  Raw ANF Sum T = " << raw.sum_T << "\n";
    std::cout << "========================================\n";

    int n_report = std::min(20, (int)results.size());
    for (int j = 0; j < n_report; j++) {
        auto& r = results[j];
        double compression = (double)raw.sum_T / std::max(int64_t(1), r.union_T);
        std::cout << "  " << (j + 1) << ". m=" << r.m << " T=" << r.total_T << " (union=" << r.union_T
                  << ", compression " << compression << "×) M=[";
        for (int row = 0; row < r.m; row++) {
            if (row > 0) std::cout << "; ";
            for (int i = 0; i < std::min(n, 8); i++)
                std::cout << ((r.M_rows[row] >> i) & 1);
            std::cout << "..";
        }
        std::cout << "]\n";
        if (r.union_T < raw.sum_T && r.m <= 32) {
            std::cout << "    Per-output T:";
            for (int oi = 0; oi < tt.n_outputs; oi++)
                std::cout << " " << circ.outputs[output_indices[oi]] << "=" << r.per_output_T[oi];
            std::cout << "\n";
        }
    }

    // Phase 6: Verify best candidate, save results
    if (!results.empty() && results[0].union_T < INT64_MAX) {
        auto& best = results[0];
        std::cout << "\nBest candidate: m=" << best.m << ", T_sum=" << best.total_T << ", T_union=" << best.union_T << "\n";
    }

    double total_time = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    // Phase 6: verification + file saving
    if (!results.empty() && results[0].union_T < INT64_MAX) {
        auto& best = results[0];

        // Print transform
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

        // Detect complement pairs dynamically (works with any row ordering)
        bool is_complement = false;
        uint64_t comp_pair_masks[16] = {0};
        int n_comp_pairs = 0;
        if (best.m > n && best.m <= 32) {
            n_comp_pairs = detect_complement_pairs(
                best.M_rows.data(), best.b, best.m, comp_pair_masks, 16);
            is_complement = (n_comp_pairs > 0);
        }

        // Re-evaluate with save_g_tt=true (with pair filtering for complement)
        bool all_verified = true;
        MbCandidate verified_cand;
        verified_cand.m = best.m;
        verified_cand.b = best.b;
        verified_cand.M_rows.resize(best.m);
        for (int i = 0; i < best.m; i++) verified_cand.M_rows[i] = best.M_rows[i];
        verified_cand.union_T = best.union_T;
        verified_cand.per_output_T = best.per_output_T;

        if (best.m == n) {
            std::cout << "Verifying best candidate (5000 random tests per output)...\n";
            verified_cand = evaluate_Mb(tt_copy, best.M_rows.data(), best.b, best.m,
                                         params.n_threads, true);
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                bool ok = verify_transform(tt_copy, verified_cand.g_tt_raw[oi], oi,
                    best.M_rows.data(), best.b, best.m, n, 5000);
                std::cout << "  " << circ.outputs[output_indices[oi]]
                          << (ok ? " ✅ Verified (5000 tests)" : " ❌ FAILED") << "\n";
                if (!ok) all_verified = false;
            }
        } else if (is_complement) {
            std::cout << "Verifying complement candidate (5000 random tests per output)...\n";
            verified_cand = evaluate_Mb(tt_copy, best.M_rows.data(), best.b, best.m,
                                         params.n_threads, true);
            // Clean g_tt_raw: Möbius → filter pairs → inverse Möbius
            int64_t n_words_g = (best.m < 6) ? 1 : (int64_t(1) << (best.m - 6));
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                auto& raw = verified_cand.g_tt_raw[oi];
                if (raw.empty()) continue;
                moebius_packed(raw.data(), best.m);
                filter_pairs(raw.data(), best.m, comp_pair_masks, n_comp_pairs);
                moebius_packed(raw.data(), best.m);
            }
            // Verify with cleaned truth table
            for (int oi = 0; oi < tt.n_outputs; oi++) {
                bool ok = verify_transform(tt_copy, verified_cand.g_tt_raw[oi], oi,
                    best.M_rows.data(), best.b, best.m, n, 5000);
                std::cout << "  " << circ.outputs[output_indices[oi]]
                          << (ok ? " ✅ Verified (5000 tests)" : " ❌ FAILED") << "\n";
                if (!ok) all_verified = false;
            }
        } else {
            std::cout << "Skipping verification (non-bijective m=" << best.m << " != n=" << n << ")\n";
        }

        if (all_verified)
            std::cout << "✅ All outputs verified!\n";

        // Save results to files
        if (!params.save_results_prefix.empty()) {
            std::string prefix = params.save_results_prefix;

            save_trans(best, circ, n, prefix + "_trans.poly");

            if (!verified_cand.g_tt_raw.empty()) {
                save_opt_expr(verified_cand.g_tt_raw, circ, output_indices, best.m,
                             prefix + "_expr.poly");
                save_opt_T(verified_cand.g_tt_raw, circ, output_indices, best.m,
                          prefix + "_T.poly");
            }

            save_verify(tt_copy, verified_cand.g_tt_raw, best, n, output_indices, circ,
                       prefix + "_verify.txt");
        }
    }

    total_time = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::cout << "\nTotal time: " << total_time << " s\n";
}
