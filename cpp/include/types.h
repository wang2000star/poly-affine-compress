#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

// ============================================================
//  Circuit types
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

// ============================================================
//  Truth table (bit-packed)
// ============================================================

struct TruthTable {
    int n;                              // number of input variables
    int n_outputs;                      // number of outputs being optimized
    int64_t n_words;                    // 2^n / 64
    std::vector<std::vector<uint64_t>> tt;   // tt[output][word]
    std::vector<int> output_indices;    // which original output indices
};

// ============================================================
//  Walsh correlation info
// ============================================================

struct WalshInfo {
    int output_idx;
    int64_t walsh_mag[32];  // |W_i| for each input bit (n ≤ 32)
    int64_t walsh_raw[32];  // raw W_i (positive = correlated with x_i=1)
};

// ============================================================
//  Single-output evaluation result (non-bijective path)
// ============================================================

struct MbResult {
    int64_t T;          // ANF term count (INT64_MAX if inconsistent)
    bool consistent;    // whether f is constant on each coset of ker(M)
};

// ============================================================
//  Multi-output candidate (one transform M,b)
// ============================================================

struct MbCandidate {
    int m;
    uint32_t M_rows[32];
    uint32_t b;
    int64_t total_T;
    std::vector<int64_t> per_output_T;
    std::vector<std::vector<uint64_t>> g_tt_raw;  // pre-Möbius per output
};

// ============================================================
//  Raw ANF statistics computed from truth table
// ============================================================

struct RawANFInfo {
    int64_t sum_T;
    std::vector<int64_t> per_output_T;
    std::vector<int> per_output_deg;
    int overall_max_deg;
};

// ============================================================
//  Search parameters
// ============================================================

struct SearchParams {
    int max_m = 12;
    int n_random = 40;
    int n_hill_climb = 10;
    int walsh_single_top = 20;
    int multi_max_rows = 10;
    int n32_random = 0;
    int n_threads = 4;
    int walsh_k = 30;
    std::string anf_out_prefix;
    std::string save_results_prefix;
};

#endif // TYPES_H
