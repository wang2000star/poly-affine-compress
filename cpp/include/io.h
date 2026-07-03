#ifndef IO_H
#define IO_H

#include "types.h"
#include <string>

// --- Raw ANF output ---
void save_raw_anf(const TruthTable& tt, const Circuit& circ,
                  const std::vector<int>& output_indices,
                  const std::string& fname);
void save_raw_T(const TruthTable& tt, const Circuit& circ,
                const std::vector<int>& output_indices,
                const std::string& fname);

// --- Optimized ANF output ---
void save_opt_expr(const std::vector<std::vector<uint64_t>>& g_tt_raw,
                   const Circuit& circ,
                   const std::vector<int>& output_indices,
                   int m, const std::string& fname);
void save_opt_T(const std::vector<std::vector<uint64_t>>& g_tt_raw,
                const Circuit& circ,
                const std::vector<int>& output_indices,
                int m, const std::string& fname);

// --- Transform output ---
void save_trans(const MbCandidate& best, const Circuit& circ,
                int n, const std::string& fname);

// --- Verification output ---
void save_verify(const TruthTable& tt,
                 const std::vector<std::vector<uint64_t>>& g_tt_raw,
                 const MbCandidate& best, int n,
                 const std::vector<int>& output_indices,
                 const Circuit& circ,
                 const std::string& fname);

#endif
