#ifndef SEARCH_H
#define SEARCH_H

#include "types.h"

// Hill climbing: iteratively flip M bits to reduce total_T
MbCandidate hill_climb(const TruthTable& tt, const MbCandidate& start,
                       int n, int n_threads);

// Full search pipeline (Phase 4-6): Walsh-guided + random + hill climbing
// Returns sorted candidates, writes results to cout, saves files if prefix set.
void run_search(const TruthTable& tt, const Circuit& circ,
                const std::vector<int>& output_indices,
                const SearchParams& params);

// ---- New theory-guided features ----

// Compute dependency set per output: which input variables affect each output.
// Returns vector of bitmasks (one per output).
// For n ≤ 20: uses truth table scanning (exact).
// For n > 20: uses random sampling.
std::vector<uint32_t> compute_dependency_set(const TruthTable& tt, int n,
                                              int n_samples = 500);

// Compute Aut(F) basis for n ≤ 20 via Walsh power spectrum.
// Returns basis vectors packed as uint64_t bitmasks.
// dim Aut(F) = size of returned vector.
std::vector<uint64_t> compute_aut_basis(const TruthTable& tt, int n,
                                         int n_threads);

// Progressive M construction (theory-guided):
// Start from identity, add rows one by one with exclusivity check.
// For each candidate row w:
//   1. Exclusivity: rank(M ∪ [w]) > rank(M)
//   2. Find Δ: M·Δ=0, w·Δ=1
//   3. Δ ∈ Aut(F)? → skip
// Returns best candidate found.
MbCandidate progressive_m_search(const TruthTable& tt,
                                  const std::vector<uint64_t>& aut_basis,
                                  int n, int max_m, int n_threads);

#endif // SEARCH_H
