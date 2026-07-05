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

// ---- d1c: product-based ANF simplification via exhaustive b-search ----

// Count degree-2+ nonlinear terms in a packed ANF (excludes constant + degree-1).
int64_t count_T_nl_packed(const uint64_t* anf, int n);

// Upward Möbius transform for a single bit position (selective complement).
// Applies the transform for bit i on packed ANF data: g(z) = f(z ⊕ e_i).
void upward_pass_bit(uint64_t* data, int n, int i);

// Exhaustive search over all 2^t b-vectors for optimal complement.
// Returns best T_nl, sets best_b to the corresponding t-bit pattern.
// Returns -1 if t > max_exhaustive_t.
int64_t exhaustive_search_best_b(const uint64_t* f_anf, int n,
                                  uint64_t support_mask, int max_exhaustive_t,
                                  uint64_t& best_b);

// Per-output info: support mask and T_nl (for all outputs).
struct OutputInfo {
    uint64_t support_mask;
    int t;
    int64_t T_nl;
};

// Compute T_nl and support for all outputs in tt_copy.
// Applies Möbius, computes info, restores truth table.
void compute_output_info(TruthTable& tt, int n, std::vector<OutputInfo>& info);

#endif // SEARCH_H
