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

#endif
