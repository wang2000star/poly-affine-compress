#ifndef ANF_H
#define ANF_H

#include "types.h"

int64_t count_T(const uint64_t* data, int n);

// Count T, skipping monomials containing both bits of any complement pair.
// pair_masks: precomputed (1<<a)|(1<<b) for each pair (a,b).
int64_t count_T_paired(const uint64_t* data, int m,
                       const uint64_t* pair_masks, int n_pairs);

// Remove monomials containing both bits of any pair (in-place).
void filter_pairs(uint64_t* data, int m,
                  const uint64_t* pair_masks, int n_pairs);

RawANFInfo compute_raw_anf_info(TruthTable& tt);

#endif
