#ifndef AFFINE_H
#define AFFINE_H

#include "types.h"

// Evaluate a single output under transform z = Mx + b (non-bijective path).
// Returns {T, consistent}. If inconsistent, T = INT64_MAX.
MbResult evaluate_Mb_single_output(
    const uint64_t* f_tt, int n, int m,
    const uint32_t* M_rows, uint32_t b,
    int64_t n_words_f, int n_threads);

// Evaluate all outputs under bijective (m=n, full-rank) transform.
// With save_g_tt=true, saves pre-Möbius g_tt for verification.
MbCandidate evaluate_Mb_bijective(
    const TruthTable& tt,
    const uint32_t* M_rows, uint32_t b,
    int m, int n, int n_threads,
    bool save_g_tt = false);

// Top-level dispatcher: bijective path for m=n, else per-output non-bijective.
MbCandidate evaluate_Mb(
    const TruthTable& tt,
    const uint32_t* M_rows, uint32_t b,
    int m, int n_threads,
    bool save_g_tt = false);

#endif
