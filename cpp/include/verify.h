#ifndef VERIFY_H
#define VERIFY_H

#include "types.h"
#include <vector>

// Check g(Mx+b) == f(x) for n_tests random inputs
bool verify_transform(const TruthTable& tt,
                      const std::vector<uint64_t>& g_tt,
                      int output_idx,
                      const uint32_t* M_rows,
                      uint64_t b, int m, int n,
                      int n_tests = 5000);

// Pack all output bits at input x into a uint64_t
uint64_t pack_output_vector(const TruthTable& tt, uint32_t x);

#endif
