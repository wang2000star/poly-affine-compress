#ifndef TRUTH_TABLE_H
#define TRUTH_TABLE_H

#include "types.h"

TruthTable compute_truth_table(const Circuit& circ,
                               const std::vector<int>& output_indices,
                               int n_threads);

#endif
