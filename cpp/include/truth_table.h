#ifndef TRUTH_TABLE_H
#define TRUTH_TABLE_H

#include "types.h"
#include "circuit.h"
#include <cstdint>
#include <vector>

// ============================================================
//  FastCircuit — Pre-computed indices for fast batch evaluation
// ============================================================

struct FastCircuit {
    int n_stmts, n_inputs;
    std::vector<uint8_t> stmt_op;
    std::vector<int> stmt_inp1, stmt_inp2;
    std::vector<int> out_idx;
};

FastCircuit make_fast(const Circuit& circ, const std::vector<int>& output_indices);

// Batch circuit evaluation using pre-computed FastCircuit indices
void eval_batch_fast(
    const FastCircuit& fc,
    const uint64_t* input_words,
    uint64_t* eval_buf,
    uint64_t* out_results, int n_outputs);

// Standard truth table (all 2^n x-values, suitable for n ≤ 20)
TruthTable compute_truth_table(const Circuit& circ,
                               const std::vector<int>& output_indices,
                               int n_threads);

#endif
