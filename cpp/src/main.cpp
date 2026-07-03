/**
 * optimize_anf — Truth-table-based ANF optimization for n ≤ 32.
 *
 * Core idea: given f(x) computed via truth table, search for affine
 * transform z = Mx + b (over GF(2)) such that g(z) = g(Mx + b) = f(x)
 * has fewer ANF terms (T(g) < T(f)).
 */

#include "circuit.h"
#include "truth_table.h"
#include "search.h"
#include <iostream>
#include <thread>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [options]\n";
        std::cerr << "  --max-m N      max z variables to search (default 12)\n";
        std::cerr << "  --walsh-k N    top K Walsh bits (default 30)\n";
        std::cerr << "  --random N     random candidates (default 40)\n";
        std::cerr << "  --n32-random N n>20 full-rank random m=n candidates (default 0)\n";
        std::cerr << "  --hill-climb N hill climb from top N candidates (default 10)\n";
        std::cerr << "  --anf-out PREF save raw ANF to PREFIX_expr.poly + PREFIX_T.poly (n<=16 only)\n";
        std::cerr << "  --save-results PREFIX  save best candidate: PREFIX_trans.poly,\n";
        std::cerr << "                          PREFIX_expr.poly, PREFIX_T.poly, PREFIX_verify.txt\n";
        return 1;
    }

    std::cout << std::unitbuf;
    std::string path = argv[1];

    std::cout << "--- Circuit ---\n";
    Circuit circ = read_circuit(path);
    std::cout << "  Inputs: " << circ.n_inputs << ", Outputs to optimize: " << circ.outputs.size() << "\n";

    std::vector<int> output_indices;
    for (int o = 0; o < (int)circ.outputs.size(); o++)
        output_indices.push_back(o);

    if (output_indices.empty()) {
        std::cerr << "No outputs to process\n";
        return 1;
    }

    SearchParams params;
    params.n_threads = std::thread::hardware_concurrency();
    if (params.n_threads < 1) params.n_threads = 1;

    for (int a = 2; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "--max-m" && a + 1 < argc)
            params.max_m = std::stoi(argv[++a]);
        else if (arg == "--walsh-k" && a + 1 < argc)
            params.walsh_single_top = std::stoi(argv[++a]);
        else if (arg == "--random" && a + 1 < argc)
            params.n_random = std::stoi(argv[++a]);
        else if (arg == "--n32-random" && a + 1 < argc)
            params.n32_random = std::stoi(argv[++a]);
        else if (arg == "--hill-climb" && a + 1 < argc)
            params.n_hill_climb = std::stoi(argv[++a]);
        else if (arg == "--anf-out" && a + 1 < argc)
            params.anf_out_prefix = argv[++a];
        else if (arg == "--save-results" && a + 1 < argc)
            params.save_results_prefix = argv[++a];
    }

    // Phase 1: Compute truth table
    std::cout << "\nPhase 1: Computing truth table...\n";
    TruthTable tt = compute_truth_table(circ, output_indices, params.n_threads);
    std::cout << "  Truth table: 2^" << circ.n_inputs << " = " << (int64_t(1) << circ.n_inputs)
              << " inputs, " << tt.n_words << " batches, " << tt.n_outputs
              << " output(s), " << params.n_threads << " thread(s)\n";

    run_search(tt, circ, output_indices, params);

    return 0;
}
