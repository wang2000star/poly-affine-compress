/**
 * optimize_anf — Truth-table-based ANF optimization for n ≤ 32.
 *
 * Core idea: given f(x) computed via truth table, search for affine
 * transform z = Mx + b (over GF(2)) such that g(z) = g(Mx + b) = f(x)
 * has fewer ANF terms (T(g) < T(f)).
 *
 * Output (with --save-results DIR):
 *   {DIR}/{inst}_d1a_opt1.affine     — affine transform matrix
 *   {DIR}/{inst}_d1a_opt1.poly       — ANF polynomial matrix
 *   {DIR}/{inst}_d1a_opt1_stats.txt  — 5-line stats
 *   {DIR}/{inst}_d1a_opt1_verify.txt — verification results
 *   {DIR}/{inst}_raw.poly            — raw ANF (with --anf-out DIR)
 *   {DIR}/{inst}_raw_stats.txt       — raw stats (with --anf-out DIR)
 */

#include "circuit.h"
#include "truth_table.h"
#include "search.h"
#include <iostream>
#include <thread>
#include <filesystem>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.txt> [options]\n";
        std::cerr << "  --max-m N      max z variables to search (default 12)\n";
        std::cerr << "  --walsh-k N    top K Walsh bits (default 30)\n";
        std::cerr << "  --random N     random candidates (default 40)\n";
        std::cerr << "  --n32-random N n>20 full-rank random m=n candidates (default 0)\n";
        std::cerr << "  --complement N complement candidates (n≤16, default 0)\n";
        std::cerr << "  --hill-climb N hill climb from top N candidates (default 10)\n";
        std::cerr << "  --use-progressive 0/1  progressive M construction (default 1)\n";
        std::cerr << "  --progressive-max-m N  max m for progressive search (default auto)\n";
        std::cerr << "  --dep-filter 0/1  dependency-set row filtering (default 1)\n";
        std::cerr << "  --time-budget N   time budget in seconds (default 0 = unlimited)\n";
        std::cerr << "  --anf-out DIR  save raw ANF to {DIR}/{inst}_raw.poly + _stats.txt\n";
        std::cerr << "  --save-results DIR  save best candidate to {DIR}/{inst}_d1a_opt1.*\n";
        return 1;
    }

    std::cout << std::unitbuf;
    std::string path;
    std::string _opt_root;  // project root for resolving relative paths

    // ---- Resolve project root from executable path, then chdir ----
    {
        namespace _fs = std::filesystem;
        auto orig_cwd = _fs::current_path();
        auto exe = _fs::weakly_canonical(_fs::absolute(_fs::path(argv[0])));
        _opt_root = exe.parent_path().string();
        for (int i = 0; i < 3 && !_fs::exists(_opt_root + "/examples"); i++)
            _opt_root = _fs::path(_opt_root).parent_path().string();
        _fs::current_path(_opt_root);

        // Resolve circuit path (relative to original CWD) to absolute
        _fs::path circ_fs_path(argv[1]);
        if (circ_fs_path.is_relative()) circ_fs_path = orig_cwd / circ_fs_path;
        path = _fs::weakly_canonical(circ_fs_path).string();
    }

    // Extract instance name: "examples/hd08.txt" -> "hd08"
    std::string inst = path;
    size_t slash = inst.find_last_of('/');
    if (slash != std::string::npos) inst = inst.substr(slash + 1);
    size_t dot = inst.find_last_of('.');
    if (dot != std::string::npos) inst = inst.substr(0, dot);

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
    params.inst_name = inst;
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
        else if (arg == "--complement" && a + 1 < argc)
            params.n_complement = std::stoi(argv[++a]);
        else if (arg == "--hill-climb" && a + 1 < argc)
            params.n_hill_climb = std::stoi(argv[++a]);
        else if (arg == "--anf-out" && a + 1 < argc)
            params.anf_out_dir = argv[++a];
        else if (arg == "--save-results" && a + 1 < argc)
            params.results_dir = argv[++a];
        else if (arg == "--use-progressive" && a + 1 < argc)
            params.use_progressive = (std::stoi(argv[++a]) != 0);
        else if (arg == "--progressive-max-m" && a + 1 < argc)
            params.progressive_max_m = std::stoi(argv[++a]);
        else if (arg == "--dep-filter" && a + 1 < argc)
            params.use_dep_filter = (std::stoi(argv[++a]) != 0);
        else if (arg == "--time-budget" && a + 1 < argc)
            params.time_budget = std::stod(argv[++a]);
    }
    // Make output paths absolute (relative to project root)
    {
        std::string& s = params.anf_out_dir;
        if (!s.empty()) {
            std::filesystem::path p(s);
            if (p.is_relative()) s = (std::filesystem::path(_opt_root) / p).string();
        }
    }
    {
        std::string& s = params.results_dir;
        if (!s.empty()) {
            std::filesystem::path p(s);
            if (p.is_relative()) s = (std::filesystem::path(_opt_root) / p).string();
        }
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
