#ifndef GATE_BUILDER_H
#define GATE_BUILDER_H

#include "sparse_anf.h"
#include "circuit.h"
#include <unordered_map>
#include <string>
#include <vector>

// ====================================================================
//  SignalInfo: (g, M, b) triple for one circuit signal
// ====================================================================

struct SignalInfo {
    SparseANF g;                       // ANF in z-space (z = Mx⊕b)
    std::vector<uint32_t> M_rows;      // M matrix rows (each a bitmask of x vars)
    uint64_t b;                        // b vector (bit j = b_j)
    int m{0};                          // number of z variables = M_rows.size()

    SignalInfo() = default;
    SignalInfo(const SparseANF& g_, const std::vector<uint32_t>& M_, uint64_t b_, int m_)
        : g(g_), M_rows(M_), b(b_), m(m_) {}
};

// ====================================================================
//  GateBuilder: circuit-guided (g, M, b) construction
// ====================================================================

class GateBuilder {
public:
    GateBuilder() : next_m_(0), threshold_(4096) {}

    // Build (g, M, b) for all signals in the circuit
    // signals_ map is populated for each statement
    void build(const Circuit& circ);

    // Access results
    const SignalInfo* get_signal(const std::string& name) const;
    const std::unordered_map<std::string, SignalInfo>& all_signals() const { return signals_; }

    // Compression controls
    void set_threshold(int t) { threshold_ = t; }

private:
    std::unordered_map<std::string, SignalInfo> signals_;
    int next_m_;   // next available z variable index
    int threshold_;

    // Get signal info, handling aliases and inputs
    const SignalInfo* _resolve(const std::string& name, const Circuit& circ);

    // Gate handlers
    void _process_input(const std::string& lhs, const std::string& rhs, const Circuit& circ);
    void _process_not(const std::string& lhs, const std::string& arg, const Circuit& circ);
    void _process_xor(const std::string& lhs,
                      const std::string& arg1, const std::string& arg2,
                      const Circuit& circ);
    void _process_and(const std::string& lhs,
                      const std::string& arg1, const std::string& arg2,
                      const Circuit& circ);

    // Compression
    void _compress(SignalInfo& si, bool verbose = false);

    // Greedy merge: iteratively merge pairs to reduce m
    // Returns (M_comp, b_comp) that transforms old z → new z'
    void _greedy_merge(SparseANF& g, int& m,
                       std::vector<uint32_t>& comp_rows,
                       uint64_t& comp_b);

    // Compose M = M_new @ M_old, b = M_new @ b_old ⊕ b_new
    // M_new has (m_new × m_old), M_old has (m_old × n)
    static std::vector<uint32_t> _compose_M(
        const std::vector<uint32_t>& M_new, int m_new,
        const std::vector<uint32_t>& M_old, int m_old, int n);
    static uint64_t _compose_b(
        const std::vector<uint32_t>& M_new, int m_new,
        uint64_t b_old, uint64_t b_new);
};

#endif // GATE_BUILDER_H
