#include "gate_builder.h"
#include <cstdio>
#include <cassert>
#include <algorithm>

// ====================================================================
//  Internal helpers
// ====================================================================

// Remove variable var_idx from ANF, shifting higher bits down
static void _remove_var_from_anf(SparseANF& g, int var_idx) {
    AnfMap new_terms;
    uint64_t low_mask = (var_idx > 0) ? ((1ULL << var_idx) - 1) : 0;
    for (auto& [mask, v] : g.terms()) {
        if (!v) continue;
        // Bits below var_idx stay, bits above shift down by 1
        uint64_t new_mask = (mask & low_mask) | ((mask >> 1) & ~low_mask);
        new_terms[new_mask] ^= 1;
    }
    AnfMap clean;
    for (auto& [m, v] : new_terms) if (v) clean[m] = 1;
    g = SparseANF(clean, g.n() - 1);
}

// ====================================================================
//  get_signal
// ====================================================================

const SignalInfo* GateBuilder::get_signal(const std::string& name) const {
    auto it = signals_.find(name);
    return (it != signals_.end()) ? &it->second : nullptr;
}

// ====================================================================
//  _resolve: find signal info for a name
// ====================================================================

const SignalInfo* GateBuilder::_resolve(const std::string& name, const Circuit& circ) {
    auto it = signals_.find(name);
    if (it != signals_.end()) return &it->second;
    // Check if it's an input
    for (int i = 0; i < (int)circ.inputs.size(); i++) {
        if (circ.inputs[i] == name) {
            AnfMap terms;
            terms[1ULL << i] = 1;
            SparseANF g(terms, circ.n_inputs);
            std::vector<uint32_t> M;
            for (int j = 0; j < circ.n_inputs; j++)
                M.push_back(1U << j);
            signals_[name] = SignalInfo(g, M, 0, circ.n_inputs);
            return &signals_[name];
        }
    }
    return nullptr;
}

// ====================================================================
//  build
// ====================================================================

void GateBuilder::build(const Circuit& circ) {
    next_m_ = circ.n_inputs;
    int n_stmts = (int)circ.stmts.size();

    for (int si = 0; si < n_stmts; si++) {
        const auto& st = circ.stmts[si];
        const std::string& lhs = st.name;

        switch (st.op) {
            case Op::CONST0: {
                SparseANF g(0);
                signals_[lhs] = SignalInfo(g, {}, 0, 0);
                break;
            }
            case Op::CONST1: {
                AnfMap terms;
                terms[0] = 1;
                SparseANF g(terms, 0);
                signals_[lhs] = SignalInfo(g, {}, 0, 0);
                break;
            }
            case Op::INPUT:
                _process_input(lhs, st.arg1, circ);
                break;
            case Op::NOT:
                _process_not(lhs, st.arg1, circ);
                break;
            case Op::XOR:
                _process_xor(lhs, st.arg1, st.arg2, circ);
                break;
            case Op::AND:
                _process_and(lhs, st.arg1, st.arg2, circ);
                break;
        }

        if ((si + 1) % 200 == 0)
            printf("    processed %d/%d stmts, next_m=%d\n", si + 1, n_stmts, next_m_);
    }
}

// ====================================================================
//  _process_input
// ====================================================================

void GateBuilder::_process_input(const std::string& lhs, const std::string& rhs,
                                  const Circuit& circ) {
    for (int i = 0; i < (int)circ.inputs.size(); i++) {
        if (circ.inputs[i] == rhs) {
            AnfMap terms;
            terms[1ULL << i] = 1;
            SparseANF g(terms, circ.n_inputs);
            std::vector<uint32_t> M;
            for (int j = 0; j < circ.n_inputs; j++)
                M.push_back(1U << j);
            signals_[lhs] = SignalInfo(g, M, 0, circ.n_inputs);
            return;
        }
    }
    auto* info = _resolve(rhs, circ);
    if (info) signals_[lhs] = *info;
    else printf("  WARNING: alias %s = %s unresolved\n", lhs.c_str(), rhs.c_str());
}

// ====================================================================
//  _process_not
// ====================================================================

void GateBuilder::_process_not(const std::string& lhs, const std::string& arg,
                                const Circuit& circ) {
    auto* info = _resolve(arg, circ);
    if (!info) { printf("  WARNING: NOT %s = !%s unresolved\n", lhs.c_str(), arg.c_str()); return; }
    SignalInfo result = *info;
    result.g.xor_const();
    signals_[lhs] = result;
}

// ====================================================================
//  _process_xor
// ====================================================================

void GateBuilder::_process_xor(const std::string& lhs,
                                const std::string& arg1, const std::string& arg2,
                                const Circuit& circ) {
    auto* info1 = _resolve(arg1, circ);
    auto* info2 = _resolve(arg2, circ);
    if (!info1 || !info2) {
        printf("  WARNING: XOR %s = %s + %s unresolved\n", lhs.c_str(), arg1.c_str(), arg2.c_str());
        return;
    }

    bool same = (info1->m == info2->m && info1->b == info2->b &&
                 info1->M_rows.size() == info2->M_rows.size());
    if (same) {
        for (int j = 0; j < (int)info1->M_rows.size(); j++)
            if (info1->M_rows[j] != info2->M_rows[j]) { same = false; break; }
    }

    SignalInfo result;
    if (same) {
        result = *info1;
        result.g.xor_op(info2->g);
    } else {
        int m_A = info1->m;
        result.M_rows = info1->M_rows;
        result.M_rows.insert(result.M_rows.end(), info2->M_rows.begin(), info2->M_rows.end());
        result.b = info1->b | (info2->b << m_A);
        result.m = m_A + info2->m;
        result.g = info1->g;
        SparseANF gB = info2->g;
        gB.shift(m_A);
        result.g.xor_op(gB);
    }

    if (result.g.T() > threshold_ || result.m > 30)
        _compress(result);

    signals_[lhs] = result;
    if (result.m > next_m_) next_m_ = result.m;
}

// ====================================================================
//  _process_and
// ====================================================================

void GateBuilder::_process_and(const std::string& lhs,
                                const std::string& arg1, const std::string& arg2,
                                const Circuit& circ) {
    auto* info1 = _resolve(arg1, circ);
    auto* info2 = _resolve(arg2, circ);
    if (!info1 || !info2) {
        printf("  WARNING: AND %s = %s * %s unresolved\n", lhs.c_str(), arg1.c_str(), arg2.c_str());
        return;
    }

    bool same = (info1->m == info2->m && info1->b == info2->b &&
                 info1->M_rows.size() == info2->M_rows.size());
    if (same) {
        for (int j = 0; j < (int)info1->M_rows.size(); j++)
            if (info1->M_rows[j] != info2->M_rows[j]) { same = false; break; }
    }

    SignalInfo result;
    if (same) {
        result = *info1;
        result.g.and_op(info2->g);
    } else {
        int m_A = info1->m;
        result.M_rows = info1->M_rows;
        result.M_rows.insert(result.M_rows.end(), info2->M_rows.begin(), info2->M_rows.end());
        result.b = info1->b | (info2->b << m_A);
        result.m = m_A + info2->m;
        result.g = info1->g;
        SparseANF gB = info2->g;
        gB.shift(m_A);
        result.g.and_op(gB);
    }

    if (result.g.T() > threshold_ || result.m > 30)
        _compress(result);

    signals_[lhs] = result;
    if (result.m > next_m_) next_m_ = result.m;
}

// ====================================================================
//  _compress: complement → greedy merge
// ====================================================================

void GateBuilder::_compress(SignalInfo& si, bool verbose) {
    int m0 = si.m;
    int T0 = si.g.T();

    if (verbose)
        printf("    Compress: m=%d T=%d\n", m0, T0);

    // Step 1: complement selection
    {
        std::vector<uint32_t> comp_M;
        uint64_t comp_b = 0;
        si.g.complement_search_greedy(comp_M, comp_b);
        if (comp_b) {
            si.b ^= comp_b;
            if (verbose)
                printf("      Complement: Δb=%lx T=%d\n", comp_b, si.g.T_nonlinear());
        }
    }

    // Step 2: greedy merge
    {
        std::vector<uint32_t> merge_rows;
        uint64_t merge_b = 0;
        int old_m = si.m;
        _greedy_merge(si.g, si.m, merge_rows, merge_b);
        if (si.m < old_m) {
            // Compose: M_new = merge_rows @ si.M_rows
            // merge_rows is m_new × old_m, si.M_rows is old_m × n
            int n_vars = (int)si.M_rows.size();
            std::vector<uint32_t> new_M_rows(si.m, 0);
            for (int r = 0; r < si.m; r++) {
                uint32_t row = merge_rows[r];
                uint32_t accum = 0;
                while (row) {
                    int j = __builtin_ctz(row);
                    row &= row - 1;
                    if (j < old_m) accum ^= si.M_rows[j];
                }
                new_M_rows[r] = accum;
            }
            si.M_rows = new_M_rows;

            // b_new[r] = XOR of merge_rows[r][j] & b_old[j], then XOR merge_b[r]
            si.b = _compose_b(merge_rows, si.m, si.b, merge_b);

            if (verbose)
                printf("      Merge: m %d→%d T=%d\n", old_m, si.m, si.g.T_nonlinear());
        }
    }

    if (verbose)
        printf("    Result: m %d→%d T %d→%d\n", m0, si.m, T0, si.g.T());
}

// ====================================================================
//  _greedy_merge
// ====================================================================

void GateBuilder::_greedy_merge(SparseANF& g, int& m,
                                 std::vector<uint32_t>& comp_rows,
                                 uint64_t& comp_b)
{
    comp_rows.clear();
    comp_b = 0;
    for (int j = 0; j < m; j++)
        comp_rows.push_back(1U << j);

    if (m <= 1) return;

    bool improved = true;
    while (improved && m > 1) {
        improved = false;
        int best_i = -1, best_j = -1;
        int best_T = g.T_nonlinear();

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (i == j) continue;
                SparseANF copy = g;
                copy.substitute_zi_eq_zi_xor_zj(i, j);
                int t = copy.T_nonlinear();
                if (t < best_T) {
                    best_T = t;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        if (best_i >= 0) {
            // Apply merge: z_i = z_i ⊕ z_j
            g.substitute_zi_eq_zi_xor_zj(best_i, best_j);

            // Update transformation: row_i ^= row_j
            comp_rows[best_i] ^= comp_rows[best_j];

            // Remove variable best_i
            _remove_var_from_anf(g, best_i);
            comp_rows.erase(comp_rows.begin() + best_i);
            m--;

            improved = true;
        }
    }
}

// ====================================================================
//  _compose_b
// ====================================================================

uint64_t GateBuilder::_compose_b(
    const std::vector<uint32_t>& M_new, int m_new,
    uint64_t b_old, uint64_t b_new)
{
    uint64_t result = 0;
    for (int r = 0; r < m_new; r++) {
        uint32_t row = M_new[r];
        int bit = 0;
        while (row) {
            int j = __builtin_ctz(row);
            row &= row - 1;
            bit ^= ((b_old >> j) & 1);
        }
        bit ^= ((b_new >> r) & 1);
        if (bit) result |= (1ULL << r);
    }
    return result;
}
