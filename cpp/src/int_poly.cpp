#include "int_poly.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>

// ====================================================================
//  IntPoly implementation
// ====================================================================

IntPoly::IntPoly(const std::unordered_map<ExpVector, int64_t, ExpHash>& t, int n)
    : n(n) {
    for (auto& [exp, c] : t) {
        if (c != 0) terms[exp] = c;
    }
}

IntPoly::IntPoly(std::unordered_map<ExpVector, int64_t, ExpHash>&& t, int n)
    : n(n) {
    for (auto& [exp, c] : t) {
        if (c != 0) terms[std::move(exp)] = c;
    }
}

int IntPoly::degree() const {
    int d = 0;
    for (auto& [exp, c] : terms) {
        int s = 0;
        for (int e : exp) s += e;
        if (s > d) d = s;
    }
    return d;
}

std::vector<int> IntPoly::deg_by_var() const {
    std::vector<int> max_deg(n, 0);
    for (auto& [exp, c] : terms) {
        for (int i = 0; i < n && i < (int)exp.size(); ++i) {
            if (exp[i] > max_deg[i]) max_deg[i] = exp[i];
        }
    }
    return max_deg;
}

IntPoly IntPoly::operator+(const IntPoly& other) const {
    assert(n == other.n);
    std::unordered_map<ExpVector, int64_t, ExpHash> result = terms;
    for (auto& [exp, c] : other.terms) {
        auto it = result.find(exp);
        if (it != result.end()) {
            it->second += c;
            if (it->second == 0) result.erase(it);
        } else {
            result[exp] = c;
        }
    }
    return IntPoly(std::move(result), n);
}

IntPoly IntPoly::operator-(const IntPoly& other) const {
    assert(n == other.n);
    std::unordered_map<ExpVector, int64_t, ExpHash> result = terms;
    for (auto& [exp, c] : other.terms) {
        auto it = result.find(exp);
        if (it != result.end()) {
            it->second -= c;
            if (it->second == 0) result.erase(it);
        } else {
            result[exp] = -c;
        }
    }
    return IntPoly(std::move(result), n);
}

IntPoly IntPoly::operator*(const IntPoly& other) const {
    assert(n == other.n);
    std::unordered_map<ExpVector, int64_t, ExpHash> result;
    for (auto& [e1, c1] : terms) {
        for (auto& [e2, c2] : other.terms) {
            ExpVector e(n, 0);
            for (int i = 0; i < n; ++i) e[i] = e1[i] + e2[i];
            auto it = result.find(e);
            if (it != result.end()) {
                it->second += c1 * c2;
                if (it->second == 0) result.erase(it);
            } else {
                result[std::move(e)] = c1 * c2;
            }
        }
    }
    return IntPoly(std::move(result), n);
}

IntPoly IntPoly::operator*(int64_t scalar) const {
    if (scalar == 0) return IntPoly({}, n);
    std::unordered_map<ExpVector, int64_t, ExpHash> result;
    for (auto& [exp, c] : terms) result[exp] = c * scalar;
    return IntPoly(std::move(result), n);
}

IntPoly IntPoly::partial_deriv(int var) const {
    std::unordered_map<ExpVector, int64_t, ExpHash> result;
    for (auto& [exp, c] : terms) {
        int e = var < (int)exp.size() ? exp[var] : 0;
        if (e > 0) {
            ExpVector new_exp = exp;
            new_exp[var] = e - 1;
            result[std::move(new_exp)] = c * e;
        }
    }
    return IntPoly(std::move(result), n);
}

std::vector<IntPoly> IntPoly::gradient() const {
    std::vector<IntPoly> grads(n);
    for (int i = 0; i < n; ++i) grads[i] = partial_deriv(i);
    return grads;
}

std::vector<int> IntPoly::variables_used() const {
    std::vector<bool> used_flag(n, false);
    for (auto& [exp, c] : terms) {
        for (int i = 0; i < n && i < (int)exp.size(); ++i) {
            if (exp[i] > 0) used_flag[i] = true;
        }
    }
    std::vector<int> used;
    for (int i = 0; i < n; ++i) {
        if (used_flag[i]) used.push_back(i);
    }
    return used;
}

IntPoly IntPoly::substitute_linear(int var, const std::vector<int64_t>& coeffs) const {
    assert(var < n);
    assert((int)coeffs.size() == n);

    IntPoly base_poly = linear_form_poly(coeffs, n);
    IntPoly result({}, n);

    for (auto& [exp, c] : terms) {
        int e = exp[var];
        if (e == 0) {
            auto it = result.terms.find(exp);
            if (it != result.terms.end()) {
                it->second += c;
                if (it->second == 0) result.terms.erase(it);
            } else {
                result.terms[exp] = c;
            }
            continue;
        }
        ExpVector rest_exp = exp;
        rest_exp[var] = 0;
        auto rest_poly = IntPoly({{rest_exp, c}}, n);

        IntPoly term_poly = base_poly;
        for (int p = 1; p < e; ++p) term_poly = term_poly * base_poly;

        IntPoly prod = term_poly * rest_poly;
        for (auto& [e2, c2] : prod.terms) {
            auto it = result.terms.find(e2);
            if (it != result.terms.end()) {
                it->second += c2;
                if (it->second == 0) result.terms.erase(it);
            } else {
                result.terms[std::move(e2)] = c2;
            }
        }
    }
    return result;
}

// ---- Forward expansion: x = Nz + c, always works for any m ----
IntPoly IntPoly::expand_affine(const std::vector<std::vector<int64_t>>& N,
                                const std::vector<int64_t>& c) const {
    if (N.empty() || n == 0) return IntPoly({}, N.empty() ? 0 : (int)N[0].size());
    int m = (int)N[0].size();

    // g(z) = f(Nz + c): expand each term Π_i (N[i]·z + c[i])^{e_i}
    std::unordered_map<ExpVector, int64_t, ExpHash> result;
    ExpVector zero_exp(m, 0);

    for (auto& [exp, coeff] : terms) {
        IntPoly cur_poly({}, m);
        cur_poly.terms[zero_exp] = 1;

        for (int i = 0; i < n && i < (int)exp.size(); ++i) {
            int e = exp[i];
            if (e == 0) continue;

            IntPoly affine = single_affine(N[i], i < (int)c.size() ? c[i] : 0, m);
            if (e == 1) {
                cur_poly = cur_poly * affine;
            } else {
                IntPoly pow_poly = affine;
                for (int p = 1; p < e; ++p) pow_poly = pow_poly * affine;
                cur_poly = cur_poly * pow_poly;
            }
        }
        if (coeff != 1) cur_poly = cur_poly * coeff;

        for (auto& [e2, c2] : cur_poly.terms) {
            auto it = result.find(e2);
            if (it != result.end()) {
                it->second += c2;
                if (it->second == 0) result.erase(it);
            } else {
                result[std::move(e2)] = c2;
            }
        }
    }
    return IntPoly(std::move(result), m);
}

// ====================================================================
//  Integer matrix helpers (for substitute_affine inverse direction)
// ====================================================================

// Invert n×n integer matrix using double-precision Gaussian elimination.
// Returns empty vector if matrix is singular or inverse not integer.
static std::vector<std::vector<int64_t>> int_mat_invert(
    const std::vector<std::vector<int64_t>>& M, int n) {
    std::vector<std::vector<double>> aug(n, std::vector<double>(2 * n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            aug[i][j] = (double)M[i][j];
    for (int i = 0; i < n; ++i)
        aug[i][n + i] = 1.0;

    for (int col = 0; col < n; ++col) {
        int pivot = -1;
        for (int r = col; r < n; ++r) {
            if (std::abs(aug[r][col]) > 1e-12) { pivot = r; break; }
        }
        if (pivot < 0) return {};
        std::swap(aug[col], aug[pivot]);

        double piv = aug[col][col];
        for (int j = 0; j < 2 * n; ++j) aug[col][j] /= piv;

        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            double factor = aug[r][col];
            for (int j = 0; j < 2 * n; ++j)
                aug[r][j] -= factor * aug[col][j];
        }
    }

    std::vector<std::vector<int64_t>> inv(n, std::vector<int64_t>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double val = aug[i][n + j];
            int64_t ival = (int64_t)std::llround(val);
            if (std::abs(val - (double)ival) > 1e-6) return {};
            inv[i][j] = ival;
        }
    }
    return inv;
}

// Extend m×n matrix (m < n, full row rank) to n×n unimodular matrix.
static std::vector<std::vector<int64_t>> int_extend_to_invertible(
    const std::vector<std::vector<int64_t>>& M, int m, int n) {
    assert(m < n);
    std::vector<std::vector<int64_t>> rows = M;
    rows.resize(n, std::vector<int64_t>(n, 0));

    // Double-precision elimination to find pivot columns
    std::vector<std::vector<double>> work(m, std::vector<double>(n, 0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            work[i][j] = (double)M[i][j];

    std::vector<bool> pivot_col(n, false);
    int rank = 0;
    for (int col = 0; col < n && rank < m; ++col) {
        int pivot = -1;
        for (int r = rank; r < m; ++r) {
            if (std::abs(work[r][col]) > 1e-12) { pivot = r; break; }
        }
        if (pivot < 0) continue;

        std::swap(work[rank], work[pivot]);
        std::swap(rows[rank], rows[pivot]);
        pivot_col[col] = true;

        double piv = work[rank][col];
        for (int j = col; j < n; ++j) work[rank][j] /= piv;

        for (int r = 0; r < m; ++r) {
            if (r == rank) continue;
            double factor = work[r][col];
            for (int j = col; j < n; ++j) work[r][j] -= factor * work[rank][j];
        }
        ++rank;
    }

    // Add unit rows for non-pivot columns
    int slot = m;
    for (int col = 0; col < n; ++col) {
        if (!pivot_col[col]) {
            rows[slot][col] = 1;
            ++slot;
        }
    }
    return rows;
}

// Extend m×n matrix (m > n, full column rank) to m×m invertible by adding columns.
static std::vector<std::vector<int64_t>> int_extend_columns_to_invertible(
    const std::vector<std::vector<int64_t>>& M, int m, int n) {
    assert(m > n);
    // Gaussian elimination to find pivot rows
    std::vector<std::vector<double>> work(m, std::vector<double>(n, 0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            work[i][j] = (double)M[i][j];

    std::vector<int> perm(m);
    for (int i = 0; i < m; ++i) perm[i] = i;

    int rank = 0;
    for (int col = 0; col < n && rank < m; ++col) {
        int pivot = -1;
        for (int r = rank; r < m; ++r) {
            if (std::abs(work[r][col]) > 1e-12) { pivot = r; break; }
        }
        if (pivot < 0) continue;
        std::swap(work[rank], work[pivot]);
        std::swap(perm[rank], perm[pivot]);
        double piv = work[rank][col];
        for (int j = col; j < n; ++j) work[rank][j] /= piv;
        for (int r = 0; r < m; ++r) {
            if (r != rank && std::abs(work[r][col]) > 1e-12) {
                double factor = work[r][col];
                for (int j = col; j < n; ++j)
                    work[r][j] -= factor * work[rank][j];
            }
        }
        ++rank;
    }

    // Build m×m extended: first n cols = M, extra cols = unit vectors at non-pivot rows
    std::vector<std::vector<int64_t>> result(m, std::vector<int64_t>(m, 0));
    for (int r = 0; r < m; ++r)
        for (int j = 0; j < n; ++j)
            result[r][j] = M[r][j];

    for (int j = 0; j < m - n; ++j) {
        int np_row = perm[n + j];
        result[np_row][n + j] = 1;
    }
    return result;
}

// Structured substitute_affine: each row of M has ≤1 non-zero.
// z_j = a_j * x_i + b_j with a_j = ±1 → x_i = ±(z_j - b_j).
// Converts to expand_affine via x = Nz + c.
IntPoly IntPoly::substitute_affine_structured(
    const std::vector<std::vector<int64_t>>& M,
    const std::vector<int64_t>& b) const {
    int m = (int)M.size();
    // Build N and c such that x_i = Σ N[i][j]*z_j + c[i]
    std::vector<std::vector<int64_t>> N(n, std::vector<int64_t>(m, 0));
    std::vector<int64_t> c_vec(n, 0);
    std::vector<bool> covered(n, false);

    for (int j = 0; j < m; ++j) {
        int nz_col = -1;
        int64_t a = 0;
        for (int i = 0; i < n; ++i) {
            if (M[j][i] != 0) {
                if (nz_col >= 0) return IntPoly({}, m); // should not reach
                nz_col = i;
                a = M[j][i];
            }
        }
        if (nz_col < 0) continue; // all-zero row (constant z)
        if (a != 1 && a != -1) continue; // non-invertible, skip
        int i = nz_col;
        if (covered[i]) continue; // already have a mapping for x_i
        N[i][j] = a;
        c_vec[i] = -a * b[j];
        covered[i] = true;
    }

    // Verify all variables appearing in f are covered
    for (auto& [exp, coeff] : terms) {
        for (int i = 0; i < n && i < (int)exp.size(); ++i) {
            if (exp[i] > 0 && !covered[i]) {
                return IntPoly({}, m);
            }
        }
    }

    return expand_affine(N, c_vec);
}

// ---- substitute_affine: z = Mx + b direction, compute g(z) s.t. f(x) = g(Mx + b) ----
// Only works for m ≤ n (needs right inverse of M).
IntPoly IntPoly::substitute_affine(const std::vector<std::vector<int64_t>>& M,
                                    const std::vector<int64_t>& b) const {
    int m = (int)M.size();
    assert(m > 0 && "substitute_affine requires m > 0");
    int cols = (int)M[0].size();
    assert(cols == n);

    // For m > n, check if M is "structured" (each row has at most one non-zero)
    if (m > n) {
        bool structured = true;
        for (auto& row : M) {
            int nz = 0;
            for (auto& v : row) if (v != 0) ++nz;
            if (nz > 1) { structured = false; break; }
        }
        if (structured) {
            return substitute_affine_structured(M, b);
        }
        // General m > n: extend columns to (m×m), invert, take first n rows as left inverse
        auto M_ext = int_extend_columns_to_invertible(M, m, n);
        auto inv_ext = int_mat_invert(M_ext, m);
        if (inv_ext.empty()) {
            std::cerr << "substitute_affine: failed to compute left inverse (m="
                      << m << ", n=" << n << ")\n";
            return IntPoly({}, m);
        }
        // Left inverse N is first n rows (n×m)
        std::vector<std::vector<int64_t>> N(inv_ext.begin(), inv_ext.begin() + n);
        std::vector<int64_t> c_vec(n, 0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                c_vec[i] -= N[i][j] * b[j];
        return expand_affine(N, c_vec);
    }

    // 1) Find right inverse N (n×m) of M (m×n): M·N = I_m
    std::vector<std::vector<int64_t>> N;
    std::vector<int64_t> c(n, 0);

    if (m == n) {
        // Direct inversion
        auto inv = int_mat_invert(M, n);
        if (inv.empty()) {
            std::cerr << "substitute_affine: M not invertible over Z\n";
            return IntPoly({}, m);
        }
        N = inv;  // n×n, serves as n×m since m=n

        // c = -N·b
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                c[i] -= N[i][j] * b[j];
    } else {
        // m < n: extend to n×n, invert, take first m columns
        auto M_ext = int_extend_to_invertible(M, m, n);
        auto inv_ext = int_mat_invert(M_ext, n);
        if (inv_ext.empty()) {
            std::cerr << "substitute_affine: failed to extend M to invertible\n";
            return IntPoly({}, m);
        }

        // N = first m columns of inv_ext (n×m)
        N.resize(n, std::vector<int64_t>(m, 0));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                N[i][j] = inv_ext[i][j];  // inv_ext[i] has m columns

        // c = -N·b  (b is m, extended b_ext has 0 for extra vars)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                c[i] -= N[i][j] * b[j];
    }

    // 2) g(z) = f(Nz + c) via forward expansion
    return expand_affine(N, c);
}

// Z = Z1 ∪ Z2: Z1 = X (n vars), Z2 = Mx+b (m vars)
// Uses left inverse path directly to avoid structured shortcut.
IntPoly IntPoly::substitute_affine_union(
    const std::vector<std::vector<int64_t>>& M,
    const std::vector<int64_t>& b) const {
    int m = (int)M.size();
    int total = m + n;
    // Build M_ext = [M; I_n] ((m+n)×n), b_ext = [b; 0]
    std::vector<std::vector<int64_t>> M_ext = M;
    for (int i = 0; i < n; ++i) {
        std::vector<int64_t> row(n, 0);
        row[i] = 1;
        M_ext.push_back(row);
    }
    std::vector<int64_t> b_ext = b;
    b_ext.resize(total, 0);
    // Left inverse: N·M_ext = I_n via column extension + inversion
    auto M_sq = int_extend_columns_to_invertible(M_ext, total, n);
    auto inv = int_mat_invert(M_sq, total);
    if (inv.empty()) {
        std::cerr << "substitute_affine_union: left inverse failed\n";
        return IntPoly({}, total);
    }
    std::vector<std::vector<int64_t>> N(inv.begin(), inv.begin() + n);
    std::vector<int64_t> c(n, 0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < total; ++j)
            c[i] -= N[i][j] * b_ext[j];
    return expand_affine(N, c);
}

IntPoly IntPoly::linear_form_poly(const std::vector<int64_t>& coeffs, int n) {
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    for (int j = 0; j < (int)coeffs.size() && j < n; ++j) {
        if (coeffs[j] != 0) {
            ExpVector exp(n, 0);
            exp[j] = 1;
            terms[std::move(exp)] = coeffs[j];
        }
    }
    return IntPoly(std::move(terms), n);
}

IntPoly IntPoly::single_affine(const std::vector<int64_t>& row,
                                int64_t const_val, int m) {
    std::unordered_map<ExpVector, int64_t, ExpHash> terms;
    if (const_val != 0) terms[ExpVector(m, 0)] = const_val;
    for (int j = 0; j < (int)row.size() && j < m; ++j) {
        if (row[j] != 0) {
            ExpVector exp(m, 0);
            exp[j] = 1;
            terms[std::move(exp)] = row[j];
        }
    }
    return IntPoly(std::move(terms), m);
}

int64_t IntPoly::eval(const std::vector<int64_t>& values) const {
    int64_t result = 0;
    for (auto& [exp, c] : terms) {
        int64_t term_val = c;
        for (int i = 0; i < n && i < (int)exp.size(); ++i) {
            for (int p = 0; p < exp[i]; ++p)
                term_val *= i < (int)values.size() ? values[i] : 0;
        }
        result += term_val;
    }
    return result;
}

std::string IntPoly::to_string() const {
    std::ostringstream os;
    os << "IntPoly(n=" << n << ", T=" << T() << "): ";
    bool first = true;
    for (auto& [exp, c] : terms) {
        if (!first) os << " + ";
        first = false;
        if (c != 1 || exp.empty()) os << c;
        for (int i = 0; i < (int)exp.size(); ++i) {
            if (exp[i] > 0) {
                os << "x" << i;
                if (exp[i] > 1) os << "^" << exp[i];
            }
        }
    }
    return os.str();
}

void IntPoly::print() const { std::cout << to_string() << "\n"; }

// ====================================================================
//  Simplification helpers
// ====================================================================

static std::vector<std::vector<int64_t>> identity_int(int n) {
    std::vector<std::vector<int64_t>> I(n, std::vector<int64_t>(n, 0));
    for (int i = 0; i < n; ++i) I[i][i] = 1;
    return I;
}

static IntSimplifyResult drop_unused_int(const IntPoly& f,
                                          const std::vector<std::vector<int64_t>>& M,
                                          const std::vector<int64_t>& b) {
    auto used = f.variables_used();
    if ((int)used.size() == f.n) return {f, M, b};

    std::unordered_map<ExpVector, int64_t, ExpHash> new_terms;
    for (auto& [exp, c] : f.terms) {
        ExpVector new_exp;
        for (int idx : used) new_exp.push_back(idx < (int)exp.size() ? exp[idx] : 0);
        new_terms[std::move(new_exp)] = c;
    }

    std::vector<std::vector<int64_t>> new_M;
    std::vector<int64_t> new_b;
    for (int idx : used) {
        if (idx < (int)M.size()) new_M.push_back(M[idx]);
        new_b.push_back(idx < (int)b.size() ? b[idx] : 0);
    }
    return {IntPoly(std::move(new_terms), (int)used.size()), new_M, new_b};
}

// ====================================================================
//  Greedy merge
// ====================================================================

IntPoly try_merge_int(const IntPoly& f, int i, int j, int64_t k) {
    std::vector<int64_t> coeffs(f.n, 0);
    coeffs[i] = 1;
    coeffs[j] = -k;
    return f.substitute_linear(i, coeffs);
}

IntSimplifyResult greedy_merge_simplify_int(const IntPoly& f,
                                             int max_iter, bool verbose) {
    if (f.terms.empty())
        return {IntPoly({}, f.n), identity_int(f.n), std::vector<int64_t>(f.n, 0)};

    IntPoly cur = f;
    auto M = identity_int(f.n);
    std::vector<int64_t> b(f.n, 0);
    int orig_T = f.T();

    if (verbose)
        std::cout << "\nGreedy merge (int): n=" << f.n << ", T₀=" << orig_T << "\n";

    std::vector<int64_t> k_values = {1, -1, 2, -2};

    for (int iter = 0; iter < max_iter; ++iter) {
        auto active = cur.variables_used();
        if ((int)active.size() <= 1) break;

        int best_delta = 0;
        int best_i = -1, best_j = -1;
        int64_t best_k = 0;

        for (int i : active) {
            for (int j : active) {
                if (i == j) continue;
                for (int64_t k : k_values) {
                    auto g_test = try_merge_int(cur, i, j, k);
                    int delta = g_test.T() - cur.T();
                    if (delta < best_delta) {
                        best_delta = delta;
                        best_i = i; best_j = j; best_k = k;
                    }
                }
            }
        }
        if (best_delta >= 0) break;

        for (int col = 0; col < f.n; ++col)
            M[best_i][col] += best_k * M[best_j][col];
        cur = try_merge_int(cur, best_i, best_j, best_k);
        b[best_i] += best_k * b[best_j];

        auto result = drop_unused_int(cur, M, b);
        cur = result.g; M = result.M; b = result.b;

        if (verbose) {
            double pct = (orig_T - cur.T()) * 100.0 / orig_T;
            std::cout << "  iter " << iter << ": x" << best_i << "→x" << best_i
                      << "+" << best_k << "x" << best_j
                      << "  T=" << cur.T() << "/" << orig_T
                      << " (" << pct << "%↓) m=" << cur.n << "\n";
        }
        if (cur.T() <= 1) break;
    }
    if (verbose)
        std::cout << "  Final: T=" << cur.T() << "/" << orig_T << "\n";
    return {cur, M, b};
}

// ====================================================================
//  Gradient-guided simplification
// ====================================================================

static bool gradient_match_int(const IntPoly& g_i, const IntPoly& g_j, int64_t k) {
    if (g_i.T() != g_j.T()) return false;
    for (auto& [exp, c] : g_i.terms) {
        auto it = g_j.terms.find(exp);
        if (it == g_j.terms.end() || it->second != k * c) return false;
    }
    return true;
}

IntSimplifyResult simplify_by_gradient_int(const IntPoly& f, bool verbose) {
    if (f.terms.empty())
        return {IntPoly({}, f.n), identity_int(f.n), std::vector<int64_t>(f.n, 0)};

    IntPoly cur = f;
    auto M = identity_int(f.n);
    std::vector<int64_t> b(f.n, 0);
    int orig_T = f.T();

    if (verbose)
        std::cout << "\nGradient-guided merge: n=" << f.n << ", T₀=" << orig_T << "\n";

    std::vector<int64_t> k_values = {1, -1, 2, -2};

    for (int iter = 0; iter < 30; ++iter) {
        auto grads = cur.gradient();
        auto active = cur.variables_used();
        if ((int)active.size() <= 1) break;

        int best_delta = 0;
        int best_i = -1, best_j = -1;
        int64_t best_k = 0;

        // Phase A: gradient match
        for (int i : active) {
            if (grads[i].T() == 0) continue;
            for (int j : active) {
                if (i == j || grads[j].T() == 0) continue;
                for (int64_t k : k_values) {
                    if (gradient_match_int(grads[i], grads[j], k)) {
                        auto g_test = try_merge_int(cur, i, j, k);
                        int delta = g_test.T() - cur.T();
                        if (delta < best_delta) {
                            best_delta = delta;
                            best_i = i; best_j = j; best_k = k;
                        }
                    }
                }
            }
        }

        // Phase B: exhaustive
        if (best_i < 0) {
            for (int i : active) {
                for (int j : active) {
                    if (i == j) continue;
                    for (int64_t k : k_values) {
                        auto g_test = try_merge_int(cur, i, j, k);
                        int delta = g_test.T() - cur.T();
                        if (delta < best_delta) {
                            best_delta = delta;
                            best_i = i; best_j = j; best_k = k;
                        }
                    }
                }
            }
        }
        if (best_i < 0 || best_delta >= 0) break;

        for (int col = 0; col < f.n; ++col)
            M[best_i][col] += best_k * M[best_j][col];
        cur = try_merge_int(cur, best_i, best_j, best_k);
        b[best_i] += best_k * b[best_j];

        auto result = drop_unused_int(cur, M, b);
        cur = result.g; M = result.M; b = result.b;

        if (verbose) {
            double pct = (orig_T - cur.T()) * 100.0 / orig_T;
            std::cout << "  iter " << iter << ": x" << best_i << "→x" << best_i
                      << "+" << best_k << "·x" << best_j
                      << "  T=" << cur.T() << "/" << orig_T
                      << " (" << pct << "%↓)  m=" << cur.n << "\n";
        }
        if (cur.T() <= 1) break;
    }
    if (verbose)
        std::cout << "  Final: T=" << cur.T() << "/" << orig_T << ", m=" << cur.n << "\n";
    return {cur, M, b};
}

// ====================================================================
//  Random M,b search
// ====================================================================

// Check if m×n integer matrix has full row rank (m ≤ n) over Q
static bool int_has_full_row_rank(const std::vector<std::vector<int64_t>>& M, int m, int n) {
    if (m > n || m == 0) return false;
    std::vector<std::vector<double>> work(m, std::vector<double>(n, 0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            work[i][j] = (double)M[i][j];

    int rank = 0;
    for (int col = 0; col < n && rank < m; ++col) {
        int pivot = -1;
        for (int r = rank; r < m; ++r) {
            if (std::abs(work[r][col]) > 1e-12) { pivot = r; break; }
        }
        if (pivot < 0) continue;
        std::swap(work[rank], work[pivot]);
        double piv = work[rank][col];
        for (int j = col; j < n; ++j) work[rank][j] /= piv;
        for (int r = 0; r < m; ++r) {
            if (r == rank) continue;
            double factor = work[r][col];
            for (int j = col; j < n; ++j) work[r][j] -= factor * work[rank][j];
        }
        ++rank;
    }
    return rank == m;
}

IntSimplifyResult search_random_int(const IntPoly& f, int max_m, int n_trials,
                                     uint64_t seed, bool verbose) {
    int n = f.n;
    int best_T = f.T();
    IntSimplifyResult best = {f, identity_int(n), std::vector<int64_t>(n, 0)};

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> entry_dist(-2, 2);
    std::uniform_int_distribution<int> m_dist(1, max_m);

    for (int trial = 0; trial < n_trials; ++trial) {
        int m = m_dist(rng);

        std::vector<std::vector<int64_t>> M(m, std::vector<int64_t>(n, 0));
        std::vector<int64_t> b(m, 0);

        if (rng() % 3 < 2) {
            // Structured: each row has at most one non-zero entry (±1)
            // Full row rank guaranteed if rows are distinct non-zero.
            std::vector<bool> used_var(n, false);
            bool ok = true;
            for (int row = 0; row < m; ++row) {
                int var = rng() % n;
                if (used_var[var]) { ok = false; break; }
                used_var[var] = true;
                M[row][var] = (rng() % 2) ? 1LL : -1LL;
                if (rng() % 2) b[row] = entry_dist(rng);
            }
            if (!ok) continue;
        } else {
            // General small-entry matrix
            for (int row = 0; row < m; ++row) {
                for (int col = 0; col < n; ++col)
                    M[row][col] = entry_dist(rng);
                b[row] = entry_dist(rng);
            }
            // Need full row rank for right inverse
            if (!int_has_full_row_rank(M, m, n)) continue;
        }

        auto g = f.substitute_affine(M, b);
        if (g.terms.empty()) continue;

        // Try union
        auto g_best = g;
        auto M_best = M;
        auto b_best = b;
        bool use_union = false;

        auto g_u = f.substitute_affine_union(M, b);
        if (g_u.T() < g.T()) {
            for (int i = 0; i < n; ++i) {
                std::vector<int64_t> row(n, 0);
                row[i] = 1;
                M_best.push_back(row);
            }
            b_best.resize(m + n, 0);
            g_best = g_u;
            use_union = true;
        }

        if (g_best.T() < best_T) {
            if ((int)M_best.size() != f.n &&
                verify_int_poly(f, g_best, M_best, b_best, 20)) {
                continue;
            }
            best_T = g_best.T();
            best = {g_best, M_best, b_best};
            if (verbose) {
                double pct = (f.T() - best_T) * 100.0 / f.T();
                std::cout << "  trial " << trial << ": T=" << best_T
                          << "/" << f.T() << " (" << pct << "%↓)"
                          << " m=" << (int)M_best.size()
                          << (use_union ? " [union]" : "") << "\n";
            }
        }
    }
    return best;
}

// ====================================================================
//  Combined pipeline
// ====================================================================

static void compose_matrices(IntSimplifyResult& acc, const IntSimplifyResult& next,
                              int orig_n) {
    int m1 = (int)acc.M.size();
    int m2 = (int)next.M.size();
    std::vector<std::vector<int64_t>> M_new(m2, std::vector<int64_t>(orig_n, 0));
    for (int j = 0; j < m2; ++j) {
        for (int k = 0; k < m1; ++k) {
            if (next.M[j][k] == 0) continue;
            for (int i = 0; i < orig_n; ++i)
                M_new[j][i] += next.M[j][k] * acc.M[k][i];
        }
    }
    std::vector<int64_t> b_new(m2, 0);
    for (int j = 0; j < m2; ++j) {
        for (int k = 0; k < m1; ++k)
            b_new[j] += next.M[j][k] * acc.b[k];
        b_new[j] += next.b[j];
    }
    acc.g = next.g;
    acc.M = M_new;
    acc.b = b_new;
}

IntSimplifyResult simplify_int(const IntPoly& f, bool verbose) {
    int orig_T = f.T();
    int orig_n = f.n;
    if (orig_T <= 1)
        return {f, identity_int(orig_n), std::vector<int64_t>(orig_n, 0)};

    if (verbose)
        std::cout << "\n=== IntPoly simplify: n=" << orig_n << ", T₀=" << orig_T << " ===\n";

    // Phase 1: gradient-guided merge
    auto r1 = simplify_by_gradient_int(f, verbose);
    IntSimplifyResult acc = {r1.g, r1.M, r1.b};

    // Phase 2: exhaustive greedy merge on remaining
    auto r2 = greedy_merge_simplify_int(acc.g, 20, verbose);
    compose_matrices(acc, r2, orig_n);

    if (verbose) {
        double pct = (orig_T - acc.g.T()) * 100.0 / orig_T;
        std::cout << "  After merge: T=" << acc.g.T() << "/" << orig_T
                  << " (" << pct << "%↓), m=" << acc.g.n << "\n";
    }

    // Phase 3: random search with union
    if (acc.g.n > 0 && acc.g.T() > 1) {
        int max_random_m = std::max(acc.g.n, 1);
        auto r3 = search_random_int(acc.g, max_random_m, 100, 2, false);
        if (r3.g.T() < acc.g.T()) {
            compose_matrices(acc, r3, orig_n);
            if (verbose) {
                double pct = (orig_T - acc.g.T()) * 100.0 / orig_T;
                std::cout << "  After random: T=" << acc.g.T() << "/" << orig_T
                          << " (" << pct << "%↓), m=" << acc.g.n << "\n";
            }
        }
    }

    if (verbose) {
        double pct = (orig_T - acc.g.T()) * 100.0 / orig_T;
        std::cout << "  Final: T=" << acc.g.T() << "/" << orig_T
                  << " (" << pct << "%↓), m=" << acc.g.n << "\n";
    }
    return acc;
}

// ====================================================================
//  Verification
// ====================================================================

int verify_int_poly(const IntPoly& f, const IntPoly& g,
                     const std::vector<std::vector<int64_t>>& M,
                     const std::vector<int64_t>& b, int n_tests) {
    int m = (int)M.size();
    std::mt19937_64 rng(12345);
    std::uniform_int_distribution<int64_t> dist(-5, 5);
    int errors = 0;

    for (int t = 0; t < n_tests; ++t) {
        std::vector<int64_t> x(f.n);
        for (int i = 0; i < f.n; ++i) x[i] = dist(rng);

        std::vector<int64_t> z(m, 0);
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < f.n; ++i)
                z[j] += M[j][i] * x[i];
            z[j] += j < (int)b.size() ? b[j] : 0;
        }

        if (f.eval(x) != g.eval(z)) ++errors;
    }
    return errors;
}
