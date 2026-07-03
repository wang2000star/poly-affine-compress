"""
IntPoly — Sparse multivariate polynomial over Z or F_p.

Core data: terms: {exp_vector: coefficient}, n variables, mod = modulus.
mod = 0 → polynomial over Z.
mod > 0 → polynomial over F_mod (coefficient arithmetic modulo mod).

API mirrors the C++ IntPoly in cpp/src/int_poly.cpp.
"""

from __future__ import annotations

import random
from typing import Optional

# Exponent vector type
ExpVector = tuple[int, ...]


def _default_exp(n: int) -> ExpVector:
    return tuple([0] * n)


class IntPoly:
    """Sparse multivariate polynomial over Z or F_p."""

    def __init__(
        self,
        terms: Optional[dict[ExpVector, int]] = None,
        n: int = 0,
        mod: int = 0,
    ):
        self.n = n
        self.mod = mod
        self.terms: dict[ExpVector, int] = {}
        if terms is not None:
            for exp, c in terms.items():
                c = self.reduce(c)
                if not self.coeff_zero(c):
                    self.terms[exp] = c

    # ---- mod-aware helpers ----

    def reduce(self, x: int) -> int:
        if self.mod == 0:
            return x
        r = x % self.mod
        return r if r >= 0 else r + self.mod

    @staticmethod
    def mul_mod(a: int, b: int, mod: int) -> int:
        if mod == 0:
            return a * b
        return (a * b) % mod

    def coeff_zero(self, c: int) -> bool:
        return c == 0 if self.mod == 0 else c % self.mod == 0

    # ---- properties ----

    def T(self) -> int:
        return len(self.terms)

    def degree(self) -> int:
        d = 0
        for exp in self.terms:
            s = sum(exp)
            if s > d:
                d = s
        return d

    def variables_used(self) -> list[int]:
        used = [False] * self.n
        for exp in self.terms:
            for i, e in enumerate(exp):
                if e > 0:
                    used[i] = True
        return [i for i, v in enumerate(used) if v]

    def copy(self) -> "IntPoly":
        return IntPoly(dict(self.terms), self.n, self.mod)

    # ---- arithmetic ----

    def __add__(self, other: "IntPoly") -> "IntPoly":
        assert self.n == other.n
        result = dict(self.terms)
        for exp, c in other.terms.items():
            if exp in result:
                result[exp] = self.reduce(result[exp] + c)
                if self.coeff_zero(result[exp]):
                    del result[exp]
            else:
                result[exp] = self.reduce(c)
        return IntPoly(result, self.n, self.mod)

    def __sub__(self, other: "IntPoly") -> "IntPoly":
        assert self.n == other.n
        result = dict(self.terms)
        for exp, c in other.terms.items():
            if exp in result:
                result[exp] = self.reduce(result[exp] - c)
                if self.coeff_zero(result[exp]):
                    del result[exp]
            else:
                result[exp] = self.reduce(-c)
        return IntPoly(result, self.n, self.mod)

    def __mul__(self, other: "IntPoly | int") -> "IntPoly":
        if isinstance(other, int):
            return self._scalar_mul(other)
        assert self.n == other.n
        result: dict[ExpVector, int] = {}
        for e1, c1 in self.terms.items():
            for e2, c2 in other.terms.items():
                e = tuple(e1[i] + e2[i] for i in range(self.n))
                prod = self.mul_mod(c1, c2, self.mod)
                if e in result:
                    result[e] = self.reduce(result[e] + prod)
                    if self.coeff_zero(result[e]):
                        del result[e]
                else:
                    result[e] = self.reduce(prod)
        return IntPoly(result, self.n, self.mod)

    def __rmul__(self, other: int) -> "IntPoly":
        return self._scalar_mul(other)

    def _scalar_mul(self, scalar: int) -> "IntPoly":
        if self.coeff_zero(scalar):
            return IntPoly({}, self.n, self.mod)
        result: dict[ExpVector, int] = {}
        for exp, c in self.terms.items():
            val = self.mul_mod(c, scalar, self.mod)
            if not self.coeff_zero(val):
                result[exp] = val
        return IntPoly(result, self.n, self.mod)

    # ---- calculus ----

    def partial_deriv(self, var: int) -> "IntPoly":
        result: dict[ExpVector, int] = {}
        for exp, c in self.terms.items():
            e = exp[var] if var < len(exp) else 0
            if e > 0:
                new_exp = list(exp)
                new_exp[var] = e - 1
                val = self.mul_mod(c, e, self.mod)
                if not self.coeff_zero(val):
                    result[tuple(new_exp)] = self.reduce(val)
        return IntPoly(result, self.n, self.mod)

    def gradient(self) -> list["IntPoly"]:
        return [self.partial_deriv(i) for i in range(self.n)]

    # ---- substitution ----

    def substitute_linear(self, var: int, coeffs: list[int]) -> "IntPoly":
        assert var < self.n
        assert len(coeffs) == self.n
        base_poly = _linear_form_poly(coeffs, self.n)
        base_poly.mod = self.mod
        result = IntPoly({}, self.n, self.mod)

        for exp, c in self.terms.items():
            e = exp[var]
            if e == 0:
                if exp in result.terms:
                    result.terms[exp] = self.reduce(result.terms[exp] + c)
                    if self.coeff_zero(result.terms[exp]):
                        del result.terms[exp]
                else:
                    result.terms[exp] = self.reduce(c)
                continue

            rest_exp = list(exp)
            rest_exp[var] = 0
            rest_poly = IntPoly({tuple(rest_exp): c}, self.n, self.mod)

            term_poly = base_poly
            for _ in range(1, e):
                term_poly = term_poly * base_poly

            prod = term_poly * rest_poly
            for e2, c2 in prod.terms.items():
                if e2 in result.terms:
                    result.terms[e2] = self.reduce(result.terms[e2] + c2)
                    if self.coeff_zero(result.terms[e2]):
                        del result.terms[e2]
                else:
                    result.terms[e2] = self.reduce(c2)

        return result

    def expand_affine(
        self, N: list[list[int]], c: list[int]
    ) -> "IntPoly":
        """g(z) = f(Nz + c), N is n×m, c is n."""
        if not N or self.n == 0:
            m = len(N[0]) if N else 0
            return IntPoly({}, m, self.mod)
        m = len(N[0])

        result: dict[ExpVector, int] = {}
        zero_exp = _default_exp(m)

        for exp, coeff in self.terms.items():
            cur_terms = {zero_exp: 1}

            for i in range(min(len(exp), self.n)):
                e = exp[i]
                if e == 0:
                    continue

                affine = _single_affine(N[i], c[i] if i < len(c) else 0, m)
                affine.mod = self.mod
                if e == 1:
                    cur_terms = _intpoly_mul_dict(cur_terms, affine.terms, m, self.mod, self)
                else:
                    pow_poly = affine
                    for _ in range(1, e):
                        pow_poly = pow_poly * affine
                    cur_terms = _intpoly_mul_dict(cur_terms, pow_poly.terms, m, self.mod, self)

            if not self.coeff_zero(coeff) and coeff != 1:
                # scalar multiply cur_terms
                scaled = {}
                for e2, c2 in cur_terms.items():
                    val = self.mul_mod(c2, coeff, self.mod)
                    if not self.coeff_zero(val):
                        scaled[e2] = self.reduce(val)
                cur_terms = scaled

            for e2, c2 in cur_terms.items():
                if e2 in result:
                    result[e2] = self.reduce(result[e2] + c2)
                    if self.coeff_zero(result[e2]):
                        del result[e2]
                else:
                    result[e2] = self.reduce(c2)

        return IntPoly(result, m, self.mod)

    def substitute_affine(
        self, M: list[list[int]], b: list[int]
    ) -> "IntPoly":
        """z = Mx + b direction, compute g(z) s.t. f(x) = g(Mx + b)."""
        m = len(M)
        assert m > 0
        n = self.n

        # Check structured m > n
        if m > n and all(sum(1 for v in row if v != 0) <= 1 for row in M):
            return self._substitute_affine_structured(M, b)

        if m > n:
            # General: left inverse via column extension
            if self.mod > 0:
                M_sq = fp_extend_columns_to_invertible(M, m, n, self.mod)
                inv_ext = fp_mat_invert(M_sq, m, self.mod)
            else:
                M_sq = int_extend_columns_to_invertible(M, m, n)
                inv_ext = int_mat_invert(M_sq, m)
            if not inv_ext:
                return IntPoly({}, m, self.mod)
            N = inv_ext[:n]
            c_vec = [0] * n
            for i in range(n):
                for j in range(m):
                    c_vec[i] -= N[i][j] * b[j]
                c_vec[i] = self.reduce(c_vec[i])
            return self.expand_affine(N, c_vec)

        if m == n:
            if self.mod > 0:
                inv = fp_mat_invert(M, n, self.mod)
            else:
                inv = int_mat_invert(M, n)
            if not inv:
                return IntPoly({}, m, self.mod)
            N = inv
            c_vec = [0] * n
            for i in range(n):
                for j in range(m):
                    c_vec[i] -= N[i][j] * b[j]
                c_vec[i] = self.reduce(c_vec[i])
            return self.expand_affine(N, c_vec)

        # m < n: extend to n×n
        if self.mod > 0:
            M_ext = fp_extend_to_invertible(M, m, n, self.mod)
            inv_ext = fp_mat_invert(M_ext, n, self.mod)
        else:
            M_ext = int_extend_to_invertible(M, m, n)
            inv_ext = int_mat_invert(M_ext, n)
        if not inv_ext:
            return IntPoly({}, m, self.mod)
        N = [[inv_ext[i][j] for j in range(m)] for i in range(n)]
        c_vec = [0] * n
        for i in range(n):
            for j in range(m):
                c_vec[i] -= N[i][j] * b[j]
            c_vec[i] = self.reduce(c_vec[i])
        return self.expand_affine(N, c_vec)

    def _substitute_affine_structured(
        self, M: list[list[int]], b: list[int]
    ) -> "IntPoly":
        m = len(M)
        n = self.n
        N = [[0] * m for _ in range(n)]
        c_vec = [0] * n
        covered = [False] * n

        for j in range(m):
            nz_col = -1
            a = 0
            for i in range(n):
                if M[j][i] != 0:
                    if nz_col >= 0:
                        return IntPoly({}, m, self.mod)
                    nz_col = i
                    a = M[j][i]
            if nz_col < 0:
                continue
            if a not in (1, -1):
                continue
            i = nz_col
            if covered[i]:
                continue
            N[i][j] = a
            c_vec[i] = self.reduce(-a * b[j])
            covered[i] = True

        for exp in self.terms:
            for i in range(min(len(exp), n)):
                if exp[i] > 0 and not covered[i]:
                    return IntPoly({}, m, self.mod)

        return self.expand_affine(N, c_vec)

    def substitute_affine_union(
        self, M: list[list[int]], b: list[int]
    ) -> "IntPoly":
        """Z = Z1 ∪ Z2 where Z1 = X (n vars), Z2 = Mx+b (m vars)."""
        m = len(M)
        total = m + self.n
        M_ext = list(M)
        for i in range(self.n):
            row = [0] * self.n
            row[i] = 1
            M_ext.append(row)
        b_ext = list(b) + [0] * self.n

        if self.mod > 0:
            M_sq = fp_extend_columns_to_invertible(M_ext, total, self.n, self.mod)
            inv = fp_mat_invert(M_sq, total, self.mod)
        else:
            M_sq = int_extend_columns_to_invertible(M_ext, total, self.n)
            inv = int_mat_invert(M_sq, total)
        if not inv:
            return IntPoly({}, total, self.mod)
        N = inv[:self.n]
        c_vec = [0] * self.n
        for i in range(self.n):
            for j in range(total):
                c_vec[i] -= N[i][j] * b_ext[j]
            c_vec[i] = self.reduce(c_vec[i])
        return self.expand_affine(N, c_vec)

    # ---- evaluation ----

    def eval(self, values: list[int]) -> int:
        result = 0
        for exp, c in self.terms.items():
            term_val = self.reduce(c)
            for i, e in enumerate(exp):
                if e == 0:
                    continue
                v = values[i] if i < len(values) else 0
                for _ in range(e):
                    term_val = self.mul_mod(term_val, v, self.mod)
            result = self.reduce(result + term_val)
        return result

    def __repr__(self) -> str:
        return f"IntPoly(n={self.n}, T={self.T()}, mod={self.mod})"


# ====================================================================
#  Internal helpers
# ====================================================================


def _linear_form_poly(coeffs: list[int], n: int) -> IntPoly:
    terms: dict[ExpVector, int] = {}
    for j, c in enumerate(coeffs[:n]):
        if c != 0:
            exp = tuple([1 if i == j else 0 for i in range(n)])
            terms[exp] = c
    return IntPoly(terms, n)


def _single_affine(row: list[int], const_val: int, m: int) -> IntPoly:
    terms: dict[ExpVector, int] = {}
    if const_val != 0:
        terms[tuple([0] * m)] = const_val
    for j, v in enumerate(row[:m]):
        if v != 0:
            exp = tuple([1 if k == j else 0 for k in range(m)])
            terms[exp] = v
    return IntPoly(terms, m, 0)


def _intpoly_mul_dict(
    d1: dict[ExpVector, int],
    d2: dict[ExpVector, int],
    m: int,
    mod: int,
    poly: IntPoly,
) -> dict[ExpVector, int]:
    result: dict[ExpVector, int] = {}
    for e1, c1 in d1.items():
        for e2, c2 in d2.items():
            e = tuple(e1[k] + e2[k] for k in range(m))
            prod = IntPoly.mul_mod(c1, c2, mod)
            if e in result:
                result[e] = poly.reduce(result[e] + prod)
                if poly.coeff_zero(result[e]):
                    del result[e]
            else:
                result[e] = poly.reduce(prod)
    return result


# ====================================================================
#  Integer matrix helpers (double-precision Gaussian elimination)
# ====================================================================


def int_mat_invert(M: list[list[int]], n: int) -> Optional[list[list[int]]]:
    """Invert n×n integer matrix using double-precision Gaussian elimination."""
    import math

    aug = [[0.0] * (2 * n) for _ in range(n)]
    for i in range(n):
        for j in range(n):
            aug[i][j] = float(M[i][j])
        aug[i][n + i] = 1.0

    for col in range(n):
        pivot = -1
        for r in range(col, n):
            if abs(aug[r][col]) > 1e-12:
                pivot = r
                break
        if pivot < 0:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]

        piv = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= piv

        for r in range(n):
            if r == col:
                continue
            factor = aug[r][col]
            for j in range(2 * n):
                aug[r][j] -= factor * aug[col][j]

    inv = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            val = aug[i][n + j]
            ival = round(val)
            if abs(val - ival) > 1e-6:
                return None
            inv[i][j] = ival
    return inv


def int_extend_to_invertible(
    M: list[list[int]], m: int, n: int
) -> list[list[int]]:
    assert m < n
    rows = [list(row) for row in M]
    rows.extend([[0] * n for _ in range(n - m)])

    # Double-precision elimination
    work = [[float(x) for x in row] for row in M]
    pivot_col = [False] * n
    rank = 0
    for col in range(n):
        if rank >= m:
            break
        pivot = -1
        for r in range(rank, m):
            if abs(work[r][col]) > 1e-12:
                pivot = r
                break
        if pivot < 0:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_col[col] = True

        piv = work[rank][col]
        for j in range(col, n):
            work[rank][j] /= piv

        for r in range(m):
            if r == rank:
                continue
            factor = work[r][col]
            for j in range(col, n):
                work[r][j] -= factor * work[rank][j]
        rank += 1

    slot = m
    for col in range(n):
        if not pivot_col[col]:
            rows[slot][col] = 1
            slot += 1
    return rows


def int_extend_columns_to_invertible(
    M: list[list[int]], m: int, n: int
) -> list[list[int]]:
    assert m > n
    work = [[float(x) for x in row] for row in M]
    perm = list(range(m))

    rank = 0
    for col in range(n):
        pivot = -1
        for r in range(rank, m):
            if abs(work[r][col]) > 1e-12:
                pivot = r
                break
        if pivot < 0:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        perm[rank], perm[pivot] = perm[pivot], perm[rank]
        piv = work[rank][col]
        for j in range(col, n):
            work[rank][j] /= piv
        for r in range(m):
            if r != rank and abs(work[r][col]) > 1e-12:
                factor = work[r][col]
                for j in range(col, n):
                    work[r][j] -= factor * work[rank][j]
        rank += 1

    result = [[0] * m for _ in range(m)]
    for r in range(m):
        for c in range(n):
            result[r][c] = M[r][c]
    for j in range(m - n):
        np_row = perm[n + j]
        result[np_row][n + j] = 1
    return result


def int_has_full_row_rank(M: list[list[int]], m: int, n: int) -> bool:
    if m > n or m == 0:
        return False
    work = [[float(x) for x in row] for row in M]
    rank = 0
    for col in range(n):
        if rank >= m:
            break
        pivot = -1
        for r in range(rank, m):
            if abs(work[r][col]) > 1e-12:
                pivot = r
                break
        if pivot < 0:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        piv = work[rank][col]
        for j in range(col, n):
            work[rank][j] /= piv
        for r in range(m):
            if r == rank:
                continue
            factor = work[r][col]
            for j in range(col, n):
                work[r][j] -= factor * work[rank][j]
        rank += 1
    return rank == m


# ====================================================================
#  F_p matrix helpers (exact modular arithmetic)
# ====================================================================


def fp_inv_mod(a: int, p: int) -> int:
    """Modular inverse via extended Euclidean algorithm."""
    a %= p
    if a < 0:
        a += p
    t, newt = 0, 1
    r, newr = p, a
    while newr != 0:
        q = r // newr
        t, newt = newt, t - q * newt
        r, newr = newr, r - q * newr
    if r > 1:
        return 0
    if t < 0:
        t += p
    return t


def fp_mat_invert(
    M: list[list[int]], n: int, p: int
) -> Optional[list[list[int]]]:
    """Invert n×n matrix over F_p (exact modular Gaussian elimination)."""
    aug = [[0] * (2 * n) for _ in range(n)]
    for i in range(n):
        for j in range(n):
            aug[i][j] = M[i][j] % p
        aug[i][n + i] = 1

    for col in range(n):
        pivot = -1
        for r in range(col, n):
            aug[r][col] %= p
            if aug[r][col] < 0:
                aug[r][col] += p
            if aug[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]

        piv_val = aug[col][col]
        inv_piv = fp_inv_mod(piv_val, p)
        if inv_piv == 0:
            return None
        for j in range(2 * n):
            aug[col][j] = (aug[col][j] * inv_piv) % p

        for r in range(n):
            if r == col:
                continue
            factor = aug[r][col]
            if factor == 0:
                continue
            for j in range(2 * n):
                aug[r][j] = (aug[r][j] - factor * aug[col][j]) % p

    inv = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            inv[i][j] = (aug[i][n + j] % p + p) % p
    return inv


def fp_extend_to_invertible(
    M: list[list[int]], m: int, n: int, p: int
) -> list[list[int]]:
    assert m < n
    rows = [list(row) for row in M]
    rows.extend([[0] * n for _ in range(n - m)])

    work = [[x % p for x in row] for row in M]
    pivot_col = [False] * n
    rank = 0
    for col in range(n):
        if rank >= m:
            break
        pivot = -1
        for r in range(rank, m):
            if work[r][col] % p != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_col[col] = True

        inv_piv = fp_inv_mod(work[rank][col], p)
        for j in range(col, n):
            work[rank][j] = (work[rank][j] * inv_piv) % p

        for r in range(m):
            if r == rank:
                continue
            factor = work[r][col]
            if factor == 0:
                continue
            for j in range(col, n):
                work[r][j] = (work[r][j] - factor * work[rank][j]) % p
        rank += 1

    slot = m
    for col in range(n):
        if not pivot_col[col]:
            rows[slot][col] = 1
            slot += 1
    return rows


def fp_extend_columns_to_invertible(
    M: list[list[int]], m: int, n: int, p: int
) -> list[list[int]]:
    assert m > n
    work = [[x % p for x in row] for row in M]
    perm = list(range(m))

    rank = 0
    for col in range(n):
        pivot = -1
        for r in range(rank, m):
            if work[r][col] % p != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        perm[rank], perm[pivot] = perm[pivot], perm[rank]
        inv_piv = fp_inv_mod(work[rank][col], p)
        for j in range(col, n):
            work[rank][j] = (work[rank][j] * inv_piv) % p
        for r in range(m):
            if r == rank:
                continue
            factor = work[r][col]
            if factor == 0:
                continue
            for j in range(col, n):
                work[r][j] = (work[r][j] - factor * work[rank][j]) % p
        rank += 1

    result = [[0] * m for _ in range(m)]
    for r in range(m):
        for c in range(n):
            result[r][c] = M[r][c]
    for j in range(m - n):
        np_row = perm[n + j]
        result[np_row][n + j] = 1
    return result


def fp_has_full_row_rank(
    M: list[list[int]], m: int, n: int, p: int
) -> bool:
    if m > n or m == 0:
        return False
    work = [[x % p for x in row] for row in M]
    rank = 0
    for col in range(n):
        if rank >= m:
            break
        pivot = -1
        for r in range(rank, m):
            if work[r][col] % p != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inv_piv = fp_inv_mod(work[rank][col], p)
        for j in range(col, n):
            work[rank][j] = (work[rank][j] * inv_piv) % p
        for r in range(m):
            if r == rank:
                continue
            factor = work[r][col]
            if factor == 0:
                continue
            for j in range(col, n):
                work[r][j] = (work[r][j] - factor * work[rank][j]) % p
        rank += 1
    return rank == m


# ====================================================================
#  Simplification helpers
# ====================================================================


def _identity_int(n: int) -> list[list[int]]:
    I = [[0] * n for _ in range(n)]
    for i in range(n):
        I[i][i] = 1
    return I


def drop_unused_int(
    f: IntPoly, M: list[list[int]], b: list[int]
) -> tuple[IntPoly, list[list[int]], list[int]]:
    used = f.variables_used()
    if len(used) == f.n:
        return f, M, b

    new_terms: dict[ExpVector, int] = {}
    for exp, c in f.terms.items():
        new_exp = tuple(exp[idx] for idx in used)
        new_terms[new_exp] = c

    new_M = [M[idx] for idx in used]
    new_b = [b[idx] if idx < len(b) else 0 for idx in used]
    return IntPoly(new_terms, len(used), f.mod), new_M, new_b


def try_merge_int(f: IntPoly, i: int, j: int, k: int) -> IntPoly:
    coeffs = [0] * f.n
    coeffs[i] = 1
    coeffs[j] = -k
    return f.substitute_linear(i, coeffs)


def gradient_match_int(g_i: IntPoly, g_j: IntPoly, k: int) -> bool:
    mod = g_i.mod
    if g_i.T() != g_j.T():
        return False
    for exp, c in g_i.terms.items():
        it = g_j.terms.get(exp)
        if it is None:
            return False
        expected = IntPoly.mul_mod(k, c, mod)
        if mod == 0:
            if it != expected:
                return False
        else:
            if (it % mod + mod) % mod != (expected % mod + mod) % mod:
                return False
    return True


# ====================================================================
#  Greedy merge (integer/F_p)
# ====================================================================


def greedy_merge_simplify_int(
    f: IntPoly, max_iter: int = 50, verbose: bool = False
) -> tuple[IntPoly, list[list[int]], list[int]]:
    if not f.terms:
        return IntPoly({}, f.n, f.mod), _identity_int(f.n), [0] * f.n

    cur = f.copy()
    M = _identity_int(f.n)
    b = [0] * f.n
    orig_T = f.T()

    def _reduce_val(x):
        return x if cur.mod == 0 else (x % cur.mod + cur.mod) % cur.mod

    k_values = [1, -1, 2, -2]

    for iteration in range(max_iter):
        active = cur.variables_used()
        if len(active) <= 1:
            break

        best_delta = 0
        best_i = best_j = -1
        best_k = 0

        for i in active:
            for j in active:
                if i == j:
                    continue
                for k in k_values:
                    g_test = try_merge_int(cur, i, j, k)
                    delta = g_test.T() - cur.T()
                    if delta < best_delta:
                        best_delta = delta
                        best_i, best_j, best_k = i, j, k

        if best_delta >= 0:
            break

        for col in range(f.n):
            M[best_i][col] = _reduce_val(M[best_i][col] + best_k * M[best_j][col])
        cur = try_merge_int(cur, best_i, best_j, best_k)
        b[best_i] = _reduce_val(b[best_i] + best_k * b[best_j])

        cur, M, b = drop_unused_int(cur, M, b)

        if verbose:
            pct = (orig_T - cur.T()) * 100.0 / orig_T
            print(f"  iter {iteration}: x{best_i}→x{best_i}+{best_k}·x{best_j}  T={cur.T()}/{orig_T} ({pct:.1f}%↓) m={cur.n}")

        if cur.T() <= 1:
            break

    return cur, M, b


# ====================================================================
#  Gradient-guided simplification
# ====================================================================


def simplify_by_gradient_int(
    f: IntPoly, verbose: bool = False
) -> tuple[IntPoly, list[list[int]], list[int]]:
    if not f.terms:
        return IntPoly({}, f.n, f.mod), _identity_int(f.n), [0] * f.n

    cur = f.copy()
    M = _identity_int(f.n)
    b = [0] * f.n
    orig_T = f.T()

    def _reduce_val(x):
        return x if cur.mod == 0 else (x % cur.mod + cur.mod) % cur.mod

    k_values = [1, -1, 2, -2]

    for iteration in range(30):
        grads = cur.gradient()
        active = cur.variables_used()
        if len(active) <= 1:
            break

        best_delta = 0
        best_i = best_j = -1
        best_k = 0

        # Phase A: gradient match
        for i in active:
            if grads[i].T() == 0:
                continue
            for j in active:
                if i == j or grads[j].T() == 0:
                    continue
                for k in k_values:
                    if gradient_match_int(grads[i], grads[j], k):
                        g_test = try_merge_int(cur, i, j, k)
                        delta = g_test.T() - cur.T()
                        if delta < best_delta:
                            best_delta = delta
                            best_i, best_j, best_k = i, j, k

        # Phase B: exhaustive
        if best_i < 0:
            for i in active:
                for j in active:
                    if i == j:
                        continue
                    for k in k_values:
                        g_test = try_merge_int(cur, i, j, k)
                        delta = g_test.T() - cur.T()
                        if delta < best_delta:
                            best_delta = delta
                            best_i, best_j, best_k = i, j, k

        if best_i < 0 or best_delta >= 0:
            break

        for col in range(f.n):
            M[best_i][col] = _reduce_val(M[best_i][col] + best_k * M[best_j][col])
        cur = try_merge_int(cur, best_i, best_j, best_k)
        b[best_i] = _reduce_val(b[best_i] + best_k * b[best_j])

        cur, M, b = drop_unused_int(cur, M, b)

        if verbose:
            pct = (orig_T - cur.T()) * 100.0 / orig_T
            print(f"  iter {iteration}: x{best_i}→x{best_i}+{best_k}·x{best_j}  T={cur.T()}/{orig_T} ({pct:.1f}%↓) m={cur.n}")

        if cur.T() <= 1:
            break

    return cur, M, b


# ====================================================================
#  Random M,b search
# ====================================================================


def search_random_int(
    f: IntPoly,
    max_m: int,
    n_trials: int,
    seed: int = 1,
    verbose: bool = False,
) -> tuple[IntPoly, list[list[int]], list[int]]:
    n = f.n
    best_T = f.T()
    best = (f, _identity_int(n), [0] * n)

    rng = random.Random(seed)
    entry_dist = lambda: rng.randint(-2, 2)

    for trial in range(n_trials):
        m = rng.randint(1, max_m)
        M = [[0] * n for _ in range(m)]
        b = [0] * m

        if rng.randint(0, 2) < 2:
            # Structured: each row has ≤1 non-zero (±1)
            used_var = [False] * n
            ok = True
            for row in range(m):
                var = rng.randint(0, n - 1)
                if used_var[var]:
                    ok = False
                    break
                used_var[var] = True
                M[row][var] = 1 if rng.randint(0, 1) else -1
                if rng.randint(0, 1):
                    b[row] = entry_dist()
            if not ok:
                continue
        else:
            for row in range(m):
                for col in range(n):
                    M[row][col] = entry_dist()
                b[row] = entry_dist()
            ok = (
                fp_has_full_row_rank(M, m, n, f.mod)
                if f.mod > 0
                else int_has_full_row_rank(M, m, n)
            )
            if not ok:
                continue

        g = f.substitute_affine(M, b)
        if not g.terms:
            continue

        # Try union
        g_best = g
        M_best = M
        b_best = b
        use_union = False

        g_u = f.substitute_affine_union(M, b)
        if g_u.T() < g.T():
            M_best = M + [[1 if i == j else 0 for i in range(n)] for j in range(n)]
            b_best = b + [0] * n
            g_best = g_u
            use_union = True

        if g_best.T() < best_T:
            if len(M_best) != f.n:
                errs = verify_int_poly(f, g_best, M_best, b_best, 20)
                if errs:
                    continue
            best_T = g_best.T()
            best = (g_best, M_best, b_best)
            if verbose:
                pct = (f.T() - best_T) * 100.0 / f.T()
                print(f"  trial {trial}: T={best_T}/{f.T()} ({pct:.1f}%↓) m={len(M_best)}"
                      f"{' [union]' if use_union else ''}")

    return best


# ====================================================================
#  Combined pipeline
# ====================================================================


def compose_matrices(
    acc: tuple[IntPoly, list[list[int]], list[int]],
    next_step: tuple[IntPoly, list[list[int]], list[int]],
    orig_n: int,
) -> None:
    g_acc, M_acc, b_acc = acc
    g_next, M_next, b_next = next_step
    mod = g_next.mod

    def _reduce_val(x):
        return x if mod == 0 else (x % mod + mod) % mod

    m1 = len(M_acc)
    m2 = len(M_next)
    M_new = [[0] * orig_n for _ in range(m2)]
    for j in range(m2):
        for k in range(m1):
            if M_next[j][k] == 0:
                continue
            for i in range(orig_n):
                M_new[j][i] = _reduce_val(M_new[j][i] + M_next[j][k] * M_acc[k][i])

    b_new = [0] * m2
    for j in range(m2):
        for k in range(m1):
            b_new[j] = _reduce_val(b_new[j] + M_next[j][k] * b_acc[k])
        b_new[j] = _reduce_val(b_new[j] + b_next[j])

    return (g_next, M_new, b_new)


def simplify_int(
    f: IntPoly, verbose: bool = False
) -> tuple[IntPoly, list[list[int]], list[int]]:
    orig_T = f.T()
    orig_n = f.n
    if orig_T <= 1:
        return f, _identity_int(orig_n), [0] * orig_n

    if verbose:
        print(f"\n=== IntPoly simplify: n={orig_n}, T₀={orig_T} ===")

    # Phase 1: gradient-guided merge
    r1 = simplify_by_gradient_int(f, verbose)
    acc_g, acc_M, acc_b = r1

    # Phase 2: exhaustive merge
    r2 = greedy_merge_simplify_int(acc_g, 20, verbose)
    acc = compose_matrices(r1, r2, orig_n)
    acc_g, acc_M, acc_b = acc

    if verbose:
        pct = (orig_T - acc_g.T()) * 100.0 / orig_T
        print(f"  After merge: T={acc_g.T()}/{orig_T} ({pct:.1f}%↓), m={acc_g.n}")

    # Phase 3: random search with union
    if acc_g.n > 0 and acc_g.T() > 1:
        max_random_m = max(acc_g.n, 1)
        r3 = search_random_int(acc_g, max_random_m, 100, 2, verbose)
        r3_g, r3_M, r3_b = r3
        if r3_g.T() < acc_g.T():
            acc = compose_matrices(acc, r3, orig_n)
            acc_g, acc_M, acc_b = acc
            if verbose:
                pct = (orig_T - acc_g.T()) * 100.0 / orig_T
                print(f"  After random: T={acc_g.T()}/{orig_T} ({pct:.1f}%↓), m={acc_g.n}")

    if verbose:
        pct = (orig_T - acc_g.T()) * 100.0 / orig_T
        print(f"  Final: T={acc_g.T()}/{orig_T} ({pct:.1f}%↓), m={acc_g.n}")

    return acc


# ====================================================================
#  Verification
# ====================================================================


def verify_int_poly(
    f: IntPoly,
    g: IntPoly,
    M: list[list[int]],
    b: list[int],
    n_tests: int = 100,
) -> int:
    m = len(M)
    mod = f.mod
    rng = random.Random(12345)
    errors = 0

    for _ in range(n_tests):
        x = [rng.randint(-5, 5) for _ in range(f.n)]
        z = [0] * m
        for j in range(m):
            for i in range(f.n):
                z[j] += M[j][i] * x[i]
            z[j] += b[j] if j < len(b) else 0

        fv = f.eval(x)
        gv = g.eval(z)
        if mod > 0:
            fv = (fv % mod + mod) % mod
            gv = (gv % mod + mod) % mod
        if fv != gv:
            errors += 1
    return errors
