"""
SparseANF — pure algebraic ANF simplification via linear-factor extraction.

Core idea: given Boolean function f(x) in ANF, find affine-linear forms
z_j = <M_j, x> ⊕ b_j such that f(x) = g(z) has fewer ANF terms.

This is a PURE ALGEBRAIC rewriting: f(x) and g(Mx⊕b) are equal as polynomials
in F₂[x]/(x_i²−x_i).  No function evaluation, no information loss, no
rank/invertibility restrictions on M.
"""

from __future__ import annotations

import itertools
import math
import random
from typing import Optional

import numpy as np

from bool_anf import (
    _gf2_rank, gf2_matrix_inv, walsh_hadamard_transform,
)


# ====================================================================
#  SparseANF — dict-based, for n up to 64
# ====================================================================

class SparseANF:
    """Boolean function ANF over F₂.  terms: {mask: coeff}.

    Masks are Python ints: bit i = 1 means x_i appears in the monomial.
    Multiplication in the Boolean ring: x_i² = x_i, so combining two
    monomials uses bitwise OR (not XOR).
    """

    def __init__(self, terms: dict[int, int] | set[int] | None = None, n: int = 0):
        if terms is None:
            self.terms: dict[int, int] = {}
        elif isinstance(terms, set):
            self.terms = {m: 1 for m in terms}
        else:
            self.terms = {m: v & 1 for m, v in terms.items() if (v & 1)}
        self.terms = {k: v for k, v in self.terms.items() if v}
        self.n = n

    # ---- constructors ----

    @classmethod
    def from_anf_vec(cls, v: np.ndarray) -> "SparseANF":
        n = int(np.log2(len(v)))
        return cls({m: 1 for m in range(len(v)) if v[m]}, n)

    @classmethod
    def random_cubic(cls, n: int, density: float = 0.15, seed: int = 0) -> "SparseANF":
        """Random degree ≤ 3 function.  Uses combinatorial sampling, never 2^n."""
        rng = random.Random(seed)
        terms: dict[int, int] = {}
        for i in range(n):
            if rng.random() < density * 3:
                terms[1 << i] = 1
        n_pairs = n * (n - 1) // 2
        n_pair_samp = int(density * n_pairs)
        if n_pair_samp > 0 and n_pairs > 0:
            pairs = rng.sample(
                [(i,j) for i in range(n) for j in range(i+1,n)],
                min(n_pair_samp, n_pairs))
            for i,j in pairs:
                terms[(1<<i)|(1<<j)] = 1
        n_triples = n * (n - 1) * (n - 2) // 6
        n_trip_samp = int(density * n_triples)
        if n_trip_samp > 0 and n_triples > 0:
            triples = rng.sample(
                [(i,j,k) for i in range(n) for j in range(i+1,n) for k in range(j+1,n)],
                min(n_trip_samp, n_triples))
        for i,j,k in triples:
            terms[(1<<i)|(1<<j)|(1<<k)] = 1
        if not terms:
            terms[1] = 1
        return cls(terms, n)

    @classmethod
    def from_linear_structure(
        cls, n: int, inner_n: int, seed: int = 0
    ) -> "SparseANF":
        """F(x) = G(Ax) where G has inner_n vars, F has n vars.
        G is sparse (deg ≤ 2, few terms), so expansion is manageable."""
        rng = random.Random(seed + 999)
        A = np.zeros((inner_n, n), dtype=np.uint8)
        while _gf2_rank(A) < inner_n:
            A = np.random.randint(0, 2, size=(inner_n, n)).astype(np.uint8)

        g_terms: dict[int, int] = {}
        for i in range(inner_n):
            if rng.random() < 0.3:
                g_terms[1 << i] = 1
        for i in range(inner_n):
            for j in range(i + 1, inner_n):
                if rng.random() < 0.2:
                    g_terms[(1 << i) | (1 << j)] = 1
        if not g_terms:
            g_terms[1] = 1

        G = cls(g_terms.copy(), inner_n)
        f_terms: dict[int, int] = {}
        for g_mask, _ in G.terms.items():
            expanded = cls._expand_via_matrix(g_mask, A, n)
            for m, v in expanded.items():
                f_terms[m] = f_terms.get(m, 0) ^ v
        f_terms = {k: v for k, v in f_terms.items() if v}
        return cls(f_terms, n)

    @staticmethod
    def _expand_via_matrix(g_mask: int, A: np.ndarray, n: int) -> dict[int, int]:
        """Expand monomial g_mask through A: ∏ (∑ A[i][j] x_j) in Boolean ring."""
        result: dict[int, int] = {0: 1}
        bits = [i for i in range(A.shape[0]) if (g_mask >> i) & 1]
        for i in bits:
            lm = sum((1 << j) for j in range(n) if A[i, j])
            new_res: dict[int, int] = {}
            for mask, val in result.items():
                t = lm
                while t:
                    j_bit = t & -t
                    j = (j_bit.bit_length() - 1)
                    t ^= j_bit
                    new_mask = mask | (1 << j)  # OR for Boolean ring
                    new_res[new_mask] = new_res.get(new_mask, 0) ^ val
            result = new_res
        return result

    @staticmethod
    def _anf_multiply(p1: dict[int, int], p2: dict[int, int]) -> dict[int, int]:
        """Multiply two SparseANF term dicts in the Boolean ring (OR for masks)."""
        result: dict[int, int] = {}
        for m1, v1 in p1.items():
            for m2, v2 in p2.items():
                m = m1 | m2  # Boolean ring: x_i² = x_i → mask OR
                result[m] = result.get(m, 0) ^ (v1 & v2)
        return {k: v for k, v in result.items() if v}

    # ---- accessors ----

    def T(self) -> int:
        return len(self.terms)

    def degree(self) -> int:
        return max((m.bit_count() for m in self.terms), default=0)

    def variables_used(self) -> set[int]:
        """Return set of variable indices appearing in any monomial."""
        used: set[int] = set()
        for m in self.terms:
            i = 0
            while m:
                if m & 1:
                    used.add(i)
                m >>= 1
                i += 1
        return used

    def copy(self) -> "SparseANF":
        return SparseANF(dict(self.terms), self.n)

    def eval_mask(self, x: int) -> int:
        """Evaluate f(x) for x given as int mask."""
        r = 0
        for m in self.terms:
            if (m & x) == m:
                r ^= 1
        return r

    # ---- core: substitute z = Mx⊕b into f(x) to get g(z) ----

    def substitute_affine(self, M: np.ndarray, b: np.ndarray,
                          verify: bool = True) -> "SparseANF":
        """PURE ALGEBRAIC substitution: define z = Mx⊕b, compute g(z) = f(x).

        M is ANY m×n matrix (any rank), b is ANY m-vector.
        This is purely syntactic: for each monomial x^s, rewrite in terms
        of z using the Boolean ring relations z_j² = z_j.

        For m=n invertible: g is unique.
        For m<n (full row rank): f(x) = g(Mx⊕b) has a unique g if f
            factors through M; otherwise no g exists.  We detect this
            via random verification when verify=True.
        For rank-deficient M: the system is underdetermined.
        """
        m, nc = M.shape
        assert nc == self.n

        if m == 0:
            const = self.eval_mask(int(b.ravel()[0]) if len(b) else 0)
            return SparseANF({0: 1} if const else {}, 0)

        rank = _gf2_rank(M)

        if m == self.n and rank == self.n:
            # Invertible: unique g via x = M⁻¹(z⊕b)
            N = gf2_matrix_inv(M)
            c = (N @ b.ravel()) % 2
            return self._expand_with(N, c, m)

        if rank < min(m, self.n):
            return self._substitute_rank_deficient(M, b, m)

        # m < n, full row rank: extend to invertible, transform, restrict
        M_aug = gf2_extend_to_invertible(M)
        b_aug = np.zeros(self.n, dtype=np.uint8)
        b_aug[:m] = b.ravel()

        N_aug = gf2_matrix_inv(M_aug)
        c_aug = (N_aug @ b_aug) % 2

        f_prime = self._expand_with(N_aug, c_aug, self.n)

        # Restrict to first m variables
        g_terms: dict[int, int] = {}
        for zm, v in f_prime.terms.items():
            g_mask = zm & ((1 << m) - 1)
            g_terms[g_mask] = g_terms.get(g_mask, 0) ^ v

        g = SparseANF({k: v for k, v in g_terms.items() if v}, m)

        if verify and m < self.n:
            self._verify_substitution(g, M, b)

        return g

    def _verify_substitution(self, g: "SparseANF", M: np.ndarray, b: np.ndarray,
                            n_tests: int = 20) -> None:
        """Verify g(Mx⊕b) = f(x) for random x.  Raises ValueError if mismatch."""
        m = M.shape[0]
        errors = 0
        for _ in range(n_tests):
            x = random.randint(0, (1 << self.n) - 1)
            z = 0
            for j in range(m):
                zj = 0
                for i in range(self.n):
                    if M[j, i] and ((x >> i) & 1):
                        zj ^= 1
                if b.ravel()[j]:
                    zj ^= 1
                z |= (zj << j)
            if self.eval_mask(x) != g.eval_mask(z):
                errors += 1
        if errors > 0:
            raise ValueError(
                f"g(Mx⊕b) ≠ f(x) on {errors}/{n_tests} tests — "
                f"f does not factor through M(m={m},n={self.n})"
            )

    def _expand_with(self, N: np.ndarray, c: np.ndarray, m: int) -> "SparseANF":
        """Expand f(N(z⊕c)) in F₂[z]/(z_j²−z_j).  N is n×m, c is n×1."""
        result: dict[int, int] = {}
        for x_mask in self.terms:
            expanded = self._expand_monomial(x_mask, N, c, m)
            for zm, v in expanded.items():
                result[zm] = result.get(zm, 0) ^ v
        return SparseANF({k: v for k, v in result.items() if v}, m)

    def _expand_monomial(self, x_mask: int, N: np.ndarray, c: np.ndarray,
                         m: int) -> dict[int, int]:
        """Expand ∏_{i: x_mask_i=1} (c_i ⊕ Σ_j N_ij z_j) in Boolean ring.

        Uses OR for combining masks since z_j² = z_j.
        """
        terms: dict[int, int] = {0: 1}
        bits = [i for i in range(self.n) if (x_mask >> i) & 1]
        for i in bits:
            lm = sum((1 << j) for j in range(m) if N[i, j])
            new_t: dict[int, int] = {}
            for zm, val in terms.items():
                if c[i]:
                    new_t[zm] = new_t.get(zm, 0) ^ val
                if lm:
                    t = lm
                    while t:
                        j_bit = t & -t
                        j = (j_bit.bit_length() - 1)
                        t ^= j_bit
                        nzm = zm | (1 << j)  # Boolean ring: OR
                        new_t[nzm] = new_t.get(nzm, 0) ^ val
            terms = {k: v for k, v in new_t.items() if v}
        return terms

    def _substitute_rank_deficient(self, M: np.ndarray, b: np.ndarray,
                                    m: int) -> "SparseANF":
        """Handle rank-deficient M by restricting to independent rows,
        then padding with unused variables."""
        rank = _gf2_rank(M)
        # Build full-row-rank submatrix
        rows = []
        for j in range(m):
            row = M[j]
            test = np.array(rows + [row], dtype=np.uint8)
            if _gf2_rank(test) == len(rows) + 1:
                rows.append(row)
            if len(rows) == rank:
                break
        M_sub = np.array(rows, dtype=np.uint8)
        b_sub = b[:rank].ravel()
        # Recurse on the sub-problem
        try:
            g = self.substitute_affine(M_sub, b_sub.reshape(-1, 1), verify=True)
        except ValueError:
            # Even with rank-deficient M, f doesn't factor
            raise
        if m > rank:
            # Pad with extra unused z variables
            g_padded = {zm: v for zm, v in g.terms.items()}
            return SparseANF(g_padded, m)
        return g

    # ---- linear system solver for sparsest g (rank-deficient or m<n) ----

    def solve_sparsest_g(self, M: np.ndarray, b: np.ndarray,
                         max_m: int = 10) -> "SparseANF":
        """Find the SPARSEST g such that f(x) = g(Mx⊕b).

        Uses linear algebra over F₂: expands each monomial z^t of g
        into x-monomials and solves V·a = f for the sparsest a.

        This is the CORRECT approach for all cases (including m<n and
        rank-deficient), but computationally limited to m ≤ max_m.
        """
        m = M.shape[0]
        if m > max_m:
            raise ValueError(f"m={m} > max_m={max_m}: exact solve infeasible")
        n = self.n

        # Step 1: Precompute z_j = M_j x ⊕ b_j as SparseANF over x
        # Each z_j is a LINEAR FORM: XOR of individual variables x_i.
        # Represent as SUM of monomials {1<<i: 1 for each i with M_ji=1}.
        z_forms = []
        for j in range(m):
            terms: dict[int, int] = {}
            for i in range(n):
                if M[j, i]:
                    # Each x_i appears as its own monomial (single-variable)
                    terms[1 << i] = terms.get(1 << i, 0) ^ 1
            if b.ravel()[j]:
                terms[0] = terms.get(0, 0) ^ 1  # constant term
            terms = {k: v for k, v in terms.items() if v}
            z_forms.append(terms)

        # Step 2: Compute column C_t = (Mx⊕b)^t for each t=0..2^m-1
        # Using DP: C_0 = {0: 1}, C_t = C_{t\lsb} * z_{lsb_position}
        cols: list[dict[int, int]] = [None] * (1 << m)
        cols[0] = {0: 1}
        for t in range(1, 1 << m):
            lsb = t & -t
            j = (lsb.bit_length() - 1)
            t_prev = t ^ lsb
            cols[t] = self._anf_multiply(cols[t_prev], z_forms[j])

        # Step 3: Build linear system V·a = f
        # Collect all x-monomials that appear
        all_masks: set[int] = set(self.terms.keys())
        for ct in cols:
            all_masks.update(ct.keys())

        mask_list = sorted(all_masks)
        mask_to_idx = {s: i for i, s in enumerate(mask_list)}
        n_eq = len(mask_list)
        n_var = 1 << m

        V = np.zeros((n_eq, n_var), dtype=np.uint8)
        for t, ct in enumerate(cols):
            for s, v in ct.items():
                V[mask_to_idx[s], t] = v

        f_vec = np.zeros(n_eq, dtype=np.uint8)
        for s, v in self.terms.items():
            f_vec[mask_to_idx.get(s, -1)] = v

        # Step 4: Gaussian elimination on augmented matrix [V | f]
        aug = np.concatenate([V, f_vec.reshape(-1, 1)], axis=1)
        n_rows, n_cols = aug.shape  # n_cols = n_var + 1
        pivot_cols: dict[int, int] = {}
        pivot_row = 0

        for col in range(n_var):
            if pivot_row >= n_rows:
                break
            found = False
            for r in range(pivot_row, n_rows):
                if aug[r, col]:
                    found = True
                    aug[[pivot_row, r]] = aug[[r, pivot_row]]
                    break
            if not found:
                continue
            for r in range(n_rows):
                if r != pivot_row and aug[r, col]:
                    aug[r] ^= aug[pivot_row]
            pivot_cols[col] = pivot_row
            pivot_row += 1

        # Check consistency
        for r in range(pivot_row, n_rows):
            if aug[r, -1]:
                raise ValueError("No solution: f cannot be expressed through M,b")

        # Particular solution: free vars = 0
        a = np.zeros(n_var, dtype=np.uint8)
        for col, row in pivot_cols.items():
            a[col] = aug[row, -1]

        # Nullspace basis (free variables)
        null_basis: list[np.ndarray] = []
        for col in range(n_var):
            if col not in pivot_cols:
                vec = np.zeros(n_var, dtype=np.uint8)
                vec[col] = 1
                for pc, pr in pivot_cols.items():
                    if aug[pr, col]:
                        vec[pc] = 1
                null_basis.append(vec)

        # Step 5: Search for minimum weight solution
        best_a = a.copy()
        best_w = int(best_a.sum())

        if len(null_basis) <= 24:  # Up to 16M combos — feasible
            n_null = len(null_basis)
            for mask in range(1 << n_null):
                cand = a.copy()
                for j in range(n_null):
                    if (mask >> j) & 1:
                        cand ^= null_basis[j]
                w = int(cand.sum())
                if w < best_w:
                    best_w = w
                    best_a = cand.copy()
                    if best_w == 0:
                        break

        g_terms = {t: 1 for t in range(n_var) if best_a[t]}
        return SparseANF(g_terms, m)

    # ---- factor extraction ----

    def monomials_by_degree(self) -> dict[int, list[int]]:
        d: dict[int, list[int]] = {}
        for m in self.terms:
            d.setdefault(m.bit_count(), []).append(m)
        return d

    def __repr__(self) -> str:
        return f"SparseANF(n={self.n}, T={self.T()}, deg={self.degree()})"


# ====================================================================
#  Affine helpers
# ====================================================================

def gf2_extend_to_invertible(M: np.ndarray) -> np.ndarray:
    """Extend m×n (m≤n, full row rank) to n×n invertible [M; P]."""
    m, n = M.shape
    assert m <= n and _gf2_rank(M) == m
    if m == n:
        return M.copy()
    P = []
    for i in range(n):
        row = np.zeros(n, dtype=np.uint8)
        row[i] = 1
        test = np.vstack([M, np.array(P + [row], dtype=np.uint8)])
        if _gf2_rank(test) == len(P) + m + 1:
            P.append(row)
        if len(P) == n - m:
            break
    return np.vstack([M, np.array(P, dtype=np.uint8)])


# ====================================================================
#  Walsh-guided dimension search
# ====================================================================

def _walsh_spectrum(f: SparseANF, n_top: int = 100) -> list[int]:
    """Walsh spectrum peak directions.  n≤16: exact FWT; n>16: estimate."""
    if f.n <= 16:
        v = np.zeros(1 << f.n, dtype=np.uint8)
        for m in f.terms:
            v[m] = 1
        w = walsh_hadamard_transform(v)
        mag = np.abs(w)
        k = min(n_top + 10, len(mag) - 1)
        idx = np.argpartition(-mag, k)[:k]
        return list(idx[np.argsort(-mag[idx])])
    else:
        return _walsh_estimate(f, n_top)


def _walsh_estimate(f: SparseANF, n_top: int) -> list[int]:
    """Estimate Walsh directions from ANF structure (n>16)."""
    n = min(f.n, 20)
    rng = random.Random(0)
    scores: dict[int, float] = {}
    for i in range(n):
        mask = 1 << i
        s = sum(1 for m in f.terms if (m & mask).bit_count() % 2 == 0) \
            - sum(1 for m in f.terms if (m & mask).bit_count() % 2 == 1)
        scores[mask] = abs(s)
    pairs = list(itertools.combinations(range(n), 2))[:500]
    for i, j in pairs:
        mask = (1 << i) | (1 << j)
        s = sum(1 for m in f.terms if (m & mask).bit_count() % 2 == 0) \
            - sum(1 for m in f.terms if (m & mask).bit_count() % 2 == 1)
        scores[mask] = abs(s)
    for _ in range(min(2000, max(100, n_top * 5))):
        mask = rng.randint(0, (1 << 16) - 1)
        s = sum(1 for m in f.terms if (m & mask).bit_count() % 2 == 0) \
            - sum(1 for m in f.terms if (m & mask).bit_count() % 2 == 1)
        scores[mask] = abs(s)
    return sorted(scores, key=lambda k: -scores[k])[:n_top + 10]


def search_affine_simplification(
    f: SparseANF,
    max_m: int = 8,
    top_k: int = 80,
    verbose: bool = False,
    n_random: int = 50,
    seed: int = 1,
) -> list[tuple[int, int, np.ndarray, np.ndarray, "SparseANF"]]:
    """Search for M,b that minimize T(g) with m ≤ max_m.

    Strategy:
    1. Walsh directions (single-row candidates)
    2. Combinations of top Walsh directions (multi-row)
    3. Random M,b candidates (dimension reduction)

    Returns list of (m, T(g), M, b, g) sorted by T(g).
    """
    n = f.n
    results = []

    if verbose:
        print(f"\nSearch dim-reduction: n={n}, T(f)={f.T()}, max_m={max_m}")

    # ---- 1. Walsh single-row directions ----
    directions = _walsh_spectrum(f, n_top=top_k + 50)
    for idx in directions[:top_k]:
        if idx == 0 or idx.bit_count() < 1 or idx.bit_count() > 16:
            continue
        M = np.zeros((1, n), dtype=np.uint8)
        M[0, :] = np.array([(idx >> j) & 1 for j in range(n)], dtype=np.uint8)
        for bc in (0, 1):
            b = np.zeros((1, 1), dtype=np.uint8)
            b[0] = bc
            try:
                g = f.substitute_affine(M, b)
                results.append((1, g.T(), M, b, g))
            except ValueError:
                pass

    # ---- 2. Multi-row from top Walsh directions ----
    top_dirs = [idx for idx in directions[:top_k]
                if idx != 0 and idx.bit_count() <= 16][:30]

    for m in range(2, max_m + 1):
        if m > len(top_dirs):
            break
        M = np.zeros((m, n), dtype=np.uint8)
        for j in range(m):
            ω = top_dirs[j]
            M[j, :] = np.array([(ω >> i) & 1 for i in range(n)], dtype=np.uint8)
        if _gf2_rank(M) < m:
            continue
        b = np.zeros((m, 1), dtype=np.uint8)
        try:
            g = f.substitute_affine(M, b)
            results.append((m, g.T(), M, b, g))
        except ValueError:
            pass

    # ---- 3. Random M,b for dimension reduction ----
    rng = random.Random(seed)
    for m in range(1, max_m + 1):
        for _ in range(n_random):
            # Random M with controlled sparsity
            M = np.zeros((m, n), dtype=np.uint8)
            for j in range(m):
                # Each row has ~ n/2 ones
                for i in range(n):
                    if rng.random() < 0.5:
                        M[j, i] = 1
            if _gf2_rank(M) < min(m, n):
                continue
            b = np.array([[rng.randint(0, 1)] for _ in range(m)], dtype=np.uint8)
            try:
                g = f.substitute_affine(M, b)
                results.append((m, g.T(), M, b, g))
            except ValueError:
                pass

    results.sort(key=lambda x: (x[1], x[0]))

    if verbose:
        for m, T, _, _, _ in results[:10]:
            print(f"  m={m}: T(g)={T}/{f.T()} ({(f.T()-T)/f.T()*100:.1f}%↓)")

    return results


# ====================================================================
#  Greedy XOR Merge — pairwise variable merging heuristic
# ====================================================================

def _score_merge(terms: dict[int, int], i: int, j: int) -> int:
    """Exact ΔT = T(after) - T(before) for merging x_i → x_i⊕x_j.

    Uses set-XOR formula: |new| = |unchanged| + |toggles| - 2|unchanged ∩ toggles|
    Derived from computing new_terms = (unchanged) ⊕ (toggles).
    """
    toggles: dict[int, int] = {}
    for m, v in terms.items():
        if (m >> i) & 1:
            toggles[m] = toggles.get(m, 0) ^ v          # toggle off
            rest = m ^ (1 << i)
            new_m = rest if ((m >> j) & 1) else rest | (1 << j)
            toggles[new_m] = toggles.get(new_m, 0) ^ v  # toggle on

    n_changed = sum(1 for m in terms if (m >> i) & 1)
    n_toggles = sum(1 for c in toggles.values() if c)
    n_overlap = sum(1 for m, c in toggles.items()
                    if c and m in terms and not ((m >> i) & 1))
    return n_toggles - 2 * n_overlap - n_changed


def _apply_merge(terms: dict[int, int], i: int, j: int) -> dict[int, int]:
    """Replace x_i with x_i⊕x_j in the ANF.  Returns new term dict."""
    new: dict[int, int] = {}
    for m, v in terms.items():
        if (m >> i) & 1:
            new[m] = new.get(m, 0) ^ v          # toggle off
            rest = m ^ (1 << i)
            new_m = rest if ((m >> j) & 1) else rest | (1 << j)
            new[new_m] = new.get(new_m, 0) ^ v  # toggle on
        else:
            new[m] = new.get(m, 0) ^ v
    return {k: v for k, v in new.items() if v}


def greedy_merge_simplify(
    f: SparseANF,
    max_iter: int = 100,
    verbose: bool = False,
) -> tuple[SparseANF, np.ndarray, np.ndarray]:
    """Greedy pairwise XOR merge to reduce T(f) and m.

    At each step, find the pair (i,j) such that replacing x_i with
    x_i⊕x_j gives the largest reduction in T.  Apply it, drop unused
    variables, and repeat.  M tracks the accumulated transformation.

    Returns (g, M, b) where f(x) = g(Mx⊕b) and g has ≤ f.n variables.
    """
    terms = dict(f.terms)  # {mask: coeff} in current variable space
    if not terms:
        return SparseANF({}, 0), np.eye(f.n, dtype=np.uint8), np.zeros((f.n, 1), dtype=np.uint8)

    m = f.n
    M = np.eye(m, dtype=np.uint8)   # rows correspond to current variables
    b = np.zeros((m, 1), dtype=np.uint8)
    orig_T = f.T()

    if verbose:
        print(f"\nGreedy merge: n={m}, T₀={orig_T}")

    for iteration in range(max_iter):
        if len(terms) <= 1:
            break

        # ---- find active variables ----
        active = set()
        for mask in terms:
            t = mask
            while t:
                lsb = t & -t
                active.add((lsb.bit_length() - 1))
                t ^= lsb
        active_list = sorted(active)

        if len(active_list) <= 1:
            break

        # ---- precompute per-variable term lists ----
        terms_by_var: dict[int, list[int]] = {var: [] for var in active_list}
        for mask in terms:
            t = mask
            while t:
                lsb = t & -t
                i = (lsb.bit_length() - 1)
                terms_by_var[i].append(mask)
                t ^= lsb

        # ---- score all pairs (i,j) with i != j ----
        best_delta = 0
        best_i = best_j = -1

        for i in active_list:
            ti = terms_by_var[i]
            if not ti:
                continue
            n_changed_i = len(ti)
            for j in active_list:
                if i == j:
                    continue

                # Compute delta using the toggle method
                toggle_counts: dict[int, int] = {}
                for mask_m in ti:
                    toggle_counts[mask_m] = toggle_counts.get(mask_m, 0) ^ 1
                    rest = mask_m ^ (1 << i)
                    new_m = rest if ((mask_m >> j) & 1) else rest | (1 << j)
                    toggle_counts[new_m] = toggle_counts.get(new_m, 0) ^ 1

                n_toggles = sum(1 for c in toggle_counts.values() if c)
                # Overlap: toggled masks that are also in unchanged terms
                n_overlap = 0
                for mask_m, c in toggle_counts.items():
                    if c and mask_m in terms and not ((mask_m >> i) & 1):
                        n_overlap += 1

                delta = n_toggles - 2 * n_overlap - n_changed_i
                if delta < best_delta:
                    best_delta = delta
                    best_i, best_j = i, j

        if best_delta >= 0:
            break

        i, j = best_i, best_j

        # ---- apply merge ----
        M[i] ^= M[j]   # row operation: new z_i = z_i ⊕ z_j
        terms = _apply_merge(terms, i, j)

        # ---- drop unused variables ----
        used = set()
        for mask in terms:
            t = mask
            while t:
                lsb = t & -t
                used.add((lsb.bit_length() - 1))
                t ^= lsb

        if len(used) < m:
            # Remap variable indices and M rows
            sorted_used = sorted(used)
            old_to_new = {old: new for new, old in enumerate(sorted_used)}

            remapped: dict[int, int] = {}
            for mask, v in terms.items():
                new_mask = 0
                for old_k in sorted_used:
                    if (mask >> old_k) & 1:
                        new_mask |= (1 << old_to_new[old_k])
                remapped[new_mask] = remapped.get(new_mask, 0) ^ v
            terms = {k: v for k, v in remapped.items() if v}

            M = M[list(sorted_used)]  # keep only used rows
            b = b[list(sorted_used)]
            m = len(sorted_used)

        if verbose and (iteration % 10 == 0 or best_delta < 0):
            pct = (orig_T - len(terms)) / orig_T * 100
            print(f"  iter {iteration}: T={len(terms)}/{orig_T} ({pct:.1f}%↓) m={m}")

    if verbose:
        pct = (orig_T - len(terms)) / orig_T * 100
        print(f"  Final: T={len(terms)}/{orig_T} ({pct:.1f}%↓) m={m}")

    return SparseANF(terms, m), M, b


# ====================================================================
#  Complement variable selection: z_i = x_i or x_i⊕1
# ====================================================================


def simplify_by_complement(
    f: SparseANF,
    verbose: bool = False,
) -> tuple[SparseANF, np.ndarray, np.ndarray]:
    """Minimize T(g) by complementing variables: z_i = x_i or z_i = x_i⊕1.

    For n ≤ 16: exhaustive over all 2ⁿ complement patterns (Gray code, fast).
    For n > 16: greedy bit-flipping search.

    Returns (g, M, b) where f(x) = g(Mx⊕b), M=I_n, b is complement vector.
    After complement, unused z variables are dropped.
    """
    n = f.n
    M = np.eye(n, dtype=np.uint8)
    b = np.zeros((n, 1), dtype=np.uint8)

    if n <= 0:
        return f, M, b

    if n <= 16:
        # ---- exhaustive: Gray code incremental, O(2ⁿ · T(f)) ----
        # g(b) = f(z⊕b).  Flipping b_i corresponds to g'(z) = g(z⊕e_i):
        #   for each monomial m with bit i, toggle m without i.
        best_T = f.T()
        best_b = np.zeros((n, 1), dtype=np.uint8)

        if verbose:
            print(f"\nComplement search (Gray, 2^{n}={1<<n}): n={n}, T₀={best_T}")

        # Start with g = f (b=0)
        cur_terms: dict[int, int] = dict(f.terms)
        cur_T = len(cur_terms)

        # Gray code: iterate over b masks in Gray order
        prev_mask = 0
        for step in range(1, 1 << n):
            # Find which bit flips in Gray code from prev to step
            gray = step ^ (step >> 1)
            change_bit = (prev_mask ^ gray).bit_length() - 1
            prev_mask = gray

            # Update g: flip effect of bit change_bit
            # g'(z) = g(z⊕e_i): toggle m\{i} for each m with bit i
            toggled = 0
            for m, v in list(cur_terms.items()):
                if (m >> change_bit) & 1:
                    m_without = m ^ (1 << change_bit)
                    cur_terms[m_without] = cur_terms.get(m_without, 0) ^ v
                    if cur_terms[m_without] == 0:
                        del cur_terms[m_without]
                        toggled -= 1
                    else:
                        toggled += 1
            cur_T += toggled

            if cur_T < best_T:
                best_T = cur_T
                best_b = np.array([(gray >> i) & 1 for i in range(n)], dtype=np.uint8).reshape(-1, 1)

        if verbose:
            pct = (f.T() - best_T) / f.T() * 100
            print(f"  best: T={best_T}/{f.T()} ({pct:.1f}%↓)")

        # Compute best g from best_b using substitute_affine
        if best_T < f.T():
            g = f.substitute_affine(M, best_b)
        else:
            g = f
        b = best_b
    else:
        # ---- greedy: flip bits one at a time ----
        best_g = f
        best_b = np.zeros((n, 1), dtype=np.uint8)
        b_vec = np.zeros((n, 1), dtype=np.uint8)
        best_T = f.T()
        improved = True
        it = 0

        if verbose:
            print(f"\nComplement search (greedy): n={n}, T₀={best_T}")

        while improved:
            improved = False
            for i in range(n):
                b_vec[i] ^= 1  # flip bit i
                g = f.substitute_affine(M, b_vec)
                if g.T() < best_T:
                    best_T = g.T()
                    best_b = b_vec.copy()
                    best_g = g
                    improved = True
                    if verbose:
                        pct = (f.T() - best_T) / f.T() * 100
                        print(f"  iter {it}: flip[{i}] → T={best_T}/{f.T()} ({pct:.1f}%↓)")
                else:
                    b_vec[i] ^= 1  # flip back
            it += 1
            if it > 2 * n:
                break

        g = best_g
        b = best_b

    # ---- drop unused variables ----
    used = set()
    for mask in g.terms:
        t = mask
        while t:
            lsb = t & -t
            used.add((lsb.bit_length() - 1))
            t ^= lsb

    if len(used) < g.n:
        sorted_used = sorted(used)
        M_compact = M[list(sorted_used)]
        b_compact = b[list(sorted_used)]
        # Remap term indices
        old_to_new = {old: new for new, old in enumerate(sorted_used)}
        remapped = {}
        for mask, v in g.terms.items():
            new_mask = 0
            for old_k in sorted_used:
                if (mask >> old_k) & 1:
                    new_mask |= (1 << old_to_new[old_k])
            remapped[new_mask] = remapped.get(new_mask, 0) ^ v
        g = SparseANF({k: v for k, v in remapped.items() if v}, len(sorted_used))
        return g, M_compact, b_compact

    return g, M, b


# ====================================================================
#  Combined simplification pipeline
# ====================================================================


def simplify(
    f: SparseANF,
    verbose: bool = False,
) -> tuple[SparseANF, np.ndarray, np.ndarray]:
    """Combined ANF simplification pipeline.

    Stages:
    1. Complement selection (greedy/exhaustive z_i = x_i or x_i⊕1)
    2. Greedy pairwise XOR merge
    3. Walsh-guided dimension reduction (from search_affine_simplification)

    Returns (g, M, b) where f(x) = g(Mx⊕b).
    """
    n, orig_T = f.n, f.T()

    if verbose:
        print(f"\n=== Simplify: n={n}, T₀={orig_T} ===")

    g, M, b = simplify_by_complement(f, verbose=verbose)

    g2, M2, b2 = greedy_merge_simplify(g, verbose=verbose)
    M = (M2 @ M) % 2
    b = ((M2 @ b) % 2).reshape(-1, 1) ^ b2
    g = g2

    if g.n > 3 and g.T() > 0:
        try:
            results = search_affine_simplification(
                g, max_m=min(g.n, 10), top_k=80, verbose=verbose, n_random=30)
        except Exception:
            results = []
        if results and results[0][1] < g.T():
            m_w, T_w, M_w, b_w, g_w = results[0]
            if verbose:
                print(f"  Walsh: T={T_w}/{g.T()} ({(g.T()-T_w)/g.T()*100:.1f}%↓)")
            M = (M_w @ M) % 2
            b = ((M_w @ b) % 2).reshape(-1, 1) ^ b_w
            g = g_w

    if verbose:
        pct = (orig_T - g.T()) / orig_T * 100
        print(f"  Final: T={g.T()}/{orig_T} ({pct:.1f}%↓), m={g.n}")

    return g, M, b


# ====================================================================
#  Demos
# ====================================================================

def demo_factor():
    """Show factorization finds the right structure."""
    print("=" * 60)
    print("Demo: factor extraction")
    print("=" * 60)

    f = SparseANF({0b011: 1, 0b101: 1}, n=3)
    print(f"f: {f}")

    M = np.array([[1, 0, 0], [0, 1, 1]], dtype=np.uint8)
    b = np.zeros((2, 1), dtype=np.uint8)
    g = f.substitute_affine(M, b)
    print(f"M(2×3):\n{M}")
    print(f"g: T={g.T()}, terms={[bin(m) for m in g.terms]}")
    print(f"T(g)/T(f) = {g.T()}/{f.T()} = {g.T()/f.T():.2f}")

    rng = random.Random(0)
    errors = 0
    for _ in range(20):
        x = rng.randint(0, 7)
        z = 0
        for j in range(2):
            z_bit = sum(M[j][i] * ((x >> i) & 1) for i in range(3)) % 2
            z |= z_bit << j
        if f.eval_mask(x) != g.eval_mask(z):
            errors += 1
    print(f"Verify: {errors}/20 errors")
    print()


def demo_auto_search():
    """Auto-search for M,b on the original 3-variable problem."""
    print("=" * 60)
    print("Demo: auto search for M,b")
    print("=" * 60)

    f = SparseANF({0b011: 1, 0b101: 1}, n=3)
    print(f"f: T={f.T()}")
    results = search_affine_simplification(f, max_m=3, top_k=50, verbose=True)
    if results:
        best = results[0]
        print(f"\nBest: m={best[0]}, T(g)={best[1]}/{f.T()} "
              f"({(f.T()-best[1])/f.T()*100:.1f}%↓)")
    print()


def demo_random_n12():
    """Test on random n=12 cubic function."""
    print("=" * 60)
    print("Demo: random n=12 cubic")
    print("=" * 60)

    f = SparseANF.random_cubic(12, density=0.15, seed=42)
    print(f"f: {f}")
    results = search_affine_simplification(f, max_m=6, top_k=100, verbose=True)
    if results:
        best = results[0]
        print(f"\nBest: m={best[0]}, T(g)={best[1]}/{f.T()} "
              f"({(f.T()-best[1])/f.T()*100:.1f}%↓)")


def demo_n16_structured():
    """Structured n=16 that factors through m=4 linear forms."""
    print("=" * 60)
    print("Demo: n=16 factoring through m=4")
    print("=" * 60)

    # g(z) = z0z1 + z2z3, m=4
    # z0 = x0+x1+x2, z1=x3+x4+x5, z2=x6+x7+x8, z3=x9+x10+x11
    n = 16

    # f(x) = (x0+x1+x2)(x3+x4+x5) + (x6+x7+x8)(x9+x10+x11)
    #      = Σ_{i∈{0,1,2},j∈{3,4,5}} x_i x_j + Σ_{i∈{6,7,8},j∈{9,10,11}} x_i x_j
    f_terms = {}
    for i in range(3):
        for j in range(3, 6):
            mask = (1 << i) | (1 << j)
            f_terms[mask] = 1
    for i in range(6, 9):
        for j in range(9, 12):
            mask = (1 << i) | (1 << j)
            f_terms[mask] = f_terms.get(mask, 0) ^ 1
    f_terms = {k: v for k, v in f_terms.items() if v}
    f = SparseANF(f_terms, n)
    print(f"f(x): {f}, T={f.T()}")

    print("\nSearch for dimension reduction:")
    results = search_affine_simplification(f, max_m=6, top_k=100, verbose=True,
                                          n_random=200)

    if results:
        best = results[0]
        g = best[4]
        print(f"\nBest: m={best[0]}, T(g)={best[1]}/{f.T()} "
              f"({(f.T()-best[1])/f.T()*100:.1f}%↓)")

        rng = random.Random(0)
        errors = 0
        for _ in range(50):
            x = rng.randint(0, (1 << n) - 1)
            z = 0
            for j in range(best[0]):
                z_bit = sum(best[2][j][i] * ((x >> i) & 1) for i in range(n)) % 2
                z |= z_bit << j
            if f.eval_mask(x) != g.eval_mask(z):
                errors += 1
        print(f"Verify: {errors}/50 errors {'OK' if errors == 0 else 'FAIL'}")
    print()


def demo_sparsest_g_via_linalg():
    """Demo: use solve_sparsest_g on a rank-deficient case."""
    print("=" * 60)
    print("Demo: solve_sparsest_g (linear algebra method)")
    print("=" * 60)

    # f(x) = x0x1 + x0x2 + x1x2 + x0
    # M = [[1,1,0]] (m=1, z0 = x0⊕x1)
    # Can we find g(z0) with ≤ 1 term?
    f = SparseANF({0b011: 1, 0b101: 1, 0b110: 1, 0b001: 1}, n=3)
    M = np.array([[1, 1, 0]], dtype=np.uint8)
    b = np.zeros((1, 1), dtype=np.uint8)

    print(f"f: {f}")
    print(f"M(1×3): [{M[0]}]")

    try:
        g = f.solve_sparsest_g(M, b)
        print(f"g (sparsest): T={g.T()}, terms={[bin(m) for m in g.terms]}")

        # Verify
        errors = 0
        rng = random.Random(0)
        for _ in range(20):
            x = rng.randint(0, 7)
            z = (M[0, 0] * ((x>>0)&1) ^ M[0,1] * ((x>>1)&1) ^ M[0,2] * ((x>>2)&1)) & 1
            if f.eval_mask(x) != g.eval_mask(z):
                errors += 1
        print(f"Verify: {errors}/20 errors")
    except ValueError as e:
        print(f"  No solution: {e}")
    print()


# ====================================================================
#  Main
# ====================================================================

if __name__ == "__main__":
    demo_factor()
    demo_auto_search()
    demo_random_n12()
    demo_n16_structured()
    demo_sparsest_g_via_linalg()
