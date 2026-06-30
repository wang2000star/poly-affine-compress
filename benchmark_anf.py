"""
Medium-Scale Boolean Function ANF Simplification Benchmark
============================================================

Tests for n = 8, 12, 16, 24, 32, 48, 64 using sparse ANF representations.

Key design for large n (≥ 24):
  - No iteration over all 2^n masks — combinatorial generation only
  - Functions have explicit sparse ANFs generated from combinatorial formulas
  - Affine transformation applied directly to ANF terms (not via truth table)
  - Heuristic Walsh-guided search
"""

from __future__ import annotations

import itertools
import random
import sys
import time
from typing import Optional

import numpy as np

from bool_anf import (
    _gf2_rank,
    gf2_matrix_inv,
    truth_table_to_anf,
    walsh_hadamard_transform,
)


# ====================================================================
#  SparseANF — dict-based ANF for n up to 64
# ====================================================================


class SparseANF:
    """Sparse Boolean function ANF: dict {mask: coeff}, mask is Python int.

    n can be up to 64. Storage ~ O(T(f)), never allocates 2^n.
    """

    def __init__(self, terms: dict[int, int] | set[int] | None = None, n: int = 0):
        if terms is None:
            self.terms: dict[int, int] = {}
        elif isinstance(terms, set):
            self.terms = {m: 1 for m in terms}
        else:
            self.terms = {m: (v & 1) for m, v in terms.items() if (v & 1)}
        self.terms = {k: v for k, v in self.terms.items() if v}
        self.n = n

    @classmethod
    def from_anf_vec(cls, v: np.ndarray) -> "SparseANF":
        n = int(np.log2(len(v)))
        return cls({m: 1 for m in range(len(v)) if v[m]}, n)

    @classmethod
    def random_cubic(cls, n: int, density: float = 0.15, seed: int = 0) -> "SparseANF":
        """Random degree ≤ 3 function. Uses combinatorial, NOT iterating 2^n."""
        rng = random.Random(seed)
        terms: dict[int, int] = {}
        # deg 1
        for i in range(n):
            if rng.random() < density * 3:
                terms[1 << i] = 1
        # deg 2
        n_pairs = min(n * (n - 1) // 2, 50000)
        pairs = rng.sample([(i, j) for i in range(n) for j in range(i + 1, n)],
                           min(n_pairs, max(100, int(density * n * (n - 1) / 2))))
        for i, j in pairs:
            terms[(1 << i) | (1 << j)] = 1
        # deg 3
        n_triples = min(n * (n - 1) * (n - 2) // 6, 20000)
        triples = rng.sample([(i, j, k) for i in range(n) for j in range(i + 1, n) for k in range(j + 1, n)],
                             min(n_triples, max(50, int(density * n * (n - 1) * (n - 2) / 12))))
        for i, j, k in triples:
            terms[(1 << i) | (1 << j) | (1 << k)] = 1
        if not terms:
            terms[1] = 1
        return cls(terms, n)

    @classmethod
    def from_linear_structure(
        cls, n: int, inner_n: int, seed: int = 0
    ) -> "SparseANF":
        """F(x) = G(Ax) where G has inner_n vars, F has n vars.
        Uses very sparse G (deg ≤ 2, few terms) so expansion is manageable."""
        rng = random.Random(seed + 999)
        # A: inner_n × n, full row rank
        A = np.zeros((inner_n, n), dtype=np.uint8)
        while _gf2_rank(A) < inner_n:
            A = np.random.randint(0, 2, size=(inner_n, n)).astype(np.uint8)

        # Sparse G: deg ≤ 2, 2-4 terms
        g_terms: dict[int, int] = {}
        # some deg 1 terms
        for i in range(inner_n):
            if rng.random() < 0.3:
                g_terms[1 << i] = 1
        # some deg 2 terms
        for i in range(inner_n):
            for j in range(i + 1, inner_n):
                if rng.random() < 0.2:
                    g_terms[(1 << i) | (1 << j)] = 1
        if not g_terms:
            g_terms[1] = 1  # fallback

        G = cls(g_terms.copy(), inner_n)

        # Expand F(x) = G(Ax)
        f_terms: dict[int, int] = {}
        for g_mask, _ in G.terms.items():
            expanded = cls._expand_via_matrix(g_mask, A, n)
            for m, v in expanded.items():
                f_terms[m] = f_terms.get(m, 0) ^ v
        f_terms = {k: v for k, v in f_terms.items() if v}
        return cls(f_terms, n)

    @staticmethod
    def _expand_via_matrix(g_mask: int, A: np.ndarray, n: int) -> dict[int, int]:
        """Expand monomial g_mask of inner_n vars through A: ∏ (∑ A[i][j] x_j)."""
        inner_n = A.shape[0]
        result: dict[int, int] = {0: 1}
        bits = [i for i in range(inner_n) if (g_mask >> i) & 1]
        for i in bits:
            # linear_mask: which x_j appear in y_i = sum_j A[i][j] x_j
            lm = sum((1 << j) for j in range(n) if A[i, j])
            if not lm:
                continue
            new_res: dict[int, int] = {}
            for mask, val in result.items():
                t = lm
                while t:
                    j_bit = t & -t
                    j = (j_bit.bit_length() - 1)
                    t ^= j_bit
                    new_mask = mask | (1 << j)
                    new_res[new_mask] = new_res.get(new_mask, 0) ^ val
            result = new_res
        return result

    def T(self) -> int:
        return len(self.terms)

    def degree(self) -> int:
        return max((m.bit_count() for m in self.terms), default=0)

    def copy(self) -> "SparseANF":
        return SparseANF(dict(self.terms), self.n)

    def _eval_at_mask(self, x: int) -> int:
        """Evaluate f(x) given x as int mask."""
        r = 0
        for m in self.terms:
            if (m & x) == m:
                r ^= 1
        return r

    def apply_affine(self, M: np.ndarray, b: np.ndarray) -> "SparseANF":
        """g(z) = f(x) where z = Mx ⊕ b, M invertible n×n."""
        m, nc = M.shape
        assert nc == self.n and m == self.n and _gf2_rank(M) == self.n
        N = gf2_matrix_inv(M)
        c = (N @ b.ravel()) % 2
        result: dict[int, int] = {}
        for x_mask in self.terms:
            exp = self._expand_affine(x_mask, N, c, m)
            for zm, v in exp.items():
                result[zm] = result.get(zm, 0) ^ v
        return SparseANF({k: v for k, v in result.items() if v}, m)

    def _expand_affine(self, x_mask: int, N: np.ndarray, c: np.ndarray, m: int) -> dict[int, int]:
        """Expand ∏_{i in x_mask} (c_i ⊕ Σ_j N_ij z_j) in F₂[z]/(z_j²−z_j)."""
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
                        nzm = zm | (1 << j)
                        new_t[nzm] = new_t.get(nzm, 0) ^ val
            terms = {k: v for k, v in new_t.items() if v}
        return terms

    def walsh_indices(self, n_top: int = 100) -> list[int]:
        """Get Walsh spectrum peak indices.
        n≤16: exact FWT (fast on 2^n entries).
        n>16: estimate from ANF terms directly."""
        if self.n <= 16:
            v = np.zeros(1 << self.n, dtype=np.uint8)
            for m in self.terms:
                v[m] = 1
            w = walsh_hadamard_transform(v)
            mag = np.abs(w)
            # argpartition is O(n) instead of O(n log n)
            k = min(n_top + 10, len(mag) - 1)
            idx = np.argpartition(-mag, k)[:k]
            return list(idx[np.argsort(-mag[idx])])  # sort only the top-k
        else:
            return self._walsh_estimate(n_top)

    def _walsh_estimate(self, n_top: int) -> list[int]:
        """Estimate good Walsh directions from ANF structure.
        Only considers directions up to n=20 (larger directions are unlikely
        to give useful signal for sparse functions)."""
        n = min(self.n, 20)
        rng = random.Random(0)
        scores: dict[int, float] = {}

        # deg-1 directions: each single variable
        for i in range(n):
            mask = 1 << i
            s = sum(
                1 for m in self.terms
                if (m & mask).bit_count() % 2 == 0
            ) - sum(
                1 for m in self.terms
                if (m & mask).bit_count() % 2 == 1
            )
            scores[mask] = abs(s)

        # deg-2 directions: pairs of variables
        pairs = list(itertools.combinations(range(n), 2))[:500]
        for i, j in pairs:
            mask = (1 << i) | (1 << j)
            s = sum(
                1 for m in self.terms
                if (m & mask).bit_count() % 2 == 0
            ) - sum(
                1 for m in self.terms
                if (m & mask).bit_count() % 2 == 1
            )
            scores[mask] = abs(s)

        # random directions (limited to n=16 max for tractability)
        for _ in range(min(2000, max(100, n_top * 5))):
            mask = rng.randint(0, (1 << 16) - 1)  # limit to 16-bit masks
            s = sum(
                1 for m in self.terms
                if (m & mask).bit_count() % 2 == 0
            ) - sum(
                1 for m in self.terms
                if (m & mask).bit_count() % 2 == 1
            )
            scores[mask] = abs(s)

        return sorted(scores, key=lambda k: -scores[k])[:n_top + 10]

    def monomials_by_degree(self) -> dict[int, list[int]]:
        dct: dict[int, list[int]] = {}
        for m in self.terms:
            d = m.bit_count()
            dct.setdefault(d, []).append(m)
        return dct

    def __repr__(self) -> str:
        return f"SparseANF(n={self.n}, T={self.T()}, deg={self.degree()})"


# ====================================================================
#  Heuristic Minimization
# ====================================================================


def repeated_walsh_minimize(
    f: SparseANF,
    max_iter: int = 8,
    top_k: int = 80,
    verbose: bool = False,
) -> tuple[SparseANF, np.ndarray, np.ndarray]:
    """Iterative Walsh-guided ANF minimization.

    At each step: compute Walsh spectrum of current function,
    test top directions as single-row replacements, apply best.
    """
    n = f.n
    current = f.copy()
    M_acc = np.eye(n, dtype=np.uint8)
    b_acc = np.zeros((n, 1), dtype=np.uint8)
    orig_T = f.T()

    if verbose:
        print(f"\n  {{Repeated Walsh}}  n={n}, T₀={orig_T}, deg={f.degree()}")

    for it in range(max_iter):
        indices = current.walsh_indices(n_top=top_k + 20)
        best_T = current.T()
        best_g = None
        best_Ms = None
        best_bs = None
        tested = 0

        for idx in indices:
            if idx == 0:
                continue
            if tested >= top_k:
                break
            ω = idx
            # skip directions with too many bits (computationally expensive)
            if ω.bit_count() > 12:
                continue
            i0 = int(ω & -ω).bit_length() - 1

            M_step = np.eye(n, dtype=np.uint8)
            M_step[i0, :] = np.array([(ω >> j) & 1 for j in range(n)], dtype=np.uint8)

            for bc in (0, 1):
                b_step = np.zeros((n, 1), dtype=np.uint8)
                b_step[i0] = bc
                tested += 1
                try:
                    g = current.apply_affine(M_step, b_step)
                    T = g.T()
                    if T < best_T:
                        best_T, best_g = T, g
                        best_Ms, best_bs = M_step.copy(), b_step.copy()
                except Exception:
                    continue

        if best_g is None:
            if verbose:
                print(f"    iter {it}: no improvement, stop")
            break

        M_acc = (best_Ms @ M_acc) % 2
        b_acc = (best_Ms @ b_acc + best_bs) % 2
        current = best_g

        red = (orig_T - best_T) / orig_T * 100
        if verbose:
            print(f"    iter {it}: T={best_T}/{orig_T} ({red:.1f}%↓)  tested={tested}")

        if best_T <= 1:
            break

    return current, M_acc, b_acc


# ====================================================================
#  Verification: n=8 with correctness check
# ====================================================================


def verify_n8():
    """n=8: detailed verification with correctness."""
    print("=" * 66)
    print("  VERIFY n=8 — correctness & detailed output")
    print("=" * 66)

    rng = random.Random(42)
    n = 8
    terms: dict[int, int] = {}

    # Mix of degrees
    for i in range(n):
        if rng.random() < 0.5:
            terms[1 << i] = 1
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.25:
                terms[(1 << i) | (1 << j)] = 1
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if rng.random() < 0.10:
                    terms[(1 << i) | (1 << j) | (1 << k)] = 1

    f = SparseANF(terms, n)
    print(f"  f: {f}")
    bd = f.monomials_by_degree()
    for d in sorted(bd):
        print(f"    deg={d}: {len(bd[d])} terms")

    t0 = time.time()
    g, M, b = repeated_walsh_minimize(f, max_iter=10, top_k=200, verbose=True)
    elapsed = time.time() - t0
    red = (f.T() - g.T()) / f.T() * 100
    print(f"\n  → T: {f.T()} → {g.T()} ({red:.1f}%↓)  deg: {f.degree()}→{g.degree()}")
    print(f"  → time: {elapsed:.3f}s")

    # correctness check
    N = gf2_matrix_inv(M)
    c = (N @ b.ravel()) % 2
    errors = 0
    for _ in range(50):
        x = rng.randint(0, (1 << n) - 1)
        z_arr = (M @ np.array([(x >> i) & 1 for i in range(n)], dtype=np.uint8)) % 2
        z_arr = (z_arr ^ b.ravel()) % 2
        z = sum(int(z_arr[i]) << i for i in range(g.n))
        fx = f._eval_at_mask(x)
        gz = g._eval_at_mask(z)
        if fx != gz:
            errors += 1

    status = "✅ PASS" if errors == 0 else f"❌ {errors} errors"
    print(f"  correctness: {status}")
    print()


# ====================================================================
#  Speed benchmarks: transformation timing
# ====================================================================


def bench_transform_time():
    """Measure SparseANF.apply_affine time for various sizes."""
    print("=" * 66)
    print("  SPEED: SparseANF.apply_affine() scaling")
    print("=" * 66)
    print(f"  {'n':<4} {'T(f)':<8} {'deg':<5} {'time(ms)':<10}")
    print("  " + "-" * 30)
    for n in [8, 12, 16, 24, 32]:
        f = SparseANF.random_cubic(n, density=0.15, seed=0)
        M = np.eye(n, dtype=np.uint8)
        # non-trivial M: random invertible
        while _gf2_rank(M) < n:
            M = np.random.randint(0, 2, size=(n, n)).astype(np.uint8)
        b = np.zeros((n, 1), dtype=np.uint8)
        t0 = time.time()
        _ = f.apply_affine(M, b)
        t = (time.time() - t0) * 1000
        print(f"  {n:<4} {f.T():<8} {f.degree():<5} {t:<9.1f}ms")
    print()


# ====================================================================
#  Full benchmark
# ====================================================================


def make_power_function(n: int, e: int) -> SparseANF:
    """F(x) = LSB(x^e) over F_{2^n}. Only for n ≤ 16."""
    k = n
    if k > 16:
        raise ValueError(f"Power function ANF only feasible for n ≤ 16, got {n}")

    primitives = {
        4: 0b10011,
        5: 0b100101,
        6: 0b1000011,
        8: 0b100011011,
        12: 0b1000000010011,
        16: 0b10000000000001011,
    }
    poly = primitives.get(k, (1 << k) | (1 << (k // 2)) | 1)
    size = 1 << k

    def gf_mult(a: int, b: int) -> int:
        result = 0
        for i in range(k):
            if (b >> i) & 1:
                result ^= a << i
        for i in range(2 * k - 2, k - 1, -1):
            if (result >> i) & 1:
                result ^= poly << (i - k)
        return result & ((1 << k) - 1)

    tt = np.zeros(size, dtype=np.uint8)
    for x in range(size):
        result = 1
        exp = e
        base = x
        while exp:
            if exp & 1:
                result = gf_mult(result, base)
            exp >>= 1
            base = gf_mult(base, base)
        tt[x] = (result >> 0) & 1

    anf_vec = truth_table_to_anf(tt)
    return SparseANF.from_anf_vec(anf_vec)


def make_hidden_structure(n: int, inner_n: int, seed: int) -> SparseANF:
    """F(x) = G(H⁻¹x) where G is sparse in its natural coordinates.
    This means there EXISTS a known good affine transformation (H) but
    F looks random in the original coordinates.

    We apply a random GL(n) matrix to a structured function, creating
    a challenging test for the search algorithm.
    """
    rng = random.Random(seed + 7777)
    # Step 1: Create a structured G over n vars with only few terms
    # Using a quadratic bent-like pattern in first inner_n coordinates
    g_terms: dict[int, int] = {}
    # deg 1: some linear terms
    for i in range(inner_n):
        if rng.random() < 0.4:
            g_terms[1 << i] = 1
    # deg 2: half of the pairs in first inner_n vars
    for i in range(inner_n):
        for j in range(i + 1, inner_n):
            if rng.random() < 0.5:
                g_terms[(1 << i) | (1 << j)] = 1
    # Also add a few terms involving vars beyond inner_n
    for i in range(inner_n, min(n, inner_n + 4)):
        if rng.random() < 0.3:
            g_terms[1 << i] = 1

    if not g_terms:
        g_terms[3] = 1  # x0x1 fallback

    G = SparseANF(g_terms.copy(), n)

    # Step 2: Hide it behind a random invertible transformation
    # F(x) = G(Hx) where H is a random invertible matrix
    H = np.eye(n, dtype=np.uint8)
    while _gf2_rank(H) < n:
        H = np.random.randint(0, 2, size=(n, n)).astype(np.uint8)

    # Compute F(x) = G(Hx)
    # x = H⁻¹y, but we want F(x) = G(Hx), so substitute y = Hx
    # F(x) = sum of monomials of G evaluated at y = Hx
    # Each monomial y_i = sum_j H[i][j] x_j
    # Expand: ∏_{i in mask} (sum_j H[i][j] x_j)
    f_terms: dict[int, int] = {}
    for g_mask, _ in G.terms.items():
        expanded = SparseANF._expand_via_matrix(g_mask, H, n)
        for m, v in expanded.items():
            f_terms[m] = f_terms.get(m, 0) ^ v
    f_terms = {k: v for k, v in f_terms.items() if v}

    return SparseANF(f_terms, n)
    """F(x) = x^e over F_{2^k} where n = k.
    Only for small n (≤ 16). Represents field multiplication as ANF."""
    k = n
    if k > 16:
        raise ValueError(f"Power function ANF only feasible for n ≤ 16, got {n}")

    # Use a primitive polynomial for F_{2^k}
    # Pre-computed primitive polys
    primitives = {
        4: 0b10011,    # x^4 + x + 1
        5: 0b100101,   # x^5 + x^2 + 1
        6: 0b1000011,  # x^6 + x + 1
        8: 0b100011011,  # x^8 + x^4 + x^3 + x^2 + 1
        12: 0b1000000010011,
        16: 0b10000000000001011,
    }
    if k not in primitives:
        # fallback: try a sparse poly
        primitives[k] = (1 << k) | (1 << (k // 2)) | 1

    poly = primitives.get(k, (1 << k) | 1)

    # Compute truth table of x → x^e
    # Represent field elements as integers 0..2^k-1
    size = 1 << k
    tt = np.zeros(size, dtype=np.uint8)

    # Build exp and log tables for fast multiplication
    # Only works if poly is primitive (generates field)
    # This is a simplified version using direct field multiplication
    def gf_mult(a: int, b: int) -> int:
        """Multiply in F_{2^k} with given primitive poly."""
        result = 0
        for i in range(k):
            if (b >> i) & 1:
                result ^= a << i
        # reduce
        for i in range(2 * k - 2, k - 1, -1):
            if (result >> i) & 1:
                result ^= poly << (i - k)
        return result & ((1 << k) - 1)

    for x in range(size):
        result = 1
        exp = e
        base = x
        while exp:
            if exp & 1:
                result = gf_mult(result, base)
            exp >>= 1
            base = gf_mult(base, base)
        tt[result] = 1  # set bit at position = x^e

    # The parity function x → x^e gives truth table of each output bit?
    # Wait, x^e is an element, we need a Boolean function.
    # Let's take f(x) = 0th bit of x^e (a Boolean function)

    # Actually, for a power function x → x^e over F_{2^k},
    # each output bit is a Boolean function of the input bits.
    # We'll use the least significant bit.
    # But this isn't quite right — let me reconsider.

    # x^e maps field element → field element. Each output bit
    # is a Boolean function of k input bits (the bits of x).
    # Let's compute the ANF of the LSB of x^e.

    # Actually, let's compute the algebraic normal form directly.
    # The function f(x) = LSB(x^e) is a Boolean function on k vars.
    # We'll compute its truth table by trying all x.

    # Re-do this properly:
    tt_func = np.zeros(size, dtype=np.uint8)
    for x in range(size):
        result = 1
        exp = e
        base = x
        while exp:
            if exp & 1:
                result = gf_mult(result, base)
            exp >>= 1
            base = gf_mult(base, base)
        tt_func[x] = (result >> 0) & 1  # LSB of x^e

    anf_vec = truth_table_to_anf(tt_func)
    return SparseANF.from_anf_vec(anf_vec)


def run_benchmark():
    """Run full medium-scale benchmark."""
    print("=" * 66)
    print("  MEDIUM-SCALE ANF SIMPLIFICATION BENCHMARK")
    print("=" * 66)

    cases = [
        # (name, SparseANF factory, (args))
    ]

    # n = 8, 12, 16: random cubic + power function
    for n in [8, 12, 16]:
        cases.append((f"Cubic n={n}", SparseANF.random_cubic, (n, 0.20, 100 + n)))
        if n <= 12:
            cases.append((f"PowerF3 n={n}", make_power_function, (n, 3)))
            cases.append((f"PowerF5 n={n}", make_power_function, (n, 5)))

    # n = 24: sparse cubic (lower density)
    cases.append((f"Cubic n=24", SparseANF.random_cubic, (24, 0.08, 200)))

    # n = 32, 48, 64: linear structure + hidden structure
    for n, inn in [(32, 5), (48, 5), (64, 5)]:
        cases.append((f"LinStr n={n}/k={inn}", SparseANF.from_linear_structure, (n, inn, 300 + n)))
    # Hidden structure: function is sparse in some coordinate system
    for n in [32, 64]:
        cases.append((f"Hidden n={n}", make_hidden_structure, (n, 8, 900 + n)))

    # timeout per case (seconds)
    TIMEOUT = 120

    results = []

    for name, maker, args in cases:
        print(f"\n{'─' * 66}")
        print(f"  [{name}]")

        t_gen = time.time()
        f = maker(*args)
        t_gen = time.time() - t_gen

        if f.T() == 0 or f.T() > 50000:
            print(f"  SKIP: T={f.T()} (empty or too large)")
            continue

        print(f"  gen: {t_gen:.2f}s  {f}")

        # show degree breakdown
        bd = f.monomials_by_degree()
        deg_info = ", ".join(f"d{d}:{len(bd[d])}" for d in sorted(bd))
        print(f"  terms: {deg_info}")

        # Run minimization with timeout
        import signal

        class TimeoutError(Exception):
            pass

        def handler(signum, frame):
            raise TimeoutError()

        signal.signal(signal.SIGALRM, handler)
        signal.alarm(TIMEOUT)

        t_start = time.time()
        try:
            if f.n <= 12:
                g, M, b = repeated_walsh_minimize(f, max_iter=12, top_k=200, verbose=True)
            elif f.n <= 24:
                g, M, b = repeated_walsh_minimize(f, max_iter=8, top_k=150, verbose=True)
            else:
                g, M, b = repeated_walsh_minimize(f, max_iter=6, top_k=100, verbose=True)
            elapsed = time.time() - t_start
            signal.alarm(0)

            red = (f.T() - g.T()) / f.T() * 100
            deg_ok = g.degree() <= f.degree()
            print(f"  ──────────────────────────────────────")
            print(f"  ✅ T: {f.T()} → {g.T()}  ({red:.1f}%↓)")
            print(f"  ✅ deg: {f.degree()} → {g.degree()}  {'✓' if deg_ok else '⚠'}")
            print(f"  ⏱  {elapsed:.2f}s")

            results.append({
                "name": name,
                "n": f.n,
                "T_in": f.T(),
                "T_out": g.T(),
                "red%": red,
                "deg_in": f.degree(),
                "deg_out": g.degree(),
                "time_s": elapsed,
            })

        except TimeoutError:
            elapsed = time.time() - t_start
            print(f"  ⏱  TIMEOUT after {elapsed:.0f}s")
            results.append({
                "name": name,
                "n": f.n,
                "T_in": f.T(),
                "T_out": -1,
                "red%": 0.0,
                "deg_in": f.degree(),
                "deg_out": -1,
                "time_s": elapsed,
                "note": "TIMEOUT",
            })
        except Exception as e:
            elapsed = time.time() - t_start
            print(f"  ❌ ERROR: {e}")
            import traceback
            traceback.print_exc()
        finally:
            signal.alarm(0)

    # Summary table
    print("\n" + "=" * 66)
    print("  SUMMARY")
    print("=" * 66)
    hdr = f"  {'Name':<20} {'n':<4} {'T_in':<7} {'T_out':<7} {'Red%':<8} {'deg':<10} {'Time':<8}"
    print(hdr)
    print("  " + "-" * 66)
    for r in results:
        note = r.get("note", "")
        tout = r.get("T_out", -1)
        if tout < 0:
            t_str = "N/A"
            deg_str = f"{r['deg_in']}→N/A"
        else:
            t_str = str(tout)
            deg_str = f"{r['deg_in']}→{r['deg_out']}{'✓' if r['deg_out'] <= r['deg_in'] else '⚠'}"
        to = f" ({r['red%']:.1f}%)" if r['T_out'] > 0 else ""
        print(f"  {r['name']:<20} {r['n']:<4} {r['T_in']:<7} {t_str:<7} {to:<8} {deg_str:<10} {r['time_s']:<7.2f}s{(' ' + note) if note else ''}")
    print()


# ====================================================================
#  Main
# ====================================================================

if __name__ == "__main__":
    bench_transform_time()
    print()
    verify_n8()
    print()
    run_benchmark()
