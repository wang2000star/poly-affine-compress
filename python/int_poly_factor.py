"""
Sparse integer multivariate polynomial simplification via affine substitution.

Find M ∈ Z^{m×n}, b ∈ Z^m such that f(x) = g(Mx + b) has fewer terms than f.
Core idea: use gradient (directional derivative) structure to detect linear forms.
"""

from __future__ import annotations

import itertools
import math
import random
from typing import Optional

import numpy as np


class IntPoly:
    """Sparse multivariate integer polynomial.

    terms: dict[tuple[int, ...], int]  exponent_tuple -> coefficient
    n: number of variables
    """

    def __init__(self, terms: dict[tuple[int, ...], int] | None = None, n: int = 0):
        self.terms = {}
        self.n = n
        if terms:
            for exp, coeff in terms.items():
                if coeff != 0:
                    self.terms[exp] = coeff

    def copy(self) -> "IntPoly":
        return IntPoly(dict(self.terms), self.n)

    def T(self) -> int:
        return len(self.terms)

    def degree(self) -> int:
        return max(sum(e) for e in self.terms) if self.terms else 0

    def deg_by_var(self) -> tuple[int, ...]:
        if not self.terms:
            return (0,) * self.n
        max_deg = [0] * self.n
        for exp in self.terms:
            for i, e in enumerate(exp):
                if e > max_deg[i]:
                    max_deg[i] = e
        return tuple(max_deg)

    def __repr__(self) -> str:
        return f"IntPoly(n={self.n}, T={self.T()}, deg={self.degree()})"

    def __add__(self, other: "IntPoly") -> "IntPoly":
        assert self.n == other.n
        new = dict(self.terms)
        for exp, coeff in other.terms.items():
            new[exp] = new.get(exp, 0) + coeff
            if new[exp] == 0:
                del new[exp]
        return IntPoly(new, self.n)

    def __sub__(self, other: "IntPoly") -> "IntPoly":
        assert self.n == other.n
        new = dict(self.terms)
        for exp, coeff in other.terms.items():
            new[exp] = new.get(exp, 0) - coeff
            if new[exp] == 0:
                del new[exp]
        return IntPoly(new, self.n)

    def __mul__(self, other: "IntPoly") -> "IntPoly":
        assert self.n == other.n
        new: dict[tuple[int, ...], int] = {}
        for e1, c1 in self.terms.items():
            for e2, c2 in other.terms.items():
                e = tuple(e1[i] + e2[i] for i in range(self.n))
                new[e] = new.get(e, 0) + c1 * c2
                if new[e] == 0:
                    del new[e]
        return IntPoly(new, self.n)

    def __rmul__(self, scalar: int) -> "IntPoly":
        new = {e: c * scalar for e, c in self.terms.items() if c * scalar != 0}
        return IntPoly(new, self.n)

    def partial_deriv(self, var: int) -> "IntPoly":
        """∂f/∂x_var"""
        new: dict[tuple[int, ...], int] = {}
        for exp, coeff in self.terms.items():
            e = exp[var]
            if e > 0:
                new_exp = list(exp)
                new_exp[var] = e - 1
                new[tuple(new_exp)] = coeff * e
        return IntPoly(new, self.n)

    def gradient(self) -> list["IntPoly"]:
        return [self.partial_deriv(i) for i in range(self.n)]

    def variables_used(self) -> set[int]:
        used = set()
        for exp in self.terms:
            for i, e in enumerate(exp):
                if e > 0:
                    used.add(i)
        return used

    def substitute_linear(self, var: int, coeffs: list[int]) -> "IntPoly":
        """Replace x_var with Σ coeffs[j] * x_j (no constant term)."""
        new: dict[tuple[int, ...], int] = {}
        for exp, coeff in self.terms.items():
            e = exp[var]
            if e == 0:
                # Copy term unchanged
                new[exp] = new.get(exp, 0) + coeff
                continue

            # For each occurrence of x_var^e, expand (Σ coeffs[j]*x_j)^e
            # as a polynomial in the remaining variables
            # Use multinomial expansion: (Σ a_j x_j)^e = Σ_{t: |t|=e} (e! / Π t_j!) Π (a_j x_j)^{t_j}
            base_poly = IntPoly._linear_form_poly(coeffs, self.n)

            # Expand x_var^e * rest
            rest_exp = list(exp)
            rest_exp[var] = 0
            rest_poly = IntPoly({tuple(rest_exp): coeff}, self.n)

            if e == 1:
                term = base_poly
            else:
                term = base_poly
                for _ in range(e - 1):
                    term = term * base_poly

            result = term * rest_poly
            for e2, c2 in result.terms.items():
                new[e2] = new.get(e2, 0) + c2
                if new[e2] == 0:
                    del new[e2]

        return IntPoly(new, self.n)

    @staticmethod
    def _linear_form_poly(coeffs: list[int], n: int) -> "IntPoly":
        """Build polynomial Σ coeffs[j] * x_j."""
        terms = {}
        for j, c in enumerate(coeffs):
            if c != 0:
                exp = tuple(1 if i == j else 0 for i in range(n))
                terms[exp] = c
        return IntPoly(terms, n)

    def substitute_affine(self, M: np.ndarray, b: np.ndarray) -> "IntPoly":
        """Substitute x = Mz + b, compute g(z).

        M: n×m integer matrix, b: n×1 integer vector.
        Returns g(z) = f(Mz + b) as IntPoly in m variables.
        """
        n, m = M.shape
        assert n == self.n
        assert b.shape == (n, 1)

        # This is a complex operation: replace each x_i with Σ M[i,j] z_j + b_i
        # We'll do it term by term
        result: dict[tuple[int, ...], int] = {tuple([0] * m): 0}

        for exp, coeff in self.terms.items():
            # Expand Π_i (Σ_j M[i,j] z_j + b_i)^{exp[i]}
            # Each factor is an IntPoly in m variables
            factor_polys = []
            for i, e in enumerate(exp):
                if e == 0:
                    continue
                # Build (M[i]·z + b[i])^e
                poly = IntPoly._single_affine(M[i], b[i, 0], m)
                if e > 1:
                    orig = poly
                    for _ in range(e - 1):
                        poly = poly * orig
                factor_polys.append(poly)

            if not factor_polys:
                # constant term
                cur = IntPoly({tuple([0] * m): coeff}, m)
            else:
                cur = factor_polys[0]
                for p in factor_polys[1:]:
                    cur = cur * p
                if coeff != 1:
                    cur = cur.__rmul__(coeff)

            for e2, c2 in cur.terms.items():
                result[e2] = result.get(e2, 0) + c2
                if result[e2] == 0:
                    del result[e2]

        return IntPoly(result, m)

    @staticmethod
    def _single_affine(row: np.ndarray, const: int, m: int) -> "IntPoly":
        """Build polynomial (Σ row[j] * z_j + const) in m variables."""
        terms = {}
        if const != 0:
            terms[tuple([0] * m)] = const
        for j in range(m):
            if row[j] != 0:
                exp = tuple([0] * j + [1] + [0] * (m - j - 1))
                terms[exp] = int(row[j])
        return IntPoly(terms, m)

    def eval(self, values: list[int]) -> int:
        """Evaluate at given integer values."""
        result = 0
        for exp, coeff in self.terms.items():
            term_val = coeff
            for i, e in enumerate(exp):
                term_val *= values[i] ** e
            result += term_val
        return result


# ====================================================================
#  Greedy merge for integer polynomials
# ====================================================================


def _try_merge(f: IntPoly, i: int, j: int, k: int) -> IntPoly:
    """Try merge: M[i] += k*M[j] ⇒ g'(z) = g(z with z_i → z_i - k*z_j)."""
    coeffs = [0] * f.n
    coeffs[i] = 1
    coeffs[j] = -k
    return f.substitute_linear(i, coeffs)


def greedy_merge_simplify_int(
    f: IntPoly,
    max_iter: int = 50,
    verbose: bool = False,
) -> tuple[IntPoly, np.ndarray, np.ndarray]:
    """Greedy pairwise merge for integer polynomials.

    Invariant: MERGE x_i → x_i + k*x_j in the BASIS means:
      M'[i] = M[i] + k*M[j]     (forward)
      g'(z) = g(..., z_i - k*z_j, ...)  (inverse in polynomial)

    Returns (g, M, b) where f(x) = g(Mx + b).
    """
    if not f.terms:
        return IntPoly({}, f.n), np.eye(f.n, dtype=np.int64), np.zeros((f.n, 1), dtype=np.int64)

    cur = f.copy()
    M = np.eye(f.n, dtype=np.int64)
    b = np.zeros((f.n, 1), dtype=np.int64)

    if verbose:
        print(f"\nGreedy merge (int): n={f.n}, T₀={f.T()}")

    for iteration in range(max_iter):
        active = cur.variables_used()
        active_list = sorted(active)

        if len(active_list) <= 1:
            break

        best_delta = 0
        best_i = best_j = -1
        best_k = 0

        for i in active_list:
            for j in active_list:
                if i == j:
                    continue
                for k in (1, -1, 2, -2):
                    g_test = _try_merge(cur, i, j, k)
                    delta = g_test.T() - cur.T()
                    if delta < best_delta:
                        best_delta = delta
                        best_i, best_j, best_k = i, j, k

        if best_delta >= 0:
            break

        i, j, k = best_i, best_j, best_k

        M[i] += k * M[j]
        cur = _try_merge(cur, i, j, k)
        cur, M, b = _drop_unused_variables(cur, M, b)

        if verbose:
            print(f"  iter {iteration}: x_{i}→x_{i}+{k}x_{j}  T={cur.T()}/{f.T()} ({(f.T()-cur.T())/f.T()*100:.1f}%↓)")

        if cur.T() <= 1:
            break

    if verbose:
        print(f"  Final: T={cur.T()}/{f.T()} ({(f.T()-cur.T())/f.T()*100:.1f}%↓)")

    return cur, M, b


def _drop_unused_variables(
    cur: IntPoly, M: np.ndarray, b: np.ndarray
) -> tuple[IntPoly, np.ndarray, np.ndarray]:
    """Drop variables not appearing in cur, compacting indices.

    Returns (compacted_cur, compacted_M, compacted_b) where compacted_cur
    has m ≤ cur.n variables and M/b have m rows.
    """
    used = cur.variables_used()
    if len(used) == cur.n:
        return cur, M, b

    sorted_used = sorted(used)

    # Compact polynomial: remap exponent indices
    new_terms = {}
    for exp, coeff in cur.terms.items():
        new_exp = tuple(exp[i] for i in sorted_used)
        new_terms[new_exp] = coeff

    # Compact M and b: keep only used rows
    new_M = M[sorted_used].copy()
    new_b = b[sorted_used].copy()

    return IntPoly(new_terms, len(sorted_used)), new_M, new_b


def _gradient_match(g_i: IntPoly, g_j: IntPoly, k: int) -> bool:
    """Check if g_i == k * g_j as polynomials."""
    all_terms = set(g_i.terms.keys()) | set(g_j.terms.keys())
    for exp in all_terms:
        if g_i.terms.get(exp, 0) != k * g_j.terms.get(exp, 0):
            return False
    return True


def simplify_by_gradient_int(
    f: IntPoly,
    verbose: bool = False,
) -> tuple[IntPoly, np.ndarray, np.ndarray]:
    """Gradient-guided integer polynomial simplification.

    Idea: if ∂f/∂x_i = k·∂f/∂x_j, then f depends on x_i, x_j only through
    x_i + k·x_j.  This suggests the merge x_i → x_i + k·x_j reduces dimension.

    Checks gradient pairs for k ∈ {1, -1, 2, -2}, applies beneficial merges.
    Falls back to exhaustive pair search when gradient is not informative.

    Returns (g, M, b) where f(x) = g(Mx + b).
    """
    if not f.terms:
        return IntPoly({}, f.n), np.eye(f.n, dtype=np.int64), np.zeros((f.n, 1), dtype=np.int64)

    cur = f.copy()
    M = np.eye(f.n, dtype=np.int64)
    b = np.zeros((f.n, 1), dtype=np.int64)

    if verbose:
        print(f"\nGradient-guided merge: n={f.n}, T₀={f.T()}")

    for iteration in range(30):
        grads = cur.gradient()
        active = cur.variables_used()
        active_list = sorted(active)

        if len(active_list) <= 1:
            break

        best_delta = 0
        best_pair: tuple[int, int, int] | None = None

        # Phase A: gradient pair search
        for i in active_list:
            if grads[i].T() == 0:
                continue
            for j in active_list:
                if i == j or grads[j].T() == 0:
                    continue
                for k in (1, -1, 2, -2):
                    if _gradient_match(grads[i], grads[j], k):
                        g_test = _try_merge(cur, i, j, k)
                        delta = g_test.T() - cur.T()
                        if delta < best_delta:
                            best_delta = delta
                            best_pair = (i, j, k)

        # Phase B: if gradient found nothing, try exhaustive
        if best_pair is None:
            for i in active_list:
                for j in active_list:
                    if i == j:
                        continue
                    for k in (1, -1, 2, -2):
                        g_test = _try_merge(cur, i, j, k)
                        delta = g_test.T() - cur.T()
                        if delta < best_delta:
                            best_delta = delta
                            best_pair = (i, j, k)

        if best_pair is None or best_delta >= 0:
            break

        i, j, k = best_pair
        M[i] += k * M[j]
        cur = _try_merge(cur, i, j, k)
        cur, M, b = _drop_unused_variables(cur, M, b)

        if verbose:
            pct = (f.T() - cur.T()) / f.T() * 100
            print(f"  iter {iteration}: x_{i}→x_{i}+{k}·x_{j}  T={cur.T()}/{f.T()} ({pct:.1f}%↓)  m={cur.n}")

        if cur.T() <= 1:
            break

    if verbose:
        print(f"  Final: T={cur.T()}/{f.T()}, m={cur.n}")

    return cur, M, b


def simplify_int(
    f: IntPoly,
    verbose: bool = False,
) -> tuple[IntPoly, np.ndarray, np.ndarray]:
    """Combined integer polynomial simplification pipeline.

    Phase 1: gradient-guided merge (fast detection of linear dependencies)
    Phase 2: exhaustive pairwise merge on remaining variables

    Returns (g, M, b) where f(x) = g(Mx + b) with g having ≤ f.n variables.
    """
    # Phase 1
    g1, M1, b1 = simplify_by_gradient_int(f, verbose=verbose)

    if g1.T() <= 1:
        return g1, M1, b1

    # Phase 2
    if verbose:
        print(f"\nPhase 2: exhaustive merge on m={g1.n}")
    g2, M2, b2 = greedy_merge_simplify_int(g1, max_iter=20, verbose=verbose)

    # Compose: z2 = M2·z1 + b2 = M2·(M1·x + b1) + b2 = (M2@M1)·x + M2·b1 + b2
    M = M2.astype(np.int64) @ M1.astype(np.int64)
    b = (M2.astype(np.int64) @ b1.astype(np.int64) + b2.astype(np.int64))

    return g2, M, b


# ====================================================================
#  Test / demo
# ====================================================================


def _verify(f: IntPoly, g: IntPoly, M: np.ndarray, b: np.ndarray, n_tests: int = 100) -> int:
    """Verify f(x) == g(Mx + b) for random integer inputs."""
    rng = random.Random(42)
    errors = 0
    for _ in range(n_tests):
        x = [rng.randint(-5, 5) for _ in range(f.n)]
        z = [sum(M[j, i] * x[i] for i in range(f.n)) + b[j, 0] for j in range(M.shape[0])]
        if f.eval(x) != g.eval(z):
            errors += 1
    return errors


def demo_int_poly():
    print("=" * 60)
    print("IntPoly: integer polynomial simplification")
    print("=" * 60)

    # f = (x + y)^3 + (u - v)^2  → has structure
    n = 4
    terms = {}
    # (x0 + x1)^3 = x0^3 + 3x0^2 x1 + 3x0 x1^2 + x1^3
    terms[(3, 0, 0, 0)] = 1
    terms[(2, 1, 0, 0)] = 3
    terms[(1, 2, 0, 0)] = 3
    terms[(0, 3, 0, 0)] = 1
    # (x2 - x3)^2 = x2^2 - 2x2x3 + x3^2
    terms[(0, 0, 2, 0)] = 1
    terms[(0, 0, 1, 1)] = -2
    terms[(0, 0, 0, 2)] = 1

    f = IntPoly(terms, n)
    print(f"f: {f}  T={f.T()}")
    print()

    print("--- simplify_int (combined pipeline) ---")
    g, M, b = simplify_int(f, verbose=True)
    print(f"\ng: {g}  T={g.T()}, m={g.n}")
    err = _verify(f, g, M, b)
    print(f"verify: {err}/100 errors")

    print()


def demo_struct_example():
    """f(x) = (x₁ + 2x₂)³ + (x₃ - x₄)² + x₅  →  m=3"""
    print("=" * 60)
    print("Structured example")
    print("=" * 60)
    n = 5
    terms = {}
    # (x0 + 2x1)^3 = x0^3 + 6x0^2 x1 + 12x0 x1^2 + 8x1^3
    for t in range(4):
        terms[(3 - t, t, 0, 0, 0)] = [1, 6, 12, 8][t]
    # (x2 - x3)^2 = x2^2 - 2x2x3 + x3^2
    terms[(0, 0, 2, 0, 0)] = 1
    terms[(0, 0, 1, 1, 0)] = -2
    terms[(0, 0, 0, 2, 0)] = 1
    # + x4
    terms[(0, 0, 0, 0, 1)] = 1

    f = IntPoly(terms, n)
    print(f"f: {f}  T={f.T()}")

    g, M, b = simplify_int(f, verbose=True)
    print(f"\ng: {g}  T={g.T()}, m={g.n}")
    err = _verify(f, g, M, b)
    print(f"verify: {err}/100 errors")
    print()


def demo_gradient_merge():
    """Show gradient-guided merge finds structure faster."""
    print("=" * 60)
    print("Gradient-guided merge demo")
    print("=" * 60)
    n = 6
    # f = (x0 + x1 + x2)^3 + (x3 - x4 + 2x5)^2
    terms = {}
    # (x0 + x1 + x2)^3 expanded
    for a in range(4):
        for b in range(4 - a):
            c = 3 - a - b
            coeff = math.factorial(3) // (math.factorial(a) * math.factorial(b) * math.factorial(c))
            terms[(a, b, c, 0, 0, 0)] = coeff
    # (x3 - x4 + 2x5)^2 expanded
    for a in range(3):
        for b in range(3 - a):
            c = 2 - a - b
            coeff = math.factorial(2) // (math.factorial(a) * math.factorial(b) * math.factorial(c))
            coeff *= (1 if a == 0 else (-1) ** a) * (2 ** c)
            if coeff != 0:
                terms[(0, 0, 0, a, b, c)] = coeff

    f = IntPoly(terms, n)
    print(f"f: {f}  T={f.T()}")

    # Show gradient structure
    grads = f.gradient()
    print("Gradients:")
    for i, g in enumerate(grads):
        print(f"  ∂f/∂x_{i}: T={g.T()}")
    for i in range(n):
        for j in range(i + 1, n):
            for k in (1, -1, 2, -2):
                if grads[i].T() > 0 and grads[j].T() > 0 and _gradient_match(grads[i], grads[j], k):
                    print(f"  → ∂f/∂x_{i} = {k}·∂f/∂x_{j}")

    print()
    g, M, b = simplify_int(f, verbose=True)
    print(f"\ng: {g}  T={g.T()}, m={g.n}")
    err = _verify(f, g, M, b)
    print(f"verify: {err}/100 errors")
    print()


def random_benchmark():
    """Benchmark on random polynomials of various sizes."""
    print("=" * 60)
    print("Random polynomial benchmark")
    print("=" * 60)
    for n, T_in, seed in [(6, 15, 0), (8, 25, 1), (10, 40, 2)]:
        rng = random.Random(seed)
        terms = {}
        for _ in range(T_in):
            exp = tuple(rng.randint(0, 3) for _ in range(n))
            if sum(exp) == 0:
                exp = (1,) + (0,) * (n - 1)
            coeff = rng.randint(-3, 3)
            if coeff == 0:
                coeff = 1
            terms[exp] = coeff
        f = IntPoly(terms, n)
        g, M, b = simplify_int(f, verbose=False)
        err = _verify(f, g, M, b, 50)
        red = (f.T() - g.T()) / f.T() * 100
        print(f"  n={n} T={f.T():>2}→{g.T():<2}  m={f.n}→{g.n}  {red:+.1f}%  verify:{'OK' if err==0 else err}")
    print()


if __name__ == "__main__":
    demo_int_poly()
    demo_struct_example()
    demo_gradient_merge()
    random_benchmark()
