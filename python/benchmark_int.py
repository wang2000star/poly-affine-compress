"""
Integer polynomial simplification benchmark.

Tests structured functions where f(x) = g(Mx) with g sparse.
Generates proper expanded forms and measures compression.
"""

from __future__ import annotations

import math
import random
import time
from typing import Callable

import numpy as np

from int_poly_factor import (
    IntPoly,
    simplify_int,
    greedy_merge_simplify_int,
    simplify_by_gradient_int,
    _gradient_match,
)


# ====================================================================
#  Test generators
# ====================================================================


def make_single_form(
    coeffs: list[int], degrees: list[int], n: int, seed: int = 0
) -> tuple[IntPoly, IntPoly, np.ndarray]:
    """f(x) = Σ g_j(L_j(x)) where L_j are linear forms and g_j(t) = t^{d_j}.

    Generates expanded polynomial in n variables.

    Returns (f, g_ideal, M_ideal) where f(x) = g_ideal(M_ideal · x).
    g_ideal is IntPoly in m = len(coeffs) variables (one per form).
    M_ideal is m×n matrix.  b=0.
    """
    rng = random.Random(seed)
    m = len(coeffs)

    # Build random linear forms (coefficients in [-3, 3])
    M_ideal = np.zeros((m, n), dtype=np.int64)
    for j in range(m):
        # each form uses 2-4 distinct variables
        used_vars = rng.sample(range(n), min(rng.randint(2, 4), n))
        for v in used_vars:
            M_ideal[j, v] = rng.randint(-3, 3)
            if M_ideal[j, v] == 0:
                M_ideal[j, v] = 1

    # Build g_ideal: sum of t^{d} for each form, with coefficient 1
    g_terms: dict[tuple[int, ...], int] = {}
    for j, d in enumerate(degrees):
        exp = tuple(1 if i == j else 0 for i in range(m))
        if d == 0:
            g_terms[tuple([0] * m)] = g_terms.get(tuple([0] * m), 0) + 1
        else:
            # Represent t^d as a monomial in the j-th variable
            # (g_ideal will be further expanded when substituted)
            # For now: just use t^d → g_ideal will have x_j^d
            g_terms[tuple(d if i == j else 0 for i in range(m))] = 1

    g_ideal = IntPoly(g_terms, m)

    # Build f by expanding g_ideal(M_ideal · x)
    f = _expand_via_matrix(g_ideal, M_ideal, n)

    return f, g_ideal, M_ideal


def _expand_via_matrix(g: IntPoly, M: np.ndarray, n: int) -> IntPoly:
    """Compute f(x) = g(M · x) by expanding each term."""
    m = M.shape[0]
    f_terms: dict[tuple[int, ...], int] = {}

    for g_exp, g_coeff in g.terms.items():
        # Term: Π_j (M[j]·x)^{g_exp[j]}
        # Each factor (M[j]·x)^{e} expands via multinomial
        factors: list[dict[tuple[int, ...], int]] = []
        for j, e in enumerate(g_exp):
            if e == 0:
                continue
            form_coeffs = [int(M[j, i]) for i in range(n)]
            factor_terms = _multinomial_expand(form_coeffs, e, n)
            factors.append(factor_terms)

        if not factors:
            # constant term
            key = tuple([0] * n)
            f_terms[key] = f_terms.get(key, 0) + g_coeff
            continue

        # multiply factors
        result: dict[tuple[int, ...], int] = {tuple([0] * n): g_coeff}
        for factor in factors:
            new = {}
            for e1, c1 in result.items():
                for e2, c2 in factor.items():
                    e = tuple(e1[i] + e2[i] for i in range(n))
                    new[e] = new.get(e, 0) + c1 * c2
                    if new[e] == 0:
                        del new[e]
            result = new

        for e, c in result.items():
            f_terms[e] = f_terms.get(e, 0) + c
            if f_terms[e] == 0:
                del f_terms[e]

    return IntPoly(f_terms, n)


def _multinomial_expand(
    coeffs: list[int], exp: int, n: int
) -> dict[tuple[int, ...], int]:
    """Expand (Σ coeffs[i]·x_i)^exp into dict of exponent → coefficient."""
    if exp == 0:
        return {tuple([0] * n): 1}
    if exp == 1:
        terms = {}
        for i, c in enumerate(coeffs):
            if c != 0:
                e = tuple([1 if j == i else 0 for j in range(n)])
                terms[e] = c
        return terms

    # Multinomial expansion via recursion / combinatorics
    terms: dict[tuple[int, ...], int] = {}

    # Generate all exponent tuples summing to exp
    def _rec(idx: int, remaining: int, current: list[int], coeff: int):
        if idx == n - 1:
            current[idx] = remaining
            # compute multinomial coefficient: exp! / Π current[i]!
            mc = math.factorial(exp)
            for ei in current:
                mc //= math.factorial(ei)
            # multiply by Π coeffs[i]^{current[i]}
            for i, ei in enumerate(current):
                if ei > 0:
                    mc *= coeffs[i] ** ei
            if mc != 0:
                terms[tuple(current)] = mc
            current[idx] = 0
            return
        for ei in range(remaining + 1):
            current[idx] = ei
            _rec(idx + 1, remaining - ei, current, coeff)
        current[idx] = 0

    _rec(0, exp, [0] * n, 1)
    return terms


def make_hidden_structure_int(
    n: int,
    inner_n: int,
    max_deg: int = 4,
    density: float = 0.3,
    seed: int = 0,
) -> tuple[IntPoly, IntPoly, np.ndarray]:
    """f(x) = g(Ax) where g is sparse in inner_n variables.

    A is a random full-rank integer matrix (inner_n × n).
    g has ~density * C(inner_n+max_deg, inner_n) random terms.

    Returns (f, g_ideal, A).
    """
    rng = random.Random(seed)

    # Build random full-rank A matrix
    A = np.zeros((inner_n, n), dtype=np.int64)
    # Each inner variable uses 1-3 distinct x variables
    for i in range(inner_n):
        n_used = rng.randint(1, min(3, n))
        used = rng.sample(range(n), n_used)
        for v in used:
            A[i, v] = rng.randint(-2, 2)
            if A[i, v] == 0:
                A[i, v] = 1

    # Build sparse random g
    g_terms: dict[tuple[int, ...], int] = {}
    # Generate possible exponent tuples up to max_deg
    all_exps: list[tuple[int, ...]] = []

    def _gen_exps(idx: int, remaining: int, current: list[int]):
        if idx == inner_n - 1:
            current[idx] = remaining
            all_exps.append(tuple(current))
            current[idx] = 0
            return
        for ei in range(remaining + 1):
            current[idx] = ei
            _gen_exps(idx + 1, remaining - ei, current)
        current[idx] = 0

    for d in range(max_deg + 1):
        _gen_exps(0, d, [0] * inner_n)

    for exp in all_exps:
        if rng.random() < density:
            coeff = rng.randint(-3, 3)
            if coeff != 0:
                g_terms[exp] = coeff

    if not g_terms:
        g_terms[tuple([1] + [0] * (inner_n - 1))] = 1

    g_ideal = IntPoly(g_terms, inner_n)

    # Expand: f = g(A · x)
    f = _expand_via_matrix(g_ideal, A, n)

    return f, g_ideal, A


def make_multi_form(
    n_forms: int,
    n_vars_per_form: int,
    deg_per_form: list[int],
    n: int,
    overlap: bool = False,
    seed: int = 0,
) -> tuple[IntPoly, IntPoly, np.ndarray]:
    """Sum of several independent linear forms raised to powers.

    f(x) = Σ_i (L_i(x))^{d_i}

    Each L_i uses n_vars_per_form distinct variables.
    If overlap=False, variable sets are disjoint.
    """
    rng = random.Random(seed)
    m = n_forms

    M_ideal = np.zeros((m, n), dtype=np.int64)
    all_vars = list(range(n))
    rng.shuffle(all_vars)
    vars_per = min(n_vars_per_form, n // m)

    for i in range(m):
        if overlap:
            used = rng.sample(range(n), vars_per)
        else:
            start = i * vars_per
            used = all_vars[start:start + vars_per]
        for v in used:
            M_ideal[i, v] = rng.randint(-2, 2)
            if M_ideal[i, v] == 0:
                M_ideal[i, v] = 1

    # Build g_ideal
    g_terms: dict[tuple[int, ...], int] = {}
    for i, d in enumerate(deg_per_form):
        exp = tuple(d if j == i else 0 for j in range(m))
        g_terms[exp] = 1

    g_ideal = IntPoly(g_terms, m)
    f = _expand_via_matrix(g_ideal, M_ideal, n)

    return f, g_ideal, M_ideal


# ====================================================================
#  Verification
# ====================================================================


def verify(f: IntPoly, g: IntPoly, M: np.ndarray, b: np.ndarray, n_tests: int = 50) -> int:
    """Check f(x) == g(Mx + b) for random inputs."""
    rng = random.Random(12345)
    errors = 0
    for _ in range(n_tests):
        x = [rng.randint(-5, 5) for _ in range(f.n)]
        z = [sum(M[j, i] * x[i] for i in range(f.n)) + b[j, 0] for j in range(M.shape[0])]
        if f.eval(x) != g.eval(z):
            errors += 1
    return errors


# ====================================================================
#  Method comparison
# ====================================================================


def compare_methods(
    f: IntPoly, name: str, results: list[dict]
):
    """Run all methods and compare."""
    print(f"  {name}")
    print(f"    f: n={f.n}, T={f.T()}, deg={f.degree()}")

    methods = [
        ("Gradient-only", lambda: simplify_by_gradient_int(f, verbose=False)),
        ("Greedy-only", lambda: greedy_merge_simplify_int(f, max_iter=50, verbose=False)),
        ("Pipeline", lambda: simplify_int(f, verbose=False)),
    ]

    for mname, mfunc in methods:
        t0 = time.time()
        g, M, b = mfunc()
        elapsed = time.time() - t0
        err = verify(f, g, M, b, 30)
        red = (f.T() - g.T()) / f.T() * 100 if f.T() > 0 else 0

        status = "✓" if err == 0 else f"✗({err})"
        label = f"    {mname:<16} T: {f.T()}→{g.T():<4} ({red:+.1f}%)  m: {f.n}→{g.n:<2}  {elapsed:>7.3f}s  {status}"
        print(label)

        results.append({
            "name": name,
            "method": mname.strip(),
            "n": f.n,
            "m_out": g.n,
            "T_in": f.T(),
            "T_out": g.T(),
            "red%": red,
            "deg": f.degree(),
            "time_s": elapsed,
            "errors": err,
        })

    print()


# ====================================================================
#  Main benchmark
# ====================================================================


def run_benchmark():
    print("=" * 70)
    print("  INTEGER POLYNOMIAL SIMPLIFICATION BENCHMARK")
    print("=" * 70)
    results: list[dict] = []

    # ── Small: single forms ──
    print("\n── Single linear form ──")
    # (x₀ + 2x₁ - x₂)⁵  → should reduce to 1 variable
    f, g_ideal, M_ideal = make_single_form([1], [5], n=8, seed=10)
    compare_methods(f, "Form⁵ n=8", results)

    # (2x₀ - 3x₁ + x₂ + x₃)⁴ + (x₄ + 2x₅ - x₆)³  → should reduce to 2 vars
    f, g_ideal, M_ideal = make_single_form([1, 1], [4, 3], n=10, seed=20)
    compare_methods(f, "TwoForm⁴⁺³ n=10", results)

    # ── Medium: multi-form ──
    print("\n── Multi-form (disjoint) ──")
    # 3 forms: (L₁)⁴ + (L₂)³ + (L₃)² with disjoint vars
    f, _, _ = make_multi_form(3, 3, [4, 3, 2], n=12, overlap=False, seed=30)
    compare_methods(f, "3Form⁴⁺³⁺² n=12", results)

    # 4 forms: (L₁)³ + ... + (L₄)³
    f, _, _ = make_multi_form(4, 2, [3, 3, 3, 3], n=12, overlap=False, seed=31)
    compare_methods(f, "4Form³ n=12", results)

    # ── Hidden structure ──
    print("\n── Hidden low-dim structure ──")
    for n, inner_n, deg, density, seed in [
        (15, 4, 3, 0.25, 100),
        (25, 5, 3, 0.20, 101),
    ]:
        f, g_ideal, A = make_hidden_structure_int(n, inner_n, deg, density, seed)
        label = f"Hidden n={n}/m={inner_n} d≤{deg}"
        compare_methods(f, label, results)

    # ── Larger: scale test ──
    print("\n── Scale tests ──")
    for n, inner_n, deg in [
        (40, 5, 3),
        (60, 5, 3),
        (80, 5, 3),
    ]:
        f, g_ideal, A = make_hidden_structure_int(n, inner_n, deg, density=0.15, seed=200 + n)
        label = f"Hidden n={n}/m={inner_n} d≤{deg}"
        TIMEOUT = 60
        t0 = time.time()
        # Only run pipeline
        t1 = time.time()
        g, M, b = simplify_int(f, verbose=False)
        elapsed = time.time() - t1
        err = verify(f, g, M, b, 20)
        red = (f.T() - g.T()) / f.T() * 100 if f.T() > 0 else 0
        total = time.time() - t0
        status = "✓" if err == 0 else f"✗({err})"
        print(f"  {label:<25} T: {f.T()}→{g.T():<4} ({red:+.1f}%)  m: {f.n}→{g.n:<2}  {elapsed:>7.3f}s  {status}")
        results.append({
            "name": label, "method": "Pipeline", "n": f.n, "m_out": g.n,
            "T_in": f.T(), "T_out": g.T(), "red%": red, "deg": f.degree(),
            "time_s": elapsed, "errors": err,
        })

    # ── Summary ──
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    hdr = f"  {'Name':<25} {'Method':<14} {'T_in':<6} {'T_out':<6} {'Red%':<8} {'n→m':<8} {'Time':<8} {'Err':<5}"
    print(hdr)
    print("  " + "-" * 70)
    for r in results:
        red_s = f"{r['red%']:+.1f}%" if r['red%'] <= 0 else f"+{r['red%']:.1f}%"
        print(f"  {r['name']:<25} {r['method']:<14} {r['T_in']:<6} {r['T_out']:<6} {red_s:<8} {r['n']}→{r['m_out']:<3} {r['time_s']:<7.3f}s {r['errors']:<5}")


# ====================================================================
#  Main
# ====================================================================


if __name__ == "__main__":
    run_benchmark()
