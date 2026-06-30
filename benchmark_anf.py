"""
Medium-Scale Boolean Function ANF Simplification Benchmark
============================================================

Tests for n = 8, 12, 16, 24, 32, 48, 64 using sparse ANF representations.

Consolidated with anf_factor.py — uses SparseANF and substitute_affine from there.
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
)

from anf_factor import (
    SparseANF,
    gf2_extend_to_invertible,
    search_affine_simplification,
    _walsh_spectrum,
)


# ====================================================================
#  Iterative Walsh-guided minimization
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
        indices = _walsh_spectrum(current, n_top=top_k + 20)
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
                    g = current.substitute_affine(M_step, b_step)
                    T = g.T()
                    if T < best_T:
                        best_T, best_g = T, g
                        best_Ms, best_bs = M_step.copy(), b_step.copy()
                except ValueError:
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

    N = gf2_matrix_inv(M)
    errors = 0
    for _ in range(50):
        x = rng.randint(0, (1 << n) - 1)
        z_arr = (M @ np.array([(x >> i) & 1 for i in range(n)], dtype=np.uint8)) % 2
        z_arr = (z_arr ^ b.ravel()) % 2
        z = sum(int(z_arr[i]) << i for i in range(g.n))
        fx = f.eval_mask(x)
        gz = g.eval_mask(z)
        if fx != gz:
            errors += 1

    status = "PASS" if errors == 0 else f"FAIL ({errors} errors)"
    print(f"  correctness: {status}")
    print()


# ====================================================================
#  Speed benchmarks: transformation timing
# ====================================================================


def bench_transform_time():
    """Measure SparseANF.substitute_affine time for various sizes."""
    print("=" * 66)
    print("  SPEED: SparseANF.substitute_affine() scaling")
    print("=" * 66)
    print(f"  {'n':<4} {'T(f)':<8} {'deg':<5} {'time(ms)':<10}")
    print("  " + "-" * 30)
    for n in [8, 12, 16, 24, 32]:
        f = SparseANF.random_cubic(n, density=0.15, seed=0)
        M = np.eye(n, dtype=np.uint8)
        while _gf2_rank(M) < n:
            M = np.random.randint(0, 2, size=(n, n)).astype(np.uint8)
        b = np.zeros((n, 1), dtype=np.uint8)
        t0 = time.time()
        _ = f.substitute_affine(M, b)
        t = (time.time() - t0) * 1000
        print(f"  {n:<4} {f.T():<8} {f.degree():<5} {t:<9.1f}ms")
    print()


# ====================================================================
#  Test function generators
# ====================================================================


def make_power_function(n: int, e: int) -> SparseANF:
    """F(x) = LSB(x^e) over F_{2^n}.  Only for n ≤ 16."""
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
        tt[x] = (result >> 0) & 1  # LSB of x^e

    anf_vec = truth_table_to_anf(tt)
    return SparseANF.from_anf_vec(anf_vec)


def make_hidden_structure(n: int, inner_n: int, seed: int) -> SparseANF:
    """F(x) = G(Hx) where G is sparse but H 'hides' it."""
    rng = random.Random(seed + 7777)

    g_terms: dict[int, int] = {}
    for i in range(inner_n):
        if rng.random() < 0.4:
            g_terms[1 << i] = 1
    for i in range(inner_n):
        for j in range(i + 1, inner_n):
            if rng.random() < 0.5:
                g_terms[(1 << i) | (1 << j)] = 1
    for i in range(inner_n, min(n, inner_n + 4)):
        if rng.random() < 0.3:
            g_terms[1 << i] = 1
    if not g_terms:
        g_terms[3] = 1

    G = SparseANF(g_terms.copy(), n)

    H = np.eye(n, dtype=np.uint8)
    while _gf2_rank(H) < n:
        H = np.random.randint(0, 2, size=(n, n)).astype(np.uint8)

    f_terms: dict[int, int] = {}
    for g_mask, _ in G.terms.items():
        expanded = SparseANF._expand_via_matrix(g_mask, H, n)
        for m, v in expanded.items():
            f_terms[m] = f_terms.get(m, 0) ^ v
    f_terms = {k: v for k, v in f_terms.items() if v}
    return SparseANF(f_terms, n)


# ====================================================================
#  Full benchmark
# ====================================================================


def run_benchmark():
    """Run full medium-scale benchmark."""
    print("=" * 66)
    print("  MEDIUM-SCALE ANF SIMPLIFICATION BENCHMARK")
    print("=" * 66)

    cases = []

    for n in [8, 12, 16]:
        cases.append((f"Cubic n={n}", SparseANF.random_cubic, (n, 0.20, 100 + n)))
        if n <= 12:
            cases.append((f"PowerF3 n={n}", make_power_function, (n, 3)))
            cases.append((f"PowerF5 n={n}", make_power_function, (n, 5)))

    cases.append((f"Cubic n=24", SparseANF.random_cubic, (24, 0.08, 200)))

    for n, inn in [(32, 5), (48, 5), (64, 5)]:
        cases.append((f"LinStr n={n}/k={inn}", SparseANF.from_linear_structure,
                      (n, inn, 300 + n)))

    for n in [32, 64]:
        cases.append((f"Hidden n={n}", make_hidden_structure, (n, 8, 900 + n)))

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

        bd = f.monomials_by_degree()
        deg_info = ", ".join(f"d{d}:{len(bd[d])}" for d in sorted(bd))
        print(f"  terms: {deg_info}")

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
                g, M, b = repeated_walsh_minimize(f, max_iter=12, top_k=200,
                                                   verbose=True)
            elif f.n <= 24:
                g, M, b = repeated_walsh_minimize(f, max_iter=8, top_k=150,
                                                   verbose=True)
            else:
                g, M, b = repeated_walsh_minimize(f, max_iter=6, top_k=100,
                                                   verbose=True)
            elapsed = time.time() - t_start
            signal.alarm(0)

            red = (f.T() - g.T()) / f.T() * 100
            deg_ok = g.degree() <= f.degree()
            print(f"  ──────────────────────────────────────")
            print(f"  T: {f.T()} → {g.T()}  ({red:.1f}%↓)")
            print(f"  deg: {f.degree()} → {g.degree()}  {'✓' if deg_ok else '⚠'}")
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
