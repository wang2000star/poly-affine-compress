"""
Vector integer multivariate polynomial simplification.

F: Zⁿ → Zᵏ,  find M,b such that the number of DISTINCT monomials
across all g_i(z) where z = Mx + b is minimized.
"""

from __future__ import annotations

import argparse
import random
import time
from typing import Callable

import numpy as np

from int_poly_factor import (
    IntPoly,
    simplify_int,
    greedy_merge_simplify_int,
    simplify_by_gradient_int,
    _try_merge,
)


class VectorIntPoly:
    """Vector integer polynomial: k components, each an IntPoly over n variables."""

    def __init__(self, components: list[IntPoly] | None = None):
        if components is None:
            self.components: list[IntPoly] = []
        else:
            self.components = components
        self.k = len(self.components)
        self.n = self.components[0].n if self.components else 0

    def union_T(self) -> int:
        """Number of distinct exponent tuples across all components."""
        seen: set[tuple[int, ...]] = set()
        for f in self.components:
            seen.update(f.terms.keys())
        return len(seen)

    def copy(self) -> "VectorIntPoly":
        return VectorIntPoly([f.copy() for f in self.components])

    def substitute_affine(self, M: np.ndarray, b: np.ndarray) -> "VectorIntPoly":
        """Apply z = Mx + b to all components."""
        new_comps = [f.substitute_affine(M, b) for f in self.components]
        return VectorIntPoly(new_comps)

    def variables_used(self) -> set[int]:
        """Union of variables used across all components."""
        used: set[int] = set()
        for f in self.components:
            used |= f.variables_used()
        return used

    def __repr__(self) -> str:
        return f"VectorIntPoly(n={self.n}, k={self.k}, T_union={self.union_T()})"


# ====================================================================
#  Greedy merge for vector integer polynomial
# ====================================================================


def _vector_drop_unused_variables(
    cur: VectorIntPoly, M: np.ndarray, b: np.ndarray
) -> tuple[VectorIntPoly, np.ndarray, np.ndarray]:
    """Drop variables not used in ANY component, compacting indices."""
    used = cur.variables_used()
    if len(used) == cur.n:
        return cur, M, b

    sorted_used = sorted(used)

    # Compact each component: remap exponent indices
    new_comps = []
    for f in cur.components:
        new_terms = {}
        for exp, coeff in f.terms.items():
            new_exp = tuple(exp[i] for i in sorted_used)
            new_terms[new_exp] = coeff
        new_comps.append(IntPoly(new_terms, len(sorted_used)))

    # Compact M and b: keep only used rows
    new_M = M[np.array(sorted_used)].copy()
    new_b = b[np.array(sorted_used)].copy()

    return VectorIntPoly(new_comps), new_M, new_b


def _vector_try_merge_int(vec: VectorIntPoly, i: int, j: int, k: int) -> tuple[VectorIntPoly, int]:
    """Try merge x_i → x_i + k·x_j.  Returns (new_vec, union_T)."""
    new_comps = []
    for f in vec.components:
        new_comps.append(_try_merge(f, i, j, k))
    new_vec = VectorIntPoly(new_comps)
    return new_vec, new_vec.union_T()


def vector_greedy_merge_int(
    vec: VectorIntPoly,
    max_iter: int = 50,
    k_values: list[int] | None = None,
    verbose: bool = False,
) -> tuple[VectorIntPoly, np.ndarray, np.ndarray]:
    """Greedy merge for vector integer polynomial.

    At each step: try all active (i,j) pairs and k values,
    pick the one that reduces union_T the most.
    Apply merge and drop unused variables.
    """
    if k_values is None:
        k_values = [1, -1]

    if vec.union_T() == 0:
        M = np.eye(vec.n, dtype=np.int64)
        b = np.zeros((vec.n, 1), dtype=np.int64)
        return vec, M, b

    cur = vec.copy()
    M = np.eye(vec.n, dtype=np.int64)
    b = np.zeros((vec.n, 1), dtype=np.int64)
    orig_T = cur.union_T()

    if verbose:
        print(f"  Vector greedy merge: n={vec.n}, k={vec.k}, T₀={orig_T}")

    for iteration in range(max_iter):
        active = cur.variables_used()
        active_list = sorted(active)

        if len(active_list) <= 1:
            break

        best_T = cur.union_T()
        best_i = best_j = -1
        best_k = 0
        best_vec = None

        for i in active_list:
            for j in active_list:
                if i == j:
                    continue
                for k in k_values:
                    new_vec, T = _vector_try_merge_int(cur, i, j, k)
                    if T < best_T:
                        best_T = T
                        best_i, best_j, best_k = i, j, k
                        best_vec = new_vec

        if best_vec is None:
            break

        # Apply merge
        M[best_i] = M[best_i] + best_k * M[best_j]
        cur = best_vec

        # Drop unused variables
        cur, M, b = _vector_drop_unused_variables(cur, M, b)

        if verbose:
            pct = (orig_T - cur.union_T()) / orig_T * 100
            print(f"    iter {iteration}: x_{best_i}→x_{best_i}+{best_k}x_{best_j}  T={cur.union_T()}/{orig_T} ({pct:.1f}%↓)  m={cur.n}")

        if cur.union_T() <= 1:
            break

    return cur, M, b


# ====================================================================
#  Combined vector integer polynomial simplification
# ====================================================================


def vector_simplify_int(
    vec: VectorIntPoly,
    verbose: bool = False,
) -> tuple[VectorIntPoly, np.ndarray, np.ndarray]:
    """Combined vector integer polynomial simplification pipeline.

    1. Greedy merge
    2. Gradient-guided merge (TODO)

    Returns (g_vec, M, b) where f_i(x) = g_i(Mx + b) for all i.
    """
    n = vec.n
    M_acc = np.eye(n, dtype=np.int64)
    b_acc = np.zeros((n, 1), dtype=np.int64)
    cur = vec.copy()

    if verbose:
        print(f"Vector simplify: n={n}, k={vec.k}, T_union={vec.union_T()}")

    # Greedy merge
    g_merge, M_m, b_m = vector_greedy_merge_int(cur, verbose=verbose)
    if g_merge.union_T() < cur.union_T():
        M_acc = M_m @ M_acc
        b_acc = (M_m @ b_acc + b_m).reshape(-1, 1)
        cur = g_merge
        if verbose:
            print(f"  After merge: T={cur.union_T()}, m={cur.n}")

    return cur, M_acc, b_acc


def vector_simplify_int_per_component(
    vec: VectorIntPoly,
    verbose: bool = False,
    return_union: bool = False,
    seed: int = 0,
) -> list[tuple[IntPoly, np.ndarray, np.ndarray]] | tuple[list[tuple[IntPoly, np.ndarray, np.ndarray]], VectorIntPoly, np.ndarray, np.ndarray]:
    """Strategy 2: simplify each component independently.

    If return_union=True, also returns a shared representation:
    all distinct rows across all M_i collected into one M.

    Returns:
    - list of (g_i, M_i, b_i) per component
    - if return_union: also (joint_vec, M_joint, b_joint)
    """
    results = []
    for idx, f in enumerate(vec.components):
        if verbose:
            print(f"  Component {idx}: n={f.n}, T={f.T()}")
        g, M, b = simplify_int(f, verbose=verbose)
        results.append((g, M, b))
        if verbose:
            print(f"    → T={g.T()}, m={g.n}")

    if not return_union:
        return results

    # Build shared M containing distinct rows from all M_i
    all_rows: list[np.ndarray] = []
    row_to_idx: dict[tuple, int] = {}
    comp_maps: list[list[int]] = []

    for g, M_i, b_i in results:
        cmap = []
        for r in range(M_i.shape[0]):
            row = tuple(M_i[r].tolist())
            if row not in row_to_idx:
                row_to_idx[row] = len(all_rows)
                all_rows.append(M_i[r])
            cmap.append(row_to_idx[row])
        comp_maps.append(cmap)

    if not all_rows:
        empty_vec = VectorIntPoly([IntPoly({}, 0)] * vec.k)
        return results, empty_vec, np.zeros((0, vec.n), dtype=np.int64), np.zeros((0, 1), dtype=np.int64)

    M_joint = np.array(all_rows, dtype=np.int64)
    m_total = M_joint.shape[0]

    b_joint = np.zeros((m_total, 1), dtype=np.int64)
    for (g, M_i, b_i), cmap in zip(results, comp_maps):
        for local_idx, global_idx in enumerate(cmap):
            b_joint[global_idx, 0] = b_i[local_idx, 0]

    # Remap each g_i to the shared variable space
    joint_comps = []
    for (g, M_i, b_i), cmap in zip(results, comp_maps):
        remapped_terms: dict[tuple[int, ...], int] = {}
        for exp, coeff in g.terms.items():
            new_exp = [0] * m_total
            for var_idx, e in enumerate(exp):
                if e:
                    new_exp[cmap[var_idx]] += e
            key = tuple(new_exp)
            remapped_terms[key] = remapped_terms.get(key, 0) + coeff
            if remapped_terms[key] == 0:
                del remapped_terms[key]
        joint_comps.append(IntPoly(remapped_terms, m_total))

    joint_vec = VectorIntPoly(joint_comps)
    return results, joint_vec, M_joint, b_joint


# ====================================================================
#  Demo
# ====================================================================


def make_vector_from_hidden_structure(
    n: int,
    inner_n: int,
    k: int,
    max_deg: int = 3,
    density: float = 0.25,
    seed: int = 0,
) -> VectorIntPoly:
    """Generate k-component vector integer polynomial.

    Each component: f_i(x) = g_i(A_i x) where g_i is sparse in inner_n variables.
    Different components may share some structure.
    """
    rng = random.Random(seed)
    comps = []
    for comp_idx in range(k):
        # Build random A matrix (inner_n × n)
        A = np.zeros((inner_n, n), dtype=np.int64)
        for i in range(inner_n):
            n_used = rng.randint(1, min(3, n))
            used = rng.sample(range(n), n_used)
            for v in used:
                A[i, v] = rng.randint(-2, 2)
                if A[i, v] == 0:
                    A[i, v] = 1

        # Build sparse random g
        g_terms: dict[tuple[int, ...], int] = {}
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
        from benchmark_int import _expand_via_matrix
        f = _expand_via_matrix(g_ideal, A, n)
        comps.append(f)

    return VectorIntPoly(comps)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vector integer polynomial simplification")
    parser.add_argument("--seed", type=int, default=0, help="random seed (default: 0)")
    parser.add_argument("--n", type=int, default=12, help="number of variables (default: 12)")
    parser.add_argument("--inner", type=int, default=4, help="inner dimension (default: 4)")
    parser.add_argument("--k", type=int, default=3, help="number of components (default: 3)")
    args = parser.parse_args()

    print("=" * 70)
    print("  VECTOR INTEGER POLYNOMIAL SIMPLIFICATION")
    print("=" * 70)

    vec = make_vector_from_hidden_structure(args.n, args.inner, args.k, seed=args.seed)
    print(f"\n  n={vec.n}, k={vec.k}")
    for i, f in enumerate(vec.components):
        print(f"  f_{i}: T={f.T()}, deg={f.degree()}")
    print(f"  UNION T: {vec.union_T()}")

    # Strategy 1: shared transformation
    print(f"\n{'='*70}")
    print("  STRATEGY 1: Shared transformation")
    print("=" * 70)
    t0 = time.time()
    g_vec, M, b = vector_simplify_int(vec, verbose=True)
    vec_elapsed = time.time() - t0
    print(f"\n  Result: T_union: {vec.union_T()} → {g_vec.union_T()}  ({(vec.union_T()-g_vec.union_T())/vec.union_T()*100:.1f}%↓)")
    print(f"  n={vec.n} → {g_vec.n},  time={vec_elapsed:.2f}s")
    for i, g in enumerate(g_vec.components):
        print(f"  g_{i}: T={g.T()} (individual)")

    # Verify
    rng_check = random.Random(999)
    errors = 0
    for _ in range(50):
        x = [rng_check.randint(-3, 3) for _ in range(vec.n)]
        z = M @ np.array(x, dtype=np.int64) + b.ravel()
        for i in range(vec.k):
            if vec.components[i].eval(x) != g_vec.components[i].eval(z.tolist()):
                errors += 1
    print(f"  Verify: {errors}/{50*vec.k} errors")

    # Strategy 2: per-component
    print(f"\n{'='*70}")
    print("  STRATEGY 2: Per-component simplification")
    print("=" * 70)
    per_comp_results, joint_vec, M_joint, b_joint = vector_simplify_int_per_component(
        vec, verbose=True, return_union=True, seed=args.seed)
    print()
    simple_sum = 0
    for idx, (g, M_i, b_i) in enumerate(per_comp_results):
        print(f"  g_{idx}: T={g.T()}, m={g.n}")
        simple_sum += g.T()
    print(f"\n  Metric 1 (sum T_i): {simple_sum}")
    print(f"  Metric 2 (union T after merging rows): {joint_vec.union_T()}  (M: {M_joint.shape[0]}×{M_joint.shape[1]})")
