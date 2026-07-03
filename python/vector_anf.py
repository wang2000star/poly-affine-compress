"""
Vector Boolean function ANF simplification.

F: {0,1}ⁿ → {0,1}ᵏ,  find M,b such that the number of DISTINCT monomials
across all g_i(z) where z = Mx⊕b is minimized.
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys
import time
from typing import Callable

import numpy as np

from anf_factor import SparseANF, simplify, greedy_merge_simplify, simplify_by_complement, search_affine_simplification


class VectorANF:
    """Vector Boolean function: k components, each an SparseANF over n variables.

    Goal: find M,b minimizing |Union_i supp(g_i)| where g_i(z) = f_i(Mx⊕b).
    """

    def __init__(self, components: list[SparseANF]):
        assert len(components) > 0
        self.n = components[0].n
        self.k = len(components)
        for f in components:
            assert f.n == self.n
        self.components = [f.copy() for f in components]

    def copy(self) -> "VectorANF":
        return VectorANF([f.copy() for f in self.components])

    def union_T(self) -> int:
        """Number of distinct monomials across all components."""
        all_masks = set()
        for f in self.components:
            all_masks.update(f.terms.keys())
        return len(all_masks)

    def T(self) -> int:
        return self.union_T()

    def union_terms(self) -> dict[int, int]:
        """Union of all component monomials (mask→1)."""
        result = {}
        for f in self.components:
            for m in f.terms:
                result[m] = 1
        return result

    def __repr__(self) -> str:
        return f"VectorANF(n={self.n}, k={self.k}, T_union={self.union_T()})"

    def substitute_affine(self, M: np.ndarray, b: np.ndarray) -> "VectorANF":
        """Apply z = Mx⊕b to all components, return new VectorANF."""
        new_comps = [f.substitute_affine(M, b) for f in self.components]
        return VectorANF(new_comps)

    def substitute_linear(self, var: int, coeffs: list[int]) -> "VectorANF":
        """Replace x_var with Σ coeffs[j]*x_j in all components.

        Boolean ring: for each j ≠ var, x_j·rest = rest with bit j SET (OR),
        because x_j² = x_j  (if x_j already in rest, x_j·rest = rest).
        """
        new_comps = []
        for f in self.components:
            new = {}
            for mask, v in f.terms.items():
                if not ((mask >> var) & 1):
                    new[mask] = new.get(mask, 0) ^ v
                    continue
                rest = mask ^ (1 << var)  # m without var
                for j, c in enumerate(coeffs):
                    if c:
                        if j == var:
                            # x_var·rest = m (original monomial)
                            new[mask] = new.get(mask, 0) ^ v
                        else:
                            # x_j·rest = rest with bit j SET (Boolean idempotence)
                            new_m = rest if ((mask >> j) & 1) else rest | (1 << j)
                            new[new_m] = new.get(new_m, 0) ^ v
            clean = {m: c for m, c in new.items() if c != 0}
            new_comps.append(SparseANF(clean, self.n))
        return VectorANF(new_comps)

    def gradient(self) -> list[list[SparseANF]]:
        """Gradients of each component: grads[comp_idx][var_idx] = ∂f_i/∂x_var."""
        return [f.gradient() for f in self.components]

    def variables_used(self) -> set[int]:
        """Union of variables used across all components."""
        used = set()
        for f in self.components:
            used |= f.variables_used()
        return used


def _vector_try_merge(
    vec: VectorANF, i: int, j: int, k: int = 1
) -> tuple[VectorANF, int]:
    """Try merge x_i → x_i⊕(k·x_j).  k ∈ {1}.  Returns (new_vec, union_T)."""
    coeffs = [0] * vec.n
    coeffs[i] = 1
    coeffs[j] = k
    new_vec = vec.substitute_linear(i, coeffs)
    return new_vec, new_vec.union_T()


def vector_greedy_merge(
    vec: VectorANF,
    max_iter: int = 50,
    verbose: bool = False,
) -> tuple[VectorANF, np.ndarray, np.ndarray]:
    """Greedy XOR merge for vector Boolean function.

    At each step: try all active (i,j) pairs, pick the one that reduces
    union_T the most.  Apply merge and drop unused variables.
    """
    if vec.union_T() == 0:
        M = np.eye(vec.n, dtype=np.uint8)
        b = np.zeros((vec.n, 1), dtype=np.uint8)
        return vec, M, b

    cur = vec.copy()
    M = np.eye(vec.n, dtype=np.uint8)
    b = np.zeros((vec.n, 1), dtype=np.uint8)
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
        best_vec = None

        for i in active_list:
            for j in active_list:
                if i == j:
                    continue
                new_vec, T = _vector_try_merge(cur, i, j, 1)
                if T < best_T:
                    best_T = T
                    best_i, best_j = i, j
                    best_vec = new_vec

        if best_vec is None:
            break

        # Apply merge
        M[best_i] ^= M[best_j]  # XOR since we use k=1
        cur = best_vec

        # Drop unused variables, trim M and b
        cur, M, b = _vector_drop_unused(cur, M, b)

        if verbose:
            pct = (orig_T - cur.union_T()) / orig_T * 100
            print(f"    iter {iteration}: x_{best_i}→x_{best_i}⊕x_{best_j}  T={cur.union_T()}/{orig_T} ({pct:.1f}%↓)  m={cur.n}")

        if cur.union_T() <= 1:
            break

    return cur, M, b


def _vector_drop_unused(
    vec: VectorANF, M: np.ndarray, b: np.ndarray
) -> tuple[VectorANF, np.ndarray, np.ndarray]:
    """Drop variables not used in ANY component, compact M and b accordingly."""
    used = vec.variables_used()
    if len(used) == vec.n:
        return vec, M, b

    sorted_used = sorted(used)
    new_comps = []
    for f in vec.components:
        new_terms = {}
        for mask, v in f.terms.items():
            new_mask = 0
            for new_idx, old_bit in enumerate(sorted_used):
                if (mask >> old_bit) & 1:
                    new_mask |= 1 << new_idx
            new_terms[new_mask] = v
        new_comps.append(SparseANF(new_terms, len(sorted_used)))

    # Trim M and b: keep only rows for variables still in use
    # M maps current z-variables (rows) → original x-variables (columns)
    new_M = M[sorted_used]
    new_b = b[sorted_used]

    return VectorANF(new_comps), new_M, new_b


def vector_simplify_by_complement(
    vec: VectorANF,
    verbose: bool = False,
) -> tuple[VectorANF, np.ndarray, np.ndarray]:
    """Find complement assignment minimizing union_T.

    For n ≤ 16: exhaustive Gray code (O(2ⁿ · T_union)).
    For n > 16: greedy per-variable.
    """
    n = vec.n
    # Build cache of component terms for faster Gray-code evaluation
    # Actually, just use substitute_affine with M=I and varying b
    best_T = vec.union_T()
    best_b = np.zeros((n, 1), dtype=np.uint8)
    best_vec = vec

    if n > 16:
        # Greedy per-variable
        if verbose:
            print(f"  Greedy complement (n={n}):")
        cur_vec = vec.copy()
        cur_b = np.zeros((n, 1), dtype=np.uint8)
        for var in range(n):
            b_test = cur_b.copy()
            b_test[var, 0] ^= 1
            M_test = np.eye(n, dtype=np.uint8)
            g_test = cur_vec.substitute_affine(M_test, b_test)
            if g_test.union_T() < cur_vec.union_T():
                cur_vec = g_test
                cur_b = b_test
        return cur_vec, np.eye(n, dtype=np.uint8), cur_b
    else:
        # Exhaustive Gray code
        if verbose:
            print(f"  Complement search (Gray, 2^{n}={1<<n}):")
            t0 = time.time()

        # For Gray code, we flip one bit at a time
        cur_vec = vec.copy()
        cur_b = np.zeros((n, 1), dtype=np.uint8)
        best_vec = cur_vec
        best_b = cur_b.copy()

        for gray_idx in range(1, 1 << n):
            # Find which bit changes in Gray code
            bit = (gray_idx ^ (gray_idx >> 1)) ^ ((gray_idx - 1) ^ ((gray_idx - 1) >> 1))
            bit = (bit & -bit).bit_length() - 1

            # Flip that bit in b and apply to components
            # Since complementing bit means: x_bit → x_bit⊕1
            # In ANF: substitute x_bit with x_bit⊕1
            # For monomials containing x_bit, toggle the monomial
            new_comps = []
            for f in cur_vec.components:
                new = {}
                for mask, v in f.terms.items():
                    new[mask] = new.get(mask, 0) ^ v
                    if (mask >> bit) & 1:
                        rest = mask ^ (1 << bit)
                        new[rest] = new.get(rest, 0) ^ v
                clean = {m: c for m, c in new.items() if c != 0}
                new_comps.append(SparseANF(clean, f.n))
            cur_vec = VectorANF(new_comps)
            cur_b[bit, 0] ^= 1

            T = cur_vec.union_T()
            if T < best_T:
                best_T = T
                best_vec = cur_vec
                best_b = cur_b.copy()

        if verbose:
            print(f"    best: T={best_T}/{vec.union_T()} ({(vec.union_T()-best_T)/vec.union_T()*100:.1f}%↓)  time={time.time()-t0:.2f}s")

        return best_vec, np.eye(n, dtype=np.uint8), best_b


def vector_simplify(
    vec: VectorANF,
    verbose: bool = False,
) -> tuple[VectorANF, np.ndarray, np.ndarray]:
    """Combined vector Boolean function simplification pipeline.

    1. Complement search
    2. Greedy XOR merge
    3. Dimension reduction (via Walsh search on each component, shared M)

    Returns (g_vec, M, b) where f_i(x) = g_i(Mx⊕b) for all i.
    """
    n = vec.n
    M_acc = np.eye(n, dtype=np.uint8)
    b_acc = np.zeros((n, 1), dtype=np.uint8)
    cur = vec.copy()

    if verbose:
        print(f"Vector simplify: n={n}, k={vec.k}, T_union={vec.union_T()}")

    # Phase 1: complement search
    g_comp, M_c, b_c = vector_simplify_by_complement(cur, verbose=verbose)
    if g_comp.union_T() < cur.union_T():
        M_acc = M_c @ M_acc
        b_acc = (M_c @ b_acc + b_c) % 2
        cur = g_comp
        if verbose:
            print(f"  After complement: T={cur.union_T()}")

    # Phase 2: greedy XOR merge
    g_merge, M_m, b_m = vector_greedy_merge(cur, verbose=verbose)
    if g_merge.union_T() < cur.union_T():
        # Compose: M_m maps from cur's variables → g_merge's variables
        # But M_m was built starting from identity of size cur.n
        M_acc = M_m @ M_acc
        b_acc = (M_m @ b_acc + b_m) % 2
        cur = g_merge
        if verbose:
            print(f"  After merge: T={cur.union_T()}, m={cur.n}")

    return cur, M_acc, b_acc


def vector_simplify_per_component(
    vec: VectorANF,
    verbose: bool = False,
    return_union: bool = False,
    n_walsh_trials: int = 10,
    seed: int = 0,
) -> list[tuple[SparseANF, np.ndarray, np.ndarray]] | tuple[list[tuple[SparseANF, np.ndarray, np.ndarray]], VectorANF, np.ndarray, np.ndarray]:
    """Strategy 2: simplify each component independently.

    Each component f_i finds its own (M_i, b_i) and g_i,
    minimizing T(g_i) individually.

    Runs complement + greedy merge once (deterministic), then
    repeats the Walsh search n_walsh_trials times with different
    seeds to find the best dimension reduction.

    If return_union=True, also returns a shared representation.

    Returns:
    - list of (g_i, M_i, b_i) per component
    - if return_union: also (joint_vec, M_joint, b_joint)
    """
    results = []
    for idx, f in enumerate(vec.components):
        if verbose:
            print(f"  Component {idx}: n={f.n}, T={f.T()}", end="")

        # Phase 1-2: complement + greedy merge (deterministic)
        g, M, b = simplify_by_complement(f, verbose=False)
        g2, M2, b2 = greedy_merge_simplify(g, verbose=False)
        M_acc = (M2 @ M) % 2
        b_acc = ((M2 @ b) % 2).reshape(-1, 1) ^ b2
        g_acc = g2

        if verbose:
            print(f", after greedy: T={g_acc.T()}, m={g_acc.n}")

        # Phase 3: multi-trial Walsh search
        best_g, best_M, best_b = g_acc, M_acc, b_acc
        best_T = g_acc.T()

        for trial in range(n_walsh_trials):
            if g_acc.n <= 3 or g_acc.T() == 0:
                break

            import random as _random
            _random.seed(seed + trial)

            max_m = min(g_acc.n, 10)
            try:
                walsh_results = search_affine_simplification(
                    g_acc, max_m=max_m, top_k=80, verbose=False, n_random=30,
                    seed=seed + trial)
            except Exception:
                walsh_results = []

            if walsh_results and walsh_results[0][1] < best_T:
                m_w, T_w, M_w, b_w, g_w = walsh_results[0]
                best_T = T_w
                best_M = (M_w @ M_acc) % 2
                best_b = ((M_w @ b_acc) % 2).reshape(-1, 1) ^ b_w
                best_g = g_w
                if verbose:
                    print(f"    trial {trial+1}: Walsh found T={T_w}, m={m_w}")

        results.append((best_g, best_M, best_b))
        if verbose:
            print(f"    → T={best_T}, m={best_g.n}")

    if not return_union:
        return results

    # Build shared M containing distinct rows from all M_i
    all_rows: list[np.ndarray] = []
    row_to_idx: dict[tuple, int] = {}
    comp_maps: list[list[int]] = []  # for each component, which row in joint M

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
        return results, VectorANF([SparseANF({}, 0)] * vec.k), np.zeros((0, vec.n), dtype=np.uint8), np.zeros((0, 1), dtype=np.uint8)

    M_joint = np.array(all_rows, dtype=np.uint8)
    m_total = M_joint.shape[0]

    # Build joint b vector: collect all distinct b entries
    b_joint = np.zeros((m_total, 1), dtype=np.uint8)
    for (g, M_i, b_i), cmap in zip(results, comp_maps):
        for local_idx, global_idx in enumerate(cmap):
            b_joint[global_idx, 0] = b_i[local_idx, 0]

    # Remap each g_i to the shared variable space
    joint_comps = []
    for (g, M_i, b_i), cmap in zip(results, comp_maps):
        remapped_terms = {}
        for mask, v in g.terms.items():
            new_mask = 0
            t = mask
            while t:
                lsb = t & -t
                local_var = (lsb.bit_length() - 1)
                global_var = cmap[local_var]
                new_mask |= 1 << global_var
                t ^= lsb
            remapped_terms[new_mask] = v
        joint_comps.append(SparseANF(remapped_terms, m_total))

    joint_vec = VectorANF(joint_comps)
    return results, joint_vec, M_joint, b_joint


# ====================================================================
#  Demo / test on the 5-component function
# ====================================================================


def make_vector_from_expressions(
    exprs: list[str], var_to_bit: dict[str, int], n: int
) -> VectorANF:
    """Parse expressions and create VectorANF."""
    comps = []
    for expr in exprs:
        terms = {}
        for term_str in expr.split('+'):
            term_str = term_str.strip()
            if not term_str:
                continue
            mask = 0
            for v in term_str.split('*'):
                mask |= 1 << var_to_bit[v.strip()]
            terms[mask] = terms.get(mask, 0) ^ 1
        comps.append(SparseANF(terms, n))
    return VectorANF(comps)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vector Boolean function ANF simplification")
    parser.add_argument("--seed", type=int, default=0, help="random seed (default: 0)")
    parser.add_argument("--walsh-trials", type=int, default=10, help="Walsh search trials per component (default: 10)")
    args = parser.parse_args()
    seed = args.seed
    walsh_trials = args.walsh_trials

    # Build variable map: i_2..i_9 → bits 0..7, i_10..i_17 → bits 8..15
    var_to_bit = {}
    for k in range(2, 18):
        var_to_bit[f'i_{k}'] = k - 2  # i_2→0, i_3→1, ..., i_17→15

    m0_expr = ("i_2 + i_10 + i_2*i_10 + i_2*i_3*i_11 + i_3*i_10*i_11 + i_2*i_3*i_4*i_12 + i_3*i_4*i_10*i_12 + i_2*i_4*i_11*i_12 + i_4*i_10*i_11*i_12 + "
               "i_2*i_3*i_4*i_5*i_13 + i_3*i_4*i_5*i_10*i_13 + i_2*i_4*i_5*i_11*i_13 + i_4*i_5*i_10*i_11*i_13 + i_2*i_3*i_5*i_12*i_13 + i_3*i_5*i_10*i_12*i_13 + "
               "i_2*i_5*i_11*i_12*i_13 + i_5*i_10*i_11*i_12*i_13 + "
               "i_2*i_3*i_4*i_5*i_6*i_14 + i_3*i_4*i_5*i_6*i_10*i_14 + i_2*i_4*i_5*i_6*i_11*i_14 + i_4*i_5*i_6*i_10*i_11*i_14 + "
               "i_2*i_3*i_5*i_6*i_12*i_14 + i_3*i_5*i_6*i_10*i_12*i_14 + i_2*i_5*i_6*i_11*i_12*i_14 + i_5*i_6*i_10*i_11*i_12*i_14 + "
               "i_2*i_3*i_4*i_6*i_13*i_14 + i_3*i_4*i_6*i_10*i_13*i_14 + i_2*i_4*i_6*i_11*i_13*i_14 + i_4*i_6*i_10*i_11*i_13*i_14 + "
               "i_2*i_3*i_6*i_12*i_13*i_14 + i_3*i_6*i_10*i_12*i_13*i_14 + i_2*i_6*i_11*i_12*i_13*i_14 + i_6*i_10*i_11*i_12*i_13*i_14 + "
               "i_2*i_3*i_4*i_5*i_6*i_7*i_15 + i_3*i_4*i_5*i_6*i_7*i_10*i_15 + i_2*i_4*i_5*i_6*i_7*i_11*i_15 + i_4*i_5*i_6*i_7*i_10*i_11*i_15 + "
               "i_2*i_3*i_5*i_6*i_7*i_12*i_15 + i_3*i_5*i_6*i_7*i_10*i_12*i_15 + i_2*i_5*i_6*i_7*i_11*i_12*i_15 + i_5*i_6*i_7*i_10*i_11*i_12*i_15 + "
               "i_2*i_3*i_4*i_6*i_7*i_13*i_15 + i_3*i_4*i_6*i_7*i_10*i_13*i_15 + i_2*i_4*i_6*i_7*i_11*i_13*i_15 + i_4*i_6*i_7*i_10*i_11*i_13*i_15 + "
               "i_2*i_3*i_6*i_7*i_12*i_13*i_15 + i_3*i_6*i_7*i_10*i_12*i_13*i_15 + i_2*i_6*i_7*i_11*i_12*i_13*i_15 + i_6*i_7*i_10*i_11*i_12*i_13*i_15 + "
               "i_2*i_3*i_4*i_5*i_7*i_14*i_15 + i_3*i_4*i_5*i_7*i_10*i_14*i_15 + i_2*i_4*i_5*i_7*i_11*i_14*i_15 + i_4*i_5*i_7*i_10*i_11*i_14*i_15 + "
               "i_2*i_3*i_5*i_7*i_12*i_14*i_15 + i_3*i_5*i_7*i_10*i_12*i_14*i_15 + i_2*i_5*i_7*i_11*i_12*i_14*i_15 + i_5*i_7*i_10*i_11*i_12*i_14*i_15 + "
               "i_2*i_3*i_4*i_7*i_13*i_14*i_15 + i_3*i_4*i_7*i_10*i_13*i_14*i_15 + i_2*i_4*i_7*i_11*i_13*i_14*i_15 + i_4*i_7*i_10*i_11*i_13*i_14*i_15 + "
               "i_2*i_3*i_7*i_12*i_13*i_14*i_15 + i_3*i_7*i_10*i_12*i_13*i_14*i_15 + i_2*i_7*i_11*i_12*i_13*i_14*i_15 + i_7*i_10*i_11*i_12*i_13*i_14*i_15 + "
               "i_2*i_3*i_4*i_5*i_6*i_7*i_8*i_16 + i_3*i_4*i_5*i_6*i_7*i_8*i_10*i_16 + i_2*i_4*i_5*i_6*i_7*i_8*i_11*i_16 + i_4*i_5*i_6*i_7*i_8*i_10*i_11*i_16 + "
               "i_2*i_3*i_5*i_6*i_7*i_8*i_12*i_16 + i_3*i_5*i_6*i_7*i_8*i_10*i_12*i_16 + i_2*i_5*i_6*i_7*i_8*i_11*i_12*i_16 + i_5*i_6*i_7*i_8*i_10*i_11*i_12*i_16 + "
               "i_2*i_3*i_4*i_6*i_7*i_8*i_13*i_16 + i_3*i_4*i_6*i_7*i_8*i_10*i_13*i_16 + i_2*i_4*i_6*i_7*i_8*i_11*i_13*i_16 + i_4*i_6*i_7*i_8*i_10*i_11*i_13*i_16 + "
               "i_2*i_3*i_6*i_7*i_8*i_12*i_13*i_16 + i_3*i_6*i_7*i_8*i_10*i_12*i_13*i_16 + i_2*i_6*i_7*i_8*i_11*i_12*i_13*i_16 + i_6*i_7*i_8*i_10*i_11*i_12*i_13*i_16 + "
               "i_2*i_3*i_4*i_5*i_7*i_8*i_14*i_16 + i_3*i_4*i_5*i_7*i_8*i_10*i_14*i_16 + i_2*i_4*i_5*i_7*i_8*i_11*i_14*i_16 + i_4*i_5*i_7*i_8*i_10*i_11*i_14*i_16 + "
               "i_2*i_3*i_5*i_7*i_8*i_12*i_14*i_16 + i_3*i_5*i_7*i_8*i_10*i_12*i_14*i_16 + i_2*i_5*i_7*i_8*i_11*i_12*i_14*i_16 + i_5*i_7*i_8*i_10*i_11*i_12*i_14*i_16 + "
               "i_2*i_3*i_4*i_7*i_8*i_13*i_14*i_16 + i_3*i_4*i_7*i_8*i_10*i_13*i_14*i_16 + i_2*i_4*i_7*i_8*i_11*i_13*i_14*i_16 + i_4*i_7*i_8*i_10*i_11*i_13*i_14*i_16 + "
               "i_2*i_3*i_7*i_8*i_12*i_13*i_14*i_16 + i_3*i_7*i_8*i_10*i_12*i_13*i_14*i_16 + i_2*i_7*i_8*i_11*i_12*i_13*i_14*i_16 + i_7*i_8*i_10*i_11*i_12*i_13*i_14*i_16 + "
               "i_2*i_3*i_4*i_5*i_6*i_8*i_15*i_16 + i_3*i_4*i_5*i_6*i_8*i_10*i_15*i_16 + i_2*i_4*i_5*i_6*i_8*i_11*i_15*i_16 + i_4*i_5*i_6*i_8*i_10*i_11*i_15*i_16 + "
               "i_2*i_3*i_5*i_6*i_8*i_12*i_15*i_16 + i_3*i_5*i_6*i_8*i_10*i_12*i_15*i_16 + i_2*i_5*i_6*i_8*i_11*i_12*i_15*i_16 + i_5*i_6*i_8*i_10*i_11*i_12*i_15*i_16 + "
               "i_2*i_3*i_4*i_6*i_8*i_13*i_15*i_16 + i_3*i_4*i_6*i_8*i_10*i_13*i_15*i_16 + i_2*i_4*i_6*i_8*i_11*i_13*i_15*i_16 + i_4*i_6*i_8*i_10*i_11*i_13*i_15*i_16 + "
               "i_2*i_3*i_6*i_8*i_12*i_13*i_15*i_16 + i_3*i_6*i_8*i_10*i_12*i_13*i_15*i_16 + i_2*i_6*i_8*i_11*i_12*i_13*i_15*i_16 + i_6*i_8*i_10*i_11*i_12*i_13*i_15*i_16 + "
               "i_2*i_3*i_4*i_5*i_8*i_14*i_15*i_16 + i_3*i_4*i_5*i_8*i_10*i_14*i_15*i_16 + i_2*i_4*i_5*i_8*i_11*i_14*i_15*i_16 + i_4*i_5*i_8*i_10*i_11*i_14*i_15*i_16 + "
               "i_2*i_3*i_5*i_8*i_12*i_14*i_15*i_16 + i_3*i_5*i_8*i_10*i_12*i_14*i_15*i_16 + i_2*i_5*i_8*i_11*i_12*i_14*i_15*i_16 + i_5*i_8*i_10*i_11*i_12*i_14*i_15*i_16 + "
               "i_2*i_3*i_4*i_8*i_13*i_14*i_15*i_16 + i_3*i_4*i_8*i_10*i_13*i_14*i_15*i_16 + i_2*i_4*i_8*i_11*i_13*i_14*i_15*i_16 + i_4*i_8*i_10*i_11*i_13*i_14*i_15*i_16 + "
               "i_2*i_3*i_8*i_12*i_13*i_14*i_15*i_16 + i_3*i_8*i_10*i_12*i_13*i_14*i_15*i_16 + i_2*i_8*i_11*i_12*i_13*i_14*i_15*i_16 + i_8*i_10*i_11*i_12*i_13*i_14*i_15*i_16 + "
               "i_2*i_3*i_4*i_5*i_6*i_7*i_8*i_9*i_17 + i_3*i_4*i_5*i_6*i_7*i_8*i_9*i_10*i_17 + i_2*i_4*i_5*i_6*i_7*i_8*i_9*i_11*i_17 + i_4*i_5*i_6*i_7*i_8*i_9*i_10*i_11*i_17 + "
               "i_2*i_3*i_5*i_6*i_7*i_8*i_9*i_12*i_17 + i_3*i_5*i_6*i_7*i_8*i_9*i_10*i_12*i_17 + i_2*i_5*i_6*i_7*i_8*i_9*i_11*i_12*i_17 + i_5*i_6*i_7*i_8*i_9*i_10*i_11*i_12*i_17 + "
               "i_2*i_3*i_4*i_6*i_7*i_8*i_9*i_13*i_17 + i_3*i_4*i_6*i_7*i_8*i_9*i_10*i_13*i_17 + i_2*i_4*i_6*i_7*i_8*i_9*i_11*i_13*i_17 + i_4*i_6*i_7*i_8*i_9*i_10*i_11*i_13*i_17 + "
               "i_2*i_3*i_6*i_7*i_8*i_9*i_12*i_13*i_17 + i_3*i_6*i_7*i_8*i_9*i_10*i_12*i_13*i_17 + i_2*i_6*i_7*i_8*i_9*i_11*i_12*i_13*i_17 + i_6*i_7*i_8*i_9*i_10*i_11*i_12*i_13*i_17 + "
               "i_2*i_3*i_4*i_5*i_7*i_8*i_9*i_14*i_17 + i_3*i_4*i_5*i_7*i_8*i_9*i_10*i_14*i_17 + i_2*i_4*i_5*i_7*i_8*i_9*i_11*i_14*i_17 + i_4*i_5*i_7*i_8*i_9*i_10*i_11*i_14*i_17 + "
               "i_2*i_3*i_5*i_7*i_8*i_9*i_12*i_14*i_17 + i_3*i_5*i_7*i_8*i_9*i_10*i_12*i_14*i_17 + i_2*i_5*i_7*i_8*i_9*i_11*i_12*i_14*i_17 + i_5*i_7*i_8*i_9*i_10*i_11*i_12*i_14*i_17 + "
               "i_2*i_3*i_4*i_7*i_8*i_9*i_13*i_14*i_17 + i_3*i_4*i_7*i_8*i_9*i_10*i_13*i_14*i_17 + i_2*i_4*i_7*i_8*i_9*i_11*i_13*i_14*i_17 + i_4*i_7*i_8*i_9*i_10*i_11*i_13*i_14*i_17 + "
               "i_2*i_3*i_7*i_8*i_9*i_12*i_13*i_14*i_17 + i_3*i_7*i_8*i_9*i_10*i_12*i_13*i_14*i_17 + i_2*i_7*i_8*i_9*i_11*i_12*i_13*i_14*i_17 + i_7*i_8*i_9*i_10*i_11*i_12*i_13*i_14*i_17 + "
               "i_2*i_3*i_4*i_5*i_6*i_8*i_9*i_15*i_17 + i_3*i_4*i_5*i_6*i_8*i_9*i_10*i_15*i_17 + i_2*i_4*i_5*i_6*i_8*i_9*i_11*i_15*i_17 + i_4*i_5*i_6*i_8*i_9*i_10*i_11*i_15*i_17 + "
               "i_2*i_3*i_5*i_6*i_8*i_9*i_12*i_15*i_17 + i_3*i_5*i_6*i_8*i_9*i_10*i_12*i_15*i_17 + i_2*i_5*i_6*i_8*i_9*i_11*i_12*i_15*i_17 + i_5*i_6*i_8*i_9*i_10*i_11*i_12*i_15*i_17 + "
               "i_2*i_3*i_4*i_6*i_8*i_9*i_13*i_15*i_17 + i_3*i_4*i_6*i_8*i_9*i_10*i_13*i_15*i_17 + i_2*i_4*i_6*i_8*i_9*i_11*i_13*i_15*i_17 + i_4*i_6*i_8*i_9*i_10*i_11*i_13*i_15*i_17 + "
               "i_2*i_3*i_6*i_8*i_9*i_12*i_13*i_15*i_17 + i_3*i_6*i_8*i_9*i_10*i_12*i_13*i_15*i_17 + i_2*i_6*i_8*i_9*i_11*i_12*i_13*i_15*i_17 + i_6*i_8*i_9*i_10*i_11*i_12*i_13*i_15*i_17 + "
               "i_2*i_3*i_4*i_5*i_8*i_9*i_14*i_15*i_17 + i_3*i_4*i_5*i_8*i_9*i_10*i_14*i_15*i_17 + i_2*i_4*i_5*i_8*i_9*i_11*i_14*i_15*i_17 + i_4*i_5*i_8*i_9*i_10*i_11*i_14*i_15*i_17 + "
               "i_2*i_3*i_5*i_8*i_9*i_12*i_14*i_15*i_17 + i_3*i_5*i_8*i_9*i_10*i_12*i_14*i_15*i_17 + i_2*i_5*i_8*i_9*i_11*i_12*i_14*i_15*i_17 + i_5*i_8*i_9*i_10*i_11*i_12*i_14*i_15*i_17 + "
               "i_2*i_3*i_4*i_8*i_9*i_13*i_14*i_15*i_17 + i_3*i_4*i_8*i_9*i_10*i_13*i_14*i_15*i_17 + i_2*i_4*i_8*i_9*i_11*i_13*i_14*i_15*i_17 + i_4*i_8*i_9*i_10*i_11*i_13*i_14*i_15*i_17 + "
               "i_2*i_3*i_8*i_9*i_12*i_13*i_14*i_15*i_17 + i_3*i_8*i_9*i_10*i_12*i_13*i_14*i_15*i_17 + i_2*i_8*i_9*i_11*i_12*i_13*i_14*i_15*i_17 + i_8*i_9*i_10*i_11*i_12*i_13*i_14*i_15*i_17 + "
               "i_2*i_3*i_4*i_5*i_6*i_7*i_9*i_16*i_17 + i_3*i_4*i_5*i_6*i_7*i_9*i_10*i_16*i_17 + i_2*i_4*i_5*i_6*i_7*i_9*i_11*i_16*i_17 + i_4*i_5*i_6*i_7*i_9*i_10*i_11*i_16*i_17 + "
               "i_2*i_3*i_5*i_6*i_7*i_9*i_12*i_16*i_17 + i_3*i_5*i_6*i_7*i_9*i_10*i_12*i_16*i_17 + i_2*i_5*i_6*i_7*i_9*i_11*i_12*i_16*i_17 + i_5*i_6*i_7*i_9*i_10*i_11*i_12*i_16*i_17 + "
               "i_2*i_3*i_4*i_6*i_7*i_9*i_13*i_16*i_17 + i_3*i_4*i_6*i_7*i_9*i_10*i_13*i_16*i_17 + i_2*i_4*i_6*i_7*i_9*i_11*i_13*i_16*i_17 + i_4*i_6*i_7*i_9*i_10*i_11*i_13*i_16*i_17 + "
               "i_2*i_3*i_6*i_7*i_9*i_12*i_13*i_16*i_17 + i_3*i_6*i_7*i_9*i_10*i_12*i_13*i_16*i_17 + i_2*i_6*i_7*i_9*i_11*i_12*i_13*i_16*i_17 + i_6*i_7*i_9*i_10*i_11*i_12*i_13*i_16*i_17 + "
               "i_2*i_3*i_4*i_5*i_7*i_9*i_14*i_16*i_17 + i_3*i_4*i_5*i_7*i_9*i_10*i_14*i_16*i_17 + i_2*i_4*i_5*i_7*i_9*i_11*i_14*i_16*i_17 + i_4*i_5*i_7*i_9*i_10*i_11*i_14*i_16*i_17 + "
               "i_2*i_3*i_5*i_7*i_9*i_12*i_14*i_16*i_17 + i_3*i_5*i_7*i_9*i_10*i_12*i_14*i_16*i_17 + i_2*i_5*i_7*i_9*i_11*i_12*i_14*i_16*i_17 + i_5*i_7*i_9*i_10*i_11*i_12*i_14*i_16*i_17 + "
               "i_2*i_3*i_4*i_7*i_9*i_13*i_14*i_16*i_17 + i_3*i_4*i_7*i_9*i_10*i_13*i_14*i_16*i_17 + i_2*i_4*i_7*i_9*i_11*i_13*i_14*i_16*i_17 + i_4*i_7*i_9*i_10*i_11*i_13*i_14*i_16*i_17 + "
               "i_2*i_3*i_7*i_9*i_12*i_13*i_14*i_16*i_17 + i_3*i_7*i_9*i_10*i_12*i_13*i_14*i_16*i_17 + i_2*i_7*i_9*i_11*i_12*i_13*i_14*i_16*i_17 + i_7*i_9*i_10*i_11*i_12*i_13*i_14*i_16*i_17 + "
               "i_2*i_3*i_4*i_5*i_6*i_9*i_15*i_16*i_17 + i_3*i_4*i_5*i_6*i_9*i_10*i_15*i_16*i_17 + i_2*i_4*i_5*i_6*i_9*i_11*i_15*i_16*i_17 + i_4*i_5*i_6*i_9*i_10*i_11*i_15*i_16*i_17 + "
               "i_2*i_3*i_5*i_6*i_9*i_12*i_15*i_16*i_17 + i_3*i_5*i_6*i_9*i_10*i_12*i_15*i_16*i_17 + i_2*i_5*i_6*i_9*i_11*i_12*i_15*i_16*i_17 + i_5*i_6*i_9*i_10*i_11*i_12*i_15*i_16*i_17 + "
               "i_2*i_3*i_4*i_6*i_9*i_13*i_15*i_16*i_17 + i_3*i_4*i_6*i_9*i_10*i_13*i_15*i_16*i_17 + i_2*i_4*i_6*i_9*i_11*i_13*i_15*i_16*i_17 + i_4*i_6*i_9*i_10*i_11*i_13*i_15*i_16*i_17 + "
               "i_2*i_3*i_6*i_9*i_12*i_13*i_15*i_16*i_17 + i_3*i_6*i_9*i_10*i_12*i_13*i_15*i_16*i_17 + i_2*i_6*i_9*i_11*i_12*i_13*i_15*i_16*i_17 + i_6*i_9*i_10*i_11*i_12*i_13*i_15*i_16*i_17 + "
               "i_2*i_3*i_4*i_5*i_9*i_14*i_15*i_16*i_17 + i_3*i_4*i_5*i_9*i_10*i_14*i_15*i_16*i_17 + i_2*i_4*i_5*i_9*i_11*i_14*i_15*i_16*i_17 + i_4*i_5*i_9*i_10*i_11*i_14*i_15*i_16*i_17 + "
               "i_2*i_3*i_5*i_9*i_12*i_14*i_15*i_16*i_17 + i_3*i_5*i_9*i_10*i_12*i_14*i_15*i_16*i_17 + i_2*i_5*i_9*i_11*i_12*i_14*i_15*i_16*i_17 + i_5*i_9*i_10*i_11*i_12*i_14*i_15*i_16*i_17 + "
               "i_2*i_3*i_4*i_9*i_13*i_14*i_15*i_16*i_17 + i_3*i_4*i_9*i_10*i_13*i_14*i_15*i_16*i_17 + i_2*i_4*i_9*i_11*i_13*i_14*i_15*i_16*i_17 + i_4*i_9*i_10*i_11*i_13*i_14*i_15*i_16*i_17 + "
               "i_2*i_3*i_9*i_12*i_13*i_14*i_15*i_16*i_17 + i_3*i_9*i_10*i_12*i_13*i_14*i_15*i_16*i_17 + i_2*i_9*i_11*i_12*i_13*i_14*i_15*i_16*i_17 + i_9*i_10*i_11*i_12*i_13*i_14*i_15*i_16*i_17")

    m1_expr = ("i_2 + i_10 + i_3*i_11 + i_3*i_4*i_12 + i_4*i_11*i_12 + i_3*i_4*i_5*i_13 + i_4*i_5*i_11*i_13 + i_3*i_5*i_12*i_13 + i_5*i_11*i_12*i_13 + "
               "i_3*i_4*i_5*i_6*i_14 + i_4*i_5*i_6*i_11*i_14 + i_3*i_5*i_6*i_12*i_14 + i_5*i_6*i_11*i_12*i_14 + i_3*i_4*i_6*i_13*i_14 + i_4*i_6*i_11*i_13*i_14 + i_3*i_6*i_12*i_13*i_14 + i_6*i_11*i_12*i_13*i_14 + "
               "i_3*i_4*i_5*i_6*i_7*i_15 + i_4*i_5*i_6*i_7*i_11*i_15 + i_3*i_5*i_6*i_7*i_12*i_15 + i_5*i_6*i_7*i_11*i_12*i_15 + i_3*i_4*i_6*i_7*i_13*i_15 + i_4*i_6*i_7*i_11*i_13*i_15 + i_3*i_6*i_7*i_12*i_13*i_15 + i_6*i_7*i_11*i_12*i_13*i_15 + "
               "i_3*i_4*i_5*i_7*i_14*i_15 + i_4*i_5*i_7*i_11*i_14*i_15 + i_3*i_5*i_7*i_12*i_14*i_15 + i_5*i_7*i_11*i_12*i_14*i_15 + i_3*i_4*i_7*i_13*i_14*i_15 + i_4*i_7*i_11*i_13*i_14*i_15 + i_3*i_7*i_12*i_13*i_14*i_15 + i_7*i_11*i_12*i_13*i_14*i_15 + "
               "i_3*i_4*i_5*i_6*i_7*i_8*i_16 + i_4*i_5*i_6*i_7*i_8*i_11*i_16 + i_3*i_5*i_6*i_7*i_8*i_12*i_16 + i_5*i_6*i_7*i_8*i_11*i_12*i_16 + "
               "i_3*i_4*i_6*i_7*i_8*i_13*i_16 + i_4*i_6*i_7*i_8*i_11*i_13*i_16 + i_3*i_6*i_7*i_8*i_12*i_13*i_16 + i_6*i_7*i_8*i_11*i_12*i_13*i_16 + "
               "i_3*i_4*i_5*i_7*i_8*i_14*i_16 + i_4*i_5*i_7*i_8*i_11*i_14*i_16 + i_3*i_5*i_7*i_8*i_12*i_14*i_16 + i_5*i_7*i_8*i_11*i_12*i_14*i_16 + "
               "i_3*i_4*i_7*i_8*i_13*i_14*i_16 + i_4*i_7*i_8*i_11*i_13*i_14*i_16 + i_3*i_7*i_8*i_12*i_13*i_14*i_16 + i_7*i_8*i_11*i_12*i_13*i_14*i_16 + "
               "i_3*i_4*i_5*i_6*i_8*i_15*i_16 + i_4*i_5*i_6*i_8*i_11*i_15*i_16 + i_3*i_5*i_6*i_8*i_12*i_15*i_16 + i_5*i_6*i_8*i_11*i_12*i_15*i_16 + "
               "i_3*i_4*i_6*i_8*i_13*i_15*i_16 + i_4*i_6*i_8*i_11*i_13*i_15*i_16 + i_3*i_6*i_8*i_12*i_13*i_15*i_16 + i_6*i_8*i_11*i_12*i_13*i_15*i_16 + "
               "i_3*i_4*i_5*i_8*i_14*i_15*i_16 + i_4*i_5*i_8*i_11*i_14*i_15*i_16 + i_3*i_5*i_8*i_12*i_14*i_15*i_16 + i_5*i_8*i_11*i_12*i_14*i_15*i_16 + "
               "i_3*i_4*i_8*i_13*i_14*i_15*i_16 + i_4*i_8*i_11*i_13*i_14*i_15*i_16 + i_3*i_8*i_12*i_13*i_14*i_15*i_16 + i_8*i_11*i_12*i_13*i_14*i_15*i_16 + "
               "i_3*i_4*i_5*i_6*i_7*i_8*i_9*i_17 + i_4*i_5*i_6*i_7*i_8*i_9*i_11*i_17 + i_3*i_5*i_6*i_7*i_8*i_9*i_12*i_17 + i_5*i_6*i_7*i_8*i_9*i_11*i_12*i_17 + "
               "i_3*i_4*i_6*i_7*i_8*i_9*i_13*i_17 + i_4*i_6*i_7*i_8*i_9*i_11*i_13*i_17 + i_3*i_6*i_7*i_8*i_9*i_12*i_13*i_17 + i_6*i_7*i_8*i_9*i_11*i_12*i_13*i_17 + "
               "i_3*i_4*i_5*i_7*i_8*i_9*i_14*i_17 + i_4*i_5*i_7*i_8*i_9*i_11*i_14*i_17 + i_3*i_5*i_7*i_8*i_9*i_12*i_14*i_17 + i_5*i_7*i_8*i_9*i_11*i_12*i_14*i_17 + "
               "i_3*i_4*i_7*i_8*i_9*i_13*i_14*i_17 + i_4*i_7*i_8*i_9*i_11*i_13*i_14*i_17 + i_3*i_7*i_8*i_9*i_12*i_13*i_14*i_17 + i_7*i_8*i_9*i_11*i_12*i_13*i_14*i_17 + "
               "i_3*i_4*i_5*i_6*i_8*i_9*i_15*i_17 + i_4*i_5*i_6*i_8*i_9*i_11*i_15*i_17 + i_3*i_5*i_6*i_8*i_9*i_12*i_15*i_17 + i_5*i_6*i_8*i_9*i_11*i_12*i_15*i_17 + "
               "i_3*i_4*i_6*i_8*i_9*i_13*i_15*i_17 + i_4*i_6*i_8*i_9*i_11*i_13*i_15*i_17 + i_3*i_6*i_8*i_9*i_12*i_13*i_15*i_17 + i_6*i_8*i_9*i_11*i_12*i_13*i_15*i_17 + "
               "i_3*i_4*i_5*i_8*i_9*i_14*i_15*i_17 + i_4*i_5*i_8*i_9*i_11*i_14*i_15*i_17 + i_3*i_5*i_8*i_9*i_12*i_14*i_15*i_17 + i_5*i_8*i_9*i_11*i_12*i_14*i_15*i_17 + "
               "i_3*i_4*i_8*i_9*i_13*i_14*i_15*i_17 + i_4*i_8*i_9*i_11*i_13*i_14*i_15*i_17 + i_3*i_8*i_9*i_12*i_13*i_14*i_15*i_17 + i_8*i_9*i_11*i_12*i_13*i_14*i_15*i_17 + "
               "i_3*i_4*i_5*i_6*i_7*i_9*i_16*i_17 + i_4*i_5*i_6*i_7*i_9*i_11*i_16*i_17 + i_3*i_5*i_6*i_7*i_9*i_12*i_16*i_17 + i_5*i_6*i_7*i_9*i_11*i_12*i_16*i_17 + "
               "i_3*i_4*i_6*i_7*i_9*i_13*i_16*i_17 + i_4*i_6*i_7*i_9*i_11*i_13*i_16*i_17 + i_3*i_6*i_7*i_9*i_12*i_13*i_16*i_17 + i_6*i_7*i_9*i_11*i_12*i_13*i_16*i_17 + "
               "i_3*i_4*i_5*i_7*i_9*i_14*i_16*i_17 + i_4*i_5*i_7*i_9*i_11*i_14*i_16*i_17 + i_3*i_5*i_7*i_9*i_12*i_14*i_16*i_17 + i_5*i_7*i_9*i_11*i_12*i_14*i_16*i_17 + "
               "i_3*i_4*i_7*i_9*i_13*i_14*i_16*i_17 + i_4*i_7*i_9*i_11*i_13*i_14*i_16*i_17 + i_3*i_7*i_9*i_12*i_13*i_14*i_16*i_17 + i_7*i_9*i_11*i_12*i_13*i_14*i_16*i_17 + "
               "i_3*i_4*i_5*i_6*i_9*i_15*i_16*i_17 + i_4*i_5*i_6*i_9*i_11*i_15*i_16*i_17 + i_3*i_5*i_6*i_9*i_12*i_15*i_16*i_17 + i_5*i_6*i_9*i_11*i_12*i_15*i_16*i_17 + "
               "i_3*i_4*i_6*i_9*i_13*i_15*i_16*i_17 + i_4*i_6*i_9*i_11*i_13*i_15*i_16*i_17 + i_3*i_6*i_9*i_12*i_13*i_15*i_16*i_17 + i_6*i_9*i_11*i_12*i_13*i_15*i_16*i_17 + "
               "i_3*i_4*i_5*i_9*i_14*i_15*i_16*i_17 + i_4*i_5*i_9*i_11*i_14*i_15*i_16*i_17 + i_3*i_5*i_9*i_12*i_14*i_15*i_16*i_17 + i_5*i_9*i_11*i_12*i_14*i_15*i_16*i_17 + "
               "i_3*i_4*i_9*i_13*i_14*i_15*i_16*i_17 + i_4*i_9*i_11*i_13*i_14*i_15*i_16*i_17 + i_3*i_9*i_12*i_13*i_14*i_15*i_16*i_17 + i_9*i_11*i_12*i_13*i_14*i_15*i_16*i_17")

    m2_expr = ("i_3 + i_11 + i_4*i_12 + i_4*i_5*i_13 + i_5*i_12*i_13 + i_4*i_5*i_6*i_14 + i_5*i_6*i_12*i_14 + i_4*i_6*i_13*i_14 + i_6*i_12*i_13*i_14 + "
               "i_4*i_5*i_6*i_7*i_15 + i_5*i_6*i_7*i_12*i_15 + i_4*i_6*i_7*i_13*i_15 + i_6*i_7*i_12*i_13*i_15 + i_4*i_5*i_7*i_14*i_15 + i_5*i_7*i_12*i_14*i_15 + i_4*i_7*i_13*i_14*i_15 + i_7*i_12*i_13*i_14*i_15 + "
               "i_4*i_5*i_6*i_7*i_8*i_16 + i_5*i_6*i_7*i_8*i_12*i_16 + i_4*i_6*i_7*i_8*i_13*i_16 + i_6*i_7*i_8*i_12*i_13*i_16 + "
               "i_4*i_5*i_7*i_8*i_14*i_16 + i_5*i_7*i_8*i_12*i_14*i_16 + i_4*i_7*i_8*i_13*i_14*i_16 + i_7*i_8*i_12*i_13*i_14*i_16 + "
               "i_4*i_5*i_6*i_8*i_15*i_16 + i_5*i_6*i_8*i_12*i_15*i_16 + i_4*i_6*i_8*i_13*i_15*i_16 + i_6*i_8*i_12*i_13*i_15*i_16 + "
               "i_4*i_5*i_8*i_14*i_15*i_16 + i_5*i_8*i_12*i_14*i_15*i_16 + i_4*i_8*i_13*i_14*i_15*i_16 + i_8*i_12*i_13*i_14*i_15*i_16 + "
               "i_4*i_5*i_6*i_7*i_8*i_9*i_17 + i_5*i_6*i_7*i_8*i_9*i_12*i_17 + i_4*i_6*i_7*i_8*i_9*i_13*i_17 + i_6*i_7*i_8*i_9*i_12*i_13*i_17 + "
               "i_4*i_5*i_7*i_8*i_9*i_14*i_17 + i_5*i_7*i_8*i_9*i_12*i_14*i_17 + i_4*i_7*i_8*i_9*i_13*i_14*i_17 + i_7*i_8*i_9*i_12*i_13*i_14*i_17 + "
               "i_4*i_5*i_6*i_8*i_9*i_15*i_17 + i_5*i_6*i_8*i_9*i_12*i_15*i_17 + i_4*i_6*i_8*i_9*i_13*i_15*i_17 + i_6*i_8*i_9*i_12*i_13*i_15*i_17 + "
               "i_4*i_5*i_8*i_9*i_14*i_15*i_17 + i_5*i_8*i_9*i_12*i_14*i_15*i_17 + i_4*i_8*i_9*i_13*i_14*i_15*i_17 + i_8*i_9*i_12*i_13*i_14*i_15*i_17 + "
               "i_4*i_5*i_6*i_7*i_9*i_16*i_17 + i_5*i_6*i_7*i_9*i_12*i_16*i_17 + i_4*i_6*i_7*i_9*i_13*i_16*i_17 + i_6*i_7*i_9*i_12*i_13*i_16*i_17 + "
               "i_4*i_5*i_7*i_9*i_14*i_16*i_17 + i_5*i_7*i_9*i_12*i_14*i_16*i_17 + i_4*i_7*i_9*i_13*i_14*i_16*i_17 + i_7*i_9*i_12*i_13*i_14*i_16*i_17 + "
               "i_4*i_5*i_6*i_9*i_15*i_16*i_17 + i_5*i_6*i_9*i_12*i_15*i_16*i_17 + i_4*i_6*i_9*i_13*i_15*i_16*i_17 + i_6*i_9*i_12*i_13*i_15*i_16*i_17 + "
               "i_4*i_5*i_9*i_14*i_15*i_16*i_17 + i_5*i_9*i_12*i_14*i_15*i_16*i_17 + i_4*i_9*i_13*i_14*i_15*i_16*i_17 + i_9*i_12*i_13*i_14*i_15*i_16*i_17")

    m3_expr = ("i_4 + i_12 + i_5*i_13 + i_5*i_6*i_14 + i_6*i_13*i_14 + i_5*i_6*i_7*i_15 + i_6*i_7*i_13*i_15 + i_5*i_7*i_14*i_15 + i_7*i_13*i_14*i_15 + "
               "i_5*i_6*i_7*i_8*i_16 + i_6*i_7*i_8*i_13*i_16 + i_5*i_7*i_8*i_14*i_16 + i_7*i_8*i_13*i_14*i_16 + i_5*i_6*i_8*i_15*i_16 + i_6*i_8*i_13*i_15*i_16 + i_5*i_8*i_14*i_15*i_16 + i_8*i_13*i_14*i_15*i_16 + "
               "i_5*i_6*i_7*i_8*i_9*i_17 + i_6*i_7*i_8*i_9*i_13*i_17 + i_5*i_7*i_8*i_9*i_14*i_17 + i_7*i_8*i_9*i_13*i_14*i_17 + "
               "i_5*i_6*i_8*i_9*i_15*i_17 + i_6*i_8*i_9*i_13*i_15*i_17 + i_5*i_8*i_9*i_14*i_15*i_17 + i_8*i_9*i_13*i_14*i_15*i_17 + "
               "i_5*i_6*i_7*i_9*i_16*i_17 + i_6*i_7*i_9*i_13*i_16*i_17 + i_5*i_7*i_9*i_14*i_16*i_17 + i_7*i_9*i_13*i_14*i_16*i_17 + "
               "i_5*i_6*i_9*i_15*i_16*i_17 + i_6*i_9*i_13*i_15*i_16*i_17 + i_5*i_9*i_14*i_15*i_16*i_17 + i_9*i_13*i_14*i_15*i_16*i_17")

    m4_expr = ("i_5 + i_13 + i_6*i_14 + i_6*i_7*i_15 + i_7*i_14*i_15 + i_6*i_7*i_8*i_16 + i_7*i_8*i_14*i_16 + i_6*i_8*i_15*i_16 + i_8*i_14*i_15*i_16 + "
               "i_6*i_7*i_8*i_9*i_17 + i_7*i_8*i_9*i_14*i_17 + i_6*i_8*i_9*i_15*i_17 + i_8*i_9*i_14*i_15*i_17 + i_6*i_7*i_9*i_16*i_17 + i_7*i_9*i_14*i_16*i_17 + "
               "i_6*i_9*i_15*i_16*i_17 + i_9*i_14*i_15*i_16*i_17")

    exprs = [m0_expr, m1_expr, m2_expr, m3_expr, m4_expr]
    vec = make_vector_from_expressions(exprs, var_to_bit, 16)

    print("=" * 70)
    print("  VECTOR BOOLEAN FUNCTION: k=5, n=16")
    print("=" * 70)
    for i, f in enumerate(vec.components):
        print(f"  m_{i}: T={f.T()}, deg={f.degree()}")
    print(f"  UNION T: {vec.union_T()}")
    print()

    # Run joint simplification
    t0 = time.time()
    g_vec, M, b = vector_simplify(vec, verbose=True)
    elapsed = time.time() - t0

    print(f"\n{'='*70}")
    print(f"  Result: T_union: {vec.union_T()} → {g_vec.union_T()}  ({(vec.union_T()-g_vec.union_T())/vec.union_T()*100:.1f}%↓)")
    print(f"  n={vec.n} → {g_vec.n},  time={elapsed:.2f}s")

    for i, g in enumerate(g_vec.components):
        print(f"  g_{i}: T={g.T()} (individual)")

    # Verify all components
    rng = __import__('random').Random(42)
    errors = 0
    n_orig = vec.n
    for _ in range(100):
        x = rng.randint(0, (1 << n_orig) - 1)
        z_arr = (M.astype(np.uint8) @ np.array([(x >> i) & 1 for i in range(n_orig)], dtype=np.uint8)) % 2
        z_arr = (z_arr ^ b.ravel()) % 2
        z = sum(int(z_arr[i]) << i for i in range(g_vec.n))
        for i in range(vec.k):
            if vec.components[i].eval_mask(x) != g_vec.components[i].eval_mask(z):
                errors += 1
    print(f"\n  Verify: {errors}/{100*vec.k} errors")
    print(f"  M ({M.shape[0]}×{M.shape[1]}):\n{M}")
    print(f"  b = {b.ravel()}")

    # Also run per-component separate simplification (strategy 2)
    print(f"\n{'='*70}")
    print("  PER-COMPONENT SIMPLIFICATION (Strategy 2)")
    print("=" * 70)
    per_component_results, joint_vec, M_joint, b_joint = vector_simplify_per_component(vec, verbose=True, return_union=True, seed=seed, n_walsh_trials=walsh_trials)
    print()
    for idx, (g, M_i, b_i) in enumerate(per_component_results):
        print(f"  g_{idx}: T={g.T()}, m={g.n}")
    print(f"\n  After merging distinct rows into shared M ({M_joint.shape[0]}×{M_joint.shape[1]}):")
    print(f"  Union T (shared variable space): {joint_vec.union_T()}")
