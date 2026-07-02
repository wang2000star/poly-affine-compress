"""
用随机可逆仿射变换 z = Mx⊕b 探索 GL(7,2) 空间 — 小规模试跑
"""
import sys
import numpy as np
from anf_factor import SparseANF
from simplify_poly import parse_poly, compute_union_T

def gf2_rank(M):
    A = np.array(M, dtype=np.uint8)
    n = len(A)
    rank = 0
    for col in range(n):
        pivot = None
        for row in range(rank, n):
            if A[row, col]:
                pivot = row
                break
        if pivot is None:
            continue
        A[[rank, pivot]] = A[[pivot, rank]]
        for row in range(n):
            if row != rank and A[row, col]:
                A[row] ^= A[rank]
        rank += 1
    return rank

def random_invertible(n, rng):
    while True:
        M = rng.integers(0, 2, (n, n))
        if gf2_rank(M) == n:
            return M

print("Loading...", flush=True)
inputs, outputs, funcs = parse_poly('/home/wangfz/bool/ctrl/ctrl.poly')
n = len(inputs)
base_funcs = [SparseANF(dict(funcs[name]), n) for name in outputs]
orig_union = compute_union_T(base_funcs)
print(f"原始: n={n}, 输出={len(outputs)}, union T={orig_union}", flush=True)

rng = np.random.default_rng(42)
best_u = orig_union

N_TRIALS = 5000
print(f"\n搜索 {N_TRIALS} 个随机可逆变换...", flush=True)
for trial in range(N_TRIALS):
    M = random_invertible(n, rng)
    b = rng.integers(0, 2, n)

    cand = []
    for f in base_funcs:
        g = f.substitute_affine(M, b.reshape(-1,1), verify=False)
        cand.append(g)
    u = compute_union_T(cand)

    if u < best_u:
        best_u = u
        print(f"  Trial {trial}: union T={best_u}", flush=True)

    if trial % 1000 == 0 and trial > 0:
        print(f"  ... {trial}/{N_TRIALS}, best={best_u}", flush=True)

print(f"\n结果: 原始 T={orig_union} → 最佳 T={best_u}", flush=True)
