"""
用随机可逆仿射变换 z = Mx⊕b 探索 GL(7,2) 空间。
M 为 n×n GF(2) 可逆矩阵，b 为 n 维向量。
这比 greedy merge 覆盖更广的变换空间。
"""
import numpy as np
from anf_factor import SparseANF
from simplify_poly import parse_poly, compute_union_T

# GF(2) 矩阵求逆（用于验证 M 可逆）
def gf2_rank(M):
    A = M.copy()
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

# 加载
inputs, outputs, funcs = parse_poly('/home/wangfz/bool/ctrl.poly')
n = len(inputs)
base_funcs = [SparseANF(dict(funcs[name]), n) for name in outputs]
orig_union = compute_union_T(base_funcs)

print(f"原始: n={n}, 输出={len(outputs)}, union T={orig_union}")
print()
print("=" * 60)
print(f"随机可逆仿射变换搜索 GL({n},2)")
print("=" * 60)

rng = np.random.default_rng(42)
best_u = orig_union
best_pair = None

# 也尝试 complement 后再随机 M
best_u2 = orig_union
best_pair2 = None

# Apply best complement first (comp=0000100 = bit 2)
comp = 0b0000100
comp_funcs = []
for f in base_funcs:
    b_comp = np.array([(comp >> j) & 1 for j in range(n)], dtype=np.int64)
    M_comp = np.eye(n, dtype=np.int64)
    g = f.substitute_affine(M_comp, b_comp.reshape(-1,1), verify=False)
    comp_funcs.append(g)

comp_union = compute_union_T(comp_funcs)
print(f"\n补码 0000100 后: union T={comp_union}\n")

N_TRIALS = 100000

for trial in range(N_TRIALS):
    M = random_invertible(n, rng)
    b = rng.integers(0, 2, n)

    # 原始函数上应用随机 (M,b)
    cand = []
    for f in base_funcs:
        g = f.substitute_affine(M, b.reshape(-1,1), verify=False)
        cand.append(g)
    u = compute_union_T(cand)

    if u < best_u:
        best_u = u
        best_pair = (M.copy(), b.copy())
        print(f"  Direct trial {trial}: union T={best_u}")

    # 补码后的函数应用随机 (M,b)
    cand2 = []
    for f in comp_funcs:
        g = f.substitute_affine(M, b.reshape(-1,1), verify=False)
        cand2.append(g)
    u2 = compute_union_T(cand2)

    if u2 < best_u2:
        best_u2 = u2
        best_pair2 = (M.copy(), b.copy())
        print(f"  Post-comp trial {trial}: union T={best_u2}")

    if trial % 10000 == 0 and trial > 0:
        print(f"  ... {trial}/{N_TRIALS}, best={best_u}, post-comp best={best_u2}")

print()
print("=" * 60)
print("结果汇总")
print("=" * 60)
print(f"  原始:                     T={orig_union}")
print(f"  随机可逆变换:             T={best_u}")
print(f"  补码 + 随机可逆变换:     T={best_u2}")
