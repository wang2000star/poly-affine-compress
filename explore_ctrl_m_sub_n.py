"""
尝试 m < n 的共享变换（降维）。条件苛刻：所有 26 个输出必须在 ker(M) 上取常值。
试 10000 个随机 M (m=3..6)，只保留通过验证的。
"""
import numpy as np
from anf_factor import SparseANF
from simplify_poly import parse_poly, compute_union_T

def gf2_rank(M):
    A = np.array(M, dtype=np.uint8)
    m, n = A.shape
    rank = 0
    for col in range(n):
        pivot = None
        for row in range(rank, m):
            if A[row, col]:
                pivot = row
                break
        if pivot is None:
            continue
        A[[rank, pivot]] = A[[pivot, rank]]
        for row in range(m):
            if row != rank and A[row, col]:
                A[row] ^= A[rank]
        rank += 1
    return rank

def factors_through(M, funcs):
    """检查所有 funcs 是否都因子化通过 M"""
    n = M.shape[1]
    # 对每个 func, 检查 Mx₁ = Mx₂ ⇒ f(x₁) = f(x₂)
    for f in funcs:
        seen = {}
        for x_mask in range(1 << n):
            x = np.array([(x_mask >> j) & 1 for j in range(n)], dtype=np.int64)
            z = (M @ x) % 2
            z_key = tuple(z)
            f_val = f.eval_mask(x_mask)
            if z_key in seen and seen[z_key] != f_val:
                return False
            seen[z_key] = f_val
    return True

print("Loading...", flush=True)
inputs, outputs, funcs = parse_poly('/home/wangfz/bool/ctrl.poly')
n = len(inputs)
base_funcs = [SparseANF(dict(funcs[name]), n) for name in outputs]
orig_union = compute_union_T(base_funcs)
print(f"原始: n={n}, 输出={len(outputs)}, union T={orig_union}", flush=True)

rng = np.random.default_rng(42)
best_u = orig_union
best_info = None
valid_count = 0
trial_count = 0

print(f"\n搜索随机 M(m<n)...", flush=True)
for trial in range(5000):
    m = rng.integers(3, 7)
    M = rng.integers(0, 2, (m, n))
    while gf2_rank(M) < m:
        M = rng.integers(0, 2, (m, n))
    b = rng.integers(0, 2, m)

    trial_count += 1

    if not factors_through(M, base_funcs):
        continue

    valid_count += 1

    cand = []
    for f in base_funcs:
        try:
            g = f.substitute_affine(M, b.reshape(-1,1), verify=False)
            cand.append(g)
        except:
            break
    else:
        u = compute_union_T(cand)
        if u < best_u:
            best_u = u
            best_info = (m, M.copy(), b.copy())
            print(f"   Trial {trial}: m={m}, union T={best_u}", flush=True)

print(f"\n验证通过: {valid_count}/{trial_count}", flush=True)
print(f"最佳: union T={best_u}", flush=True)
