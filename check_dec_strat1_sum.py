"""检查策略 1 前后 ΣT 的变化"""
import numpy as np
from anf_factor import SparseANF
from simplify_poly import parse_poly, compute_union_T, _apply_merge_to_all

inputs, outputs, funcs = parse_poly('/home/wangfz/bool/dec.poly')
n = len(inputs)
func_list = [SparseANF(dict(funcs[name]), n) for name in outputs]

orig_sum_T = sum(f.T() for f in func_list)
orig_union = compute_union_T(func_list)
print(f"原始: ΣT={orig_sum_T}, union T={orig_union}")

# Strategy 1: best complement (exhaustive)
best_comp, best_union = 0, orig_union
for comp in range(1 << n):
    cand = []
    for f in func_list:
        M = np.eye(n, dtype=np.int64)
        b = np.array([(comp >> j) & 1 for j in range(n)], dtype=np.int64).reshape(-1, 1)
        cand.append(f.substitute_affine(M, b, verify=False))
    u = compute_union_T(cand)
    if u < best_union:
        best_union, best_comp = u, comp

print(f"\n最佳 complement: {best_comp:08b}")
print(f"  union T = {best_union}")

# ΣT after complement
b_arr = np.array([(best_comp >> j) & 1 for j in range(n)], dtype=np.int64).reshape(-1, 1)
comp_funcs = [f.substitute_affine(np.eye(n, dtype=np.int64), b_arr, verify=False) for f in func_list]
comp_sum = sum(f.T() for f in comp_funcs)
# Also count per-output T distribution
from collections import Counter
t_dist = Counter(f.T() for f in comp_funcs)
print(f"  ΣT after complement = {comp_sum}")
print(f"  T distribution: {dict(sorted(t_dist.items()))}")

# Greedy merge
cur_funcs = comp_funcs
M_s = np.eye(n, dtype=np.int64)
improved = True
iterations = 0
while improved and iterations < 50:
    improved = False
    cur_u = compute_union_T(cur_funcs)
    best_i = best_j = -1
    best_mu = cur_u
    active = sorted(set(
        i for f in cur_funcs for mask in f.terms
        for i in range(n) if mask & (1 << i)
    ))
    for i in active:
        for j in active:
            if i == j: continue
            cand = _apply_merge_to_all(cur_funcs, i, j, n)
            u = compute_union_T(cand)
            if u < best_mu:
                best_mu, best_i, best_j = u, i, j
    if best_mu < cur_u:
        cur_funcs = _apply_merge_to_all(cur_funcs, best_i, best_j, n)
        M_s[best_i] = (M_s[best_i] + M_s[best_j]) % 2
        improved = True
        iterations += 1

merge_sum = sum(f.T() for f in cur_funcs)
merge_union = compute_union_T(cur_funcs)
print(f"\nAfter greedy merge ({iterations} merges):")
print(f"  ΣT = {merge_sum}")
print(f"  union T = {merge_union}")

# Also compute what ΣT would be for each output under optimal individual complement (strategy 2)
print(f"\n策略 2 对比:")
print(f"  ΣT = {sum(f.T() for f in func_list)} → {sum(1 for _ in outputs)} (all T=1)")
