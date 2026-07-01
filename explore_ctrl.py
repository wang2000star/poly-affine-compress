"""
探索 ctrl 共享变换能否做得更好。
策略：随机 M,b 搜索 + 穷举所有可能的单行变换组合。
"""
import re, sys, os, time
import numpy as np
from anf_factor import SparseANF, _apply_merge
from simplify_poly import parse_poly, compute_union_T, anf_dict_to_str, _apply_complement_to_all, _apply_merge_to_all

# 加载
inputs, outputs, funcs = parse_poly('/home/wangfz/bool/ctrl.poly')
n = len(inputs)
func_list = [SparseANF(dict(funcs[name]), n) for name in outputs]
orig_union = compute_union_T(func_list)
print(f"原始: n={n}, 输出={len(outputs)}, union T={orig_union}")
print()

# 1. 穷举 complement (2^7 = 128)
print("=" * 60)
print("1. 穷举 Complement 搜索")
print("=" * 60)
best_comp = 0
best_u = orig_union
for comp in range(1 << n):
    cand = _apply_complement_to_all(func_list, comp, n)
    u = compute_union_T(cand)
    if u < best_u:
        best_u = u
        best_comp = comp
        print(f"  complement {comp:07b}: union T={u}/{orig_union}")
print(f"  最佳 complement: {best_comp:07b}, union T={best_u}")

# 2. 对最佳 complement 之后做 greedy merge（模拟 simplify_shared_all 行为）
print()
print("=" * 60)
print("2. Best complement + Greedy merge")
print("=" * 60)
cur_funcs = _apply_complement_to_all(func_list, best_comp, n)
M = np.eye(n, dtype=np.int64)
b = np.array([(best_comp >> j) & 1 for j in range(n)], dtype=np.int64)

# Greedy merge
improved = True
while improved:
    improved = False
    cur_u = compute_union_T(cur_funcs)
    best_i = best_j = -1
    best_merge_u = cur_u

    active = set()
    for f in cur_funcs:
        for mask in f.terms:
            for i in range(n):
                if mask & (1 << i):
                    active.add(i)
    active = sorted(active)

    for i in active:
        for j in active:
            if i == j:
                continue
            cand = _apply_merge_to_all(cur_funcs, i, j, n)
            u = compute_union_T(cand)
            if u <= best_merge_u:
                best_merge_u = u
                best_i, best_j = i, j

    if best_merge_u < cur_u:
        cur_funcs = _apply_merge_to_all(cur_funcs, best_i, best_j, n)
        for col in range(n):
            M[best_i][col] = (M[best_i][col] + M[best_j][col]) % 2
        b[best_i] = (b[best_i] + b[best_j]) % 2
        improved = True
        print(f"  merge z_{best_i}→z_{best_i}+z_{best_j}: union T={best_merge_u}")

print(f"  Greedy merge 后 union T={compute_union_T(cur_funcs)}")

drop_u = compute_union_T(cur_funcs)

# 3. 随机 M,b 搜索
print()
print("=" * 60)
print("3. 随机 M,b 搜索")
print("=" * 60)

rng = np.random.default_rng(42)

# 先试试 complement + 随机行组合
best_random_u = orig_union
best_random = None

# 用 complement 后的函数做 baseline
base_funcs = func_list

for trial in range(20000):
    m = rng.integers(3, 8)
    # 随机 M: 每个元素 0/1
    M_rnd = rng.integers(0, 2, (m, n))
    # 确保行不全为 0
    valid = False
    for row in range(m):
        if np.any(M_rnd[row]):
            valid = True
            break
    if not valid:
        continue
    b_rnd = rng.integers(0, 2, m)

    # 对每个输出做 substitute_affine
    cand_funcs = []
    for f in base_funcs:
        g = f.substitute_affine(M_rnd, b_rnd, verify=False)
        cand_funcs.append(g)
    u = compute_union_T(cand_funcs)

    if u < best_random_u:
        best_random_u = u
        best_random = (M_rnd.copy(), b_rnd.copy())
        print(f"  试验 {trial}: m={m}, union T={u}")

    if u <= 1:
        break

print(f"\n  最佳随机: m={len(best_random[0])}, union T={best_random_u}")
if best_random:
    M_best, b_best = best_random
    print(f"  M:")
    for row in M_best:
        print(f"    {''.join(str(int(x)) for x in row)}")
    print(f"  b: {''.join(str(int(x)) for x in b_best)}")

# 4. 最佳 complement + 随机 M,b
print()
print("=" * 60)
print("4. 最佳 complement + 随机 M,b")
print("=" * 60)

best_comp_random_u = orig_union
best_comp_random = None

for trial in range(20000):
    m = rng.integers(3, 8)
    M_rnd = rng.integers(0, 2, (m, n))
    valid = False
    for row in range(m):
        if np.any(M_rnd[row]):
            valid = True
            break
    if not valid:
        continue
    b_rnd = rng.integers(0, 2, m)

    cand_funcs = []
    for f in base_funcs:
        g = f.substitute_affine(M_rnd, b_rnd, verify=False)
        cand_funcs.append(g)
    u = compute_union_T(cand_funcs)

    if u < best_comp_random_u:
        best_comp_random_u = u
        best_comp_random = (M_rnd.copy(), b_rnd.copy())
        print(f"  试验 {trial}: m={m}, union T={u}")

    if u <= 1:
        break

# 5. 在 complement 后的空间上做随机变换
print()
print("=" * 60)
print("5. Complement 后空间 + 随机 M,b")
print("=" * 60)

comp_funcs = _apply_complement_to_all(func_list, best_comp, n)

best_post_u = compute_union_T(comp_funcs)
best_post = None

for trial in range(20000):
    m = rng.integers(3, 8)
    M_rnd = rng.integers(0, 2, (m, n))
    valid = False
    for row in range(m):
        if np.any(M_rnd[row]):
            valid = True
            break
    if not valid:
        continue
    b_rnd = rng.integers(0, 2, m)

    cand_funcs = []
    for f in comp_funcs:
        g = f.substitute_affine(M_rnd, b_rnd, verify=False)
        cand_funcs.append(g)
    u = compute_union_T(cand_funcs)

    if u < best_post_u:
        best_post_u = u
        best_post = (M_rnd.copy(), b_rnd.copy())
        print(f"  试验 {trial}: m={m}, union T={u}")

    if u <= 1:
        break

print()
print("=" * 60)
print("汇总")
print("=" * 60)
print(f"  Original:               {orig_union}")
print(f"  Best complement:        {best_u}")
print(f"  Complement + greedy:    {drop_u}")
print(f"  Random M,b:             {best_random_u}")
print(f"  Complement + random:    {best_comp_random_u}")
print(f"  Post-complement random: {best_post_u}")
