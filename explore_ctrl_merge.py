"""
探索 ctrl 共享变换多步 merge 搜索。

有效变换空间: 补码 (z_i → z_i ⊕ c_i) + merge (z_i → z_i ⊕ z_j)
两者都是可逆的，保证所有函数因子存在。
"""
import re, sys, os, time
import numpy as np
from anf_factor import SparseANF, _apply_merge
from simplify_poly import parse_poly, compute_union_T, _apply_complement_to_all, _apply_merge_to_all


def all_merge_pairs(n):
    pairs = []
    for i in range(n):
        for j in range(n):
            if i != j:
                pairs.append((i, j))
    return pairs


def main():
    inputs, outputs, funcs = parse_poly('/home/wangfz/bool/ctrl.poly')
    n = len(inputs)
    base_funcs = [SparseANF(dict(funcs[name]), n) for name in outputs]
    orig_union = compute_union_T(base_funcs)
    print(f"原始: n={n}, 输出={len(outputs)}, union T={orig_union}")

    all_merges = all_merge_pairs(n)

    # ---- 多步 merge 搜索 ----
    # 对每个 complement:
    #   1. 尝试 0 步 merge (baseline)
    #   2. 尝试 1 步 merge
    #   3. 尝试 2 步 merge (best_1 + all 2nd)
    # 最终取最小

    print()
    print("=" * 60)
    print("多步 Merge 搜索（所有 complement × 最多 2 步 merge）")
    print("=" * 60)

    best_overall_u = orig_union
    best_overall = None

    for comp in range(1 << n):
        # Apply complement
        cur = _apply_complement_to_all(base_funcs, comp, n)
        cur_u = compute_union_T(cur)
        info = f"comp={comp:07b}"

        # Track best for this complement
        comp_best_u = cur_u
        comp_best_seq = []

        # Try 1-step merges
        for i, j in all_merges:
            cand1 = _apply_merge_to_all(cur, i, j, n)
            u1 = compute_union_T(cand1)
            if u1 < comp_best_u:
                comp_best_u = u1
                comp_best_seq = [(i, j)]

        # Try 2-step merges: for each (i,j) that improves (or matches), try again
        # To keep it fast, only try best K 1-step merges per complement
        candidates_1 = [(i, j) for i, j in all_merges]
        u1_vals = []
        for i, j in candidates_1:
            cand1 = _apply_merge_to_all(cur, i, j, n)
            u1_vals.append((i, j, compute_union_T(cand1)))
        # Sort by union T, take top 5
        u1_vals.sort(key=lambda x: x[2])
        for i, j, u1 in u1_vals[:8]:
            cand1 = _apply_merge_to_all(cur, i, j, n)
            for k, l in all_merges:
                if k == l:
                    continue
                cand2 = _apply_merge_to_all(cand1, k, l, n)
                u2 = compute_union_T(cand2)
                if u2 < comp_best_u:
                    comp_best_u = u2
                    comp_best_seq = [(i, j), (k, l)]

        if comp_best_u < best_overall_u:
            best_overall_u = comp_best_u
            best_overall = (comp, comp_best_seq)
            print(f"  新最佳: {info}, seq={comp_best_seq}, union T={comp_best_u}")

    print()
    print("=" * 60)
    print(f"最佳结果: union T = {best_overall_u}")
    if best_overall:
        comp, seq = best_overall
        print(f"  complement: {comp:07b}")
        print(f"  merge 序列: {seq}")

    # ---- 验证最佳结果 ----
    comp, seq = best_overall
    cur = _apply_complement_to_all(base_funcs, comp, n)
    M = np.eye(n, dtype=np.int64)
    b = np.array([(comp >> j) & 1 for j in range(n)], dtype=np.int64)
    for i, j in seq:
        cur = _apply_merge_to_all(cur, i, j, n)
        M[i] = (M[i] + M[j]) % 2
        b[i] = (b[i] + b[j]) % 2

    # Drop unused
    used = set()
    for f in cur:
        for mask in f.terms:
            for i in range(n):
                if mask & (1 << i):
                    used.add(i)
    used = sorted(used)

    print(f"\n  最终 m={len(used)} (原 n={n})")
    print(f"  共享线性变换:")
    for idx in used:
        terms = [inputs[j] for j in range(n) if M[idx][j]]
        const = b[idx] % 2
        if const:
            terms.append("1")
        print(f"    z_{idx} = {' ⊕ '.join(terms) if terms else '0'}")

    # 各输出
    for name, f in zip(outputs, cur):
        g_str = anf_dict_to_str(f.terms, [f"z_{i}" for i in range(n)])
        print(f"  {name}(z) = {g_str}")

    # 验证
    errors = 0
    rng = np.random.default_rng(42)
    M_arr = np.array(M, dtype=np.int64)
    b_arr = b.flatten()
    for f_orig, f_new in zip(base_funcs, cur):
        for _ in range(100):
            x_vec = rng.integers(0, 2, n)
            x_mask = sum(int(x_vec[j]) << j for j in range(n))
            z_vec = (M_arr @ x_vec + b_arr) % 2
            z_mask = sum(int(z_vec[j]) << j for j in range(len(used)))
            if f_orig.eval_mask(x_mask) != f_new.eval_mask(z_mask):
                errors += 1
    print(f"  验证: {errors}/100 错误")


def anf_dict_to_str(terms, var_names):
    if not terms:
        return "0"
    parts = []
    for mask, coeff in sorted(terms.items(), key=lambda x: (bin(x[0]).count('1'), x[0])):
        if coeff == 0:
            continue
        if mask == 0:
            parts.append("1")
        else:
            monom = [var_names[j] for j in range(len(var_names)) if mask & (1 << j)]
            parts.append("·".join(monom) if monom else "1")
    return " ⊕ ".join(parts) if parts else "0"


if __name__ == '__main__':
    main()
