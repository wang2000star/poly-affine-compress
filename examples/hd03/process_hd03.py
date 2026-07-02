"""
hd03.txt → ANF → 优化
16 输入 8 输出（hd03 电路）
"""
import sys
sys.path.insert(0, '/home/wangfz/bool')
import numpy as np
import time, os, re
from circuit_simplify import CircuitSimplify
from anf_factor import SparseANF, simplify as anf_simplify
from simplify_poly import compute_union_T, _apply_merge_to_all, anf_dict_to_str
from collections import Counter


def extract_anf(node, n):
    M = np.array(node.M, dtype=np.int64)
    b = np.array(node.b, dtype=np.int64)
    g = SparseANF(node.g, node.m)
    tt_len = 1 << n
    tt = np.zeros(tt_len, dtype=np.int64)
    for x_mask in range(tt_len):
        x_arr = np.array([(x_mask >> i) & 1 for i in range(n)], dtype=np.int64)
        z = (M @ x_arr + b) % 2
        z_mask = sum(int(z[j]) << j for j in range(node.m))
        tt[x_mask] = g.eval_mask(z_mask)
    anf_mask = tt.copy()
    for i in range(n):
        step = 1 << i
        for j in range(0, tt_len, step << 1):
            for k in range(step):
                anf_mask[j + k + step] ^= anf_mask[j + k]
    return {mask: 1 for mask, val in enumerate(anf_mask) if val}


def main():
    out_dir = '/home/wangfz/bool/hd03'

    print("=" * 60)
    print("hd03.txt → 优化")
    print("=" * 60)

    # ---- Preprocess XOR expressions ----
    print("预处理 XOR 模式...")
    text = open(f'{out_dir}/hd03.txt').read()
    text = re.sub(
        r'\(\s*(\w+)\s*\*\s*!(\w+)\s*\)\s*\+\s*\(\s*!\1\s*\*\s*\2\s*\)',
        r'\1 + \2',
        text
    )

    # ---- Parse ----
    print("解析网表 (xor_semantics=True)...")
    cs = CircuitSimplify(threshold=4096, verbose=False, xor_semantics=True)
    cs.parse(text)
    n, inputs, outputs = cs.n, cs.inputs, cs.outputs
    print(f"  输入: {len(inputs)}, 输出: {len(outputs)}")
    print(f"  输入: {', '.join(inputs)}")
    print(f"  输出: {', '.join(outputs)}")

    # ---- Extract ANF ----
    print("提取 ANF...")
    funcs = {}
    total_T0 = 0
    t0 = time.time()
    for name in outputs:
        anf = extract_anf(cs.nodes[name], n)
        funcs[name] = anf
        total_T0 += len(anf)
    dt = time.time() - t0
    print(f"  ΣT₀ = {total_T0}, 耗时 {dt:.1f}s")
    for name in outputs:
        print(f"    {name}: T={len(funcs[name])}")

    func_list = [SparseANF(dict(funcs[name]), n) for name in outputs]
    orig_union = compute_union_T(func_list)
    print(f"  union T₀ = {orig_union}")

    # ==============================
    # 策略 1: 共享变换
    # ==============================
    print(f"\n策略 1（共享变换）...")

    # Phase 1: greedy complement search
    best_comp, best_union = 0, orig_union
    cur = list(func_list)
    improved = True
    while improved:
        improved = False
        best_var = -1
        for var in range(n):
            if best_comp & (1 << var):
                continue  # already flipped
            test_comp = best_comp | (1 << var)
            bc = np.array([(test_comp >> j) & 1 for j in range(n)], dtype=np.int64).reshape(-1, 1)
            cand = [f.substitute_affine(np.eye(n), bc, verify=False) for f in func_list]
            u = compute_union_T(cand)
            if u < best_union:
                best_union, best_comp, best_var = u, test_comp, var
                improved = True
        if improved:
            print(f"    flip var {best_var}: union T = {best_union}")
    print(f"  最佳 complement: {best_comp:016b}, union T: {orig_union} → {best_union}")

    # Build accumulated M (m×n), b (m)
    M1 = np.eye(n, dtype=np.int64)
    b1 = np.array([(best_comp >> j) & 1 for j in range(n)], dtype=np.int64)

    cur_funcs = [SparseANF(dict(funcs[name]), n) for name in outputs]
    if best_comp:
        cur_funcs = [f.substitute_affine(np.eye(n), b1.reshape(-1, 1), verify=False) for f in cur_funcs]

    # Phase 2: greedy merge
    improved = True
    merge_iter = 0
    while improved and merge_iter < 50:
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
            for col in range(n):
                M1[best_i][col] = (M1[best_i][col] + M1[best_j][col]) % 2
            b1[best_i] = (b1[best_i] + b1[best_j]) % 2
            improved = True
            merge_iter += 1

    # Drop unused variables
    used = sorted(set(
        i for f in cur_funcs for mask in f.terms
        for i in range(n) if mask & (1 << i)
    ))
    m1 = len(used)
    if m1 < n:
        M1 = M1[used]
        b1 = b1[used]
        new_funcs = []
        for f in cur_funcs:
            new_terms = {}
            for mask, coeff in f.terms.items():
                new_mask = 0
                for dst, src in enumerate(used):
                    if mask & (1 << src):
                        new_mask |= 1 << dst
                new_terms[new_mask] = coeff
            new_funcs.append(SparseANF(new_terms, m1))
        cur_funcs = new_funcs

    final_union_1 = compute_union_T(cur_funcs)
    sum_T_1 = sum(f.T() for f in cur_funcs)
    print(f"  greedy merge ({merge_iter} iter): m={n}→{m1}, union T = {final_union_1}")
    print(f"  ΣT after strategy 1 = {sum_T_1}")

    # ==============================
    # 策略 2: 各自变换
    # ==============================
    print(f"\n策略 2（各自变换）...")
    total_T1 = 0
    simplified = {}
    t_dist = Counter()
    t0 = time.time()
    for name in outputs:
        anf = funcs[name]
        if not anf:
            simplified[name] = None
            continue
        f = SparseANF(dict(anf), n)
        g, M_arr, b_arr = anf_simplify(f, verbose=False)
        simplified[name] = (g, M_arr, b_arr)
        total_T1 += g.T()
        t_dist[g.T()] += 1
    dt = time.time() - t0
    pct2 = (total_T0 - total_T1) / total_T0 * 100

    print(f"  ΣT: {total_T0} → {total_T1} ({pct2:.1f}%↓)")
    print(f"  耗时: {dt:.1f}s")
    print(f"  T 分布（简化后）:")
    for k in sorted(t_dist):
        print(f"    T={k}: {t_dist[k]} 个输出")

    # Verify strategy 2
    errors = 0
    rng = np.random.default_rng(42)
    for name in outputs:
        item = simplified.get(name)
        if item is None: continue
        g, M_arr, b_arr = item
        f = SparseANF(dict(funcs[name]), n)
        m = g.n
        for _ in range(20):
            x = rng.integers(0, 2, n)
            xm = sum(int(x[j]) << j for j in range(n))
            z = (M_arr @ x + b_arr.flatten()) % 2
            zm = sum(int(z[j]) << j for j in range(m))
            if f.eval_mask(xm) != g.eval_mask(zm):
                errors += 1
    print(f"  验证: {errors}/{len(outputs)*20} 错误")

    # ==============================
    # 输出文件: 策略 1（共享变换）
    # ==============================
    z_names_1 = [f"z_{i}" for i in range(m1)]

    lines1 = []
    lines1.append("=" * 70)
    lines1.append("  策略 1（共享变换）: hd03.txt")
    lines1.append("=" * 70)
    lines1.append("")
    lines1.append(f"  原始输入: {', '.join(inputs)}")
    lines1.append(f"  输出函数: {len(outputs)} 个")
    lines1.append("")
    lines1.append(f"  union T: {orig_union} → {final_union_1}")
    lines1.append(f"  ΣT: {total_T0} → {sum_T_1}")
    lines1.append("")
    lines1.append("  简化过程:")
    lines1.append(f"    最佳 complement: {best_comp:016b}")
    lines1.append(f"    greedy merge: {merge_iter} 次")
    lines1.append(f"    最终变量数: m={m1}")
    lines1.append("")
    lines1.append(f"  共享线性变换 z = Mx ⊕ b (m={m1}):")
    for i in range(m1):
        terms = [inputs[j] for j in range(n) if M1[i][j] % 2]
        const = b1[i] % 2
        if const:
            terms.append("1")
        lines1.append(f"    {z_names_1[i]} = {' ⊕ '.join(terms) if terms else '0'}")
    lines1.append("")
    lines1.append("  各输出 g(z)（策略 1 后）:")
    lines1.append("")
    for name, f in zip(outputs, cur_funcs):
        lines1.append(f"  {name}: T={len(funcs[name])} → {f.T()}")
        g_str = anf_dict_to_str(f.terms, z_names_1)
        lines1.append(f"    {g_str}")
        lines1.append("")
    lines1.append("=" * 70)

    # ==============================
    # 输出文件: 策略 2
    # ==============================
    lines2 = []
    lines2.append("=" * 70)
    lines2.append("  策略 2（各自变换）: hd03.txt")
    lines2.append("=" * 70)
    lines2.append("")
    lines2.append(f"  原始输入: {', '.join(inputs)}")
    lines2.append(f"  输出函数: {len(outputs)} 个")
    lines2.append("")
    lines2.append(f"  ΣT: {total_T0} → {total_T1} ({pct2:.1f}%↓)")
    lines2.append("")
    lines2.append("  T 分布（简化后）:")
    for k in sorted(t_dist):
        lines2.append(f"    T={k}: {t_dist[k]} 个输出")
    lines2.append("")
    lines2.append(f"  验证: {errors}/{len(outputs)*20} 错误")
    lines2.append("")
    lines2.append("  各输出详情:")
    lines2.append("")

    for name in outputs:
        item = simplified.get(name)
        if item is None:
            lines2.append(f"  {name}(z) = 0")
            lines2.append("")
            continue
        g, M_arr, b_arr = item
        m = g.n
        z_names = [f"z_{name}_{i}" for i in range(m)]
        lines2.append(f"  {name}: T={len(funcs[name])} → {g.T()}, m={n}→{m}")
        lines2.append(f"    线性变换 z = Mx ⊕ b:")
        for i in range(m):
            terms = [inputs[j] for j in range(n) if int(M_arr[i][j]) % 2]
            const = int(b_arr[i]) % 2 if i < len(b_arr) else 0
            if const:
                terms.append("1")
            lines2.append(f"      {z_names[i]} = {' ⊕ '.join(terms) if terms else '0'}")
        lines2.append(f"    g(z) = {anf_dict_to_str(g.terms, z_names)}")
        lines2.append("")

    lines2.append("=" * 70)

    # Write files
    with open(f"{out_dir}/hd03_opt1.poly", 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines1))
    print(f"\n  → hd03/hd03_opt1.poly")

    with open(f"{out_dir}/hd03_opt2.poly", 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines2))
    print(f"  → hd03/hd03_opt2.poly")

    print(f"\n完成: ΣT {total_T0} → {total_T1} ({pct2:.1f}%↓)")
    print(f"      union T {orig_union} → {final_union_1}")
    print("=" * 60)


if __name__ == '__main__':
    main()
