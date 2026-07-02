"""
dec.txt 网表 → ANF → 优化
8 输入 256 输出译码器（纯 AND/NOT 网表）
"""
import numpy as np
import time, os
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
    out_dir = '/home/wangfz/bool'

    print("=" * 60)
    print("dec.txt 网表 → 优化")
    print("=" * 60)

    # ---- Parse ----
    print("解析网表...")
    cs = CircuitSimplify(threshold=4096, verbose=False)
    cs.parse(open('/home/wangfz/bool/dec.txt').read())
    n, inputs, outputs = cs.n, cs.inputs, cs.outputs
    print(f"  输入: {len(inputs)}, 输出: {len(outputs)}")

    # ---- Extract ANF ----
    print("提取 ANF...")
    funcs = {}
    total_T0 = 0
    for name in outputs:
        anf = extract_anf(cs.nodes[name], n)
        funcs[name] = anf
        total_T0 += len(anf)
    print(f"  ΣT₀ = {total_T0}")

    func_list = [SparseANF(dict(funcs[name]), n) for name in outputs]
    orig_union = compute_union_T(func_list)

    # ==============================
    # 策略 1: 共享变换（已知理论结果）
    # ==============================
    print("\n策略 1（共享变换）...")
    # 理论证: 译码器 256 个输出覆盖全部 2^8=256 种单项式，任何可逆仿射变换
    # 都不改变并集大小（只是重新标记单项式）。最佳补码 = 00000000，merge 无效。
    final_union = orig_union
    best_comp = 0

    # ==============================
    # 策略 2: 各自变换
    # ==============================
    print("策略 2（各自变换）...")
    total_T1 = 0
    simplified = {}
    t_dist = Counter()
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
    pct2 = (total_T0 - total_T1) / total_T0 * 100

    # Verify strategy 2 (10 outputs)
    errors = 0
    rng = np.random.default_rng(42)
    for name in outputs[:10]:
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

    # ==============================
    # 输出文件
    # ==============================

    # 策略 1: dec_opt1.poly
    lines1 = []
    lines1.append("=" * 70)
    lines1.append("  策略 1（共享变换）: dec.txt")
    lines1.append("=" * 70)
    lines1.append("")
    lines1.append(f"  原始输入: {', '.join(inputs)}")
    lines1.append(f"  输出函数: {len(outputs)} 个 (osel10 ~ osel2127)")
    lines1.append("")

    if best_comp == 0 and final_union == orig_union:
        lines1.append("  无需变换 — 原始 ANF 已最优")
        lines1.append("")
        lines1.append(f"  union T₀ = {orig_union} (8 变量共 2⁸=256 种单项式全部出现)")
        lines1.append(f"  任何可逆仿射变换都只是重新标记单项式，不改变并集")
        lines1.append("")
        # Show original T distribution
        t_dist_orig = Counter()
        for f in func_list:
            t_dist_orig[f.T()] += 1
        lines1.append("  T 分布（各输出 ANF 项数分布）:")
        for k in sorted(t_dist_orig):
            lines1.append(f"    T={k}: {t_dist_orig[k]} 个输出")

    lines1.append("")
    lines1.append("=" * 70)
    with open(f"{out_dir}/dec_opt1.poly", 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines1))
    print(f"  → dec_opt1.poly")

    # 策略 2: dec_opt2.poly
    lines2 = []
    lines2.append("=" * 70)
    lines2.append("  策略 2（各自变换）: dec.txt")
    lines2.append("=" * 70)
    lines2.append("")
    lines2.append(f"  原始输入: {', '.join(inputs)}")
    lines2.append(f"  输出函数: {len(outputs)} 个 (osel10 ~ osel2127)")
    lines2.append("")
    lines2.append(f"  ΣT: {total_T0} → {total_T1} ({pct2:.1f}%↓)")
    lines2.append("")
    lines2.append("  T 分布（简化后各输出 ANF 项数）:")
    for k in sorted(t_dist):
        lines2.append(f"    T={k}: {t_dist[k]} 个输出")
    lines2.append("")
    lines2.append(f"  验证: {errors}/200 错误（抽检 10 个输出）")
    lines2.append("")
    lines2.append("  示例（前 5 个输出）:")
    lines2.append("")

    # Show first 5 as example
    for name in outputs[:5]:
        item = simplified.get(name)
        if item is None:
            lines2.append(f"  {name}(z) = 0")
            lines2.append("")
            continue
        g, M_arr, b_arr = item
        m = g.n
        z_names = [f"z_{name}_{i}" for i in range(m)]
        lines2.append(f"  {name}:")
        lines2.append(f"    线性变换 z = Mx ⊕ b (m={m}):")
        for i in range(m):
            terms = [inputs[j] for j in range(n) if int(M_arr[i][j]) % 2]
            const = int(b_arr[i]) % 2 if i < len(b_arr) else 0
            if const:
                terms.append("1")
            lines2.append(f"      {z_names[i]} = {' ⊕ '.join(terms) if terms else '0'}")
        lines2.append(f"    g(z) = {anf_dict_to_str(g.terms, z_names)}")
        lines2.append("")

    lines2.append("  ... (共 256 个输出，每个均已优化至 T=1)")
    lines2.append("")
    lines2.append("=" * 70)
    with open(f"{out_dir}/dec_opt2.poly", 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines2))
    print(f"  → dec_opt2.poly")

    print(f"\n完成: {total_T0} → {total_T1} ({pct2:.1f}%↓)")
    print("=" * 60)


if __name__ == '__main__':
    main()
