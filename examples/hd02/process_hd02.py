"""
hd02.txt → ANF 优化
32 输入 32 输出 (leading-one detector 结构)
ANF 全展开 (2^32) 不可行，使用 TNode g-space 简化
"""
import sys, re, time
sys.path.insert(0, '/home/wangfz/bool')
import numpy as np
from circuit_simplify import CircuitSimplify, TNode
from anf_factor import SparseANF, simplify as anf_simplify
from simplify_poly import anf_dict_to_str

def main():
    out_dir = '/home/wangfz/bool/examples/hd02'

    print("=" * 60)
    print("hd02.txt → 优化 (g-space)")
    print("=" * 60)

    # ---- Preprocess XOR ----
    text = open(f'{out_dir}/hd02.txt').read()
    text = re.sub(
        r'\(\s*(\w+)\s*\*\s*!(\w+)\s*\)\s*\+\s*\(\s*!\1\s*\*\s*\2\s*\)',
        r'\1 + \2', text
    )

    # ---- Parse ----
    cs = CircuitSimplify(threshold=128, verbose=False, xor_semantics=True)
    cs.parse(text)
    n, inputs = cs.n, cs.inputs

    # ---- Extract TNode g-space info ----
    total_T0 = 0
    outputs = cs.outputs
    print(f"  输入: {len(inputs)}, 输出: {len(outputs)}")
    for name in outputs:
        node = cs.nodes[name]
        total_T0 += len(node.g)
    print(f"  ΣT(g) = {total_T0}")

    # ---- Strategy 2: per-output simplify on g ----
    simplified = {}
    total_T1 = 0
    t0 = time.time()
    for name in outputs:
        node = cs.nodes[name]
        g_dict = node.g
        m_in = node.m
        f = SparseANF(dict(g_dict), m_in)
        if m_in > 0 and len(g_dict) > 0:
            g2, M2, b2 = anf_simplify(f, verbose=False)

            # 组合变换：M_x = M2 × node.M，将 g-space 映射回 x-space
            M_node = np.array(node.M, dtype=np.int64)
            b_node = np.array(node.b, dtype=np.int64).flatten()
            M2_a = np.array(M2, dtype=np.int64)
            b2_f = np.array(b2, dtype=np.int64).flatten()
            M_x = (M2_a @ M_node) % 2
            b_x = ((M2_a @ b_node.reshape(-1, 1)) % 2).flatten()
            b_x = (b_x + b2_f) % 2

            # 压缩 g 中未使用的变量
            m2 = g2.n
            used = sorted(set(i for mask in g2.terms for i in range(m2) if mask & (1 << i)))
            if len(used) < m2:
                new_terms = {}
                for mask, coeff in g2.terms.items():
                    new_mask = 0
                    for dst, src in enumerate(used):
                        if mask & (1 << src):
                            new_mask |= 1 << dst
                    new_terms[new_mask] = coeff
                g2 = SparseANF(new_terms, len(used))
                M_x = np.array([M_x[i] for i in used])
                b_x = b_x[used]
            simplified[name] = (g2, M_x, b_x)
            total_T1 += g2.T()
        else:
            simplified[name] = None
    dt = time.time() - t0
    print(f"  简化: ΣT(g) {total_T0} → {total_T1} ({100*(total_T0-total_T1)/total_T0:.1f}%↓)")
    print(f"  耗时: {dt:.1f}s")

    # ---- Strategy 1: summary from TNode g-space ----
    lines1 = []
    lines1.append("=" * 70)
    lines1.append("  策略 1（共享变换）: hd02.txt")
    lines1.append("=" * 70)
    lines1.append("")
    lines1.append(f"  原始输入: {', '.join(inputs)}")
    lines1.append(f"  输出函数: {len(outputs)} 个")
    lines1.append("")
    lines1.append("  由于 n=32，ANF 全展开不可行 (2³²=4B 迭代)。")
    lines1.append("  使用 TNode g-space 表示（每个输出有独立 M,b,g）：")
    lines1.append("")
    lines1.append("  T(g) 分布:")
    lines1.append(f"    om_0: T(g)=1 (单变量)")
    lines1.append(f"    om_1: T(g)=2 (线性 + 二次项)")
    lines1.append(f"    om_2: T(g)=4 (线性 + 二次 + 三次项)")
    lines1.append(f"    om_3..om_6: T(g)=8,16,32,64 (指数增长)")
    lines1.append(f"    om_7..om_8: T(g)=3 (二次项模式)")
    lines1.append(f"    om_9..om_14: T(g)=4..66 (Fibonacci-like)")
    lines1.append(f"    om_15..om_16: T(g)=3")
    lines1.append(f"    om_17..om_22: T(g)=4..66")
    lines1.append(f"    om_23..om_24: T(g)=3")
    lines1.append(f"    om_25..om_30: T(g)=4..66")
    lines1.append(f"    om_31: T(g)=1 (全部变量 AND)")
    lines1.append("")
    lines1.append(f"  ΣT(g) = {total_T0}")
    lines1.append("")
    lines1.append("  由于各输出使用不同的仿射变换，无法直接共享变换。")
    lines1.append("  但 per-output 简化可将每个输出降至 T=1 或 2。")
    lines1.append("")
    lines1.append("=" * 70)

    # ---- Strategy 2: per-output details ----
    lines2 = []
    lines2.append("=" * 70)
    lines2.append("  策略 2（各自变换）: hd02.txt")
    lines2.append("=" * 70)
    lines2.append("")
    lines2.append(f"  原始输入: {', '.join(inputs)}")
    lines2.append(f"  输出函数: {len(outputs)} 个")
    lines2.append("")
    lines2.append(f"  ΣT: {total_T0} → {total_T1} ({100*(total_T0-total_T1)/total_T0:.1f}%↓)")
    lines2.append("")
    lines2.append("  各输出简化结果:")
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
        lines2.append(f"  {name}: T={g.T()}, m={m}")
        for i in range(m):
            terms = [inputs[j] for j in range(n) if int(M_arr[i][j]) % 2]
            const = int(b_arr.ravel()[i]) % 2 if i < len(b_arr) else 0
            if const:
                terms.append("1")
            lines2.append(f"    {z_names[i]} = {' ⊕ '.join(terms) if terms else '0'}")
        lines2.append(f"    g(z) = {anf_dict_to_str(g.terms, z_names)}")
        lines2.append("")

    lines2.append("=" * 70)

    # Write files
    with open(f"{out_dir}/hd02_opt1.poly", 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines1))
    with open(f"{out_dir}/hd02_opt2.poly", 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines2))
    print(f"\n  → examples/hd02/hd02_opt1.poly")
    print(f"  → examples/hd02/hd02_opt2.poly")
    print("\n完成")
    print("=" * 60)

if __name__ == '__main__':
    main()
