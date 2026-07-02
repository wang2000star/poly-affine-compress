"""
ctrl.txt → ctrl.poly → ctrl_opt1/2.poly

ctrl.txt 是 AND/NOT/XOR 网表（+ 表示 XOR）。
流程：
1. 用 CircuitSimplify(xor_semantics=True) 解析网表
2. 对每个输出提取 ANF（真值表 → Mobius 变换）
3. 写出 ctrl.poly
4. 运行 simplify_poly 获得两种策略的优化结果
"""
import numpy as np
from circuit_simplify import CircuitSimplify
from anf_factor import SparseANF


def extract_anf(node, n):
    """从 TNode 提取 f(x) ANF: 计算 2^n 真值表后做 Mobius 变换"""
    M = np.array(node.M, dtype=np.int64)
    b = np.array(node.b, dtype=np.int64)
    g = SparseANF(node.g, node.m)

    tt_len = 1 << n
    tt = np.zeros(tt_len, dtype=np.int64)

    # 真值表
    for x_mask in range(tt_len):
        x_arr = np.array([(x_mask >> i) & 1 for i in range(n)], dtype=np.int64)
        z = (M @ x_arr + b) % 2
        z_mask = sum(int(z[j]) << j for j in range(node.m))
        tt[x_mask] = g.eval_mask(z_mask)

    # 快速 Mobius 变换 → ANF 系数
    anf_mask = tt.copy()
    for i in range(n):
        step = 1 << i
        for j in range(0, tt_len, step << 1):
            for k in range(step):
                anf_mask[j + k + step] ^= anf_mask[j + k]

    anf = {}
    for mask, val in enumerate(anf_mask):
        if val:
            anf[mask] = 1
    return anf


def main():
    path = '/home/wangfz/bool/ctrl/ctrl.txt'
    out_dir = '/home/wangfz/bool'

    print("=" * 60)
    print("ctrl.txt 网表解析及 ANF 简化")
    print("=" * 60)
    print()

    text = open(path).read()
    cs = CircuitSimplify(threshold=4096, verbose=True, xor_semantics=True)
    cs.parse(text)

    n = cs.n
    outputs = cs.outputs
    inputs = cs.inputs

    print(f"\n输入: {len(inputs)} 变量")
    print(f"输出: {len(outputs)} 个")
    print()

    # 提取每个输出的 ANF
    funcs = {}
    for name in outputs:
        node = cs.nodes[name]
        anf = extract_anf(node, n)
        funcs[name] = anf
        print(f"  {name}: T={len(anf):>4}", end="")
        if not anf:
            print(" (常数 0)")
        elif anf == {0: 1}:
            print(" (常数 1)")
        else:
            print()

    # 写出 ctrl.poly
    poly_path = f"{out_dir}/ctrl/ctrl.poly"
    with open(poly_path, 'w') as f:
        f.write(f"INORDER = {' '.join(inputs)};\n")
        f.write(f"OUTORDER = {' '.join(outputs)};\n")
        for name in outputs:
            anf = funcs[name]
            if not anf:
                f.write(f"{name} = 0;\n")
            else:
                terms = []
                for mask, coeff in sorted(anf.items(), key=lambda x: (bin(x[0]).count('1'), x[0])):
                    if coeff == 0:
                        continue
                    if mask == 0:
                        terms.append("1")
                    else:
                        monom = [inputs[j] for j in range(n) if mask & (1 << j)]
                        terms.append("*".join(monom))
                f.write(f"{name} = {' + '.join(terms)};\n")

    print(f"\n输出: {poly_path}")

    # 运行 simplify_poly
    from simplify_poly import simplify_poly
    print("\n" + "=" * 60)
    print("运行 simplify_poly...")
    print("=" * 60)
    simplify_poly(poly_path, verbose=True)


if __name__ == '__main__':
    main()
