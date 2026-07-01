"""
.poly 文件 ANF 简化工具

输入: hdNN.poly  — ANF 表达式（XOR, AND, NOT）
输出: hdNN_opt.poly — 简化结果 (g, M, b)

.poly 格式:
  INORDER = i0 i1 ...;
  OUTORDER = om_0 ...;
  om_0 = monom1 + monom2 + ...;   # XOR = +, AND = *, NOT = !
"""
import re, sys, os, time
import numpy as np
from anf_factor import SparseANF, simplify, search_affine_simplification


def parse_poly(path):
    """解析 .poly 文件，返回 (inputs, outputs, {name: anf_dict})"""
    with open(path) as f:
        text = f.read()

    # 去掉注释和换行
    text = re.sub(r'#.*', '', text)
    text = text.replace('\n', ' ').replace('\r', '')
    text = re.sub(r'\s+', ' ', text).strip()

    minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
    mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
    if not minp or not mout:
        raise ValueError("Missing INORDER or OUTORDER")
    inputs = minp.group(1).strip().split()
    outputs = mout.group(1).strip().split()
    n = len(inputs)

    var_index = {v: i for i, v in enumerate(inputs)}
    body_start = max(text.find(';', text.find('INORDER')),
                     text.find(';', text.find('OUTORDER'))) + 1

    funcs = {}
    for stmt in text[body_start:].split(';'):
        stmt = stmt.strip()
        if not stmt:
            continue
        m = re.match(r'(\w+)\s*=\s*(.+)$', stmt)
        if not m:
            continue
        name, expr = m.group(1), m.group(2).strip()
        funcs[name] = _parse_anf_expr(expr, var_index)
    return inputs, outputs, funcs


def _parse_anf_expr(expr, var_index):
    """解析 ANF 表达式（XOR of AND monomials）为 dict {mask: 1}

    支持:
      i0                    → 单变量
      i0*i1*i3              → AND 单项式
      i0 + i1 + i0*i1       → XOR 组合
      !i0                   → 取反（等价于 1 ⊕ i0）
    """
    from anf_factor import SparseANF
    # 按 '+' 分割为 term 列表
    terms = re.split(r'\s*\+\s*', expr)
    result = {}
    for term in terms:
        term = term.strip()
        if not term:
            continue
        mask = _parse_monomial(term, var_index)
        if mask == -1:  # 常数 0
            continue
        result[mask] = result.get(mask, 0) ^ 1
    return {k: v for k, v in result.items() if v}


def _parse_monomial(term, var_index):
    """解析单项式（AND），返回 mask。支持 !var 取反。"""
    # 常数 0
    if term == '0':
        return -1  # sentinel: will be skipped in _parse_anf_expr
    if term == '1':
        return 0   # constant term (mask=0)

    # NOT: !var = 1 ⊕ var → 需要特殊处理
    if term.startswith('!'):
        var = term[1:]
        if var in var_index:
            # !var 在 ANF 中 = 1 ⊕ var
            return 0 | (1 << var_index[var])
        raise ValueError(f"Unknown variable: {var}")

    # AND: i0*i1*i3
    vars_in_term = term.split('*')
    mask = 0
    for v in vars_in_term:
        v = v.strip()
        if not v:
            continue
        if v in var_index:
            mask |= 1 << var_index[v]
        elif v.startswith('!'):
            # !v in a monomial: not directly representable as mask
            # But we handle it as (1 ⊕ v) embedded in the XOR
            # This shouldn't happen since we split on + before
            raise ValueError(f"NOT inside monomial: {term}")
        else:
            raise ValueError(f"Unknown variable: {v} in term {term}")
    return mask


def anf_dict_to_str(terms, var_names):
    """ANF dict → 可读表达式字符串"""
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


def simplify_poly(input_path, output_path, threshold=4096, verbose=True):
    """简化 .poly 文件中的函数"""
    inputs, outputs, funcs = parse_poly(input_path)
    n = len(inputs)

    if verbose:
        print(f"读取: {input_path}")
        print(f"  输入: {n} 变量, 输出: {outputs}")
        for name, terms in funcs.items():
            print(f"  {name}: {len(terms)} 项")

    lines = []
    lines.append("=" * 70)
    lines.append(f"  ANF 简化结果: {os.path.basename(input_path)} → {os.path.basename(output_path)}")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"  原始输入: {', '.join(inputs)}")
    lines.append(f"  输出函数: {', '.join(outputs)}")
    lines.append("")

    total_T0 = 0
    total_T1 = 0
    for name in outputs:
        if name not in funcs:
            print(f"  ⚠ 输出 {name} 未定义，跳过")
            continue

        f = SparseANF(funcs[name], n)
        T0 = f.T()
        total_T0 += T0

        lines.append(f"--- {name} ---")
        lines.append(f"  原始: T={T0}, deg={f.degree()}, vars={f.n}")

        # 处理常数 0（无项）
        if T0 == 0:
            t0 = time.time()
            import numpy as np
            g = SparseANF({}, n)
            M = np.eye(n, dtype=np.int64)
            b = np.zeros(n, dtype=np.int64)
            t1 = time.time()
            T1 = 0
            total_T1 += T1
            m = n
            lines.append(f"  常数 0，无需简化")
            lines.append(f"  耗时: {t1-t0:.3f}s")
            lines.append("")
            z_names = [f"z_{i}" for i in range(m)]
            lines.append(f"  g(z) = 0")
            lines.append("")
            lines.append(f"  验证: 0/100 错误")
            lines.append("")
            lines.append("-" * 50)
            lines.append("")
            continue

        # 简化管线
        t0 = time.time()
        g, M, b = simplify(f, verbose=verbose)
        t1 = time.time()
        T1 = g.T()
        total_T1 += T1
        m = g.n

        ratio = (T0 - T1) / T0 * 100 if T0 else 0
        lines.append(f"  简化: T={T0}→{T1} ({ratio:.1f}%↓), m={n}→{m}")
        lines.append(f"  耗时: {t1-t0:.3f}s")
        lines.append("")

        # 线性变换 z = Mx ⊕ b
        z_names = [f"z_{i}" for i in range(m)]
        lines.append(f"  线性变换 z = Mx ⊕ b (m={m}):")
        for i in range(m):
            terms_z = []
            for j in range(n):
                if int(M[i][j]) % 2:
                    terms_z.append(inputs[j])
            const = int(b[i]) % 2 if i < len(b) else 0
            if const:
                terms_z.append("1")
            lines.append(f"    {z_names[i]} = {' ⊕ '.join(terms_z) if terms_z else '0'}")
        lines.append("")

        # g(z) 表达式
        lines.append(f"  g(z) = g({', '.join(z_names)}) =")
        g_str = anf_dict_to_str(g.terms, z_names)
        lines.append(f"    {g_str}")
        lines.append("")

        # 验证
        errors = 0
        rng = np.random.default_rng(42)
        M_arr = np.array(M, dtype=np.int64)
        b_arr = np.array(b, dtype=np.int64).flatten()
        for _ in range(100):
            x_vec = rng.integers(0, 2, n)
            x_mask = sum(int(x_vec[j]) << j for j in range(n))
            z_vec = (M_arr @ x_vec + b_arr) % 2
            z_mask = sum(int(z_vec[j]) << j for j in range(m))
            f_val = f.eval_mask(x_mask)
            g_val = g.eval_mask(z_mask)
            if f_val != g_val:
                errors += 1
        lines.append(f"  验证: {errors}/100 错误")
        lines.append("")
        lines.append("-" * 50)
        lines.append("")

    # 汇总
    if len(outputs) > 1:
        pct = (total_T0 - total_T1) / total_T0 * 100 if total_T0 else 0
        lines.append(f"  汇总: T₀={total_T0} → T₁={total_T1} ({pct:.1f}%↓)")

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f"\n输出: {output_path}")
    return lines


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python simplify_poly.py <input.poly> [output.poly]")
        sys.exit(1)
    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else re.sub(r'\.poly$', '_opt.poly', input_path)
    simplify_poly(input_path, output_path, verbose=True)
