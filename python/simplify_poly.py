"""
.poly 文件 ANF 简化工具

输入: hdNN.poly  — ANF 表达式（XOR, AND, NOT）
输出: hdNN_opt1.poly — 策略 1（共享变换）
      hdNN_opt2.poly — 策略 2（各自变换）

.poly 格式:
  INORDER = i0 i1 ...;
  OUTORDER = om_0 ...;
  om_0 = monom1 + monom2 + ...;   # XOR = +, AND = *, NOT = !

策略 1（共享变换）：所有输出共用同一组 M,b，极小化 union T = |⋃_k terms(g_k)|
策略 2（各自变换）：每个输出独立优化自己的 M_k,b_k，极小化各自 T(g_k)
"""
import re, sys, os, time
import numpy as np
from anf_factor import SparseANF, _apply_merge


# ====================================================================
#  .poly 解析
# ====================================================================

def parse_poly(path):
    with open(path) as f:
        text = f.read()
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
        funcs[m.group(1)] = _parse_anf_expr(m.group(2).strip(), var_index)
    return inputs, outputs, funcs


def _parse_anf_expr(expr, var_index):
    terms = re.split(r'\s*\+\s*', expr)
    result = {}
    for term in terms:
        term = term.strip()
        if not term:
            continue
        mask = _parse_monomial(term, var_index)
        if mask == -1:
            continue
        result[mask] = result.get(mask, 0) ^ 1
    return {k: v for k, v in result.items() if v}


def _parse_monomial(term, var_index):
    if term == '0':
        return -1
    if term == '1':
        return 0
    if term.startswith('!'):
        var = term[1:]
        if var in var_index:
            return 0 | (1 << var_index[var])
        raise ValueError(f"Unknown variable: {var}")
    vars_in_term = term.split('*')
    mask = 0
    for v in vars_in_term:
        v = v.strip()
        if not v:
            continue
        if v in var_index:
            mask |= 1 << var_index[v]
        else:
            raise ValueError(f"Unknown variable: {v} in term {term}")
    return mask


# ====================================================================
#  格式化
# ====================================================================

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


# ====================================================================
#  工具函数
# ====================================================================

def compute_union_T(funcs):
    """计算多个 SparseANF 的 union T"""
    seen = set()
    for f in funcs:
        seen.update(f.terms.keys())
    return len(seen)


def compose(M_sub, b_sub, M_acc, b_acc):
    """组合变换: M_new, b_new = M_sub ∘ (M_acc, b_acc)"""
    M_new = (np.array(M_sub, dtype=np.int64) @ np.array(M_acc, dtype=np.int64)) % 2
    b_sub_a = np.array(b_sub, dtype=np.int64).reshape(-1, 1)
    b_acc_a = np.array(b_acc, dtype=np.int64)
    b_new = (b_sub_a + np.array(M_sub, dtype=np.int64) @ b_acc_a.reshape(-1, 1)) % 2
    return M_new.tolist(), b_new.flatten().tolist()


def mat_to_lists(M):
    return M.tolist() if hasattr(M, 'tolist') else list(M)


# ====================================================================
#  策略 1：共享变换 — 所有分量共用同一组 M,b
# ====================================================================

def _apply_complement_to_all(funcs, complement, n, verbose=False):
    """对所有分量应用 x_i → x_i ⊕ c_i（使用正确的代数代换）"""
    M = np.eye(n, dtype=np.int64)
    b = np.array([(complement >> j) & 1 for j in range(n)], dtype=np.int64).reshape(-1, 1)
    new_funcs = []
    for f in funcs:
        g = f.substitute_affine(M, b, verify=False)
        new_funcs.append(g)
    return new_funcs


def _apply_merge_to_all(funcs, i, j, n):
    """对所有分量应用 z_i → z_i + z_j (GF(2) merge)"""
    new_funcs = []
    for f in funcs:
        new_terms = _apply_merge(f.terms, i, j)
        new_funcs.append(SparseANF(new_terms, n))
    return new_funcs


def simplify_shared_all(funcs, n, max_iter=50, verbose=False):
    """策略 1：共享变换

    Returns (M, b, [g_k]) 其中 M 是 m×n, b 是 m 维, 各 g_k 共享 z = Mx⊕b
    """
    M = np.eye(n, dtype=np.int64)
    b = np.zeros(n, dtype=np.int64)
    cur_funcs = [SparseANF(dict(f.terms), n) for f in funcs]

    orig_union = compute_union_T(cur_funcs)

    # ---- Phase 1: Shared complement search ----
    if verbose:
        print(f"\n  共享 Complement 搜索: n={n}, union T₀={orig_union}")

    best_comp = 0
    best_union = orig_union

    if n <= 16:
        for comp in range(1 << n):
            cand = _apply_complement_to_all(cur_funcs, comp, n)
            u = compute_union_T(cand)
            if u < best_union:
                best_union = u
                best_comp = comp
                if verbose:
                    pct = (orig_union - u) / orig_union * 100
                    print(f"    complement {comp:0{n}b}: union T={u}/{orig_union} ({pct:.1f}%↓)")
    else:
        # Greedy for n > 16
        best_comp = 0
        improved = True
        while improved:
            improved = False
            for i in range(n):
                comp_try = best_comp ^ (1 << i)
                cand = _apply_complement_to_all(cur_funcs, comp_try, n)
                u = compute_union_T(cand)
                if u < best_union:
                    best_union = u
                    best_comp = comp_try
                    improved = True
                    if verbose:
                        pct = (orig_union - u) / orig_union * 100
                        print(f"    greedy flip x{i}: union T={u}/{orig_union} ({pct:.1f}%↓)")

    if best_comp:
        cur_funcs = _apply_complement_to_all(cur_funcs, best_comp, n)
        # Update M, b: z_i = x_i ⊕ c_i → M stays I, b = complement
        b = (b + np.array([(best_comp >> j) & 1 for j in range(n)], dtype=np.int64)) % 2
        if verbose:
            print(f"  Complement 后: union T={best_union}/{orig_union}")

    # ---- Phase 2: Shared greedy merge ----
    if verbose:
        print(f"\n  共享 Greedy merge: n={n}, union T={compute_union_T(cur_funcs)}")

    if orig_union > 1:
        for _ in range(max_iter):
            active = set()
            for f in cur_funcs:
                for mask in f.terms:
                    for i in range(n):
                        if mask & (1 << i):
                            active.add(i)
            active = sorted(active)
            if len(active) <= 1:
                break

            cur_union = compute_union_T(cur_funcs)
            best_i = best_j = -1
            best_u = cur_union

            for i in active:
                for j in active:
                    if i == j:
                        continue
                    cand = _apply_merge_to_all(cur_funcs, i, j, n)
                    u = compute_union_T(cand)
                    if u < best_u:
                        best_u = u
                        best_i, best_j = i, j

            if best_i < 0:
                break

            cur_funcs = _apply_merge_to_all(cur_funcs, best_i, best_j, n)
            # Update M: new_row_i = row_i + row_j (GF: M[best_i][:] += M[best_j][:])
            for col in range(n):
                M[best_i][col] = (M[best_i][col] + M[best_j][col]) % 2
            b[best_i] = (b[best_i] + b[best_j]) % 2

            if verbose:
                pct = (orig_union - best_u) / orig_union * 100
                print(f"    merge x{best_i}→x{best_i}+x{best_j}: union T={best_u}/{orig_union} ({pct:.1f}%↓)")

            if best_u <= 1:
                break

    # ---- Drop unused variables ----
    used = set()
    for f in cur_funcs:
        for mask in f.terms:
            for i in range(n):
                if mask & (1 << i):
                    used.add(i)
    used = sorted(used)

    if len(used) < n and used:
        new_funcs = []
        for f in cur_funcs:
            new_terms = {}
            for mask, coeff in f.terms.items():
                new_mask = 0
                for dst, src in enumerate(used):
                    if mask & (1 << src):
                        new_mask |= 1 << dst
                new_terms[new_mask] = coeff
            new_funcs.append(SparseANF(new_terms, len(used)))
        cur_funcs = new_funcs
        M = M[used]
        n_new = len(used)
    else:
        n_new = n

    if verbose:
        print(f"  最终: m={n_new}, union T={compute_union_T(cur_funcs)}/{orig_union}")

    return cur_funcs, M.tolist(), b.tolist()


# ====================================================================
#  输出生成
# ====================================================================

def _fmt_single(name, f, inputs, n, verbose=False):
    """策略 2 单个输出简化 + 格式化"""
    from anf_factor import simplify as _simplify
    lines = []
    f_sp = SparseANF(dict(f), n)
    T0 = f_sp.T()

    lines.append(f"--- {name} ---")
    lines.append(f"  原始: T={T0}, deg={f_sp.degree()}, vars={f_sp.n}")

    if T0 == 0:
        lines.append(f"  常数 0，无需简化")
        lines.append(f"")
        lines.append(f"  g(z) = 0")
        lines.append(f"")
        lines.append(f"  验证: 0/100 错误")
        lines.append("")
        lines.append("-" * 50)
        lines.append("")
        return lines, T0, 0

    t0 = time.time()
    g, M, b = _simplify(f_sp, verbose=False)
    t1 = time.time()
    T1 = g.T()
    m = g.n

    ratio = (T0 - T1) / T0 * 100 if T0 else 0
    lines.append(f"  简化: T={T0}→{T1} ({ratio:.1f}%↓), m={n}→{m}")
    lines.append(f"  耗时: {t1-t0:.3f}s")
    lines.append("")

    z_names = [f"z_{name}_{i}" for i in range(m)]
    lines.append(f"  线性变换 z = Mx ⊕ b (m={m}):")
    for i in range(m):
        terms_z = [inputs[j] for j in range(n) if int(M[i][j]) % 2]
        const = int(b[i]) % 2 if i < len(b) else 0
        if const:
            terms_z.append("1")
        lines.append(f"    {z_names[i]} = {' ⊕ '.join(terms_z) if terms_z else '0'}")
    lines.append("")

    lines.append(f"  g(z) =")
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
        if f_sp.eval_mask(x_mask) != g.eval_mask(z_mask):
            errors += 1
    lines.append(f"  验证: {errors}/100 错误")
    lines.append("")
    lines.append("-" * 50)
    lines.append("")
    return lines, T0, T1


def _fmt_shared(name, funcs, names, inputs, n, verbose=False):
    """策略 1 共享变换 + 格式化"""
    lines = []
    func_list = [SparseANF(dict(f), n) for f in funcs]
    orig_union = compute_union_T(func_list)

    lines.append(f"--- 共享变换（所有 {len(names)} 个输出）---")
    for nm, f in zip(names, func_list):
        lines.append(f"  {nm}: T={f.T()}, deg={f.degree()}")
    lines.append(f"  union T₀={orig_union}")
    lines.append("")

    t0 = time.time()
    gs, M, b = simplify_shared_all(func_list, n, verbose=False)
    t1 = time.time()
    m = len(M)
    new_union = compute_union_T(gs)

    ratio = (orig_union - new_union) / orig_union * 100 if orig_union else 0
    lines.append(f"  简化: union T={orig_union}→{new_union} ({ratio:.1f}%↓), m={n}→{m}")
    lines.append(f"  耗时: {t1-t0:.3f}s")
    lines.append("")

    z_names = [f"z_{i}" for i in range(m)]
    lines.append(f"  共享线性变换 z = Mx ⊕ b (m={m}):")
    for i in range(m):
        terms_z = [inputs[j] for j in range(n) if int(M[i][j]) % 2]
        const = int(b[i]) % 2 if i < len(b) else 0
        if const:
            terms_z.append("1")
        lines.append(f"    {z_names[i]} = {' ⊕ '.join(terms_z) if terms_z else '0'}")
    lines.append("")

    # 各分量 g_k(z)
    for nm, g in zip(names, gs):
        lines.append(f"  {nm}(z) =")
        g_str = anf_dict_to_str(g.terms, z_names)
        lines.append(f"    {g_str}")
        lines.append("")

    # 验证
    errors = 0
    rng = np.random.default_rng(42)
    M_arr = np.array(M, dtype=np.int64)
    b_arr = np.array(b, dtype=np.int64).flatten()
    for f_orig, g in zip(func_list, gs):
        for _ in range(100):
            x_vec = rng.integers(0, 2, n)
            x_mask = sum(int(x_vec[j]) << j for j in range(n))
            z_vec = (M_arr @ x_vec + b_arr) % 2
            z_mask = sum(int(z_vec[j]) << j for j in range(m))
            if f_orig.eval_mask(x_mask) != g.eval_mask(z_mask):
                errors += 1
    lines.append(f"  验证: {errors}/100 错误（各分量各 100 随机测试）")
    lines.append("")
    lines.append("=" * 50)
    lines.append("")
    return lines, orig_union, new_union


# ====================================================================
#  主入口
# ====================================================================

def simplify_poly(input_path, verbose=True):
    """简化 .poly 文件，输出策略 1 和策略 2"""
    inputs, outputs, funcs = parse_poly(input_path)
    n = len(inputs)

    if verbose:
        print(f"读取: {input_path}")
        print(f"  输入: {n} 变量, 输出: {outputs}")
        for name, terms in funcs.items():
            print(f"  {name}: {len(terms)} 项")

    base = re.sub(r'\.poly$', '', os.path.basename(input_path))

    if len(outputs) == 1:
        # 单输出：只输出一份
        lines, header_lines = [], []
        header_lines.append("=" * 70)
        header_lines.append(f"  ANF 简化结果: {base}.poly")
        header_lines.append("=" * 70)
        header_lines.append("")
        header_lines.append(f"  原始输入: {', '.join(inputs)}")
        header_lines.append(f"  输出函数: {outputs[0]}")
        header_lines.append("")

        l, T0, T1 = _fmt_single(outputs[0], funcs[outputs[0]], inputs, n, verbose)
        lines = header_lines + l

        out_path = f"{base}_opt.poly"
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"\n输出: {out_path}")

    else:
        # 多输出：策略 1 + 策略 2

        # ---- 策略 2（各自变换）----
        lines2 = []
        lines2.append("=" * 70)
        lines2.append(f"  策略 2（各自变换）: {base}.poly")
        lines2.append("=" * 70)
        lines2.append("")
        lines2.append(f"  原始输入: {', '.join(inputs)}")
        lines2.append(f"  输出函数: {', '.join(outputs)}")
        lines2.append("")

        total_T0 = total_T1 = 0
        for name in outputs:
            if name not in funcs:
                continue
            l, T0, T1 = _fmt_single(name, funcs[name], inputs, n, verbose)
            lines2.extend(l)
            total_T0 += T0
            total_T1 += T1

        pct = (total_T0 - total_T1) / total_T0 * 100 if total_T0 else 0
        lines2.append(f"  汇总: T₀={total_T0} → T₁={total_T1} ({pct:.1f}%↓)")

        out2 = f"{base}_opt2.poly"
        with open(out2, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines2))
        print(f"\n输出: {out2}")

        # ---- 策略 1（共享变换）----
        func_list = [funcs[name] for name in outputs if name in funcs]
        out_names = [name for name in outputs if name in funcs]

        lines1 = []
        lines1.append("=" * 70)
        lines1.append(f"  策略 1（共享变换）: {base}.poly")
        lines1.append("=" * 70)
        lines1.append("")
        lines1.append(f"  原始输入: {', '.join(inputs)}")
        lines1.append(f"  输出函数: {', '.join(out_names)}")
        lines1.append("")

        l1, u0, u1 = _fmt_shared(base, func_list, out_names, inputs, n, verbose)
        lines1.extend(l1)

        pct = (u0 - u1) / u0 * 100 if u0 else 0
        lines1.append(f"  汇总: union T₀={u0} → union T₁={u1} ({pct:.1f}%↓)")

        out1 = f"{base}_opt1.poly"
        with open(out1, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines1))
        print(f"\n输出: {out1}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python simplify_poly.py <input.poly>")
        sys.exit(1)
    simplify_poly(sys.argv[1], verbose=True)
