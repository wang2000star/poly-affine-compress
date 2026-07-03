"""将简化结果输出为可读的 txt 文件"""
import numpy as np
from anf_factor import SparseANF


def masks_to_anf_str(terms: dict, var_names: list[str]) -> str:
    """ANF dict → 可读表达式字符串"""
    if not terms:
        return "0"
    parts = []
    for mask, coeff in sorted(terms.items(), key=lambda x: bin(x[0])):
        if coeff == 0:
            continue
        if mask == 0:
            parts.append("1")
        else:
            monom = []
            for j in range(len(var_names)):
                if mask & (1 << j):
                    monom.append(var_names[j])
            if monom:
                parts.append("·".join(monom))
            else:
                parts.append("1")
    if not parts:
        return "0"
    return " ⊕ ".join(parts)


def matrix_to_str(M, row_names=None, col_names=None) -> str:
    """矩阵格式化输出"""
    if not M or not M[0]:
        return "(空)"
    lines = []
    for i, row in enumerate(M):
        label = f"{row_names[i]:>6}" if row_names and i < len(row_names) else ""
        bits = " ".join(str(int(x)) for x in row)
        lines.append(f"  {label}  [{bits}]")
    return "\n".join(lines)


def vector_to_str(v, names=None) -> str:
    if not v:
        return "(空)"
    lines = []
    for i, val in enumerate(v):
        label = f"{names[i]:>6}" if names and i < len(names) else ""
        lines.append(f"  {label}   {int(val)}")
    return "\n".join(lines)


def fmt_transform(M, b, orig_names: list[str], new_names: list[str]) -> list[str]:
    """z = Mx ⊕ b 的公式表示"""
    lines = [f"z = Mx ⊕ b, 其中 z ∈ F₂^{len(M)}, x ∈ F₂^{len(b)}"]
    lines.append("")
    lines.append(f"M ({len(M)}×{len(M[0])}):")
    for row_name, row in zip(new_names, M):
        terms = []
        for j, val in enumerate(row):
            if val % 2:
                terms.append(orig_names[j])
        if not terms:
            terms.append("0")
        lines.append(f"  {row_name} = {' ⊕ '.join(terms)}")
    lines.append("")
    lines.append(f"b ({len(b)}):")
    for row_name, bval in zip(new_names, b):
        if int(bval) % 2:
            lines.append(f"  {row_name} = 1")
        else:
            lines.append(f"  {row_name} = 0")
    return lines


def write_simplify_result(
    output_path: str,
    g: SparseANF,
    M,
    b,
    orig_names: list[str] = None,
    output_name: str = "f",
    orig_T: int = None,
    extra: dict = None,
):
    """将简化结果写入 txt 文件

    Args:
        output_path: 输出文件路径
        g: 简化后的 SparseANF (新变量 z)
        M: m×n 变换矩阵
        b: m×1 偏移向量
        orig_names: 原始变量名列表 (默认 z0..zm-1)
        output_name: 输出函数名
        orig_T: 原始 ANF 项数 (可选)
        extra: 额外信息 dict
    """
    m, n = g.n, len(M[0]) if M else 0
    z_names = [f"z{i}" for i in range(m)]
    x_names = orig_names or [f"x{i}" for i in range(n)]

    lines = []
    lines.append("=" * 70)
    lines.append(f"  ANF 简化结果: {output_name}")
    lines.append("=" * 70)
    lines.append("")

    # 基本信息
    lines.append("--- 基本信息 ---")
    lines.append(f"  原始变量数 (n): {n}")
    lines.append(f"  变换后变量数 (m): {m}")
    if orig_T is not None:
        lines.append(f"  原始 ANF 项数: {orig_T}")
    lines.append(f"  简化后 ANF 项数: {g.T()}")
    if orig_T and orig_T > 0:
        pct = (orig_T - g.T()) / orig_T * 100
        lines.append(f"  压缩率: {pct:.1f}%")
    lines.append("")

    # 线性变换
    lines.append("--- 线性变换: z = Mx ⊕ b ---")
    lines.append("")
    for row_name, row in zip(z_names, M):
        terms = []
        for j, val in enumerate(row):
            if int(val) % 2:
                terms.append(x_names[j])
        constant = int(b[len(lines) - len(z_names) - 3]) if len(lines) - 3 < len(b) else 0
        # Actually let me redo this more cleanly
    lines.append("")  # placeholder, will be replaced below

    # Redo cleanly
    del lines[-3:]  # remove the bad attempt

    for i in range(m):
        terms = []
        for j in range(n):
            if int(M[i][j]) % 2:
                terms.append(x_names[j])
        const = int(b[i]) % 2 if i < len(b) else 0
        if const:
            terms.append("1")
        expr = " ⊕ ".join(terms) if terms else "0"
        lines.append(f"  {z_names[i]} = {expr}")

    lines.append("")

    # 简化后的 ANF (在新变量 z 中)
    lines.append(f"--- 简化后 g(z) ---")
    g_str = masks_to_anf_str(g.terms, z_names)
    lines.append(f"  g({', '.join(z_names)}) = {g_str}")
    lines.append("")

    # 展开回原始变量 (如果 m 较小)
    if m <= 10:
        lines.append("--- 展开到原始变量: g(Mx⊕b) ---")
        fx = g.substitute_affine(
            np.array(M, dtype=np.int64),
            np.array(b, dtype=np.int64).reshape(-1, 1),
        )
        fx_str = masks_to_anf_str(fx.terms, x_names)
        lines.append(f"  {output_name}({', '.join(x_names)}) = {fx_str}")
        lines.append("")

    # 验证: 随机测试
    lines.append("--- 正确性验证 ---")
    rng = np.random.default_rng(42)
    errors = 0
    for _ in range(20):
        x_bits = rng.integers(0, 2, n)
        x_mask = sum(int(x_bits[j]) << j for j in range(n))

        # Compute g(z) where z = Mx⊕b
        z_bits = ((np.array(M, dtype=np.int64) @ x_bits + np.array(b, dtype=np.int64)) % 2).flatten()
        z_mask = sum(int(z_bits[j]) << j for j in range(m))

        f_val = (x_mask,)  # placeholder — can't eval without original ANF
        g_val = g.eval_mask(z_mask)

        lines.append(f"  x={''.join(str(int(b)) for b in x_bits)} → z={''.join(str(int(b)) for b in z_bits)} → g(z)={g_val}")
    lines.append(f"  随机测试完成 (如需验证请传入原始函数)")
    lines.append("")

    # 额外信息
    if extra:
        lines.append("--- 额外信息 ---")
        for k, v in extra.items():
            lines.append(f"  {k}: {v}")
        lines.append("")

    lines.append("=" * 70)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return "\n".join(lines)
