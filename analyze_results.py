"""Analyze optimization results: full 6-metric breakdown."""

import os, re

def parse_expr_poly(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        lines = f.read().strip().splitlines()
    result = {}
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        m = re.match(r'(\w+)\s*=\s*(.+)', line)
        if not m:
            continue
        name, expr = m.group(1), m.group(2).strip()
        if expr == '0':
            result[name] = set()
            continue
        masks = set()
        monomials = expr.split('+')
        for mon in monomials:
            mon = mon.strip()
            if not mon:
                continue
            if mon == '1':
                masks.add(0)
                continue
            mask = 0
            for var in mon.split('*'):
                var = var.strip()
                if var.startswith('z_'):
                    mask |= (1 << int(var[2:]))
            masks.add(mask)
        result[name] = masks
    return result


def get_raw_per_output_T(path):
    """Get per-output T values from opt1_T.poly."""
    result = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'(\w+)\s+T=(\d+)', line.strip())
            if m:
                result[m.group(1)] = int(m.group(2))
    return result


instances = ['hd08', 'hd07', 'hd03', 'hd04', 'ctrl', 'dec',
             'int2float', 'hd01', 'hd02', 'cavlc']
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'examples')

print("=" * 140)
header = (f"{'实例':>10} | {'原始累加T':>10} | {'原始并集T':>10} | "
          f"{'opt1累加T':>10} | {'opt1并集T':>10} | "
          f"{'opt2累加T':>10} | {'opt2并集T':>10} | {'压缩比':>8} | {'策略':>6}")
print(header)
print("-" * 140)

for inst in instances:
    d = os.path.join(base, inst)
    raw_T_file = os.path.join(d, f"{inst}_opt1_T.poly")

    # 原始累加T = sum of per-output T from raw parser (stored in opt1 T file)
    raw_per_T = get_raw_per_output_T(raw_T_file)
    raw_sum = sum(raw_per_T.values())

    # 原始并集T: raw parser has independent z per output → union = sum
    raw_union = raw_sum

    # Parse opt1 and opt2 expressions (already in merged z-space)
    opt1_expr = parse_expr_poly(os.path.join(d, f"{inst}_opt1_expr.poly"))
    opt2_expr = parse_expr_poly(os.path.join(d, f"{inst}_opt2_expr.poly"))

    opt1_sum = 0
    opt1_union = 0
    if opt1_expr:
        all_masks = set()
        for masks in opt1_expr.values():
            opt1_sum += len(masks)
            all_masks |= masks
        opt1_union = len(all_masks)

    opt2_sum = 0
    opt2_union = 0
    if opt2_expr:
        all_masks = set()
        for masks in opt2_expr.values():
            opt2_sum += len(masks)
            all_masks |= masks
        opt2_union = len(all_masks)

    # Best strategy: min(opt1并集T, opt2并集T)
    has_opt2 = opt2_expr is not None
    if has_opt2 and opt2_union < opt1_union:
        best_union = opt2_union
        strategy = "opt2"
    else:
        best_union = opt1_union
        strategy = "opt1" if not has_opt2 or opt2_union >= opt1_union else "opt2"

    ratio = raw_union / best_union if best_union > 0 else 1.0

    o2_sum_str = f"{opt2_sum}" if has_opt2 else "—"
    o2_union_str = f"{opt2_union}" if has_opt2 else "—"

    print(f"{inst:>10} | {raw_sum:>10} | {raw_union:>10} | "
          f"{opt1_sum:>10} | {opt1_union:>10} | "
          f"{o2_sum_str:>10} | {o2_union_str:>10} | "
          f"{ratio:>7.2f}× | {strategy:>6}")

print("=" * 140)
print()
print("说明：")
print("  原始累加T = Σ T(f_i)，解析器各输出 ANF 项数之和")
print("  原始并集T = 各输出独立 z 空间下的并集 = 累加T（因独立，无共享）")
print("  opt1累加T = Σ T(opt1 g_i)，解析器输出经 z 合并后的项数和")
print("  opt1并集T = z 合并后所有输出的唯一单项式数（即实际 AND 数）")
print("  opt2累加T = Σ T(opt2 g_i)，逐个简化输出再 z 合并的项数和")
print("  opt2并集T = 同上，唯一单项式数")
print("  压缩比 = 原始并集T / min(opt1并集T, opt2并集T)")
print()
print("注意：原表格的「Opt1 T」列实为累加T而非并集T，现已修正。")
print("      原表格 opt1 和原始并集T相同是因为它们都是累加T（非并集T）。")
