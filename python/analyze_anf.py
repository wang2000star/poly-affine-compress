"""正确计算原始 x 空间 ANF T 值，并与优化后 z 空间 T 对比。"""

import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from verify_opt import make_circuit_func, make_transform_func


def build_truth_table(circuit_func, inputs, outputs):
    n = len(inputs)
    tt = {o: [] for o in outputs}
    for x_bits in range(1 << n):
        x_dict = {inputs[i]: (x_bits >> i) & 1 for i in range(n)}
        for o, v in zip(outputs, circuit_func(x_dict)):
            tt[o].append(v)
    return tt


def truth_table_to_anf(tt_vec):
    """Moebius 变换：真值表 → ANF 系数 {mask: coeff}."""
    n = (len(tt_vec).bit_length() - 1)
    a = list(tt_vec)
    for i in range(n):
        step = 1 << i
        for j in range(0, len(a), step * 2):
            for k in range(j, j + step):
                a[k + step] ^= a[k]
    return {k: v for k, v in enumerate(a) if v}


def parse_expr_masks(expr_path, outputs):
    """从 _expr.poly 解析每个输出的 z 空间 mask 集合。"""
    if not os.path.exists(expr_path):
        return None
    with open(expr_path) as f:
        lines = f.read().strip().splitlines()
    result = {}
    for line in lines:
        line = line.strip()
        if not line or line[0] == '#':
            continue
        m = re.match(r'(\w+)\s*=\s*(.+)', line)
        if not m:
            continue
        oname, expr = m.group(1), m.group(2).strip()
        if oname not in outputs:
            continue
        if expr == '0':
            result[oname] = set()
            continue
        masks = set()
        for mon in expr.split('+'):
            mon = mon.strip()
            if mon == '1':
                masks.add(0)
                continue
            mask = 0
            for var in mon.split('*'):
                var = var.strip()
                mask |= (1 << int(var[2:]))
            masks.add(mask)
        result[oname] = masks
    return result


def analyze_instance(inst_name, base_path):
    txt_path = os.path.join(base_path, inst_name, f"{inst_name}.txt")
    if not os.path.exists(txt_path):
        return None

    inputs, outputs, circuit_f = make_circuit_func(txt_path)
    n = len(inputs)

    if n > 16:
        print(f"  {inst_name}: n={n}>16, skip exhaustive ANF")
        return None

    print(f"  Truth table and ANF (2^{n} = {1<<n})...")
    tt = build_truth_table(circuit_f, inputs, outputs)

    raw_sum = 0
    raw_all = set()
    for o in outputs:
        anf = truth_table_to_anf(tt[o])
        raw_sum += len(anf)
        raw_all |= set(anf.keys())

    print(f"    raw sum T = {raw_sum}")
    print(f"    raw union T = {len(raw_all)} (<= 2^{n} = {1<<n})")

    for strategy in ['opt1', 'opt2']:
        expr_path = os.path.join(base_path, inst_name, f"{inst_name}_{strategy}_expr.poly")
        if not os.path.exists(expr_path):
            continue
        masks = parse_expr_masks(expr_path, outputs)
        if masks is None:
            continue
        opt_sum = sum(len(m) for m in masks.values())
        opt_union = len(set().union(*masks.values())) if masks else 0
        print(f"    {strategy}: sum T = {opt_sum}, union T = {opt_union}")

    return raw_sum, len(raw_all)


def main():
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'examples')
    for inst in ['hd08', 'hd07', 'hd03', 'hd04', 'ctrl', 'dec', 'int2float', 'cavlc']:
        print(f"\n--- {inst} ---")
        analyze_instance(inst, base)


if __name__ == '__main__':
    main()
