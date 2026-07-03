"""Compute max ANF degree from expr.poly files for completed instances."""

import os, re

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'examples')

def parse_expr_degree(path):
    """Read expr.poly, return {output_name: max_degree} and overall max."""
    with open(path) as f:
        lines = f.read().strip().splitlines()
    degs = {}
    overall = 0
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        m = re.match(r'(\w+)\s*=\s*(.+)', line)
        if not m:
            continue
        oname, expr = m.group(1), m.group(2).strip()
        if expr == '0':
            degs[oname] = 0
            continue
        if 'not linear-form representable' in expr:
            continue
        max_d = 0
        for mon in expr.split('+'):
            mon = mon.strip()
            if not mon:
                continue
            if mon == '1':
                d = 0
            else:
                d = mon.count('*') + 1
            if d > max_d:
                max_d = d
        degs[oname] = max_d
        if max_d > overall:
            overall = max_d
    return degs, overall


INSTANCES = [
    'hd08', 'hd07', 'hd10', 'hd03', 'hd04', 'ctrl', 'dec',
    'int2float', 'hd01', 'hd02', 'cavlc'
]

print(f"{'实例':>10} | {'opt1最高次':>10} | {'opt2最高次':>10} | {'原始最高次':>10}")
print("-" * 55)

results = {}
for inst in INSTANCES:
    d = os.path.join(BASE, inst)
    opt1_file = os.path.join(d, f"{inst}_opt1_expr.poly")
    opt2_file = os.path.join(d, f"{inst}_opt2_expr.poly")

    _, opt1_overall = parse_expr_degree(opt1_file) if os.path.exists(opt1_file) else ({}, 0)
    _, opt2_overall = parse_expr_degree(opt2_file) if os.path.exists(opt2_file) else ({}, 0)

    # raw degree = opt degree (preserved under affine transform)
    raw_deg = opt1_overall if opt1_overall > opt2_overall else opt2_overall
    if raw_deg == 0:
        raw_deg = opt1_overall or opt2_overall

    results[inst] = {'raw': raw_deg, 'opt1': opt1_overall, 'opt2': opt2_overall}

    print(f"{inst:>10} | {str(opt1_overall):>10} | {str(opt2_overall):>10} | {str(raw_deg):>10}")

print()
print("Copy-paste for result.md:")
print("-" * 55)
for inst in INSTANCES:
    r = results[inst]
    print(f"{inst}: 原始最高次={r['raw']}, opt1最高次={r['opt1']}, opt2最高次={r['opt2']}")
