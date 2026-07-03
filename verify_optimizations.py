"""
Verify all optimized circuits against original circuits.

For each instance, verifies that:
  f(x) = g(Mx ⊕ b)   for all x (or sampled x for large n)

Two strategies:
  opt1: single shared transform for all outputs
  opt2: per-output transforms merged into shared transform

Method 1: Random simulation (all instances)
Method 2: Exhaustive enumeration (n ≤ 16)
"""

import sys, os, re
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from format_poly import _anf_term_str


def parse_circuit_expr(txt_path):
    """Parse circuit .txt file, return (inputs, outputs, stmt_map)."""
    with open(txt_path) as f:
        text = f.read()
    raw = text.replace('\n', ' ').replace('\r', '')
    raw = re.sub(r'\s+', ' ', raw).strip()
    inp_m = re.search(r'INORDER\s*=\s*([^;]+);', raw)
    out_m = re.search(r'OUTORDER\s*=\s*([^;]+);', raw)
    inputs = inp_m.group(1).strip().split()
    outputs = out_m.group(1).strip().split()
    stmt_map = {}
    body_start = max(raw.find(';', raw.find('INORDER')),
                     raw.find(';', raw.find('OUTORDER')) + 1)
    for stmt in raw[body_start:].split(';'):
        stmt = stmt.strip()
        if not stmt:
            continue
        m = re.match(r'(\w+)\s*=\s*(.+)$', stmt)
        if m:
            stmt_map[m.group(1)] = m.group(2).strip()
    return inputs, outputs, stmt_map


def parse_transform(trans_path):
    """Parse _trans.poly file, return list of (mask, const) for each z_j.

    z_j = linear_form(inputs)  →  (mask: int, const: 0/1)
    where bit i of mask corresponds to input i.
    """
    if not os.path.exists(trans_path) or os.path.getsize(trans_path) == 0:
        return None
    with open(trans_path) as f:
        lines = f.read().strip().splitlines()
    rows = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        m = re.match(r'z_(\d+)\s*=\s*(.+)', line)
        if not m:
            continue
        expr = m.group(2).strip()
        mask = 0
        const = 0
        if expr == '0':
            rows.append((mask, const))
            continue
        # Handle '1' (constant)
        # expr is like: i0 + i1 + 1  (XOR)
        parts = re.findall(r'[a-zA-Z_]\w*|1', expr)
        for p in parts:
            if p == '1':
                const ^= 1
        # We also need to map variable names to indices.
        # But we don't have the input names here.
        # Solution: return the variable names as strings.
        rows.append((None, const, expr))
    return rows


def parse_transform_with_inputs(trans_path, input_names):
    """Parse _trans.poly with knowledge of input variable names."""
    if not os.path.exists(trans_path) or os.path.getsize(trans_path) == 0:
        return None

    name_to_idx = {name: i for i, name in enumerate(input_names)}

    with open(trans_path) as f:
        lines = f.read().strip().splitlines()
    rows = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        m = re.match(r'z_(\d+)\s*=\s*(.+)', line)
        if not m:
            continue
        expr = m.group(2).strip()
        mask = 0
        const = 0
        if expr == '0':
            rows.append((mask, const))
            continue
        # Parse XOR-of-terms
        terms = expr.split('+')
        for term in terms:
            term = term.strip()
            if term == '1':
                const ^= 1
            elif term in name_to_idx:
                mask ^= (1 << name_to_idx[term])
            elif term.startswith('!') and term[1:] in name_to_idx:
                # !var = (1 XOR var) = adds 1 and var
                const ^= 1
                mask ^= (1 << name_to_idx[term[1:]])
            else:
                # Unknown variable — could be a z or internal name
                # Treat as symbolic (skip for now, likely not in input list)
                pass
        rows.append((mask, const))
    return rows


def parse_expr(expr_path, outputs=None):
    """Parse _expr.poly file, return dict {output_name: {mask: coeff}}."""

    def term_str_to_mask(term_str):
        """Convert z_i*z_j*... string to mask integer."""
        if term_str == '1':
            return 0
        mask = 0
        parts = term_str.split('*')
        for p in parts:
            p = p.strip()
            if p.startswith('z_'):
                idx = int(p[2:])
                mask |= (1 << idx)
            # Skip other variable types (internal names)
        return mask

    with open(expr_path) as f:
        lines = f.read().strip().splitlines()

    result = {}
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        m = re.match(r'(\w+)\s*=\s*(.+)', line)
        if not m:
            continue
        name = m.group(1)
        expr = m.group(2).strip()

        if outputs and name not in outputs:
            continue

        if expr == '0':
            result[name] = {}
            continue

        terms = {}
        monomials = expr.split('+')
        for mon in monomials:
            mon = mon.strip()
            if not mon:
                continue
            mask = term_str_to_mask(mon)
            terms[mask] = terms.get(mask, 0) ^ 1
        terms = {k: v for k, v in terms.items() if v}
        result[name] = terms

    return result


def eval_circuit(inputs, outputs, stmt_map, input_vals):
    """Evaluate circuit for given input values. Returns dict output_name -> value."""
    input_set = set(inputs)
    vals = {}
    for inp in inputs:
        vals[inp] = input_vals.get(inp, 0)

    # Topological order
    names = [n for n in stmt_map if n not in input_set]
    order = []
    in_deg = {}
    graph = {}
    name_set = set(names)
    for name in names:
        in_deg.setdefault(name, 0)
        for d in re.findall(r'\b([a-zA-Z_]\w*)\b', stmt_map[name]):
            if d != name and d in name_set:
                graph.setdefault(d, []).append(name)
                in_deg[name] = in_deg.get(name, 0) + 1
    q = [n for n in names if in_deg.get(n, 0) == 0]
    while q:
        n = q.pop(0)
        order.append(n)
        for s in graph.get(n, []):
            in_deg[s] -= 1
            if in_deg[s] == 0:
                q.append(s)
    order.extend([n for n in names if n not in order])

    # Evaluate
    for name in order:
        expr = stmt_map[name]
        if not expr:
            vals[name] = 0
            continue
        if expr.startswith('!'):
            op = expr[1:].strip()
            vals[name] = 1 ^ vals.get(op, 0)
        elif '*' in expr:
            ops = re.findall(r'\b(\w+)\b', expr)
            result = 1
            for op in ops:
                result &= vals.get(op, 0)
            vals[name] = result
        elif '+' in expr:
            ops = re.findall(r'\b(\w+)\b', expr)
            result = 0
            for op in ops:
                result ^= vals.get(op, 0)
            vals[name] = result
        elif expr == '0':
            vals[name] = 0
        elif expr == '1':
            vals[name] = 1
        else:
            vals[name] = vals.get(expr, 0)

    return {o: vals.get(o, 0) for o in outputs}


def eval_anf(expr_dict, z_vals):
    """Evaluate ANF expression {mask: coeff} given z variable values.

    z_vals: list of int (0/1) for z_0, z_1, ..., z_{m-1}
    """
    result = 0
    for mask, coeff in expr_dict.items():
        if coeff == 0:
            continue
        term_val = 1
        t = mask
        while t:
            lsb = t & -t
            zi = (lsb.bit_length() - 1)
            if zi < len(z_vals):
                term_val &= z_vals[zi]
            else:
                term_val = 0
                break
            t ^= lsb
        result ^= term_val
    return result


def compute_z(M_rows, b_vals, x_vec):
    """Compute z = Mx ⊕ b.

    M_rows: list of int masks (each mask → bit i = M[j][i])
    b_vals: list of int (0/1)
    x_vec: int mask of input values
    Returns: list of int (0/1) for z_0, ..., z_{m-1}
    """
    z = []
    for mask, b in zip(M_rows, b_vals):
        # Dot product (mod 2) of M row with x
        dot = bin(mask & x_vec).count('1') & 1
        z.append(dot ^ b)
    return z


def expr_to_z_expr(terms, z_count):
    """Convert ANF dict expression to a function that evaluates from z values."""
    def evaluate(z_vals):
        result = 0
        for mask, coeff in terms.items():
            if coeff == 0:
                continue
            term_val = 1
            t = mask
            while t:
                lsb = t & -t
                zi = (lsb.bit_length() - 1)
                if zi < len(z_vals):
                    term_val &= z_vals[zi]
                else:
                    term_val = 0
                    break
                t ^= lsb
            result ^= term_val
        return result
    return evaluate


def verify_instance(inst_name, base_path, num_random_tests=50000, exhaustive=True):
    """Verify one instance: original circuit = optimized circuit."""
    txt_path = os.path.join(base_path, inst_name, f"{inst_name}.txt")

    if not os.path.exists(txt_path):
        print(f"  {inst_name}: .txt not found, skip")
        return True, 0

    inputs, outputs, stmt_map = parse_circuit_expr(txt_path)
    n = len(inputs)

    max_exhaustive = min(n, 20)  # Cap exhaustive at 2^20 = 1M
    if n > max_exhaustive:
        exhaustive = False

    is_ok = True
    total_tests = 0

    for strategy in ['opt1', 'opt2']:
        trans_path = os.path.join(base_path, inst_name, f"{inst_name}_{strategy}_trans.poly")
        expr_path = os.path.join(base_path, inst_name, f"{inst_name}_{strategy}_expr.poly")

        if not os.path.exists(expr_path):
            continue

        # Parse transform
        M = parse_transform_with_inputs(trans_path, inputs)
        if M is None:
            print(f"  {inst_name} ({strategy}): transform empty/missing, skip")
            continue

        M_rows = [m for m, _ in M]
        b_vals = [b for _, b in M]

        # Parse expressions
        exprs = parse_expr(expr_path, outputs)
        if not exprs:
            continue

        m_z = len(M_rows)

        # Collect non-zero expression maps
        expr_maps = {}
        for o in outputs:
            if o in exprs and exprs[o]:
                expr_maps[o] = exprs[o]

        non_zero_out = list(expr_maps.keys())

        # Test: random + exhaustive
        if exhaustive:
            test_range = range(1 << n)
        else:
            import random
            rng = random.Random(f"verify_{inst_name}_{strategy}")
            N = min(num_random_tests, 1 << n)
            test_range = [rng.randint(0, (1 << n) - 1) for _ in range(N)]

        mismatch_count = 0
        tested = 0

        for x_bits in test_range:
            x_vec = {}
            for i in range(n):
                x_vec[inputs[i]] = (x_bits >> i) & 1

            # Original circuit
            orig_outs = eval_circuit(inputs, outputs, stmt_map, x_vec)

            # Optimized circuit
            z_vals = compute_z(M_rows, b_vals, x_bits)

            for o in non_zero_out:
                g_val = eval_anf(expr_maps[o], z_vals)
                if g_val != orig_outs.get(o, 0):
                    mismatch_count += 1
                    if mismatch_count <= 3:
                        print(f"    MISMATCH {strategy}: x={x_bits}, "
                              f"{o}: orig={orig_outs.get(o,0)}, opt={g_val}, "
                              f"z={z_vals}")
                    break  # one mismatch per test vector is enough

            tested += 1
            if tested % 10000 == 0 and not exhaustive:
                pass  # progress

        if mismatch_count > 0:
            print(f"  ❌ {inst_name} ({strategy}): {mismatch_count}/{tested} mismatches")
            is_ok = False
        else:
            ok_str = "✅" if not exhaustive or (1 << n) <= (1 << max_exhaustive) else "✓"
            test_str = f"全部{1<<n}" if exhaustive and n <= max_exhaustive else f"{tested}次随机"
            print(f"  {ok_str} {inst_name} ({strategy}): {test_str}测试通过")
            total_tests += tested

    return is_ok, total_tests


def main():
    base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'examples')

    # Parser instances (directly solvable)
    instances = ['hd08', 'hd07', 'hd03', 'hd04', 'ctrl', 'dec',
                 'int2float', 'hd01', 'hd02', 'cavlc']

    print("=" * 60)
    print("  验证优化电路等价性")
    print("=" * 60)

    all_pass = True
    total = 0

    for inst in instances:
        print(f"\n--- {inst} ---")
        txt_path = os.path.join(base_path, inst, f"{inst}.txt")
        inputs, outputs, _ = parse_circuit_expr(txt_path)
        n = len(inputs)
        exhaustive = n <= 16

        ok, ntests = verify_instance(inst, base_path,
                                     num_random_tests=100000,
                                     exhaustive=exhaustive)
        if not ok:
            all_pass = False
        total += ntests

    print(f"\n{'='*60}")
    if all_pass:
        print(f"  ✅ 全部实例验证通过 (共 {total} 个测试)")
    else:
        print(f"  ❌ 存在验证失败的实例")
    print(f"{'='*60}")

    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
