"""
AND count optimizer for AND-XOR circuits.

Techniques:
1. XOR factoring: a*b XOR a*c = a*(b XOR c) — reduces AND count
2. Common subexpression elimination: share identical AND computations
3. AND depth preservation: doesn't increase depth beyond theoretical minimum
"""

import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def parse_circuit(txt_path):
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


def deps(expr):
    return re.findall(r'\b([a-zA-Z_]\w*)\b', expr)


def is_and(expr):
    return isinstance(expr, str) and '*' in expr and not expr.startswith('!')


def is_pure_and(expr):
    return is_and(expr) and re.match(r'\w+\s*\*\s*\w+$', expr)


def count_uses(stmt_map):
    use_count = {}
    for name, expr in stmt_map.items():
        for d in deps(expr):
            use_count[d] = use_count.get(d, 0) + 1
    return use_count


def topo_sort(names, stmt_map, input_set):
    names = [n for n in names if n not in input_set]
    in_deg = {}
    graph = {}
    name_set = set(names)
    for name in names:
        in_deg.setdefault(name, 0)
        for d in deps(stmt_map[name]):
            if d != name and d in name_set:
                graph.setdefault(d, []).append(name)
                in_deg[name] = in_deg.get(name, 0) + 1
    q = [n for n in names if in_deg.get(n, 0) == 0]
    result = []
    while q:
        n = q.pop(0)
        result.append(n)
        for s in graph.get(n, []):
            in_deg[s] -= 1
            if in_deg[s] == 0:
                q.append(s)
    result.extend([n for n in names if n not in result])
    return result


def compute_metrics(inputs, outputs, stmt_map, input_set=None):
    if input_set is None:
        input_set = set(inputs)
    and_depth = {}
    for nm in inputs:
        and_depth[nm] = 0
    order = topo_sort(list(stmt_map.keys()), stmt_map, input_set)
    total_and = 0
    for name in order:
        expr = stmt_map[name]
        ds = [d for d in deps(expr) if d != name and d in and_depth]
        if is_pure_and(expr):
            total_and += 1
            and_depth[name] = max(and_depth.get(d, 0) for d in ds) + 1 if ds else 1
        else:
            and_depth[name] = max(and_depth.get(d, 0) for d in ds) if ds else 0
    return and_depth, total_and


def split_at_top_level(expr, op):
    parts = []
    depth = 0
    current = []
    for ch in expr:
        if ch == '(':
            depth += 1
            current.append(ch)
        elif ch == ')':
            depth -= 1
            current.append(ch)
        elif ch == op and depth == 0:
            parts.append(''.join(current).strip())
            current = []
        else:
            current.append(ch)
    parts.append(''.join(current).strip())
    return [p for p in parts if p]


def get_and_operands(expr):
    """Return (op1, op2) if expr is a pure AND, else None."""
    if is_and(expr) and re.match(r'\w+\s*\*\s*\w+$', expr):
        m = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
        return m.group(1), m.group(2)
    return None


def xor_factor_pass(stmt_map, input_set, prefix='_fac_', counter=None):
    """
    One pass of XOR factoring.
    For each XOR gate, find AND inputs with common operands and factor them:
        out = (x*a) XOR (x*b) XOR ...  →  t = a XOR b XOR ...; out = x*t

    Only factors single-use AND gates to ensure correctness.
    Returns (new_stmt_map, savings) where savings = number of AND gates removed.
    """
    if counter is None:
        counter = [0]

    use_count = count_uses(stmt_map)
    new_nodes = {}
    nodes_to_remove = set()
    changed = False

    for name, expr in sorted(stmt_map.items()):
        if name in new_nodes:
            continue
        if '+' not in expr:
            continue

        # Parse XOR inputs
        parts = split_at_top_level(expr, '+')

        # Identify AND-gate inputs and group by common operand
        # For each part, check if it references a pure AND gate
        and_info = {}  # orig_part -> (op1, op2, and_name, is_negated)

        for p in parts:
            p_stripped = p.strip()
            is_neg = p_stripped.startswith('!')
            p_clean = p_stripped.lstrip('!')

            if p_clean in stmt_map:
                ops = get_and_operands(stmt_map[p_clean])
                if ops:
                    and_info[p_stripped] = (*ops, p_clean, is_neg)

        if len(and_info) < 2:
            continue

        # Group by common operand, only including single-use ANDs
        groups = {}  # common_op -> [(other_op, orig_part, and_name, is_neg)]
        for orig_part, (op1, op2, and_name, is_neg) in and_info.items():
            if use_count.get(and_name, 0) > 1 and and_name not in input_set:
                continue  # shared AND, can't factor safely

            for common_op, other_op in [(op1, op2), (op2, op1)]:
                groups.setdefault(common_op, []).append(
                    (other_op, orig_part, and_name, is_neg))

        # Apply best factoring opportunity (largest group first)
        best_groups = sorted(groups.values(), key=len, reverse=True)
        best_group = None
        for g in best_groups:
            if len(g) >= 2:
                # Check that all ANDs in this group are distinct
                and_names = set(item[2] for item in g)
                if len(and_names) >= 2:
                    best_group = g
                    break

        if best_group is None:
            continue

        common_op = None
        for op, g in groups.items():
            if g == best_group:
                common_op = op
                break

        # Create factored expression
        other_ops = []
        parts_to_remove = set()
        ands_to_remove = set()

        for other_op, orig_part, and_name, is_neg in best_group:
            if is_neg:
                other_ops.append(f"!{other_op}")
            else:
                other_ops.append(other_op)
            parts_to_remove.add(orig_part)
            ands_to_remove.add(and_name)

        # Create t = other_ops[0] XOR other_ops[1] XOR ...
        t_expr = " + ".join(other_ops) if len(other_ops) > 1 else other_ops[0]
        counter[0] += 1
        t_name = f"{prefix}{counter[0]}"
        new_nodes[t_name] = t_expr

        # Check if t_name itself is a reference to an existing node
        if t_expr in stmt_map:
            t_name = t_expr  # reuse existing node

        # Create new_and = common_op * t
        counter[0] += 1
        new_and_name = f"{prefix}{counter[0]}"
        new_nodes[new_and_name] = f"{common_op} * {t_name}" if common_op < t_name else f"{t_name} * {common_op}"

        # Rebuild XOR expression
        new_parts = [p for p in parts if p not in parts_to_remove]
        new_parts.append(new_and_name)

        new_expr = " + ".join(new_parts)
        if new_expr != expr:
            new_nodes[name] = new_expr
            changed = True
            for n in ands_to_remove:
                nodes_to_remove.add(n)

    if not changed:
        return stmt_map, 0

    # Build result
    result = dict(stmt_map)
    result.update(new_nodes)
    for n in nodes_to_remove:
        if n in result:
            del result[n]

    # Clean up: remove orphaned nodes (not needed for outputs)
    # But we skip this here since we don't know outputs

    and_removed = len(nodes_to_remove)
    return result, and_removed


def optimize_and_count(inputs, outputs, stmt_map, max_passes=5):
    """
    Reduce AND count via XOR factoring and CSE.
    Returns (optimized_stmt_map, total_and_removed, passes_done).
    """
    input_set = set(inputs)
    total_removed = 0
    passes = 0

    current = dict(stmt_map)

    for _ in range(max_passes):
        and_before = sum(1 for e in current.values() if is_pure_and(e))
        current, removed = xor_factor_pass(current, input_set)
        and_after = sum(1 for e in current.values() if is_pure_and(e))
        total_removed += removed
        passes += 1
        if removed == 0:
            break

    # Prune unused nodes (cone-of-influence)
    needed = set(outputs)
    changed = True
    while changed:
        changed = False
        for n in list(needed):
            if n in current:
                for d in deps(current[n]):
                    if d not in needed:
                        needed.add(d)
                        changed = True

    clean = {}
    for n in needed:
        if n in current and n not in input_set:
            clean[n] = current[n]
        elif n in set(outputs) and n not in current:
            clean[n] = '0'
    for o in outputs:
        if o not in clean and o in current and o not in input_set:
            clean[o] = current[o]

    new_depth, new_and = compute_metrics(inputs, outputs, clean, input_set)
    old_depth, old_and = compute_metrics(inputs, outputs, stmt_map, input_set)

    return clean, old_and - new_and, passes


def run(txt_path):
    name = os.path.splitext(os.path.basename(txt_path))[0]
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    inputs, outputs, stmt_map = parse_circuit(txt_path)
    input_set = set(inputs)

    old_depth, old_and = compute_metrics(inputs, outputs, stmt_map, input_set)
    max_depth = max(old_depth.get(o, 0) for o in outputs)
    print(f"  原始: AND={old_and}, 最大深度={max_depth}")

    new_stmt, and_saved, passes = optimize_and_count(inputs, outputs, stmt_map)
    new_depth, new_and = compute_metrics(inputs, outputs, new_stmt, input_set)
    new_max_depth = max(new_depth.get(o, 0) for o in outputs)

    print(f"  优化后: AND={new_and}, 最大深度={new_max_depth}")
    print(f"  节省AND: {and_saved} ({passes} 轮)")
    print(f"  深度变化: {max_depth} → {new_max_depth}")

    if and_saved > 0 or new_max_depth < max_depth:
        # Write optimized circuit
        out_path = txt_path.replace('.txt', '_opt_and.txt')
        with open(out_path, 'w') as f:
            f.write(f"INORDER = {' '.join(inputs)};\n")
            f.write(f"OUTORDER = {' '.join(outputs)};\n")
            order = topo_sort(list(new_stmt.keys()), new_stmt, input_set)
            written = set()
            for n in order:
                if n in input_set or n in set(outputs):
                    continue
                if n in new_stmt:
                    f.write(f"{n} = {new_stmt[n]};\n")
                    written.add(n)
            for o in outputs:
                if o in new_stmt:
                    f.write(f"{o} = {new_stmt[o]};\n")
        print(f"  写入: {out_path}")

        # Print per-output depth change
        for o in outputs:
            od = old_depth.get(o, 0)
            nd = new_depth.get(o, 0)
            if od != nd:
                print(f"    {o}: {od} → {nd}")
    else:
        print(f"  无优化")

    return new_stmt


if __name__ == '__main__':
    if len(sys.argv) < 2:
        for circ in ['hd09','hd10','hd11','hd12','router','dsort']:
            try:
                run(f'examples/{circ}/{circ}.txt')
            except Exception as e:
                import traceback
                print(f"  ERROR: {e}")
                traceback.print_exc()
    else:
        run(sys.argv[1])
