"""
AND depth optimizer v2: aggressive rebalancing with selective duplication.

Strategy:
1. Find outputs with AND depth significantly above theoretical minimum
2. For each deep output, identify the AND chain on the critical path
3. For shared AND nodes on the chain, create duplicates (trades AND count for depth)
4. Collect all leaf factors of the AND chain (including through shared nodes)
5. Rebuild the AND tree as balanced binary tree
6. Remap the output to use the new balanced tree

Only applies when the depth improvement justifies the AND increase.
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
    and_cnt = {}
    for nm in inputs:
        and_depth[nm] = 0
        and_cnt[nm] = 0
    order = topo_sort(list(stmt_map.keys()), stmt_map, input_set)
    total_and = 0
    for name in order:
        expr = stmt_map[name]
        ds = [d for d in deps(expr) if d != name and d in and_depth]
        if is_pure_and(expr):
            total_and += 1
            depth = max(and_depth.get(d, 0) for d in ds) + 1 if ds else 1
            cnt = max(and_cnt.get(d, 0) for d in ds) + 1 if ds else 1
        else:
            depth = max(and_depth.get(d, 0) for d in ds) if ds else 0
            cnt = max(and_cnt.get(d, 0) for d in ds) if ds else 0
        and_depth[name] = depth
        and_cnt[name] = cnt
    return and_depth, and_cnt, total_and


def trace_critical_path(output, and_depth, stmt_map, input_set):
    path = []
    name = output
    while name not in input_set:
        if name not in stmt_map:
            break
        expr = stmt_map[name]
        ds = [d for d in deps(expr) if d != name and d in and_depth]
        if not ds:
            break
        next_name = max(ds, key=lambda d: and_depth.get(d, 0))
        path.append((name, and_depth.get(name, 0)))
        name = next_name
    return path


def get_and_operands(expr):
    if is_pure_and(expr):
        m = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
        if m:
            return m.group(1), m.group(2)
    return None


def collect_all_and_factors(name, stmt_map, input_set, and_depth, visited=None):
    """Recursively expand AND nodes to collect all non-AND leaf factors.
    Transparent through NOT gates. Inlines through ALL ANDs (including shared)."""
    if visited is None:
        visited = set()

    name = name.strip()
    if name in input_set or name in visited:
        return {name}
    visited.add(name)

    expr = stmt_map.get(name, '')
    if not expr:
        return {name}

    # NOT gate: transparent, recurse into operand
    if expr.startswith('!'):
        op = expr[1:].strip()
        return collect_all_and_factors(op, stmt_map, input_set, and_depth, visited)

    # Pure AND gate: expand both operands
    ops = get_and_operands(expr)
    if ops:
        result = set()
        for op in ops:
            result.update(collect_all_and_factors(op, stmt_map, input_set, and_depth, visited))
        return result

    # XOR or other: can't expand through, it's a leaf
    return {name}


def collect_all_and_factors_safe(name, stmt_map, input_set, stop_set, visited=None):
    """Expand through ANDs, but don't re-enter nodes in stop_set (chain nodes).
    Transparent through NOT gates."""
    if visited is None:
        visited = set()

    name = name.strip()
    if name in input_set or name in visited:
        return {name}
    if name in stop_set:
        return {name}  # stop at chain nodes
    visited.add(name)

    expr = stmt_map.get(name, '')
    if not expr:
        return {name}

    if expr.startswith('!'):
        op = expr[1:].strip()
        return collect_all_and_factors_safe(op, stmt_map, input_set, stop_set, visited)

    ops = get_and_operands(expr)
    if ops:
        result = set()
        for op in ops:
            result.update(collect_all_and_factors_safe(op, stmt_map, input_set, stop_set, visited))
        return result

    return {name}


def build_balanced(factors, prefix, counter):
    """Build balanced AND tree from factors. Returns (nodes, root_name)."""
    nodes = {}
    queue = list(factors)
    while len(queue) > 1:
        new_q = []
        for i in range(0, len(queue), 2):
            if i + 1 < len(queue):
                a, b = queue[i], queue[i + 1]
                if a == b:
                    new_q.append(a)
                else:
                    counter[0] += 1
                    n = f"{prefix}{counter[0]}"
                    nodes[n] = f"{a} * {b}" if a < b else f"{b} * {a}"
                    new_q.append(n)
            else:
                new_q.append(queue[i])
        queue = new_q
    return nodes, queue[0]


def optimize_depth_aggressive(inputs, outputs, stmt_map, min_chain=3, max_and_blowup=2.0):
    """
    Aggressive AND depth optimization.
    Duplicates shared nodes on critical paths and rebalances AND chains.

    max_and_blowup: maximum allowed AND count increase as a ratio (2.0 = 2x).
    """
    input_set = set(inputs)
    output_set = set(outputs)
    and_depth, and_cnt, total_and = compute_metrics(inputs, outputs, stmt_map, input_set)
    use_count = count_uses(stmt_map)

    # Find outputs with depth > theoretical minimum
    n_inputs = len(inputs)
    theoretical_min = (n_inputs - 1).bit_length()

    output_depths = [(o, and_depth.get(o, 0)) for o in outputs]
    output_depths.sort(key=lambda x: -x[1])

    new_nodes = {}
    remap = {}
    counter = [0]
    and_added = 0

    for out_name, depth in output_depths:
        if depth < min_chain:
            continue
        if depth <= theoretical_min + 1:
            continue  # already close to minimum

        # Check budget
        if total_and + and_added > total_and * max_and_blowup:
            continue

        # Trace critical path
        path = trace_critical_path(out_name, and_depth, stmt_map, input_set)
        # path is [(name, depth), ...] deepest first

        # Collect ALL AND nodes on the path
        all_chain_ands = []
        for name, d in path:
            if name in input_set:
                break
            expr = stmt_map.get(name, '')
            if is_pure_and(expr):
                all_chain_ands.append(name)
        # Reverse to deepest-first
        all_chain_ands.reverse()

        if len(all_chain_ands) < min_chain:
            continue

        chain_set = set(all_chain_ands)
        path_names = set(n for n, _ in path)

        # For each AND in the chain, find the side leaf (operand NOT on the
        # critical path, tracing through NOT gates).
        # The predecessor is the operand that connects to the next chain node
        # (possibly through NOT gates).
        def find_predecessor(and_name, stmt_map, path_names):
            """Find which operand of an AND gate is on the critical path
            (tracing through NOT gates). Returns the predecessor name or None."""
            expr = stmt_map[and_name]
            ops = get_and_operands(expr)
            if not ops:
                return None
            lhs, rhs = ops
            for op in [lhs, rhs]:
                cur = op
                # Walk through NOT gates
                while cur in stmt_map and stmt_map[cur].startswith('!'):
                    cur = stmt_map[cur][1:].strip()
                if cur in chain_set or cur in path_names:
                    return op
            return None

        all_leaves = set()

        for and_name in all_chain_ands:
            expr = stmt_map[and_name]
            ops = get_and_operands(expr)
            if not ops:
                continue
            lhs, rhs = ops
            pred = find_predecessor(and_name, stmt_map, path_names)
            for op in [lhs, rhs]:
                if op == pred:
                    continue  # chain predecessor
                # Side leaf — expand through ANDs to get true leaves,
                # but stop at chain nodes
                op_leaves = collect_all_and_factors_safe(
                    op, stmt_map, input_set, chain_set)
                all_leaves.update(op_leaves)

        if len(all_leaves) < min_chain:
            continue

        # Compute new depth
        max_leaf_depth = max(and_depth.get(l, 0) for l in all_leaves)
        tree_height = (len(all_leaves) - 1).bit_length()
        new_depth_est = max_leaf_depth + tree_height

        if new_depth_est >= depth:
            continue

        # Build balanced tree
        balanced, root = build_balanced(all_leaves, '_d2_', counter)
        if not balanced:
            continue

        # AND change
        new_and_count = len(balanced)
        old_and_count = len(all_chain_ands)
        and_delta = new_and_count - old_and_count

        depth_improvement = depth - new_depth_est
        if depth_improvement <= 0:
            continue

        # Trade-off: don't more than double ANDs per depth level
        if and_delta > 0 and depth_improvement * 10 < and_delta:
            continue

        new_nodes.update(balanced)
        # Remap: replace the output-side AND (last in deepest-first list)
        shallowest_and = all_chain_ands[0]  # first in deepest-first = closest to output
        # Wait: all_chain_ands is deepest-first, so last element is closest to output
        shallowest_and = all_chain_ands[-1]
        remap[shallowest_and] = root
        and_added += and_delta

    if not remap:
        return stmt_map, 0, 0

    # Build new statement map with remapping
    new_stmt = dict(stmt_map)
    new_stmt.update(new_nodes)

    # Remap references to the deepest node in other expressions
    for name in list(new_stmt.keys()):
        if name in input_set:
            continue
        expr = new_stmt[name]
        if isinstance(expr, str):
            for old, new_ in remap.items():
                if old in deps(expr):
                    expr = re.sub(r'\b' + re.escape(old) + r'\b', new_, expr)
            new_stmt[name] = expr

    # Cone-of-influence pruning
    needed = set(outputs)
    changed = True
    while changed:
        changed = False
        for n in list(needed):
            if n in new_stmt:
                for d in deps(new_stmt[n]):
                    if d not in needed:
                        needed.add(d)
                        changed = True

    clean = {}
    for n in needed:
        if n in new_stmt and n not in input_set:
            clean[n] = new_stmt[n]
        elif n in output_set and n not in new_stmt:
            clean[n] = '0'
    for n in outputs:
        if n not in clean and n in new_stmt and n not in input_set:
            clean[n] = new_stmt[n]

    m0 = compute_metrics(inputs, outputs, stmt_map, input_set)
    m1 = compute_metrics(inputs, outputs, clean, input_set)
    orig_max = max(m0[0].get(o, 0) for o in outputs)
    new_max = max(m1[0].get(o, 0) for o in outputs)
    depth_improv = orig_max - new_max
    and_improv = m0[2] - m1[2]  # negative means AND count increased

    return clean, depth_improv, and_improv


def run(txt_path):
    name = os.path.splitext(os.path.basename(txt_path))[0]
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    inputs, outputs, stmt_map = parse_circuit(txt_path)
    input_set = set(inputs)
    m0 = compute_metrics(inputs, outputs, stmt_map, input_set)
    orig_max_depth = max(m0[0].get(o, 0) for o in outputs)
    print(f"  原始: AND={m0[2]}, 最大AND深度={orig_max_depth}")

    new_stmt, depth_improv, and_improv = optimize_depth_aggressive(
        inputs, outputs, stmt_map, min_chain=3, max_and_blowup=3.0)

    if depth_improv > 0:
        m1 = compute_metrics(inputs, outputs, new_stmt, input_set)
        new_max_depth = max(m1[0].get(o, 0) for o in outputs)
        print(f"  优化后: AND={m1[2]}, 最大AND深度={new_max_depth}")
        print(f"  AND深度改善: {orig_max_depth} → {new_max_depth} (-{depth_improv})")
        print(f"  AND数: {m0[2]} → {m1[2]} ({m1[2]-m0[2]})")

        for o in outputs:
            od = m0[0].get(o, 0)
            nd = m1[0].get(o, 0)
            if od != nd:
                print(f"    {o}: {od} → {nd} {'↓' if nd < od else '↑'}")

        out_path = txt_path.replace('.txt', '_opt_d2.txt')
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
    else:
        print(f"  无有效优化")

    return new_stmt if depth_improv > 0 else stmt_map


if __name__ == '__main__':
    if len(sys.argv) < 2:
        for circ in ['hd09', 'hd10', 'hd11', 'hd12', 'router', 'dsort']:
            try:
                run(f'examples/{circ}/{circ}.txt')
            except Exception as e:
                import traceback
                print(f"  ERROR: {e}")
                traceback.print_exc()
    else:
        run(sys.argv[1])
