"""
Circuit optimizer: reduce AND depth by rebalancing AND chains on the critical path.

Finds the deepest AND path to each output, then rebalances sub-chains.
"""

import sys, os, re, time
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


def is_xor(expr):
    return isinstance(expr, str) and '+' in expr


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
    total_and = 0
    total_xor = 0
    for nm in inputs:
        and_depth[nm] = 0
        and_cnt[nm] = 0
    order = topo_sort(list(stmt_map.keys()), stmt_map, input_set)
    for name in order:
        expr = stmt_map[name]
        ds = [d for d in deps(expr) if d != name and d in and_depth]
        is_and_gate = is_and(expr) and re.match(r'\w+\s*\*\s*\w+$', expr)
        if is_and_gate:
            total_and += 1
            depth = max(and_depth.get(d, 0) for d in ds) + 1
            cnt = max(and_cnt.get(d, 0) for d in ds) + 1
        else:
            depth = max(and_depth.get(d, 0) for d in ds) if ds else 0
            cnt = max(and_cnt.get(d, 0) for d in ds) if ds else 0
            if is_xor(expr):
                total_xor += 1
        and_depth[name] = depth
        and_cnt[name] = cnt
    return and_depth, and_cnt, total_and, total_xor


def count_uses(stmt_map):
    """Count how many times each node is used."""
    use_count = {}
    for name, expr in stmt_map.items():
        for d in deps(expr):
            use_count[d] = use_count.get(d, 0) + 1
    return use_count


def trace_critical_path(output, and_depth, stmt_map, input_set):
    """Trace the AND-critical path from output to inputs.
    Returns list of (name, depth) on the path, in reverse order
    (deepest first)."""
    path = []
    name = output
    while name not in input_set:
        if name not in stmt_map:
            break
        expr = stmt_map[name]
        ds = [d for d in deps(expr) if d != name and d in and_depth]
        if not ds:
            break
        # Pick the dep with max depth
        next_name = max(ds, key=lambda d: and_depth.get(d, 0))
        path.append((name, and_depth.get(name, 0)))
        name = next_name
    return path


def collect_chain_leaves(start_name, stmt_map, and_depth, input_set):
    """Starting from a deep AND node, collect leaf factors.
    Only includes direct AND operands that are on the critical path."""
    leaves = []
    stack = [start_name]
    visited = set()
    while stack:
        name = stack.pop()
        if name in visited or name in input_set:
            leaves.append(name) if name not in leaves else None
            continue
        visited.add(name)
        expr = stmt_map.get(name)
        if not expr:
            leaves.append(name)
            continue
        if is_and(expr) and re.match(r'\w+\s*\*\s*\w+$', expr):
            ma = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
            lhs, rhs = ma.group(1), ma.group(2)
            # Add operands
            for op in [lhs, rhs]:
                op_expr = stmt_map.get(op, '')
                if is_and(op_expr) and re.match(r'\w+\s*\*\s*\w+$', op_expr):
                    stack.append(op)
                else:
                    if op not in leaves:
                        leaves.append(op)
        else:
            if name not in leaves:
                leaves.append(name)
    return leaves


def build_balanced(factors, prefix, counter):
    """Build balanced AND tree from factors. Returns (nodes, root_name)."""
    nodes = {}
    queue = list(factors)
    while len(queue) > 1:
        new_q = []
        for i in range(0, len(queue), 2):
            if i + 1 < len(queue):
                a, b = queue[i], queue[i+1]
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


def optimize_depth(inputs, outputs, stmt_map, min_chain=3):
    """Reduce AND depth by rebalancing chains on critical paths."""
    input_set = set(inputs)
    output_set = set(outputs)
    and_depth, and_cnt, _, _ = compute_metrics(inputs, outputs, stmt_map, input_set)
    use_count = count_uses(stmt_map)

    # Find deepest outputs
    output_depths = [(o, and_depth.get(o, 0)) for o in outputs]
    output_depths.sort(key=lambda x: -x[1])

    new_nodes = {}  # name → expr
    nodes_to_remove = set()
    remap = {}      # old_name → new_root_name
    counter = [0]

    for out_name, depth in output_depths:
        if depth < min_chain:
            continue

        # Trace critical path
        path = trace_critical_path(out_name, and_depth, stmt_map, input_set)
        # Path is [(node, depth), ...] from deepest to output

        # Find consecutive AND chain on the path
        # Look for subchains of at least `min_chain` consecutive ANDs
        chain_starts = []
        i = 0
        while i < len(path):
            name, d = path[i]
            # Check if this is an AND node
            expr = stmt_map.get(name, '')
            if is_and(expr) and re.match(r'\w+\s*\*\s*\w+$', expr):
                # Start collecting AND chain
                chain = [name]
                j = i + 1
                while j < len(path):
                    name2, d2 = path[j]
                    expr2 = stmt_map.get(name2, '')
                    if is_and(expr2) and re.match(r'\w+\s*\*\s*\w+$', expr2):
                        chain.append(name2)
                        j += 1
                    else:
                        break

                if len(chain) >= min_chain:
                    # Check if chain nodes are predominantly single-use
                    shared_in_chain = sum(1 for n in chain
                                          if use_count.get(n, 0) > 1
                                          and n not in output_set)
                    # If too many shared nodes, only optimize from first single-use
                    if shared_in_chain <= len(chain) // 2:
                        chain_starts.append(chain)
                    i = j
                else:
                    i = j
            else:
                i += 1

        if not chain_starts:
            continue

        # For each AND chain, collect leaves and rebalance
        for chain in chain_starts:
            # Chain is from deepest to shallowest: [deep_node, ..., shallow_node]
            # shallow_node is closer to output
            # We want to rebalance from the deepest node up to where shared nodes start
            # Find the deepest point where sharing starts
            split_idx = 0
            for idx, n in enumerate(chain):
                if use_count.get(n, 0) > 1 and n not in output_set:
                    # This node is shared, don't include in rebalancing
                    break
                split_idx = idx + 1

            if split_idx < min_chain:
                continue  # Not enough single-use nodes to rebalance

            # The nodes to rebalance: chain[:split_idx]
            rebal_nodes = chain[:split_idx]
            # The deepest node in the rebalancing group (furthest from output)
            deepest = rebal_nodes[-1]  # last in list, deepest

            # Collect leaves feeding into the deepest node
            leaves = collect_chain_leaves(deepest, stmt_map, and_depth, input_set)
            if len(leaves) < min_chain:
                continue

            # Build balanced tree
            balanced, root = build_balanced(leaves, '_opt_', counter)
            if not balanced:
                continue

            # Map the deepest node to the new root
            new_nodes.update(balanced)
            remap[deepest] = root

    if not remap:
        return stmt_map, 0, 0

    # Build new statement map
    new_stmt = dict(stmt_map)
    new_stmt.update(new_nodes)

    # Explicitly add new nodes for outputs that use remapped chains
    for name in list(new_stmt.keys()):
        if name in input_set:
            continue
        expr = new_stmt[name]
        if isinstance(expr, str):
            for old, new_ in remap.items():
                if old in deps(expr):
                    expr = re.sub(r'\b' + re.escape(old) + r'\b', new_, expr)
            new_stmt[name] = expr

    # Cone-of-influence: keep only what's needed
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
        if n in clean:
            pass  # already added
        elif n in new_stmt and n not in input_set:
            clean[n] = new_stmt[n]

    m0 = compute_metrics(inputs, outputs, stmt_map, input_set)
    m1 = compute_metrics(inputs, outputs, clean, input_set)
    orig_max = max(m0[0].get(o, 0) for o in outputs)
    new_max = max(m1[0].get(o, 0) for o in outputs)
    depth_improv = orig_max - new_max
    and_improv = m0[2] - m1[2]

    return clean, depth_improv, and_improv


def run(txt_path):
    name = os.path.splitext(os.path.basename(txt_path))[0]
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    inputs, outputs, stmt_map = parse_circuit(txt_path)
    input_set = set(inputs)
    output_set = set(outputs)
    m0 = compute_metrics(inputs, outputs, stmt_map, input_set)
    orig_max_depth = max(m0[0].get(o, 0) for o in outputs)
    print(f"  原始: AND={m0[2]}, XOR={m0[3]}, 最大AND深度={orig_max_depth}")

    new_stmt, depth_improv, and_improv = optimize_depth(
        inputs, outputs, stmt_map, min_chain=3)

    if depth_improv > 0 or and_improv > 0:
        m1 = compute_metrics(inputs, outputs, new_stmt, input_set)
        new_max_depth = max(m1[0].get(o, 0) for o in outputs)
        print(f"  优化后: AND={m1[2]}, 最大AND深度={new_max_depth}")
        print(f"  AND深度改善: {orig_max_depth} → {new_max_depth} (-{depth_improv})")
        print(f"  AND数: {m0[2]} → {m1[2]} ({-and_improv})")

        for o in outputs:
            od = m0[0].get(o, 0)
            nd = m1[0].get(o, 0)
            if od != nd:
                arrow = f" ↓{od-nd}" if nd < od else ""
                print(f"    {o}: {od} → {nd}{arrow}")

        # Write optimized circuit
        out_path = txt_path.replace('.txt', '_opt_dp.txt')
        with open(out_path, 'w') as f:
            f.write(f"INORDER = {' '.join(inputs)};\n")
            f.write(f"OUTORDER = {' '.join(outputs)};\n")
            # Write internal nodes in topo order
            order = topo_sort(list(new_stmt.keys()), new_stmt, input_set)
            written = set()
            for n in order:
                if n in input_set or n in output_set:
                    continue
                if n in new_stmt and isinstance(new_stmt[n], str):
                    f.write(f"{n} = {new_stmt[n]};\n")
                    written.add(n)
            for n in outputs:
                if n in new_stmt:
                    expr = new_stmt[n]
                    if isinstance(expr, str):
                        f.write(f"{n} = {expr};\n")
                    else:
                        f.write(f"{n} = 0;\n")
        print(f"  写入: {out_path}")
    else:
        print(f"  无有效优化 (AND树已较平衡或共享度过高)")

    return new_stmt if depth_improv > 0 else stmt_map


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
