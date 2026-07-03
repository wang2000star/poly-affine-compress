"""
AND depth optimizer v3: aggressive expansion + rebalancing.

Strategy:
1. Find deep outputs with critical AND depth gap
2. Collect ALL side leaves from the AND chain (expand through AND to minimize leaf depth)
3. Keep expanding the deepest AND leaves recursively until no improvement
4. Build balanced tree from the expanded leaf set
5. Remap output to use new tree (duplicates computation for this output)
"""

import sys, os, re, math
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


def expand_leaf(name, stmt_map, input_set, chain_set, and_depth, depth_limit=2):
    """
    Expand a leaf if it's an AND node with depth > 0 and non-chain.
    Returns set of replacement leaves (or {name} if no expansion).
    depth_limit: only expand if leaf.depth >= depth_limit
    """
    name = name.strip()
    if name in input_set or name in chain_set:
        return {name}

    expr = stmt_map.get(name, '')
    if not expr:
        return {name}

    # NOT: transparent
    if expr.startswith('!'):
        op = expr[1:].strip()
        return expand_leaf(op, stmt_map, input_set, chain_set, and_depth, depth_limit)

    # AND: expand only if deep enough
    ops = get_and_operands(expr)
    if ops and and_depth.get(name, 0) >= depth_limit:
        result = set()
        for op in ops:
            result.update(expand_leaf(op, stmt_map, input_set, chain_set, and_depth, depth_limit))
        return result

    return {name}


def collect_side_leaves_expanded(chain_ands, stmt_map, input_set, chain_set, and_depth):
    """
    Collect side leaves with recursive expansion of deep AND leaves.
    Finds the optimal expansion depth that minimizes max_leaf_depth + tree_height.
    """
    # First identify all side leaves (without expansion beyond chain)
    def find_predecessor(and_name):
        expr = stmt_map[and_name]
        ops = get_and_operands(expr)
        if not ops:
            return None
        lhs, rhs = ops
        for op in [lhs, rhs]:
            cur = op
            while cur in stmt_map and stmt_map[cur].startswith('!'):
                cur = stmt_map[cur][1:].strip()
            if cur in chain_set:
                return op
        return None

    # Initial side leaves (no expansion through ANDs)
    side_leaves = {}
    for and_name in chain_ands:
        expr = stmt_map[and_name]
        ops = get_and_operands(expr)
        if not ops:
            continue
        lhs, rhs = ops
        pred = find_predecessor(and_name)
        for op in [lhs, rhs]:
            if op == pred:
                continue
            # Expand through NOT only
            leaf_name = op.strip()
            while leaf_name in stmt_map and stmt_map[leaf_name].startswith('!'):
                leaf_name = stmt_map[leaf_name][1:].strip()
            if leaf_name not in side_leaves:
                side_leaves[leaf_name] = 0
            side_leaves[leaf_name] += 1

    # Now iteratively expand deep AND leaves
    # Try different depth limits
    best_leaves = set(side_leaves.keys())
    best_score = max(and_depth.get(l, 0) for l in best_leaves) + \
        max(0, (len(best_leaves) - 1)).bit_length()

    for depth_limit in range(2, 0, -1):  # try limit=2, then 1
        leaves = set(side_leaves.keys())
        changed = True
        while changed:
            changed = False
            new_leaves = set()
            for l in sorted(leaves, key=lambda x: -and_depth.get(x, 0)):
                expanded = expand_leaf(l, stmt_map, input_set, chain_set,
                                       and_depth, depth_limit)
                new_leaves.update(expanded)
            if new_leaves != leaves:
                leaves = new_leaves
                changed = True

        if len(leaves) < 2:
            continue

        max_d = max(and_depth.get(l, 0) for l in leaves)
        tree_h = max(0, (len(leaves) - 1)).bit_length()
        score = max_d + tree_h

        if score < best_score:
            best_score = score
            best_leaves = leaves

    return best_leaves


def build_balanced(factors, prefix, counter):
    nodes = {}
    queue = sorted(set(factors))
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


def optimize_depth_v3(inputs, outputs, stmt_map, min_chain=3,
                      max_and_blowup=4.0, min_depth_improvement=1):
    """
    Aggressive AND depth optimization with leaf expansion.
    For each deep output: collect side leaves, expand deep AND leaves,
    build balanced tree.
    """
    input_set = set(inputs)
    output_set = set(outputs)
    and_depth, total_and = compute_metrics(inputs, outputs, stmt_map, input_set)

    n_inputs = len(inputs)
    theoretical_min = max(1, (n_inputs - 1)).bit_length()

    output_depths = [(o, and_depth.get(o, 0)) for o in outputs
                     if and_depth.get(o, 0) >= min_chain]
    output_depths.sort(key=lambda x: -x[1])

    new_nodes = {}
    remap = {}
    counter = [0]
    and_added = 0
    and_budget = int(total_and * max_and_blowup) - total_and

    for out_name, depth in output_depths:
        if depth <= theoretical_min + 1:
            continue
        if len(new_nodes) > and_budget:
            break

        path = trace_critical_path(out_name, and_depth, stmt_map, input_set)
        path_names = set(n for n, _ in path)

        # Collect AND chain
        chain_ands = []
        for name, d in path:
            if name in input_set:
                break
            if is_pure_and(stmt_map.get(name, '')):
                chain_ands.append(name)
        chain_ands.reverse()  # deepest first

        if len(chain_ands) < min_chain:
            continue

        chain_set = set(chain_ands)

        # Collect side leaves with expansion
        leaves = collect_side_leaves_expanded(
            chain_ands, stmt_map, input_set, chain_set, and_depth)

        if len(leaves) < min_chain:
            continue

        max_leaf_depth = max(and_depth.get(l, 0) for l in leaves)
        tree_height = max(0, (len(leaves) - 1)).bit_length()
        new_depth_est = max_leaf_depth + tree_height
        depth_improvement = depth - new_depth_est

        if depth_improvement < min_depth_improvement:
            continue

        # Build balanced tree
        balanced, root = build_balanced(leaves, '_d3_', counter)
        if not balanced:
            continue

        new_and_count = len(balanced)
        old_chain_and_count = len(chain_ands)
        and_delta = new_and_count - old_chain_and_count

        # Trade-off check
        if and_delta > 0 and depth_improvement * 5 < and_delta:
            continue

        # Check budget
        if and_added + and_delta > and_budget:
            continue

        new_nodes.update(balanced)
        # Remap the output-side AND to the new tree root
        shallowest_and = chain_ands[-1]  # closest to output
        remap[shallowest_and] = root
        and_added += and_delta

    if not remap:
        return stmt_map, 0, 0

    # Build result
    new_stmt = dict(stmt_map)
    new_stmt.update(new_nodes)

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
    and_improv = m0[1] - m1[1]

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
    print(f"  原始: AND={m0[1]}, 最大AND深度={orig_max_depth}")

    new_stmt, depth_improv, and_improv = optimize_depth_v3(
        inputs, outputs, stmt_map, min_chain=3, max_and_blowup=4.0)

    if depth_improv > 0:
        m1 = compute_metrics(inputs, outputs, new_stmt, input_set)
        new_max_depth = max(m1[0].get(o, 0) for o in outputs)
        print(f"  优化后: AND={m1[1]}, 最大AND深度={new_max_depth}")
        print(f"  深度: {orig_max_depth} → {new_max_depth} (-{depth_improv})")
        print(f"  AND: {m0[1]} → {m1[1]} ({m1[1]-m0[1]})")

        for o in outputs:
            od = m0[0].get(o, 0)
            nd = m1[0].get(o, 0)
            if od != nd:
                print(f"    {o}: {od} → {nd}")

        out_path = txt_path.replace('.txt', '_opt_d3.txt')
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
