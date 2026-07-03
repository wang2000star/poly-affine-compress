"""
AND-XOR circuit optimizer v3 — targeted bottom-up with cost control.

Problems with v2:
1. Optimized ALL deep ANDs → AND count explosion
2. No targeting of critical-path-only nodes
3. Cumulative cost from many overlapping optimizations

Fixes in v3:
1. Only optimize nodes on critical paths of max-depth outputs
2. Stricter trade-off: each optimization must have depth/AND ratio ≤ 0.5
3. Multiple strategies per node: try different levels of AND expansion
4. Accept only the best trade-off per node
"""

import sys, os, re, heapq


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


def get_and_operands(expr):
    if is_pure_and(expr):
        m = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
        if m:
            return m.group(1), m.group(2)
    return None


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


def count_ands(stmt_map):
    return sum(1 for e in stmt_map.values() if is_pure_and(e))


class Dag:
    """Mutable DAG for AND-XOR circuits."""

    def __init__(self, inputs, outputs, stmt_map):
        self.inputs = inputs
        self.input_set = set(inputs)
        self.outputs = outputs
        self.output_set = set(outputs)
        self.stmt = dict(stmt_map)
        self._rebuild()

    def _rebuild(self):
        self.order = topo_sort(list(self.stmt.keys()), self.stmt, self.input_set)
        self._compute_depths()

    def _compute_depths(self):
        self.and_depth = {}
        for inp in self.inputs:
            self.and_depth[inp] = 0
        for name in self.order:
            if name in self.input_set:
                continue
            expr = self.stmt.get(name, '')
            if not expr:
                self.and_depth[name] = 0
                continue
            ds = [d for d in deps(expr) if d != name and d in self.and_depth]
            if is_pure_and(expr):
                self.and_depth[name] = max(self.and_depth.get(d, 0) for d in ds) + 1 if ds else 1
            else:
                self.and_depth[name] = max(self.and_depth.get(d, 0) for d in ds) if ds else 0

    def max_output_depth(self):
        return max(self.and_depth.get(o, 0) for o in self.outputs)

    def trace_critical_path(self, output):
        path = []
        name = output
        while name not in self.input_set:
            if name not in self.stmt:
                break
            expr = self.stmt[name]
            ds = [d for d in deps(expr) if d != name and d in self.and_depth]
            if not ds:
                break
            next_name = max(ds, key=lambda d: self.and_depth.get(d, 0))
            path.append((name, self.and_depth.get(name, 0)))
            name = next_name
        return path

    def get_and_chain(self, node_name):
        """Trace back to find AND chain on critical path.
        Returns (chain_ands_deepest_first, side_leaves_dict)."""
        path = self.trace_critical_path(node_name)
        chain_ands = []
        for name, d in path:
            if name in self.input_set:
                break
            expr = self.stmt.get(name, '')
            if is_pure_and(expr):
                chain_ands.append(name)
        chain_ands.reverse()

        if len(chain_ands) < 2:
            return [], {}

        chain_set = set(chain_ands)
        side_leaves = {}
        for and_name in chain_ands:
            expr = self.stmt.get(and_name, '')
            ops = get_and_operands(expr)
            if not ops:
                continue
            lhs, rhs = ops
            pred = None
            for op in [lhs, rhs]:
                cur = op
                while cur in self.stmt and self.stmt[cur].startswith('!'):
                    cur = self.stmt[cur][1:].strip()
                if cur in chain_set:
                    pred = op
                    break
            for op in [lhs, rhs]:
                if op == pred:
                    continue
                leaf = op.strip()
                while leaf in self.stmt and self.stmt[leaf].startswith('!'):
                    leaf = self.stmt[leaf][1:].strip()
                side_leaves[leaf] = self.and_depth.get(leaf, 0)

        return chain_ands, side_leaves


def expand_ands(leaves, dag, chain_set, max_depth=20):
    """
    Fully expand AND leaves down to atomic nodes (inputs or XOR nodes).
    Uses iterative deepening: expand deepest ANDs first.
    """
    result = {}
    for name, ld in leaves.items():
        sub = _expand(name, dag, chain_set, set(), 0, max_depth)
        for k, v in sub.items():
            if k not in result or v < result[k]:
                result[k] = v
    return result


def _expand(name, dag, chain_set, visited, depth, max_depth):
    if name in visited or name in chain_set or depth > max_depth:
        return {name: dag.and_depth.get(name, 0)}
    visited.add(name)

    if name in dag.input_set:
        return {name: 0}

    expr = dag.stmt.get(name, '')
    if not expr:
        return {name: dag.and_depth.get(name, 0)}

    # NOT: transparent
    if expr.startswith('!'):
        op = expr[1:].strip()
        return _expand(op, dag, chain_set, visited, depth + 1, max_depth)

    # AND: expand
    ops = get_and_operands(expr)
    if ops and name not in chain_set:
        result = {}
        for op in ops:
            sub = _expand(op, dag, chain_set, visited, depth + 1, max_depth)
            for k, v in sub.items():
                if k not in result or v < result[k]:
                    result[k] = v
        return result

    # XOR or other: atomic
    return {name: dag.and_depth.get(name, 0)}


def huffman_depth(leaves):
    """Compute Huffman tree depth. leaves: {name: weight}. Returns (depth, num_nodes)."""
    if not leaves:
        return 0, 0
    if len(leaves) == 1:
        return list(leaves.values())[0], 0

    heap = []
    tie = [0]
    for name, weight in leaves.items():
        tie[0] += 1
        heapq.heappush(heap, (weight, tie[0], name))

    nodes = 0
    while len(heap) >= 2:
        w1, _, a = heapq.heappop(heap)
        w2, _, b = heapq.heappop(heap)
        new_weight = max(w1, w2) + 1
        nodes += 1
        tie[0] += 1
        heapq.heappush(heap, (new_weight, tie[0], (a, b)))

    return heap[0][0], nodes


def build_tree(leaves, prefix, counter):
    """Build Huffman tree. Returns (nodes_dict, root_name)."""
    if not leaves:
        return {}, None
    if len(leaves) == 1:
        return {}, list(leaves.keys())[0]

    heap = []
    tie = [0]
    for name, weight in leaves.items():
        tie[0] += 1
        heapq.heappush(heap, (weight, tie[0], name))

    nodes = {}
    while len(heap) >= 2:
        w1, t1, a = heapq.heappop(heap)
        w2, t2, b = heapq.heappop(heap)
        new_weight = max(w1, w2) + 1
        a_name = a[2] if isinstance(a, list) else a
        b_name = b[2] if isinstance(b, list) else b

        counter[0] += 1
        new_name = f"{prefix}{counter[0]}"
        nodes[new_name] = f"{a_name} * {b_name}" if a_name < b_name else f"{b_name} * {a_name}"

        tie[0] += 1
        heapq.heappush(heap, (new_weight, tie[0], [a, b, new_name]))

    _, _, final = heapq.heappop(heap)
    root_name = final[2] if isinstance(final, list) else final
    return nodes, root_name


def try_optimize_node(dag, node_name, and_budget):
    """
    Try to optimize the AND chain feeding into node_name.
    Returns opt dict or None.
    """
    chain_ands, side_leaves = dag.get_and_chain(node_name)

    if len(chain_ands) < 3 or len(side_leaves) < 2:
        return None

    # Check that this node is near the top of the chain
    top_idx = max(0, len(chain_ands) - 5)
    top_chain = chain_ands[top_idx:]
    top_set = set(top_chain)
    if node_name not in top_set:
        return None

    chain_set = set(chain_ands)
    original_depth = dag.and_depth.get(node_name, 0)
    original_chain_len = len(chain_ands)

    # The deepest unchanged node in the chain
    bottom_node = chain_ands[0]
    bottom_depth = dag.and_depth.get(bottom_node, 0)

    best_result = None
    best_score = float('inf')

    # Strategy 1: No expansion, use side leaves + bottom anchor
    all_leaves = dict(side_leaves)
    all_leaves[bottom_node] = bottom_depth
    hd, hn = huffman_depth(all_leaves)
    impr = original_depth - hd
    if impr > 0 and hn - original_chain_len <= and_budget:
        score = (hn - original_chain_len) / impr
        if score < best_score and (hn - original_chain_len) <= impr * 5:
            best_score = score
            best_result = {
                'leaves': dict(all_leaves),
                'new_depth': hd,
                'improvement': impr,
                'and_delta': hn - original_chain_len,
                'strategy': 'no_expand',
            }

    # Strategy 2: Expand deep ANDs (threshold 3)
    expanded = expand_ands(side_leaves, dag, chain_set, max_depth=20)
    all_leaves2 = dict(expanded)
    all_leaves2[bottom_node] = bottom_depth
    hd2, hn2 = huffman_depth(all_leaves2)
    impr2 = original_depth - hd2
    if impr2 > 0 and hn2 - original_chain_len <= and_budget and len(all_leaves2) < 100:
        score = (hn2 - original_chain_len) / impr2
        if score < best_score and (hn2 - original_chain_len) <= impr2 * 8:
            best_score = score
            best_result = {
                'leaves': dict(all_leaves2),
                'new_depth': hd2,
                'improvement': impr2,
                'and_delta': hn2 - original_chain_len,
                'strategy': 'expand',
            }

    # Strategy 3: Partial chain (only top 3-4 levels)
    for window in [3, 4, 5]:
        if window >= len(chain_ands):
            continue
        cut_idx = len(chain_ands) - window  # 0 = deepest, end = output
        # Ensure bottom_node is at or below the window
        if cut_idx <= 0:
            continue
        window_bottom = chain_ands[cut_idx - 1]
        window_leaves = {}
        for and_name in chain_ands[cut_idx:]:
            expr = dag.stmt.get(and_name, '')
            ops = get_and_operands(expr)
            if not ops:
                continue
            lhs, rhs = ops
            pred = None
            for op in [lhs, rhs]:
                cur = op
                while cur in dag.stmt and dag.stmt[cur].startswith('!'):
                    cur = dag.stmt[cur][1:].strip()
                if cur in chain_set:
                    pred = op
                    break
            for op in [lhs, rhs]:
                if op == pred:
                    continue
                leaf = op.strip()
                while leaf in dag.stmt and dag.stmt[leaf].startswith('!'):
                    leaf = dag.stmt[leaf][1:].strip()
                window_leaves[leaf] = dag.and_depth.get(leaf, 0)

        if len(window_leaves) < 2:
            continue

        # Include bottom of window
        window_leaves[window_bottom] = dag.and_depth.get(window_bottom, 0)

        w_hd, w_hn = huffman_depth(window_leaves)
        # Original depth contribution from window levels:
        window_old = dag.and_depth.get(chain_ands[-1], 0) - dag.and_depth.get(window_bottom, 0)
        if w_hd < window_old + dag.and_depth.get(window_bottom, 0):
            effective_impr = (window_old + dag.and_depth.get(window_bottom, 0)) - w_hd
            old_ands = window
            and_delta_w = w_hn - old_ands
            if effective_impr > 0 and and_delta_w <= and_budget:
                score = and_delta_w / effective_impr
                if score < best_score and and_delta_w <= effective_impr * 10:
                    best_score = score
                    best_result = {
                        'leaves': dict(window_leaves),
                        'new_depth': w_hd,
                        'improvement': effective_impr,
                        'and_delta': and_delta_w,
                        'strategy': f'window_{window}',
                        'window_bottom': window_bottom,
                    }

    if best_result is None:
        return None

    # Build tree
    counter = [0]
    prefix = f'_o{node_name}_'
    nodes, root = build_tree(best_result['leaves'], prefix, counter)
    if nodes is None or root is None:
        return None

    best_result['nodes'] = nodes
    best_result['root'] = root
    best_result['chain_ands'] = chain_ands
    best_result['original_depth'] = original_depth
    return best_result


def optimize_circuit(inputs, outputs, stmt_map, max_and_blowup=1.5, max_passes=5):
    """
    Targeted AND depth optimization.
    Only optimizes nodes on critical paths of max-depth outputs.
    """
    dag = Dag(inputs, outputs, stmt_map)
    orig_and = count_ands(stmt_map)
    and_budget_total = int(orig_and * max_and_blowup) - orig_and
    total_and_added = 0

    print(f"  Original AND: {orig_and}, depth: {dag.max_output_depth()}")

    for iteration in range(max_passes):
        max_depth = dag.max_output_depth()
        if max_depth <= 3:
            break

        # Find all outputs at max depth
        max_outputs = [o for o in dag.outputs
                       if dag.and_depth.get(o, 0) >= max_depth]

        # Collect all AND nodes on critical paths of max-depth outputs
        # that are deep enough to matter (within 2 of max)
        critical_and_nodes = set()
        for out in max_outputs:
            path = dag.trace_critical_path(out)
            for name, d in path:
                if name in dag.input_set:
                    break
                if is_pure_and(dag.stmt.get(name, '')):
                    if d >= max_depth - 1:  # only top 2 levels
                        critical_and_nodes.add(name)

        if not critical_and_nodes:
            break

        # Process in topological order (bottom-up, but only critical nodes)
        # Sort by depth (shallowest of the critical nodes first)
        crit_sorted = sorted(critical_and_nodes,
                             key=lambda n: dag.and_depth.get(n, 0))

        improvements_in_pass = 0

        for node_name in crit_sorted:
            if total_and_added >= and_budget_total:
                break

            node_depth = dag.and_depth.get(node_name, 0)
            if node_depth < max_depth - 2:
                continue  # not critical enough

            # Check remaining budget for this node
            remaining = and_budget_total - total_and_added
            remaining_outputs = max(1,
                                    sum(1 for o in max_outputs
                                        if dag.and_depth.get(o, 0) >= max_depth))
            per_node_budget = max(5, remaining // remaining_outputs)

            opt = try_optimize_node(dag, node_name, per_node_budget)

            if opt is None:
                continue

            # Apply optimization
            for n, e in opt['nodes'].items():
                dag.stmt[n] = e
                if n not in dag.order:
                    dag.order.append(n)

            # Remap: replace the top chain node's expr
            top_node = opt['chain_ands'][-1]
            if top_node in dag.stmt:
                dag.stmt[top_node] = opt['root']

            total_and_added += opt['and_delta']

            # Rebuild
            dag._rebuild()

            impr = opt['improvement']
            and_delta = opt['and_delta']
            chain_len = len(opt['chain_ands'])
            n_leaves = len(opt['leaves'])
            strat = opt['strategy']

            print(f"  [{iteration}] {node_name}: depth {opt['original_depth']}→{opt['new_depth']} "
                  f"(-{impr}), AND {chain_len}→{chain_len+and_delta} ({and_delta:+d}), "
                  f"leaves={n_leaves}, {strat}")
            improvements_in_pass += 1

            # Recompute max depth
            new_max = dag.max_output_depth()
            if new_max < max_depth:
                print(f"  >>> Max depth reduced: {max_depth} → {new_max}")

        if improvements_in_pass == 0:
            if iteration == 0:
                print(f"  No improvements found in pass {iteration}")
            break

    # Cone-of-influence pruning
    needed = set(dag.outputs)
    changed = True
    while changed:
        changed = False
        for n in list(needed):
            if n in dag.stmt:
                for d in deps(dag.stmt[n]):
                    if d not in needed:
                        needed.add(d)
                        changed = True

    clean = {}
    for n in needed:
        if n in dag.stmt and n not in dag.input_set:
            clean[n] = dag.stmt[n]
        elif n in dag.output_set and n not in dag.stmt:
            clean[n] = '0'
    for o in dag.outputs:
        if o not in clean and o in dag.stmt and o not in dag.input_set:
            clean[o] = dag.stmt[o]

    new_and = count_ands(clean)
    m0 = compute_metrics(dag.inputs, dag.outputs, stmt_map, dag.input_set)
    m1 = compute_metrics(dag.inputs, dag.outputs, clean, dag.input_set)
    orig_max = max(m0[0].get(o, 0) for o in dag.outputs)
    new_max = max(m1[0].get(o, 0) for o in dag.outputs)

    return clean, orig_max - new_max, orig_and - new_and


def run(txt_path):
    name = os.path.splitext(os.path.basename(txt_path))[0]
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    inputs, outputs, stmt_map = parse_circuit(txt_path)
    input_set = set(inputs)

    m0 = compute_metrics(inputs, outputs, stmt_map, input_set)
    orig_max_depth = max(m0[0].get(o, 0) for o in outputs)
    orig_and = count_ands(stmt_map)
    print(f"  原始: AND={orig_and}, 最大深度={orig_max_depth}")

    clean, depth_improv, and_improv = optimize_circuit(
        inputs, outputs, stmt_map, max_and_blowup=1.5, max_passes=20)

    m1 = compute_metrics(inputs, outputs, clean, input_set)
    new_max_depth = max(m1[0].get(o, 0) for o in outputs)
    new_and = count_ands(clean)

    print(f"\n  最终: AND={orig_and} → {new_and} ({new_and - orig_and:+d})")
    print(f"  深度: {orig_max_depth} → {new_max_depth} ({new_max_depth - orig_max_depth:+d}")

    if depth_improv > 0:
        for o in outputs:
            od = m0[0].get(o, 0)
            nd = m1[0].get(o, 0)
            if od != nd:
                print(f"    {o}: {od} → {nd}")

        out_path = txt_path.replace('.txt', '_opt_s3.txt')
        with open(out_path, 'w') as f:
            f.write(f"INORDER = {' '.join(inputs)};\n")
            f.write(f"OUTORDER = {' '.join(outputs)};\n")
            order = topo_sort(list(clean.keys()), clean, input_set)
            for n in order:
                if n in input_set or n in set(outputs):
                    continue
                if n in clean:
                    f.write(f"{n} = {clean[n]};\n")
            for o in outputs:
                if o in clean:
                    f.write(f"{o} = {clean[o]};\n")
        print(f"  写入: {out_path}")
    else:
        print(f"  无有效优化")

    return clean, depth_improv, and_improv


if __name__ == '__main__':
    if len(sys.argv) < 2:
        for circ in ['hd09', 'hd10', 'hd11', 'hd12', 'router', 'dsort']:
            try:
                run(f'../examples/{circ}/{circ}.txt')
            except Exception as e:
                import traceback
                print(f"  ERROR: {e}")
                traceback.print_exc()
    else:
        run(sys.argv[1])
