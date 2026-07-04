"""
AND-XOR circuit optimizer v2 — bottom-up AND chain rebalancing.

Key insight: AND chains in these circuits are sequential (unbalanced).
A balanced AND tree can reduce depth when side leaves are shallow.
Processing AND nodes bottom-up ensures side leaves are already optimized.

Strategy:
1. Topological sort AND nodes by depth (shallowest first)
2. For each AND node, identify the AND chain feeding into it
3. Collect side leaves from the chain, expand deep AND leaves
4. Build Huffman-optimal AND tree
5. If depth improves and AND cost is acceptable, replace the chain
6. Recompute depths and continue
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
    """Mutable DAG representing the AND-XOR circuit."""

    def __init__(self, inputs, outputs, stmt_map):
        self.inputs = inputs
        self.input_set = set(inputs)
        self.outputs = outputs
        self.output_set = set(outputs)
        self.stmt = dict(stmt_map)
        self._rebuild()

    def _rebuild(self):
        """Recompute topological order and AND depths."""
        self.order = topo_sort(list(self.stmt.keys()), self.stmt, self.input_set)
        self._compute_depths()

    def _compute_depths(self):
        self.and_depth = {}
        for inp in self.inputs:
            self.and_depth[inp] = 0
        # Also initialize for all stmt nodes
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

    def trace_critical_path(self, output):
        """Return [(name, depth), ...] deepest first."""
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
        """Trace back from a node to find the AND chain on its critical path.
        Returns (chain_ands, side_leaves_dict) where chain_ands is deepest-first
        and side_leaves maps each chain AND to its non-predecessor operand."""
        path = self.trace_critical_path(node_name)
        chain_ands = []
        for name, d in path:
            if name in self.input_set:
                break
            expr = self.stmt.get(name, '')
            if is_pure_and(expr):
                chain_ands.append(name)
        chain_ands.reverse()  # deepest first

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
            # Find predecessor (operand on critical path)
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
                # Strip NOTs
                leaf = op.strip()
                while leaf in self.stmt and self.stmt[leaf].startswith('!'):
                    leaf = self.stmt[leaf][1:].strip()
                side_leaves[leaf] = self.and_depth.get(leaf, 0)

        return chain_ands, side_leaves


def expand_deep_ands(leaves, dag, chain_set, depth_threshold=2):
    """
    Recursively expand AND leaves with depth >= depth_threshold.
    Returns expanded {leaf: depth} dict.
    """
    result = {}
    for name, ld in leaves.items():
        sub = _expand_one(name, ld, dag, chain_set, depth_threshold, set(), 0)
        for k, v in sub.items():
            if k not in result or v < result[k]:
                result[k] = v
    return result


def _expand_one(name, ld, dag, chain_set, depth_threshold, visited, depth):
    if name in visited or name in chain_set or depth > 20:
        return {name: ld}
    visited.add(name)

    if name in dag.input_set:
        return {name: 0}

    expr = dag.stmt.get(name, '')
    if not expr:
        return {name: ld}

    # NOT: transparent
    if expr.startswith('!'):
        op = expr[1:].strip()
        return _expand_one(op, dag.and_depth.get(op, ld), dag, chain_set,
                           depth_threshold, visited, depth + 1)

    # AND: expand if depth >= threshold
    ops = get_and_operands(expr)
    if ops and ld >= depth_threshold and name not in chain_set:
        result = {}
        for op in ops:
            op_ld = dag.and_depth.get(op, 0)
            sub = _expand_one(op, op_ld, dag, chain_set,
                              depth_threshold, visited, depth + 1)
            for k, v in sub.items():
                if k not in result or v < result[k]:
                    result[k] = v
        return result

    # XOR or other: atomic
    return {name: ld}


def huffman_tree_depth(leaves):
    """Compute the optimal AND tree depth using Huffman's algorithm.
    leaves: {name: weight}. Returns (depth, num_nodes)."""
    if not leaves:
        return 0, 0
    if len(leaves) == 1:
        return list(leaves.values())[0], 0

    heap = []
    tie = [0]
    for name, weight in leaves.items():
        tie[0] += 1
        heapq.heappush(heap, (weight, tie[0], name))

    nodes_created = 0
    while len(heap) >= 2:
        w1, _, a = heapq.heappop(heap)
        w2, _, b = heapq.heappop(heap)
        new_weight = max(w1, w2) + 1
        nodes_created += 1
        tie[0] += 1
        heapq.heappush(heap, (new_weight, tie[0], (a, b)))

    return heap[0][0], nodes_created


def build_huffman_tree(leaves, prefix, counter):
    """Build Huffman-optimal AND tree.
    Returns (nodes_dict, root_name)."""
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


def compute_huffman_depth_for_leaves(leaves):
    """Quick Huffman depth computation for a set of leaves."""
    return huffman_tree_depth(leaves)[0]


def optimize_chain(dag, chain_ands, side_leaves, max_and_budget=50):
    """
    Try to optimize an AND chain by replacing with balanced tree.
    Returns (new_nodes, new_root_name, depth_improvement, and_delta) or None.
    """
    if len(chain_ands) < 3:
        return None

    chain_set = set(chain_ands)

    # Original chain contribution
    # The chain depth from chain_ands[0] to chain_ands[-1]
    bottom_node = chain_ands[0]  # deepest
    top_node = chain_ands[-1]  # closest to output
    bottom_depth = dag.and_depth.get(bottom_node, 0)
    top_depth = dag.and_depth.get(top_node, 0)
    chain_contribution = top_depth - bottom_depth if top_depth > bottom_depth else len(chain_ands)

    # Try different expansion thresholds
    best_result = None
    best_score = float('inf')

    for threshold in [1, 2, 3, 4]:
        expanded = expand_deep_ands(side_leaves, dag, chain_set,
                                    depth_threshold=threshold)

        if len(expanded) < 2:
            continue

        # Compute Huffman depth
        huff_depth, huff_nodes = huffman_tree_depth(expanded)

        # The new top depth = max(bottom_depth, huff_depth of side leaves)
        # But the bottom node is ALSO in the leaf set (chain bottom)
        # Actually: the bottom node's predecessor is below the chain bottom.
        # We include the bottom node as a leaf in the tree.
        all_leaves = dict(expanded)
        all_leaves[bottom_node] = dag.and_depth.get(bottom_node, 0)

        if len(all_leaves) < 2:
            continue

        huff_depth_full, huff_nodes_full = huffman_tree_depth(all_leaves)

        # The new depth = huff_depth_full
        old_depth = dag.and_depth.get(top_node, 0)
        new_depth = huff_depth_full
        improvement = old_depth - new_depth

        if improvement <= 0:
            continue

        old_and_count = len(chain_ands)
        new_and_count = huff_nodes_full
        and_delta = new_and_count - old_and_count

        if and_delta > max_and_budget:
            continue

        # Check if depth improvement justifies AND cost
        if and_delta > 0 and improvement * 5 < and_delta:
            continue

        score = and_delta / max(1, improvement)
        if score < best_score:
            best_score = score
            best_result = {
                'nodes': None,  # will build later
                'root': None,
                'all_leaves': all_leaves,
                'new_depth': new_depth,
                'old_depth': old_depth,
                'improvement': improvement,
                'and_delta': and_delta,
                'threshold': threshold,
                'top_node': top_node,
            }

    if best_result is None:
        return None

    # Build the tree
    counter = [0]
    prefix = f'_opt_{top_node}_'
    nodes, root = build_huffman_tree(best_result['all_leaves'], prefix, counter)
    if nodes is None:
        return None

    best_result['nodes'] = nodes
    best_result['root'] = root
    return best_result


def optimize_circuit(inputs, outputs, stmt_map, max_and_blowup=2.0, max_passes=10):
    """
    Bottom-up AND chain optimization with iterative refinement.

    Key differences from solver.py:
    1. Process all AND nodes (not just outputs) in topological order
    2. Build balanced tree for each chain level
    3. Apply optimizations immediately (recompute depths)
    4. Repeat until convergence
    """
    dag = Dag(inputs, outputs, stmt_map)
    input_set = set(inputs)
    orig_and_count = count_ands(stmt_map)
    and_budget = int(orig_and_count * max_and_blowup) - orig_and_count

    total_and_added = 0
    all_new_nodes = {}

    for iteration in range(max_passes):
        # Get all AND nodes sorted by depth
        and_nodes = [(name, dag.and_depth.get(name, 0))
                     for name in dag.order
                     if is_pure_and(dag.stmt.get(name, ''))]
        and_nodes.sort(key=lambda x: (x[1], x[0]))  # shallowest first

        improvements_found = False

        for name, depth in and_nodes:
            if depth < 3:
                continue  # skip shallow nodes

            if total_and_added >= and_budget:
                break

            # Get chain feeding into this AND
            chain_ands, side_leaves = dag.get_and_chain(name)

            if len(chain_ands) < 3 or len(side_leaves) < 2:
                continue

            # The top of this chain IS this AND node (or close to it)
            # Make sure the chain ends at or near this node
            if name not in chain_ands[-3:]:
                continue  # skip if this AND isn't near the chain top

            opt = optimize_chain(dag, chain_ands, side_leaves,
                                 max_and_budget=and_budget - total_and_added)

            if opt is None:
                continue

            improvements_found = True

            # Apply optimization
            new_nodes = opt['nodes']
            root = opt['root']
            top_node = opt['top_node']

            # Add nodes to DAG
            for n, e in new_nodes.items():
                dag.stmt[n] = e
                if n not in dag.order:
                    dag.order.append(n)

            all_new_nodes.update(new_nodes)
            total_and_added += opt['and_delta']

            # Remap: replace the top chain node's expr with the new tree root
            if top_node in dag.stmt:
                dag.stmt[top_node] = root

            # Rebuild depths
            dag._rebuild()

            improvement = opt['improvement']
            and_delta = opt['and_delta']
            chain_len = len(chain_ands)

            print(f"  [{iteration}] {name}: depth {opt['old_depth']}→{opt['new_depth']} "
                  f"(-{improvement}), AND {chain_len}→{chain_len + and_delta} ({and_delta:+d}), "
                  f"chain={chain_len}, leaves={len(opt['all_leaves'])}")

        if not improvements_found:
            print(f"  No improvements found in pass {iteration}")
            break

    # Apply cone-of-influence pruning
    needed = set(outputs)
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
        if n in dag.stmt and n not in input_set:
            clean[n] = dag.stmt[n]
        elif n in dag.output_set and n not in dag.stmt:
            clean[n] = '0'
    for o in outputs:
        if o not in clean and o in dag.stmt and o not in input_set:
            clean[o] = dag.stmt[o]

    # Final metrics
    new_and = count_ands(clean)
    m0 = compute_metrics(inputs, outputs, stmt_map, input_set)
    m1 = compute_metrics(inputs, outputs, clean, input_set)
    orig_max = max(m0[0].get(o, 0) for o in outputs)
    new_max = max(m1[0].get(o, 0) for o in outputs)

    return clean, orig_max - new_max, orig_and_count - new_and


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
        inputs, outputs, stmt_map, max_and_blowup=2.0, max_passes=10)

    m1 = compute_metrics(inputs, outputs, clean, input_set)
    new_max_depth = max(m1[0].get(o, 0) for o in outputs)
    new_and = count_ands(clean)

    if depth_improv > 0:
        print(f"\n  优化后: AND={new_and}, 最大深度={new_max_depth}")
        print(f"  深度: {orig_max_depth} → {new_max_depth} (-{depth_improv})")
        print(f"  AND: {orig_and} → {new_and} ({new_and - orig_and:+d})")

        for o in outputs:
            od = m0[0].get(o, 0)
            nd = m1[0].get(o, 0)
            if od != nd:
                print(f"    {o}: {od} → {nd}")

        out_path = txt_path.replace('.txt', '_opt_s2.txt')
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
