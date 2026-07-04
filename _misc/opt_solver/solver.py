"""
AND-XOR circuit optimizer — optimal AND tree construction with weighted leaves.

Core algorithm: Given leaves with different AND depths, build the optimal
binary tree that minimizes max(leaf_depth + distance_to_root).

This is the Huffman-like algorithm where combining nodes a,b gives
weight = max(w(a), w(b)) + 1, and we always combine the two smallest-weight nodes.
"""

import sys, os, sys, re, heapq, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')


class AndXorDag:
    """AND-XOR circuit DAG."""

    def __init__(self, inputs, outputs, stmt_map):
        self.inputs = inputs
        self.input_set = set(inputs)
        self.outputs = outputs
        self.output_set = set(outputs)
        self.stmt = dict(stmt_map)

        # Compute topological order
        self._topo_sort()

        # Compute AND depth and support
        self._compute_metrics()

    def _deps(self, expr):
        return re.findall(r'\b([a-zA-Z_]\w*)\b', expr)

    def _is_pure_and(self, expr):
        return (isinstance(expr, str) and '*' in expr
                and not expr.startswith('!')
                and re.match(r'\w+\s*\*\s*\w+$', expr))

    def _topo_sort(self):
        names = [n for n in self.stmt if n not in self.input_set]
        in_deg = {}
        graph = {}
        name_set = set(names)
        for name in names:
            in_deg.setdefault(name, 0)
            for d in self._deps(self.stmt[name]):
                if d != name and d in name_set:
                    graph.setdefault(d, []).append(name)
                    in_deg[name] = in_deg.get(name, 0) + 1
        q = [n for n in names if in_deg.get(n, 0) == 0]
        self.order = []
        while q:
            n = q.pop(0)
            self.order.append(n)
            for s in graph.get(n, []):
                in_deg[s] -= 1
                if in_deg[s] == 0:
                    q.append(s)
        self.order.extend([n for n in names if n not in self.order])

    def _compute_metrics(self):
        self.and_depth = {}
        self.degree = {}
        for inp in self.inputs:
            self.and_depth[inp] = 0
            self.degree[inp] = 1

        for name in self.order:
            expr = self.stmt[name]
            ds = [d for d in self._deps(expr)
                  if d != name and d in self.and_depth]
            if not ds:
                self.and_depth[name] = 0
                self.degree[name] = 1
                continue

            if expr.startswith('!'):
                op = expr[1:].strip()
                self.and_depth[name] = self.and_depth.get(op, 0)
                self.degree[name] = self.degree.get(op, 1)
            elif self._is_pure_and(expr):
                m = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
                if m:
                    a, b = m.group(1), m.group(2)
                    da = self.and_depth.get(a, 0) if a != name else 0
                    db = self.and_depth.get(b, 0) if b != name else 0
                    self.and_depth[name] = max(da, db) + 1
                    self.degree[name] = (self.degree.get(a, 1) +
                                         self.degree.get(b, 1))
                else:
                    self.and_depth[name] = 1
                    self.degree[name] = 2
            elif '+' in expr:
                max_ad = max(self.and_depth.get(d, 0) for d in ds)
                self.and_depth[name] = max_ad
                self.degree[name] = max(self.degree.get(d, 1) for d in ds)
            else:
                self.and_depth[name] = max(self.and_depth.get(d, 0)
                                           for d in ds)
                self.degree[name] = max(self.degree.get(d, 1) for d in ds)

    def trace_critical_path(self, output):
        """Return critical path [(name, depth), ...] deepest first."""
        path = []
        name = output
        while name not in self.input_set:
            if name not in self.stmt:
                break
            expr = self.stmt[name]
            ds = [d for d in self._deps(expr)
                  if d != name and d in self.and_depth]
            if not ds:
                break
            next_name = max(ds, key=lambda d: self.and_depth.get(d, 0))
            path.append((name, self.and_depth.get(name, 0)))
            name = next_name
        return path

    def count_and(self):
        return sum(1 for e in self.stmt.values() if self._is_pure_and(e))

    def get_and_operands(self, expr):
        if self._is_pure_and(expr):
            m = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
            if m:
                return m.group(1), m.group(2)
        return None


def optimal_and_tree(leaves, prefix='_opt_', counter=None):
    """
    Build optimal AND tree from leaves with weights.

    Uses Huffman-like algorithm: always combine two smallest-weight nodes.
    Combining: weight = max(w1, w2) + 1.

    Args:
        leaves: dict of {name: weight} where weight = AND depth
        prefix: prefix for new node names
        counter: shared counter for unique naming

    Returns:
        (nodes_dict, root_name, root_weight)
    """
    if counter is None:
        counter = [0]

    if not leaves:
        return {}, None, -1

    if len(leaves) == 1:
        name = list(leaves.keys())[0]
        return {}, name, list(leaves.values())[0]

    # Priority queue: (weight, tiebreaker, item)
    # item is either a leaf name string or a list [left, right, name]
    # Using tiebreaker ensures stable ordering even when weights are equal
    tie = [0]

    heap = []
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

    # Last remaining node
    _, _, final = heapq.heappop(heap)
    root_name = final[2] if isinstance(final, list) else final
    root_weight = 0

    max_leaf_depth = max(leaves.values())
    # Root weight is the max leaf depth + tree height
    root_weight = max_leaf_depth + max(0, (len(leaves) - 1)).bit_length()
    # The Huffman-like algorithm gives us the actual root weight
    # Let's just compute it
    temp_nodes = dict(nodes)
    for n in leaves:
        temp_nodes[n] = None
    # Recompute actual weight from the structure
    actual_weight = _compute_node_weight(root_name, temp_nodes, leaves)

    return nodes, root_name, actual_weight


def _compute_node_weight(name, nodes, leaf_weights):
    """Compute effective weight of a node in the tree."""
    if name in leaf_weights:
        return leaf_weights[name]
    if name not in nodes or nodes[name] is None:
        return leaf_weights.get(name, 0)
    expr = nodes[name]
    m = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
    if m:
        a, b = m.group(1), m.group(2)
        wa = _compute_node_weight(a, nodes, leaf_weights)
        wb = _compute_node_weight(b, nodes, leaf_weights)
        return max(wa, wb) + 1
    return 0


def compute_tree_depth(tree_root, nodes, and_depth_cache):
    """Compute the exact AND depth of the new tree."""
    if tree_root in and_depth_cache:
        return and_depth_cache[tree_root]

    if tree_root not in nodes:
        return and_depth_cache.get(tree_root, 0)

    expr = nodes[tree_root]
    if '*' in expr and not expr.startswith('!'):
        m = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
        if m:
            a, b = m.group(1), m.group(2)
            da = compute_tree_depth(a, nodes, and_depth_cache)
            db = compute_tree_depth(b, nodes, and_depth_cache)
            return max(da, db) + 1

    return and_depth_cache.get(tree_root, 0)


def extract_chain_side_leaves(dag, output_name, chain_set, expand_and=False):
    """
    Extract side leaves from the AND chain feeding an output.

    Returns dict of {leaf_name: and_depth} for each side leaf.
    expand_and: if True, recursively expand AND side leaves
    """
    path = dag.trace_critical_path(output_name)
    path_names = set(n for n, _ in path)

    # Collect AND chain
    chain_ands = []
    for name, d in path:
        if name in dag.input_set:
            break
        if dag._is_pure_and(dag.stmt.get(name, '')):
            chain_ands.append(name)
    chain_ands.reverse()

    # Find predecessors
    def find_pred(and_name):
        expr = dag.stmt[and_name]
        ops = dag.get_and_operands(expr)
        if not ops:
            return None
        lhs, rhs = ops
        for op in [lhs, rhs]:
            cur = op
            while cur in dag.stmt and dag.stmt[cur].startswith('!'):
                cur = dag.stmt[cur][1:].strip()
            if cur in chain_set:
                return op
        return None

    # Collect side leaves
    side_leaves = {}
    for and_name in chain_ands:
        expr = dag.stmt[and_name]
        ops = dag.get_and_operands(expr)
        if not ops:
            continue
        lhs, rhs = ops
        pred = find_pred(and_name)
        for op in [lhs, rhs]:
            if op == pred:
                continue
            # Get the leaf name (through NOTs)
            leaf = op.strip()
            while leaf in dag.stmt and dag.stmt[leaf].startswith('!'):
                leaf = dag.stmt[leaf][1:].strip()
            if leaf not in side_leaves:
                side_leaves[leaf] = dag.and_depth.get(leaf, 0)

    if not expand_and:
        return side_leaves

    # Recursively expand AND leaves
    changed = True
    while changed:
        changed = False
        new_leaves = {}
        for leaf, weight in side_leaves.items():
            if leaf in dag.input_set or leaf in chain_set:
                new_leaves[leaf] = weight
                continue
            expr = dag.stmt.get(leaf, '')
            if not expr:
                new_leaves[leaf] = weight
                continue
            if expr.startswith('!'):
                new_leaves[leaf] = weight
                continue
            ops = dag.get_and_operands(expr)
            if ops and weight > 0:
                # Expand this AND
                for op in ops:
                    l = op.strip()
                    while l in dag.stmt and dag.stmt[l].startswith('!'):
                        l = dag.stmt[l][1:].strip()
                    new_leaves[l] = dag.and_depth.get(l, 0)
                changed = True
            else:
                new_leaves[leaf] = weight
        side_leaves = new_leaves

    return side_leaves


def optimize_output(dag, output_name, expand_and=True, min_chain=3):
    """
    Optimize a single output's AND depth using optimal tree construction.

    Returns (new_nodes, remap_from, remap_to, depth_improvement) or None.
    """
    depth = dag.and_depth.get(output_name, 0)
    if depth < min_chain:
        return None

    # Get critical path
    path = dag.trace_critical_path(output_name)
    path_names = set(n for n, _ in path)

    # Collect AND chain
    chain_ands = []
    for name, d in path:
        if name in dag.input_set:
            break
        if dag._is_pure_and(dag.stmt.get(name, '')):
            chain_ands.append(name)
    chain_ands.reverse()

    if len(chain_ands) < min_chain:
        return None

    chain_set = set(chain_ands)

    # Extract side leaves with depths
    side_leaves = extract_chain_side_leaves(
        dag, output_name, chain_set, expand_and=expand_and)

    if len(side_leaves) < 2:
        return None

    # Build optimal tree
    counter = [0]
    nodes, root, root_weight = optimal_and_tree(
        side_leaves, prefix=f'_opt_{output_name}_', counter=counter)

    if not nodes:
        return None

    # Compute new depth
    all_nodes = dict(nodes)
    for name in side_leaves:
        all_nodes[name] = None  # marker

    new_depth = root_weight

    if new_depth >= depth:
        return None

    # The shallowest AND (closest to output) gets remapped
    shallowest_and = chain_ands[-1]
    depth_improvement = depth - new_depth

    # Calculate AND change
    new_and_count = len([n for n in nodes
                         if dag._is_pure_and(nodes[n])])
    old_chain_and_count = len(chain_ands)
    and_delta = new_and_count - old_chain_and_count

    return {
        'nodes': nodes,
        'root': root,
        'shallowest_and': shallowest_and,
        'chain_ands': chain_ands,
        'new_depth': new_depth,
        'old_depth': depth,
        'depth_improvement': depth_improvement,
        'and_delta': and_delta,
    }


def apply_optimization(dag, optimizations):
    """
    Apply multiple output optimizations to the DAG.
    Returns (new_stmt_map, total_and_change, depth_changes).
    """
    new_stmt = dict(dag.stmt)
    remap = {}
    depth_changes = {}

    for opt in optimizations:
        if opt is None:
            continue
        # Add new nodes
        for name, expr in opt['nodes'].items():
            new_stmt[name] = expr
        # Remap
        remap[opt['shallowest_and']] = opt['root']
        depth_changes[opt['shallowest_and']] = opt['depth_improvement']

    # Apply remapping
    for name in list(new_stmt.keys()):
        if name in dag.input_set:
            continue
        expr = new_stmt[name]
        if isinstance(expr, str):
            for old, new_ in remap.items():
                if old in dag._deps(expr):
                    expr = re.sub(r'\b' + re.escape(old) + r'\b',
                                  new_, expr)
            new_stmt[name] = expr

    return new_stmt, depth_changes


def run_circuit(txt_path, expand_and=True, max_depth_target=None):
    """Run the solver on a single circuit."""
    import importlib
    from optimize_depth import parse_circuit, compute_metrics, topo_sort

    name = os.path.splitext(os.path.basename(txt_path))[0]
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    inputs, outputs, stmt_map = parse_circuit(txt_path)
    input_set = set(inputs)

    dag = AndXorDag(inputs, outputs, stmt_map)

    m0 = compute_metrics(inputs, outputs, stmt_map, input_set)
    orig_max_depth = max(m0[0].get(o, 0) for o in outputs)
    orig_and = m0[2]
    print(f"  原始: AND={orig_and}, 最大深度={orig_max_depth}")

    # Find deep outputs and try to optimize each
    output_depths = [(o, dag.and_depth.get(o, 0))
                     for o in outputs if dag.and_depth.get(o, 0) >= 3]
    output_depths.sort(key=lambda x: -x[1])

    all_optimizations = []
    total_and_added = 0

    for out_name, depth in output_depths:
        opt = optimize_output(dag, out_name, expand_and=expand_and)

        if opt is None or opt['depth_improvement'] <= 0:
            continue

        # Accept optimization
        all_optimizations.append(opt)
        total_and_added += opt['and_delta']

        print(f"  {out_name}: 深度 {opt['old_depth']} → {opt['new_depth']} "
              f"(-{opt['depth_improvement']}), AND {opt['and_delta']:+d}")

    if not all_optimizations:
        print(f"  无有效优化")
        return stmt_map, 0, 0

    # Apply all optimizations
    from optimize_depth import compute_metrics, topo_sort
    new_stmt = dict(stmt_map)
    remap = {}

    for opt in all_optimizations:
        for n, e in opt['nodes'].items():
            new_stmt[n] = e
        remap[opt['shallowest_and']] = opt['root']

    for name in list(new_stmt.keys()):
        if name in input_set:
            continue
        expr = new_stmt[name]
        if isinstance(expr, str):
            for old, new_ in remap.items():
                if old in re.findall(r'\b([a-zA-Z_]\w*)\b', expr):
                    expr = re.sub(r'\b' + re.escape(old) + r'\b', new_, expr)
            new_stmt[name] = expr

    # Cone-of-influence pruning
    needed = set(outputs)
    changed = True
    while changed:
        changed = False
        for n in list(needed):
            if n in new_stmt:
                for d in re.findall(r'\b([a-zA-Z_]\w*)\b', new_stmt[n]):
                    if d not in needed:
                        needed.add(d)
                        changed = True

    clean = {}
    for n in needed:
        if n in new_stmt and n not in input_set:
            clean[n] = new_stmt[n]
        elif n in set(outputs) and n not in new_stmt:
            clean[n] = '0'
    for o in outputs:
        if o not in clean and o in new_stmt and o not in input_set:
            clean[o] = new_stmt[o]

    m1 = compute_metrics(inputs, outputs, clean, input_set)
    new_max_depth = max(m1[0].get(o, 0) for o in outputs)
    new_and = m1[2]

    print(f"\n  最终: AND={orig_and} → {new_and} ({new_and-orig_and:+d})")
    print(f"  深度: {orig_max_depth} → {new_max_depth} ({new_max_depth-orig_max_depth:+d})")

    # Write output
    out_path = txt_path.replace('.txt', '_opt_solved.txt')
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

    return clean, orig_max_depth - new_max_depth, orig_and - new_and


if __name__ == '__main__':
    for circ in ['hd09', 'hd10', 'hd11', 'hd12', 'router', 'dsort']:
        try:
            run_circuit(f'examples/{circ}/{circ}.txt', expand_and=True)
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()
