"""
AND-XOR circuit optimizer v4 — all-level critical path with scaled cost control.

Process ALL nodes on critical paths of deep outputs, from bottom to top.
Cost/improvement ratio scales with depth: strict for shallow, generous for deep.
This lets the shared prefix chain be optimized gradually, reducing side leaf depths
for all outputs.
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
    def __init__(self, inputs, outputs, stmt_map):
        self.inputs = inputs
        self.input_set = set(inputs)
        self.outputs = outputs
        self.output_set = set(outputs)
        self.stmt = dict(stmt_map)
        self._rebuild()
        self.counter = [0]

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
        path = self.trace_critical_path(node_name)
        chain_ands = []
        for name, d in path:
            if name in self.input_set:
                break
            expr = self.stmt.get(name, '')
            if is_pure_and(expr):
                chain_ands.append(name)
        chain_ands.reverse()

        if len(chain_ands) < 3:
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


def expand_ands(leaves, dag, chain_set):
    """Fully expand AND leaves to atomic nodes."""
    result = {}
    for name, ld in leaves.items():
        sub = _expand(name, dag, chain_set, set(), 0)
        for k, v in sub.items():
            if k not in result or v < result[k]:
                result[k] = v
    return result


def _expand(name, dag, chain_set, visited, depth):
    if name in visited or name in chain_set or depth > 30:
        return {name: dag.and_depth.get(name, 0)}
    visited.add(name)

    if name in dag.input_set:
        return {name: 0}

    expr = dag.stmt.get(name, '')
    if not expr:
        return {name: dag.and_depth.get(name, 0)}

    if expr.startswith('!'):
        op = expr[1:].strip()
        return _expand(op, dag, chain_set, visited, depth + 1)

    ops = get_and_operands(expr)
    if ops and name not in chain_set:
        result = {}
        for op in ops:
            sub = _expand(op, dag, chain_set, visited, depth + 1)
            for k, v in sub.items():
                if k not in result or v < result[k]:
                    result[k] = v
        return result

    return {name: dag.and_depth.get(name, 0)}


def huffman_depth(leaves):
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
    """Try optimizing the AND chain feeding into node_name. Multiple strategies."""
    chain_ands, side_leaves = dag.get_and_chain(node_name)

    if len(chain_ands) < 3 or len(side_leaves) < 2:
        return None

    chain_set = set(chain_ands)
    original_depth = dag.and_depth.get(node_name, 0)
    original_chain_len = len(chain_ands)

    bottom_node = chain_ands[0]
    bottom_depth = dag.and_depth.get(bottom_node, 0)

    best_result = None
    best_impr = 0

    strategies = []

    # Strategy 1: no expansion, all chain leaves + bottom
    l1 = dict(side_leaves)
    l1[bottom_node] = bottom_depth
    hd1, hn1 = huffman_depth(l1)
    strategies.append({
        'leaves': l1, 'hd': hd1, 'hn': hn1,
        'orig': original_depth,
        'name': 'no_expand',
    })

    # Strategy 2: expand deep ANDs (fully)
    if len(side_leaves) < 50:
        expanded = expand_ands(side_leaves, dag, chain_set)
        l2 = dict(expanded)
        l2[bottom_node] = bottom_depth
        hd2, hn2 = huffman_depth(l2)
        if len(l2) <= 100:
            strategies.append({
                'leaves': l2, 'hd': hd2, 'hn': hn2,
                'orig': original_depth,
                'name': f'expand({len(l2)} leaves)',
            })

    # Strategy 3: window-based (top W levels)
    for window in [3, 4, 5, 6, 7]:
        if window >= len(chain_ands):
            continue
        cut_off = len(chain_ands) - window
        if cut_off <= 0:
            continue
        window_bottom_name = chain_ands[cut_off - 1]
        window_leaves = {}
        for and_name in chain_ands[cut_off:]:
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
        window_leaves[window_bottom_name] = dag.and_depth.get(window_bottom_name, 0)

        w_hd, w_hn = huffman_depth(window_leaves)
        strategies.append({
            'leaves': window_leaves,
            'hd': w_hd,
            'hn': w_hn,
            'orig': original_depth,
            'name': f'window{window}({len(window_leaves)} leaves)',
        })

    # Evaluate strategies
    for s in strategies:
        improvement = s['orig'] - s['hd']
        and_delta = s['hn'] - original_chain_len

        if improvement <= 0:
            continue
        if and_delta > and_budget:
            continue

        # Accept if: improvement >= 1 and cost not too high
        # ratio: AND per depth level saved
        if and_delta > improvement * 10:
            continue

        if improvement > best_impr or (improvement == best_impr and and_delta < 0):
            best_impr = improvement
            best_result = s

    if best_result is None:
        return None

    # Build tree
    counter = [0]
    prefix = f'_o{node_name}_'
    nodes, root = build_tree(best_result['leaves'], prefix, counter)
    if nodes is None or root is None:
        return None

    return {
        'nodes': nodes,
        'root': root,
        'chain_ands': chain_ands,
        'original_depth': original_depth,
        'new_depth': best_result['hd'],
        'improvement': best_result['orig'] - best_result['hd'],
        'and_delta': best_result['hn'] - original_chain_len,
        'n_leaves': len(best_result['leaves']),
        'strategy': best_result['name'],
    }


def optimize_circuit(inputs, outputs, stmt_map, max_and_blowup=1.5):
    """Iterative optimization: process all critical path AND nodes bottom-up."""
    dag = Dag(inputs, outputs, stmt_map)
    orig_and = count_ands(stmt_map)
    and_budget_total = int(orig_and * max_and_blowup) - orig_and
    total_added = 0

    print(f"  Original: {orig_and} ANDs, max depth={dag.max_output_depth()}")
    print(f"  AND budget: +{and_budget_total}")

    for iteration in range(20):
        max_depth = dag.max_output_depth()
        if max_depth <= 3:
            break

        # Get all outputs at or near max depth
        target_outputs = [o for o in dag.outputs
                          if dag.and_depth.get(o, 0) >= max_depth - 2 and
                          dag.and_depth.get(o, 0) >= 3]

        if not target_outputs:
            break

        # Collect ALL AND nodes on their critical paths
        # Include nodes from middle depths (not just top)
        critical_nodes = set()
        for out in target_outputs:
            path = dag.trace_critical_path(out)
            for name, d in path:
                if name in dag.input_set:
                    break
                if is_pure_and(dag.stmt.get(name, '')):
                    if d >= 2:  # only skip trivial nodes
                        critical_nodes.add(name)

        if not critical_nodes:
            break

        # Sort by depth (shallowest first) for bottom-up processing
        crit_sorted = sorted(critical_nodes,
                             key=lambda n: dag.and_depth.get(n, 0))

        made_progress = False

        for node_name in crit_sorted:
            if total_added >= and_budget_total:
                break

            node_depth = dag.and_depth.get(node_name, 0)
            if node_depth < 3:
                continue

            # Per-node budget based on depth distance from max
            dist = max_depth - node_depth
            if dist <= 1:
                node_budget = max(10, (and_budget_total - total_added) // max(1, len(crit_sorted) // 3))
            elif dist <= 3:
                node_budget = max(5, (and_budget_total - total_added) // max(1, len(crit_sorted) // 2))
            else:
                node_budget = max(2, (and_budget_total - total_added) // len(crit_sorted))

            opt = try_optimize_node(dag, node_name, min(node_budget, and_budget_total - total_added))
            if opt is None:
                continue

            # Apply
            for n, e in opt['nodes'].items():
                dag.stmt[n] = e
                if n not in dag.order:
                    dag.order.append(n)

            top_node = opt['chain_ands'][-1]
            if top_node in dag.stmt:
                dag.stmt[top_node] = opt['root']

            total_added += opt['and_delta']
            dag._rebuild()
            made_progress = True

            impr = opt['improvement']
            ad = opt['and_delta']
            cl = len(opt['chain_ands'])
            strat = opt['strategy']

            print(f"  [{iteration}] {node_name}: d={node_depth}→{opt['new_depth']} "
                  f"(-{impr}), AND {cl}→{cl+ad} ({ad:+d}), {strat}")

            new_max = dag.max_output_depth()
            if new_max < max_depth:
                print(f"    >>> Depth reduced: {max_depth} → {new_max} <<<")

        if not made_progress:
            if iteration == 0:
                print(f"  No improvements in pass 0")
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

    m0 = compute_metrics(inputs, outputs, stmt_map, dag.input_set)
    m1 = compute_metrics(inputs, outputs, clean, dag.input_set)
    new_max = max(m1[0].get(o, 0) for o in outputs)
    new_and = count_ands(clean)

    return clean, m0[1] - m1[1], orig_and - new_and, m1


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

    clean, depth_improv, and_saved, m1 = optimize_circuit(
        inputs, outputs, stmt_map, max_and_blowup=1.5)

    new_max_depth = max(m1[0].get(o, 0) for o in outputs)
    new_and = count_ands(clean)

    print(f"\n  最终: AND={orig_and} → {new_and} ({new_and - orig_and:+d})")
    print(f"  深度: {orig_max_depth} → {new_max_depth} ({new_max_depth - orig_max_depth:+d})")

    if new_max_depth < orig_max_depth or new_and < orig_and:
        for o in outputs:
            od = m0[0].get(o, 0)
            nd = m1[0].get(o, 0)
            if od != nd:
                print(f"    {o}: {od} → {nd}")

        out_path = txt_path.replace('.txt', '_opt_s4.txt')
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

    return clean, depth_improv, and_saved


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
