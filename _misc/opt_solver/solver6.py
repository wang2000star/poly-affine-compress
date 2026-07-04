"""
AND-XOR circuit optimizer v6 — iterative depth-aware Huffman tree.

Key insight: deep XOR side leaves reference chain nodes. Their effective AND depth
in the NEW tree depends on where those chain nodes are placed. We use iterative
weight relaxation: build Huffman tree → update leaf weights based on new chain
node positions → rebuild until convergence.

Also: expand ALL ANDs in side leaves (no chain_set stop) to find atomic factors.
"""

import sys, os, re, heapq
from collections import defaultdict


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


def trace_critical_path(output, stmt_map, and_depth, input_set):
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


def get_and_chain(stmt_map, output, and_depth, input_set):
    path = trace_critical_path(output, stmt_map, and_depth, input_set)
    chain_ands = []
    for name, d in path:
        if name in input_set:
            break
        if is_pure_and(stmt_map.get(name, '')):
            chain_ands.append(name)
    chain_ands.reverse()
    if len(chain_ands) < 2:
        return [], {}
    chain_set = set(chain_ands)
    side_leaves = {}
    for and_name in chain_ands:
        expr = stmt_map.get(and_name, '')
        ops = get_and_operands(expr)
        if not ops:
            continue
        lhs, rhs = ops
        pred = None
        for op in [lhs, rhs]:
            cur = op
            while cur in stmt_map and stmt_map[cur].startswith('!'):
                cur = stmt_map[cur][1:].strip()
            if cur in chain_set:
                pred = op
                break
        for op in [lhs, rhs]:
            if op == pred:
                continue
            leaf = op.strip()
            while leaf in stmt_map and stmt_map[leaf].startswith('!'):
                leaf = stmt_map[leaf][1:].strip()
            side_leaves[leaf] = and_depth.get(leaf, 0)
    return chain_ands, side_leaves


def expand_xor_leaf(name, stmt_map, input_set, and_depth, expand_chain_set, visited=None):
    """
    Expand a side leaf XOR expression by recursively expanding its AND sub-terms.
    For AND operands within XOR expressions, recursively expand to atomic factors.
    Returns dict of {name: effective_depth} where effective_depth accounts for
    chain node recomputation in the new tree.
    """
    if visited is None:
        visited = set()
    if name in visited:
        return {}
    visited.add(name)

    if name in input_set:
        return {name: 0}

    expr = stmt_map.get(name, '')
    if not expr:
        return {}

    # NOT: follow the operand
    if expr.startswith('!'):
        op = expr[1:].strip()
        return expand_xor_leaf(op, stmt_map, input_set, and_depth, expand_chain_set, visited)

    # AND within XOR: further expand
    ops = get_and_operands(expr)
    if ops:
        a, b = ops
        result = {}
        for op in [a, b]:
            # Always expand, even if in chain_set
            sub = expand_xor_leaf(op, stmt_map, input_set, and_depth, expand_chain_set, visited)
            for k, v in sub.items():
                if k not in result or v < result[k]:
                    result[k] = v
        return result

    # XOR or terminal: keep as leaf with its depth
    return {name: and_depth.get(name, 0)}


def expand_side_leaf(leaf, stmt_map, input_set, and_depth, expand_chain_set):
    """
    Expand a side leaf: if it's an AND, expand to atomic factors.
    If it's an XOR, keep it as a leaf but record its depth.
    """
    if leaf in input_set:
        return {leaf: 0}

    expr = stmt_map.get(leaf, '')
    if not expr:
        return {leaf: and_depth.get(leaf, 0)}

    if expr.startswith('!'):
        op = expr[1:].strip()
        return expand_side_leaf(op, stmt_map, input_set, and_depth, expand_chain_set)

    ops = get_and_operands(expr)
    if ops:
        # AND leaf: expand through chain set too
        result = {}
        for op in ops:
            sub = expand_xor_leaf(op, stmt_map, input_set, and_depth, expand_chain_set)
            for k, v in sub.items():
                if k not in result or v < result[k]:
                    result[k] = v
        return result

    # XOR or other: keep as is
    return {leaf: and_depth.get(leaf, 0)}


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


def compute_node_weight_in_tree(stmt, node, depth_from_bottom):
    """Compute the depth of a node if it were computed in a simple tree structure."""
    if node not in stmt:
        return 0
    expr = stmt[node]
    ops = get_and_operands(expr)
    if not ops:
        return 0
    a, b = ops
    return max(compute_node_weight_in_tree(stmt, a, depth_from_bottom),
               compute_node_weight_in_tree(stmt, b, depth_from_bottom)) + 1


def trace_chain_deps(leaf, stmt_map, chain_set, input_set, max_depth=5):
    """Trace chain node dependencies in an expression up to max_depth."""
    deps_set = set()

    def trace(name, depth):
        if depth > max_depth or name in input_set:
            return
        if name in chain_set:
            deps_set.add(name)
            return
        if name not in stmt_map:
            return
        expr = stmt_map[name]
        for d in deps(expr):
            if d != name:
                trace(d, depth + 1)

    trace(leaf, 0)
    return deps_set


def optimize_circuit(inputs, outputs, stmt_map, max_and_blowup=2.0):
    """Per-output optimization with iterative depth relaxation."""
    dag = Dag(inputs, outputs, stmt_map)
    orig_and = count_ands(stmt_map)
    and_budget = int(orig_and * max_and_blowup) - orig_and
    total_added = 0

    print(f"  Original: {orig_and} ANDs, max depth={dag.max_output_depth()}")
    print(f"  AND budget: +{and_budget}")

    counter = [0]

    # Phase 1: Iterative output optimization
    for iteration in range(8):
        max_depth = dag.max_output_depth()
        if max_depth <= 3:
            break

        candidates = [(o, dag.and_depth.get(o, 0))
                       for o in outputs
                       if dag.and_depth.get(o, 0) >= max(4, max_depth - 2)]
        candidates.sort(key=lambda x: -x[1])

        made_progress = False

        for out_name, depth in candidates:
            if total_added >= and_budget:
                break

            chain_ands, side_leaves = get_and_chain(
                dag.stmt, out_name, dag.and_depth, dag.input_set)

            if len(chain_ands) < 3:
                continue

            chain_set = set(chain_ands)
            original_depth = dag.and_depth.get(out_name, 0)

            # --- New approach: expand side leaves with chain set awareness ---
            # For each side leaf:
            #   - If it's a chain node: expand through it (it'll be recomputed)
            #   - If it's an XOR: expand its AND sub-terms through chain nodes
            #   - Track which leaves are "fixed" (no chain deps) vs "adjustable"

            expanded = {}
            adjustable = {}  # leaves whose weight might change with tree structure

            for leaf, ld in side_leaves.items():
                # Expand the side leaf
                sub = expand_side_leaf(leaf, dag.stmt, dag.input_set, dag.and_depth, chain_set)
                for k, v in sub.items():
                    if k not in expanded or v < expanded[k]:
                        expanded[k] = v
                    # Track if this leaf depends on chain nodes
                    chain_deps = trace_chain_deps(k, dag.stmt, chain_set, dag.input_set, 5)
                    if chain_deps:
                        adjustable[k] = chain_deps

            # Bottom node
            bottom = chain_ands[0]
            bottom_d = dag.and_depth.get(bottom, 0)

            # --- Iterative weight relaxation ---
            # Step 1: initial weights
            current_weights = dict(expanded)
            current_weights[bottom] = bottom_d

            best_weights = dict(current_weights)
            best_hd, best_hn = huffman_depth(current_weights)
            best_impr = original_depth - best_hd

            # Step 2: try weight relaxation (reduce weights of chain-dependent leaves)
            # For each adjustable leaf, try reducing its weight by 1-3 levels
            # This simulates the effect of chain node recomputation
            for reduction in [1, 2, 3]:
                trial_weights = dict(current_weights)
                for leaf, chain_deps in adjustable.items():
                    # Only reduce if the leaf has deep chain deps
                    # The reduction is based on the number of chain-dep levels saved
                    max_chain_dep = max(dag.and_depth.get(cd, 0) for cd in chain_deps)
                    chain_pos = min(dag.and_depth.get(cd, 0) for cd in chain_deps)
                    if max_chain_dep >= reduction:
                        trial_weights[leaf] = max(0, trial_weights[leaf] - reduction)

                hd, hn = huffman_depth(trial_weights)
                impr = original_depth - hd
                if impr > best_impr:
                    best_impr = impr
                    best_hd = hd
                    best_hn = hn
                    best_weights = dict(trial_weights)

            # Also try: remove chain-dependent leaves that will be recomputed
            # and instead add their chain sub-dependencies directly
            if False:  # disabled for now
                pass

            if best_impr <= 0:
                continue

            and_delta = best_hn - len(chain_ands)
            if and_delta > and_budget - total_added:
                continue
            if best_impr < 2 and and_delta > best_impr * 8:
                continue

            # Build the tree
            prefix = f'_s{out_name}_'
            nodes, root = build_tree(best_weights, prefix, counter)
            if nodes is None or root is None:
                continue

            top_and = chain_ands[-1]
            for n, e in nodes.items():
                if n not in dag.stmt:
                    dag.stmt[n] = e
                    if n not in dag.order:
                        dag.order.append(n)
            dag.stmt[top_and] = root

            dag._rebuild()
            total_added += and_delta
            made_progress = True

            new_d = dag.and_depth.get(out_name, 0)
            print(f"  [P1:{iteration}] {out_name}: d={original_depth}→{new_d} "
                  f"(-{original_depth - new_d}), "
                  f"AND {len(chain_ands)}→{len(chain_ands)+and_delta} ({and_delta:+d}), "
                  f"expand({len(expanded)} leaves)")

            if new_d < max_depth:
                new_max = dag.max_output_depth()
                print(f"    >>> Max depth: {max_depth} → {new_max} <<<")

        if not made_progress:
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
    return clean, m0[1] - m1[1], orig_and - count_ands(clean), m1


def run(txt_path, max_and_blowup=2.0):
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

    clean, depth_improv, and_improv, m1 = optimize_circuit(
        inputs, outputs, stmt_map, max_and_blowup=max_and_blowup)

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

        out_path = txt_path.replace('.txt', '_opt_s6.txt')
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
        for circ in ['hd09', 'hd11', 'hd12', 'router', 'dsort', 'hd10']:
            try:
                run(f'../examples/{circ}/{circ}.txt', max_and_blowup=2.0)
            except Exception as e:
                import traceback
                print(f"  ERROR: {e}")
                traceback.print_exc()
    else:
        run(sys.argv[1])
