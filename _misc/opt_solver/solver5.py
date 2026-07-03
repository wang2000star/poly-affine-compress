"""
AND-XOR circuit optimizer v5 — global per-output optimal trees.

For each deep output:
1. Expand ALL AND side leaves down to atomic factors (inputs/XORs)
2. Build ONE Huffman-optimal AND tree
3. Replace the entire AND chain with this tree
4. Only if depth improves significantly (>1 level)

This is aggressive but gives the BEST possible depth for each output.
Process outputs in dependency order so later outputs benefit.
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
    """Get AND chain and side leaves for an output."""
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


def expand_fully(name, stmt_map, input_set, chain_set, and_depth, visited=None):
    """Recursively expand AND gates to atomic leaves."""
    if visited is None:
        visited = set()
    if name in visited or name in chain_set:
        return {name: and_depth.get(name, 0)}
    visited.add(name)

    if name in input_set:
        return {name: 0}

    expr = stmt_map.get(name, '')
    if not expr:
        return {name: and_depth.get(name, 0)}

    if expr.startswith('!'):
        op = expr[1:].strip()
        return expand_fully(op, stmt_map, input_set, chain_set, and_depth, visited)

    ops = get_and_operands(expr)
    if ops:
        result = {}
        for op in ops:
            sub = expand_fully(op, stmt_map, input_set, chain_set, and_depth, visited)
            for k, v in sub.items():
                if k not in result or v < result[k]:
                    result[k] = v
        return result

    return {name: and_depth.get(name, 0)}


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


def optimize_circuit(inputs, outputs, stmt_map, max_and_blowup=1.5):
    """
    Per-output full-expansion optimization.
    For each deep output, expand ALL AND side leaves to atomic factors
    and build ONE optimal tree.
    """
    dag = Dag(inputs, outputs, stmt_map)
    orig_and = count_ands(stmt_map)
    and_budget = int(orig_and * max_and_blowup) - orig_and
    total_added = 0

    print(f"  Original: {orig_and} ANDs, max depth={dag.max_output_depth()}")
    print(f"  AND budget: +{and_budget}")

    # Order outputs by depth (deepest first)
    output_order = sorted(
        [(o, dag.and_depth.get(o, 0)) for o in outputs if dag.and_depth.get(o, 0) >= 3],
        key=lambda x: (-x[1], x[0]))

    counter = [0]

    # Phase 1: optimize deep outputs
    for iteration in range(5):
        max_depth = dag.max_output_depth()
        if max_depth <= 3:
            break

        # Recompute candidate outputs
        current_outputs = [(o, dag.and_depth.get(o, 0))
                           for o in outputs
                           if dag.and_depth.get(o, 0) >= max_depth - 1
                           and dag.and_depth.get(o, 0) >= 3]
        current_outputs.sort(key=lambda x: (-x[1], x[0]))

        made_progress = False

        for out_name, depth in current_outputs:
            if total_added >= and_budget:
                break

            # Get AND chain
            chain_ands, side_leaves = get_and_chain(
                dag.stmt, out_name, dag.and_depth, dag.input_set)

            if len(chain_ands) < 3:
                continue

            chain_set = set(chain_ands)
            original_depth = dag.and_depth.get(out_name, 0)

            # Strategy 1: no expansion
            bottom = chain_ands[0]
            bottom_d = dag.and_depth.get(bottom, 0)
            all_leaves = dict(side_leaves)
            all_leaves[bottom] = bottom_d
            hd1, hn1 = huffman_depth(all_leaves)
            impr1 = original_depth - hd1

            # Strategy 2: full expansion of AND leaves
            expanded = {}
            for leaf, ld in side_leaves.items():
                sub = expand_fully(leaf, dag.stmt, dag.input_set, chain_set, dag.and_depth)
                for k, v in sub.items():
                    if k not in expanded or v < expanded[k]:
                        expanded[k] = v

            all_expanded = dict(expanded)
            all_expanded[bottom] = bottom_d
            hd2, hn2 = huffman_depth(all_expanded)
            impr2 = original_depth - hd2

            # Pick best strategy
            best_hd, best_hn, best_leaves, strat = hd1, hn1, all_leaves, 'direct'
            best_impr = impr1

            if impr2 > impr1 or (impr2 == impr1 and (hn2 - len(chain_ands)) < (hn1 - len(chain_ands))):
                best_hd, best_hn, best_leaves = hd2, hn2, all_expanded
                best_impr = impr2
                strat = f'expand({len(all_expanded)} leaves)'

            if best_impr <= 0:
                continue

            and_delta = best_hn - len(chain_ands)
            if and_delta > and_budget - total_added:
                continue

            # Accept if improvement/AND ratio is OK
            if best_impr < 2 and and_delta > best_impr * 5:
                continue

            # Build tree
            prefix = f'_s{out_name}_'
            nodes, root = build_tree(best_leaves, prefix, counter)
            if nodes is None or root is None:
                continue

            # Apply: add new nodes, remap the output
            # The output's current top chain AND gets replaced by the new tree root
            top_and = chain_ands[-1]

            for n, e in nodes.items():
                if n not in dag.stmt:
                    dag.stmt[n] = e
                    if n not in dag.order:
                        dag.order.append(n)

            # Remap: the top chain AND now computes the new tree root
            # (it already exists in dag.stmt, so we modify its expression)
            dag.stmt[top_and] = root

            # Rebuild depths
            dag._rebuild()
            total_added += and_delta
            made_progress = True

            new_out_d = dag.and_depth.get(out_name, 0)
            print(f"  [{iteration}] {out_name}: d={original_depth}→{new_out_d} "
                  f"(-{original_depth - new_out_d}), "
                  f"AND {len(chain_ands)}→{len(chain_ands)+and_delta} ({and_delta:+d}), "
                  f"{strat}")

            if new_out_d < max_depth:
                new_max = dag.max_output_depth()
                print(f"    >>> Max depth: {max_depth} → {new_max} <<<")

        if not made_progress:
            break

    # Phase 2: optimize deep internal AND nodes (shared by multiple outputs)
    if dag.max_output_depth() >= 3:
        for iteration in range(5):
            max_depth = dag.max_output_depth()
            if max_depth <= 3:
                break

            # Collect ALL deep AND nodes (not just outputs)
            all_ands = [(n, dag.and_depth.get(n, 0))
                        for n in dag.order
                        if is_pure_and(dag.stmt.get(n, ''))
                        and dag.and_depth.get(n, 0) >= max(3, max_depth - 4)
                        and n not in dag.output_set]
            all_ands.sort(key=lambda x: -x[1])  # deepest first (critical nodes get priority)

            if not all_ands:
                break

            made_progress = False
            for node_name, depth in all_ands:
                if total_added >= and_budget:
                    break

                if dag.and_depth.get(node_name, 0) < max(3, max_depth - 4):
                    continue

                chain_ands, side_leaves = get_and_chain(
                    dag.stmt, node_name, dag.and_depth, dag.input_set)

                if len(chain_ands) < 3 or len(side_leaves) < 2:
                    continue

                if node_name not in chain_ands[-3:]:
                    continue

                chain_set = set(chain_ands)
                original_depth = dag.and_depth.get(node_name, 0)

                # Full expansion
                expanded = {}
                for leaf, ld in side_leaves.items():
                    sub = expand_fully(leaf, dag.stmt, dag.input_set,
                                       chain_set, dag.and_depth)
                    for k, v in sub.items():
                        if k not in expanded or v < expanded[k]:
                            expanded[k] = v

                bottom = chain_ands[0]
                bottom_d = dag.and_depth.get(bottom, 0)
                all_leaves = dict(expanded)
                all_leaves[bottom] = bottom_d
                hd, hn = huffman_depth(all_leaves)
                impr = original_depth - hd
                and_delta = hn - len(chain_ands)

                if impr <= 0 or and_delta > and_budget - total_added:
                    continue
                if impr < 2 and and_delta > impr * 12:
                    continue

                prefix = f'_in{node_name}_'
                nodes, root = build_tree(all_leaves, prefix, counter)
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

                new_d = dag.and_depth.get(node_name, 0)
                print(f"  [P2:{iteration}] {node_name}: d={original_depth}→{new_d} "
                      f"(-{original_depth - new_d}), "
                      f"AND {len(chain_ands)}→{len(chain_ands)+and_delta} ({and_delta:+d}), "
                      f"leaves={len(all_leaves)}")

                new_max = dag.max_output_depth()
                if new_max < max_depth:
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

    clean, depth_improv, and_improv, m1 = optimize_circuit(
        inputs, outputs, stmt_map, max_and_blowup=2.5)

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

        out_path = txt_path.replace('.txt', '_opt_s5.txt')
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
