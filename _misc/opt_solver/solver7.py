"""
AND-XOR circuit optimizer v7 — AND count reduction while preserving depth.

Strategies:
1. Boolean simplification: !a*b ⊕ a*b = b
2. XOR factoring: a*b ⊕ a*c = a*(b⊕c)

All optimizations verified against original circuit by random simulation.
"""

import sys, os, re
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
    return set(re.findall(r'\b([a-zA-Z_]\w*)\b', expr))


def is_and(expr):
    return isinstance(expr, str) and '*' in expr and not expr.startswith('!')


def is_pure_and(expr):
    return is_and(expr) and re.match(r'\w+\s*\*\s*\w+$', expr)


def get_and_ops(expr):
    m = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
    if m:
        return m.group(1), m.group(2)
    return None


def count_ands(stmt_map):
    return sum(1 for e in stmt_map.values() if is_pure_and(e))


def topo_sort(names, stmt_map, input_set):
    names = [n for n in names if n not in input_set]
    in_deg = {}
    graph = {}
    name_set = set(names)
    for name in names:
        in_deg.setdefault(name, 0)
        for d in deps(stmt_map.get(name, '')):
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


def compute_depth(inputs, outputs, stmt_map, input_set=None):
    if input_set is None:
        input_set = set(inputs)
    and_depth = {}
    for nm in inputs:
        and_depth[nm] = 0
    order = topo_sort(list(stmt_map.keys()), stmt_map, input_set)
    for name in order:
        expr = stmt_map.get(name, '')
        if not expr:
            and_depth[name] = 0
            continue
        ds = [d for d in deps(expr) if d != name and d in and_depth]
        if is_pure_and(expr):
            and_depth[name] = max(and_depth.get(d, 0) for d in ds) + 1 if ds else 1
        else:
            and_depth[name] = max(and_depth.get(d, 0) for d in ds) if ds else 0
    return and_depth


def build_usage(stmt_map, inputs, outputs):
    """Count how many times each node is used as a dependency."""
    usage = defaultdict(int)
    for o in outputs:
        usage[o] += 1  # outputs are always "used"
    input_set = set(inputs)
    for name, expr in stmt_map.items():
        if name in input_set:
            continue
        for d in deps(expr):
            usage[d] += 1
    return usage


def verify_circuit(inputs, outputs, stmt_map_orig, stmt_map_new, num_tests=2000):
    """Verify that two circuits produce identical outputs for random inputs."""
    import random
    rng = random.Random(12345)

    input_set = set(inputs)
    n_inputs = len(inputs)

    # Compute topological order for original
    orig_order = topo_sort(list(stmt_map_orig.keys()), stmt_map_orig, input_set)
    new_order = topo_sort(list(stmt_map_new.keys()), stmt_map_new, input_set)

    for test_idx in range(num_tests):
        vals = {}
        for inp in inputs:
            vals[inp] = rng.randint(0, 1)

        # Eval original
        for name in orig_order:
            expr = stmt_map_orig.get(name, '')
            if not expr:
                vals[name] = 0
            elif expr.startswith('!'):
                op = expr[1:].strip()
                vals[name] = 1 ^ vals.get(op, 0)
            elif '*' in expr:
                ops = re.findall(r'\b(\w+)\b', expr)
                result = 1
                for op in ops:
                    result &= vals.get(op, 0)
                vals[name] = result
            elif '+' in expr:
                ops = re.findall(r'\b(\w+)\b', expr)
                result = 0
                for op in ops:
                    result ^= vals.get(op, 0)
                vals[name] = result
            else:
                vals[name] = vals.get(expr, 0)

        orig_out = tuple(vals.get(o, 0) for o in outputs)

        # Eval new
        for name in new_order:
            expr = stmt_map_new.get(name, '')
            if not expr:
                vals[name] = 0
            elif name in input_set:
                continue
            elif expr.startswith('!'):
                op = expr[1:].strip()
                vals[name] = 1 ^ vals.get(op, 0)
            elif '*' in expr:
                ops = re.findall(r'\b(\w+)\b', expr)
                result = 1
                for op in ops:
                    result &= vals.get(op, 0)
                vals[name] = result
            elif '+' in expr:
                ops = re.findall(r'\b(\w+)\b', expr)
                result = 0
                for op in ops:
                    result ^= vals.get(op, 0)
                vals[name] = result
            else:
                vals[name] = vals.get(expr, 0)

        new_out = tuple(vals.get(o, 0) for o in outputs)

        if orig_out != new_out:
            if test_idx < 20:
                print(f"    MISMATCH at test {test_idx}: inputs={[vals[i] for i in inputs]}")
                print(f"      orig={orig_out}, new={new_out}")
            return False, test_idx

    return True, num_tests


def boole_simplify(inputs, outputs, stmt_map):
    """
    Apply Boolean simplification: !a*b ⊕ a*b = b.

    Scans XOR expressions where both operands are ANDs that share an operand,
    and the non-shared operands are complementary (one is NOT of the other).

    Returns new statement map with AND count reduced, or None if no simplification found.
    """
    input_set = set(inputs)
    output_set = set(outputs)

    stmt = dict(stmt_map)
    usage = build_usage(stmt, inputs, outputs)
    changes = 0

    # Build NOT alias map: node_name -> base_name if node computes !base
    not_map = {}
    for n, e in stmt.items():
        if e.startswith('!'):
            not_map[n] = e[1:].strip()

    def resolve_with_not(node):
        """Return (base_name, is_not) respecting NOT aliases and direct ! prefix."""
        if node.startswith('!'):
            return node[1:], True
        if node in not_map:
            return not_map[node], True
        return node, False

    for name in list(stmt.keys()):
        if name in input_set:
            continue
        expr = stmt.get(name, '')
        if '+' not in expr:
            continue

        ops = sorted(deps(expr))
        if len(ops) != 2:
            continue

        a, b = ops

        # Both operands must be ANDs
        expr_a = stmt.get(a, '')
        expr_b = stmt.get(b, '')
        if not is_pure_and(expr_a) or not is_pure_and(expr_b):
            continue

        ops_a = get_and_ops(expr_a)
        ops_b = get_and_ops(expr_b)
        if not ops_a or not ops_b:
            continue

        a1, a2 = ops_a
        b1, b2 = ops_b

        # Resolve all operands through NOT aliases
        # This converts (n268, n305) where n268=!m_30 to (!m_30, n305)
        def resolved_set(a1, a2):
            r1 = resolve_with_not(a1)
            r2 = resolve_with_not(a2)
            return {r1, r2}

        set_a = resolved_set(a1, a2)
        set_b = resolved_set(b1, b2)

        # We need two names total: X and Y
        # One set has (!X, Y), the other has (X, Y)
        # Combined: {(!X,Y), (X,Y)} => names {X, Y}
        # X appears as both NOT and non-NOT, Y appears as non-NOT in both

        all_names = set()
        for n, notted in set_a:
            all_names.add(n)
        for n, notted in set_b:
            all_names.add(n)

        if len(all_names) != 2:
            continue

        names = list(all_names)
        n1, n2 = names[0], names[1]

        def check_complement_pattern(set_a, set_b, names):
            """Check if set_a and set_b form (!x, y) and (x, y) pattern."""
            counts = {}
            for n in names:
                counts[n] = {'not': 0, 'nonot': 0}
            # Iterate both sets separately (not union) to count multiplicity
            for n, notted in list(set_a) + list(set_b):
                if notted:
                    counts[n]['not'] += 1
                else:
                    counts[n]['nonot'] += 1

            comp_name = None
            shared_name = None
            for n in names:
                if counts[n]['not'] >= 1 and counts[n]['nonot'] >= 1:
                    comp_name = n
                elif counts[n]['not'] == 0 and counts[n]['nonot'] >= 2:
                    shared_name = n

            if comp_name is not None and shared_name is not None:
                return shared_name
            return None

        replacement = check_complement_pattern(set_a, set_b, names)

        if replacement is not None:
            # Only apply if both ANDs have usage == 1 (only used here)
            if usage.get(a, 0) > 1 or usage.get(b, 0) > 1:
                continue

            old_expr = stmt.get(name, '')
            stmt[name] = replacement

            changes += 1
            print(f"    Boolean simplify: {name} = {old_expr}  →  {replacement}")

    if changes == 0:
        return None

    # Cone-of-influence pruning
    needed = set(outputs)
    changed = True
    while changed:
        changed = False
        for n in list(needed):
            if n in stmt:
                for d in deps(stmt.get(n, '')):
                    if d not in needed:
                        needed.add(d)
                        changed = True

    clean = {}
    for n in needed:
        if n not in input_set:
            clean[n] = stmt.get(n, '')

    # Verify
    print(f"    Verifying Boolean simplification ({changes} changes)...")
    ok, ntests = verify_circuit(inputs, outputs, stmt_map, clean)
    if not ok:
        print(f"    *** VERIFICATION FAILED! Reverting. ***")
        return None

    return clean


def xor_factor(inputs, outputs, stmt_map):
    """
    Apply XOR factoring: a*b ⊕ a*c = a*(b⊕c).

    Scans for XOR nodes where both operands are ANDs that share a common operand.
    Creates a factored form and rewires.

    Does NOT fire when the non-shared operands are complementary (that's Boolean
    simplification territory: !a*b + a*b = b).
    """
    input_set = set(inputs)
    output_set = set(outputs)
    stmt = dict(stmt_map)
    usage = build_usage(stmt, inputs, outputs)
    not_map = {}
    for n, e in stmt.items():
        if e.startswith('!'):
            not_map[n] = e[1:].strip()
    counter = [0]
    changes = 0

    def resolve_with_not(node):
        if node.startswith('!'):
            return node[1:], True
        if node in not_map:
            return not_map[node], True
        return node, False

    def are_complementary(a, b):
        """Check if a and b are complementary (one is NOT of the other)."""
        ra, na = resolve_with_not(a)
        rb, nb = resolve_with_not(b)
        return ra == rb and na != nb

    # Collect XOR nodes where both operands are ANDs with shared operand
    candidates = []
    for name, expr in stmt.items():
        if name in input_set:
            continue
        if '+' not in expr:
            continue
        ops = sorted(deps(expr))
        if len(ops) != 2:
            continue
        a, b = ops

        # Both must be ANDs
        expr_a = stmt.get(a, '')
        expr_b = stmt.get(b, '')
        if not is_pure_and(expr_a) or not is_pure_and(expr_b):
            continue

        ops_a = get_and_ops(expr_a)
        ops_b = get_and_ops(expr_b)
        if not ops_a or not ops_b:
            continue

        a1, a2 = ops_a
        b1, b2 = ops_b

        # Check shared operand
        if a1 == b1 or a1 == b2 or a2 == b1 or a2 == b2:
            candidates.append((name, a, b, ops_a, ops_b))

    if not candidates:
        return None

    # Process candidates
    for name, a, b, ops_a, ops_b in candidates:
        a1, a2 = ops_a
        b1, b2 = ops_b

        # Find shared operand
        shared = None
        other_a = None
        other_b = None
        if a1 == b1:
            shared, other_a, other_b = a1, a2, b2
        elif a1 == b2:
            shared, other_a, other_b = a1, a2, b1
        elif a2 == b1:
            shared, other_a, other_b = a2, a1, b2
        elif a2 == b2:
            shared, other_a, other_b = a2, a1, b1
        else:
            continue

        # Must ensure both other_a and other_b are distinct from shared
        if other_a == shared or other_b == shared or other_a == other_b:
            continue

        # Skip if the non-shared operands are complementary
        # (Boolean simplification: !a*b + a*b = b — handles this case better)
        if are_complementary(other_a, other_b):
            continue

        # Only proceed if both a and b have usage == 1 (only used here)
        if usage.get(a, 0) > 1 or usage.get(b, 0) > 1:
            continue

        # Create new XOR: other_a + other_b
        counter[0] += 1
        new_xor_name = f'_fx{counter[0]}'
        stmt[new_xor_name] = f'{other_a} + {other_b}' if other_a < other_b else f'{other_b} + {other_a}'

        # Create new AND: shared * new_xor_name
        counter[0] += 1
        new_and_name = f'_fa{counter[0]}'
        stmt[new_and_name] = f'{shared} * {new_xor_name}' if shared < new_xor_name else f'{new_xor_name} * {shared}'

        # Replace the original XOR's expression
        old_expr = stmt.get(name, '')
        stmt[name] = new_and_name

        changes += 1
        print(f"    XOR factor: {name} = {old_expr}")
        print(f"      Factored: {shared} * ({other_a} + {other_b})")
        print(f"      Created: {new_xor_name} = {stmt[new_xor_name]}")
        print(f"      Created: {new_and_name} = {stmt[new_and_name]}")

    if changes == 0:
        return None

    # Cone-of-influence pruning
    needed = set(outputs)
    changed = True
    while changed:
        changed = False
        for n in list(needed):
            if n in stmt:
                for d in deps(stmt.get(n, '')):
                    if d not in needed:
                        needed.add(d)
                        changed = True

    clean = {}
    for n in needed:
        if n not in input_set:
            clean[n] = stmt.get(n, '')

    # Verify
    print(f"    Verifying XOR factoring ({changes} changes)...")
    ok, ntests = verify_circuit(inputs, outputs, stmt_map, clean)
    if not ok:
        print(f"    *** VERIFICATION FAILED! Reverting. ***")
        return None

    return clean


def xor_factor_deep(inputs, outputs, stmt_map):
    """
    Deep XOR factoring: expand XOR trees, collect AND leaf terms, find shared factors.

    Collects the full XOR tree rooted at each XOR expression, finds AND leaf terms,
    and factors out common sub-expressions.
    """
    input_set = set(inputs)
    output_set = set(outputs)
    stmt = dict(stmt_map)
    usage = build_usage(stmt, inputs, outputs)
    counter = [0]
    changes = 0

    # Helper: check if a node is an internal XOR node
    def is_xor(name):
        if name in input_set:
            return False
        e = stmt.get(name, '')
        return e and '+' in e

    # Collect XOR tree leaves: for an XOR node, recursively traverse
    # the XOR tree to find "base" terms (ANDs or inputs)
    def collect_xor_leaves(name, visited=None):
        if visited is None:
            visited = set()
        if name in visited or name in input_set:
            return [name]
        visited.add(name)
        expr = stmt.get(name, '')
        if not expr or expr.startswith('!'):
            return [name]
        if is_xor(name):
            ops = sorted(deps(expr))
            terms = []
            for op in ops:
                if op == name:
                    continue
                terms.extend(collect_xor_leaves(op, visited))
            return terms
        # AND or other — return as leaf
        return [name]

    # For each XOR node, try deep factoring
    xor_nodes = [n for n in stmt if is_xor(n) and n not in output_set]
    # Sort by usage count (ascending) so we factor nodes with fewer users first
    xor_nodes.sort(key=lambda n: usage.get(n, 0))

    processed = set()

    for name in xor_nodes:
        if name in processed:
            continue
        if usage.get(name, 0) == 0:
            continue

        # Collect XOR leaves
        leaves_raw = collect_xor_leaves(name)
        # Filter to AND leaves only and get their operands
        and_leaves = {}  # name -> (op1, op2)
        for leaf in leaves_raw:
            expr = stmt.get(leaf, '')
            ops_pair = get_and_ops(expr)
            if ops_pair:
                and_leaves[leaf] = ops_pair

        # Need at least 2 AND leaves to factor
        if len(and_leaves) < 2:
            continue

        # Count operand frequency
        op_freq = defaultdict(int)
        leaf_operand_map = defaultdict(list)  # operand -> [leaf_name]
        for leaf, (op1, op2) in and_leaves.items():
            op_freq[op1] += 1
            op_freq[op2] += 1
            leaf_operand_map[op1].append(leaf)
            leaf_operand_map[op2].append(leaf)

        # Find operands that appear in >= 2 leaves
        shared_ops = {op for op, freq in op_freq.items() if freq >= 2
                      and op not in input_set and not op.startswith('!')}

        if not shared_ops:
            continue

        # For each shared operand, try factoring
        for sh_op in sorted(shared_ops, key=lambda o: -op_freq[o]):
            affected_leaves = leaf_operand_map[sh_op]
            if len(affected_leaves) < 2:
                continue

            # Check that all affected leaves have usage == 1 (only used here)
            if any(usage.get(l, 0) > 1 for l in affected_leaves):
                continue

            # Check all affected leaves use sh_op in same position (or any, doesn't matter)
            other_ops = []
            can_factor = True
            for leaf in affected_leaves:
                op1, op2 = and_leaves[leaf]
                other = op2 if op1 == sh_op else (op1 if op2 == sh_op else None)
                if other is None:
                    can_factor = False
                    break
                other_ops.append(other)

            if not can_factor or len(other_ops) < 2:
                continue

            # Build XOR tree of the "other" operands
            # We replace the XOR tree that feeds into these ANDs

            # This is complex because the leaves feed into the XOR tree through
            # intermediate AND nodes. Need to restructure:
            #   Old: leaf1 ⊕ leaf2 ⊕ ... = (sh_op * other1) ⊕ (sh_op * other2) ⊕ ...
            #   New: sh_op * (other1 ⊕ other2 ⊕ ...)

            # But the XOR tree is binary and the leaves might be at different levels.
            # For now, only handle the simplest case: name itself is the XOR of two ANDs
            # that share an operand. This is already handled by xor_factor() above.

            # More advanced: trace up from affected_leaves to see if they all feed
            # into name through only XOR gates (no intermediate ANDs).
            # If they do, we can restructure.

            # For now, skip deep factoring (too complex)
            pass

    if changes == 0:
        return None

    return stmt


def xor_factor_deep(inputs, outputs, stmt_map):
    """
    Deep XOR factoring across full XOR trees.

    For each XOR node, collect all AND leaf terms in the XOR tree,
    find shared operands across AND leaves, and factor them out:
        a*X ⊕ a*Y ⊕ ... ⊕ other_terms  =  a*(X⊕Y⊕...) ⊕ other_terms

    This saves (group_size - 1) ANDs per factoring.
    """
    input_set = set(inputs)
    output_set = set(outputs)
    stmt = dict(stmt_map)
    counter = [0]
    changes = 0

    def is_xor(name):
        if name in input_set: return False
        e = stmt.get(name, '')
        return bool(e and '+' in e)

    def is_and_gate(name):
        if name in input_set: return False
        e = stmt.get(name, '')
        return bool(e and '*' in e and not e.startswith('!'))

    # Trace XOR tree leaves: follow XOR-only paths to find AND/input leaves.
    # Also returns the set of internal XOR nodes in the tree (for sharing checks).
    def xor_leaves_full(name, visited=None, internal_nodes=None):
        if internal_nodes is None:
            internal_nodes = set()
        if visited is None:
            visited = set()
        if name in visited or name in input_set:
            return {name}, internal_nodes
        visited.add(name)
        e = stmt.get(name, '')
        if not e or e.startswith('!'):
            return {name}, internal_nodes
        if is_xor(name):
            internal_nodes.add(name)
            terms = set()
            for d in deps(e):
                if d != name:
                    sub, _ = xor_leaves_full(d, visited, internal_nodes)
                    terms |= sub
            return terms, internal_nodes
        return {name}, internal_nodes

    # Build usage counts
    usage = build_usage(stmt, inputs, outputs)

    # Find XOR nodes with shared-operand AND leaves
    # Process XOR nodes by usage (ascending), prioritizing nodes with fewer users
    xor_nodes = [n for n in stmt if is_xor(n) and n not in output_set]
    xor_nodes.sort(key=lambda n: usage.get(n, 0))

    processed = set()

    for name in xor_nodes:
        if name in processed or usage.get(name, 0) == 0:
            continue

        leaves, internal = xor_leaves_full(name)
        # Check: all internal XOR nodes must have usage == 1 (this tree is private)
        if any(usage.get(n, 0) > 1 for n in internal):
            continue
        # Only AND leaves qualify for factoring
        and_leaves = {}
        for leaf in leaves:
            ops_pair = get_and_ops(stmt.get(leaf, ''))
            if ops_pair:
                and_leaves[leaf] = ops_pair

        if len(and_leaves) < 2:
            continue

        # Count operand frequencies
        op_freq = defaultdict(int)
        op_and_groups = defaultdict(list)  # operand -> [(leaf_name, other_operand)]
        for leaf, (op1, op2) in and_leaves.items():
            op_freq[op1] += 1
            op_freq[op2] += 1
            op_and_groups[op1].append((leaf, op2))
            op_and_groups[op2].append((leaf, op1))

        # Find most profitable shared operand (used by most ANDs)
        best_shared = None
        best_group_size = 0
        best_group_and_leaves = []

        for op, freq in sorted(op_freq.items(), key=lambda x: -x[1]):
            if freq < 2 or op in input_set:
                continue
            group = op_and_groups[op]
            # Filter: AND leaves must have exactly one user (only in THIS tree)
            valid = [(leaf, other) for leaf, other in group if usage.get(leaf, 0) == 1]
            if len(valid) >= 2 and len(valid) > best_group_size:
                # Also check: the "other" operands shouldn't all be the same
                others = set(other for _, other in valid)
                if len(others) >= 2:
                    best_shared = op
                    best_group_size = len(valid)
                    best_group_and_leaves = valid

        if best_shared is None or best_group_size < 2:
            continue

        # Found a group to factor!
        remaining_and = [l for l, _ in and_leaves.items() if l not in [x for x, _ in best_group_and_leaves]]
        remaining_non_and = [l for l in leaves if l not in and_leaves]

        # Build new XOR of "other" operands: other_X ⊕ other_Y ⊕ ...
        others = [other for _, other in best_group_and_leaves]
        # Lexicographic sort for determinism
        others.sort()
        # Chain XOR: iteratively build XOR tree
        if len(others) == 1:
            sum_name = others[0]
        else:
            sum_name = others[0]
            for i in range(1, len(others)):
                counter[0] += 1
                new_xor_name = f'_dx{counter[0]}'
                o1, o2 = sum_name, others[i]
                stmt[new_xor_name] = f'{o1} + {o2}' if o1 < o2 else f'{o2} + {o1}'
                sum_name = new_xor_name

        # Create new AND: best_shared * sum_name
        counter[0] += 1
        new_and_name = f'_da{counter[0]}'
        stmt[new_and_name] = f'{best_shared} * {sum_name}' if best_shared < sum_name else f'{sum_name} * {best_shared}'

        # Now the new XOR tree = new_and_name ⊕ remaining_leaves ⊕ remaining_non_and
        all_remaining = remaining_and + remaining_non_and
        all_remaining.sort()

        # Build XOR chain: new_and_name ⊕ rem1 ⊕ rem2 ⊕ ...
        xor_parts = [new_and_name] + all_remaining
        current = xor_parts[0]
        for i in range(1, len(xor_parts)):
            counter[0] += 1
            new_node = f'_dx{counter[0]}'
            o1, o2 = current, xor_parts[i]
            stmt[new_node] = f'{o1} + {o2}' if o1 < o2 else f'{o2} + {o1}'
            current = new_node

        # Replace the original XOR node's expression AND mark as processed
        stmt[name] = current

        saved_ands = best_group_size - 1
        changes += 1
        print(f"    Deep factor: {name} (shared={best_shared}, group={best_group_size})")
        print(f"      AND leaves: {len(and_leaves)} → {len(and_leaves) - saved_ands} (+{len(all_remaining)} remaining)")
        print(f"      New AND: {new_and_name} = {stmt[new_and_name]}")
        print(f"      Saved ANDs: {saved_ands}")

        # Mark affected nodes as processed
        for leaf, _ in best_group_and_leaves:
            processed.add(leaf)

    if changes == 0:
        return None

    # Cone-of-influence pruning
    needed = set(outputs)
    changed = True
    while changed:
        changed = False
        for n in list(needed):
            if n in stmt:
                for d in deps(stmt.get(n, '')):
                    if d not in needed:
                        needed.add(d)
                        changed = True

    clean = {}
    for n in needed:
        if n not in input_set:
            clean[n] = stmt.get(n, '')

    and_before_prune = count_ands(stmt)
    and_after_prune = count_ands(clean)
    print(f"    AND: {and_before_prune} (before prune) → {and_after_prune} (after prune), "
          f"orig was {count_ands(stmt_map)}")
    print(f"    Groups: {changes}, expected net: {-changes} (1 AND saved per group)")

    # Verify
    print(f"    Verifying deep XOR factoring ({changes} groups)...")
    ok, ntests = verify_circuit(inputs, outputs, stmt_map, clean)
    if not ok:
        print(f"    *** VERIFICATION FAILED! Reverting. ***")
        return None

    return clean


def apply_and_prune(stmt_map, inputs, outputs):
    """Apply cone-of-influence pruning."""
    input_set = set(inputs)
    needed = set(outputs)
    changed = True
    while changed:
        changed = False
        for n in list(needed):
            if n in stmt_map:
                for d in deps(stmt_map.get(n, '')):
                    if d not in needed:
                        needed.add(d)
                        changed = True

    clean = {}
    for n in needed:
        if n in stmt_map and n not in input_set:
            clean[n] = stmt_map[n]
        elif n in set(outputs) and n not in stmt_map:
            clean[n] = '0'
    for o in outputs:
        if o not in clean and o in stmt_map and o not in input_set:
            clean[o] = stmt_map[o]
    return clean


def run(txt_path):
    name = os.path.splitext(os.path.basename(txt_path))[0]
    input_set = set()
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    inputs, outputs, stmt_map = parse_circuit(txt_path)
    input_set = set(inputs)
    orig_and = count_ands(stmt_map)
    orig_depth = compute_depth(inputs, outputs, stmt_map, input_set)
    orig_max_depth = max(orig_depth.get(o, 0) for o in outputs)

    print(f"  原始: AND={orig_and}, 最大深度={orig_max_depth}")

    current = dict(stmt_map)
    total_and_saved = 0

    # Phase 1: Boolean simplification (!a*b ⊕ a*b = b)
    print(f"  --- Boolean simplification ---")
    result = boole_simplify(inputs, outputs, current)
    if result is not None:
        and_before = count_ands(current)
        and_after = count_ands(result)
        saved = and_before - and_after
        total_and_saved += saved
        current = result
        print(f"    AND: {and_before} → {and_after} ({-saved:+d})")

    # Phase 2: XOR factoring (a*b ⊕ a*c = a*(b⊕c))
    print(f"  --- XOR factoring ---")
    result = xor_factor(inputs, outputs, current)
    if result is not None:
        and_before = count_ands(current)
        and_after = count_ands(result)
        saved = and_before - and_after
        total_and_saved += saved
        current = result
        print(f"    AND: {and_before} → {and_after} ({-saved:+d})")

    # Phase 3: Deep XOR factoring (across full XOR trees)
    print(f"  --- Deep XOR factoring ---")
    result = xor_factor_deep(inputs, outputs, current)
    if result is not None:
        and_before = count_ands(current)
        and_after = count_ands(result)
        saved = and_before - and_after
        total_and_saved += saved
        current = result
        print(f"    AND: {and_before} → {and_after} ({-saved:+d})")

    # Iterate until no more improvements
    for iteration in range(5):
        made_progress = False
        print(f"  --- Iteration {iteration + 1} ---")

        r1 = boole_simplify(inputs, outputs, current)
        if r1 is not None:
            saved = count_ands(current) - count_ands(r1)
            total_and_saved += saved
            current = r1
            made_progress = True
            print(f"    Boolean: AND {count_ands(current) + saved} → {count_ands(current)} ({-saved:+d})")

        r2 = xor_factor(inputs, outputs, current)
        if r2 is not None:
            saved = count_ands(current) - count_ands(r2)
            total_and_saved += saved
            current = r2
            made_progress = True
            print(f"    Factor: AND {count_ands(current) + saved} → {count_ands(current)} ({-saved:+d})")

        r3 = xor_factor_deep(inputs, outputs, current)
        if r3 is not None:
            saved = count_ands(current) - count_ands(r3)
            total_and_saved += saved
            current = r3
            made_progress = True
            print(f"    Deep factor: AND {count_ands(current) + saved} → {count_ands(current)} ({-saved:+d})")

        if not made_progress:
            break

    # Final pruning
    current = apply_and_prune(current, inputs, outputs)

    new_and = count_ands(current)
    new_depth = compute_depth(inputs, outputs, current, input_set)
    new_max_depth = max(new_depth.get(o, 0) for o in outputs)

    print(f"\n  最终: AND={orig_and} → {new_and} ({new_and-orig_and:+d})")
    print(f"  深度: {orig_max_depth} → {new_max_depth} ({new_max_depth-orig_max_depth:+d})")

    if new_and < orig_and and new_max_depth <= orig_max_depth:
        for o in outputs:
            od = orig_depth.get(o, 0)
            nd = new_depth.get(o, 0)
            if od != nd:
                print(f"    {o}: {od} → {nd}")

        out_path = txt_path.replace('.txt', '_opt_s7.txt')
        with open(out_path, 'w') as f:
            f.write(f"INORDER = {' '.join(inputs)};\n")
            f.write(f"OUTORDER = {' '.join(outputs)};\n")
            order = topo_sort(list(current.keys()), current, input_set)
            for n in order:
                if n in input_set or n in set(outputs):
                    continue
                if n in current:
                    f.write(f"{n} = {current[n]};\n")
            for o in outputs:
                if o in current:
                    f.write(f"{o} = {current[o]};\n")
        print(f"  写入: {out_path}")
    else:
        if new_max_depth > orig_max_depth:
            print(f"  深度增加，不写入")
        else:
            print(f"  无有效优化")

    return current, total_and_saved


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
