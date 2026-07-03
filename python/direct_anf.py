"""
Direct SparseANF circuit parser — prunes irrelevant nodes, simplifies when
term count exceeds threshold, and integrates with simplification from anf_factor.

Usage:
    python3 direct_anf.py <circuit.txt> [threshold]

Outputs files compatibly with process_all.py:
    {name}_direct_trans.poly
    {name}_direct_expr.poly
    {name}_direct_T.poly
"""

import sys, os, re, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from anf_factor import SparseANF, simplify as anf_simplify
from circuit_simplify import anf_and, anf_xor, anf_not
from format_poly import merge_transformations, write_transformation, write_expressions

import numpy as np


class DirectANFParser:
    """Parse circuit → SparseANF for each node, no M,b framework."""

    def __init__(self, threshold=200000, verbose=False):
        self.threshold = threshold
        self.verbose = verbose
        self.inputs = []
        self.outputs = []
        self.n = 0
        self.nodes = {}       # name → SparseANF
        self._stmt_map = {}   # name → expr string
        self.simplify_count = 0
        self.max_T = 0
        # Stats for merge
        self._g_list = []
        self._M_list = []
        self._b_list = []

    def parse(self, text):
        raw = text.replace('\n', ' ').replace('\r', '')
        raw = re.sub(r'\s+', ' ', raw).strip()

        minp = re.search(r'INORDER\s*=\s*([^;]+);', raw)
        if minp:
            self.inputs = minp.group(1).strip().split()
            self.n = len(self.inputs)

        mout = re.search(r'OUTORDER\s*=\s*([^;]+);', raw)
        if mout:
            self.outputs = mout.group(1).strip().split()

        # Parse all statements
        body_start = max(
            raw.find(';', raw.find('INORDER')),
            raw.find(';', raw.find('OUTORDER')) + 1
        )
        for stmt in raw[body_start:].split(';'):
            stmt = stmt.strip()
            if not stmt:
                continue
            m = re.match(r'(\w+)\s*=\s*(.+)$', stmt)
            if m:
                name, expr = m.group(1), m.group(2).strip()
                self._stmt_map[name] = expr

        # Cone-of-influence pruning: only keep nodes needed for outputs
        needed = set(self.outputs)
        changed = True
        while changed:
            changed = False
            for name in list(needed):
                if name in self._stmt_map:
                    deps = re.findall(r'\b([a-zA-Z_]\w*)\b', self._stmt_map[name])
                    for d in deps:
                        if d == name: continue
                        if d not in needed:
                            needed.add(d)
                            changed = True
                elif name in self.inputs:
                    pass
        # Remove inputs and constants from needed set for topological sort
        stmt_names = [name for name in needed
                      if name in self._stmt_map and name not in self.inputs]

        # Topological sort the needed statements
        self._eval_order = self._topological_sort(stmt_names)
        if self.verbose:
            print(f"  DAG: {len(self._stmt_map)} stmts → {len(stmt_names)} needed nodes",
                  flush=True)

        # Initialize input nodes
        for j, nm in enumerate(self.inputs):
            self.nodes[nm] = SparseANF({1 << j: 1}, self.n)

        # Evaluate in topological order
        for name in self._eval_order:
            expr = self._stmt_map[name]
            self._eval(name, expr)

        return self

    def _topological_sort(self, stmt_names):
        """Kahn's algorithm for topological sort. Returns ordered name list."""
        in_degree = {}
        graph = {}
        name_set = set(stmt_names)

        for name in stmt_names:
            in_degree.setdefault(name, 0)
            graph.setdefault(name, [])
            deps = re.findall(r'\b([a-zA-Z_]\w*)\b', self._stmt_map[name])
            for d in deps:
                if d == name: continue
                if d in name_set:
                    graph.setdefault(d, []).append(name)
                    in_degree[name] = in_degree.get(name, 0) + 1

        queue = [n for n in stmt_names if in_degree.get(n, 0) == 0]
        result = []
        while queue:
            node = queue.pop(0)
            result.append(node)
            for succ in graph.get(node, []):
                in_degree[succ] -= 1
                if in_degree[succ] == 0:
                    queue.append(succ)

        if len(result) < len(stmt_names):
            missing = set(stmt_names) - set(result)
            if self.verbose:
                print(f"  ⚠ {len(missing)} nodes not sorted, appending", flush=True)
            result.extend(missing)
        return result

    def _eval(self, name, expr):
        while expr.startswith('(') and expr.endswith(')'):
            expr = expr[1:-1].strip()

        if expr in ('0', 'false'):
            self.nodes[name] = SparseANF({}, self.n)
            return
        if expr in ('1', 'true'):
            self.nodes[name] = SparseANF({0: 1}, self.n)
            return
        if expr in self.nodes:
            self.nodes[name] = self.nodes[expr]
            return

        result = self._parse_expr(expr)
        self.nodes[name] = result
        T = result.T()
        if T > self.max_T:
            self.max_T = T

    def _parse_expr(self, expr):
        expr = expr.strip()
        while expr.startswith('(') and expr.endswith(')'):
            expr = expr[1:-1].strip()

        if expr in self.nodes:
            return self.nodes[expr]

        # NOT
        if expr.startswith('!'):
            operand = expr[1:].strip()
            while operand.startswith('(') and operand.endswith(')'):
                operand = operand[1:-1].strip()
            a = self._parse_expr(operand)
            return SparseANF(anf_xor({0: 1}, a.terms), self.n)

        # XOR (top-level +)
        plus_parts = self._split_at_top_level(expr, '+')
        if len(plus_parts) > 1:
            result = self._parse_expr(plus_parts[0])
            for part in plus_parts[1:]:
                b = self._parse_expr(part)
                result = SparseANF(anf_xor(result.terms, b.terms), self.n)
                if result.T() > self.threshold:
                    result = self._try_simplify(result)
            return result

        # AND (top-level *)
        mult_parts = self._split_at_top_level(expr, '*')
        if len(mult_parts) > 1:
            result = self._parse_expr(mult_parts[0])
            for part in mult_parts[1:]:
                b = self._parse_expr(part)
                result = self._anf_and_cached(result, b)
                if result.T() > self.threshold:
                    result = self._try_simplify(result)
            return result

        return SparseANF({}, self.n)

    def _anf_and_cached(self, a, b):
        """AND two SparseANFs, optimizing for disjoint variable support."""
        # If supports are small, direct AND is fine
        return SparseANF(anf_and(a.terms, b.terms), self.n)

    def _split_at_top_level(self, expr, op):
        parts = []
        depth = 0
        current = []
        for ch in expr:
            if ch == '(':
                depth += 1
                current.append(ch)
            elif ch == ')':
                depth -= 1
                current.append(ch)
            elif ch == op and depth == 0:
                parts.append(''.join(current).strip())
                current = []
            else:
                current.append(ch)
        parts.append(''.join(current).strip())
        return [p for p in parts if p]

    def _try_simplify(self, f):
        """Try to reduce ANF term count via simplification."""
        if f.T() == 0 or f.n == 0:
            return f
        self.simplify_count += 1
        T_before = f.T()
        if self.verbose:
            print(f"  ⚡ simplify #{self.simplify_count}: T={T_before}, n={f.n}",
                  flush=True)
        try:
            g_s, M, b = anf_simplify(f, verbose=False)
            T_after = g_s.T()
            if T_after < T_before and T_after > 0:
                if self.verbose:
                    print(f"     T reduced: {T_before}→{T_after} (m={g_s.n})",
                          flush=True)
                # Convert back: g_s is in z-space, substitute to x-space
                from anf_factor import substitute_affine
                f_new = substitute_affine(g_s.terms, g_s.n, M, b)
                if f_new.T() <= T_before:
                    return SparseANF(f_new.terms, self.n)
        except Exception as e:
            if self.verbose:
                print(f"     simplify failed: {e}", flush=True)
        return f

    def get_output(self, name):
        if name in self.nodes:
            return self.nodes[name]
        return SparseANF({}, self.n)


def process_circuit(txt_path, threshold=200000, verbose=True):
    name = os.path.splitext(os.path.basename(txt_path))[0]
    if verbose:
        print(f"\n{'='*60}")
        print(f"  {name}")
        print(f"{'='*60}")

    t0 = time.time()
    with open(txt_path) as f:
        text = f.read()

    parser = DirectANFParser(threshold=threshold, verbose=verbose)
    parser.parse(text)
    elapsed = time.time() - t0

    if verbose:
        print(f"\n  输入: {parser.n}, 输出: {len(parser.outputs)}")
        print(f"  最大T: {parser.max_T}")
        print(f"  简化次数: {parser.simplify_count}")
        print(f"  耗时: {elapsed:.1f}s")
        for o in parser.outputs:
            if o in parser.nodes:
                T = parser.nodes[o].T()
                print(f"  {o}: T={T}")
            else:
                print(f"  {o}: (missing)")

    return parser, elapsed


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 direct_anf.py <circuit.txt> [threshold]")
        sys.exit(1)
    txt_path = sys.argv[1]
    threshold = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    process_circuit(txt_path, threshold)
