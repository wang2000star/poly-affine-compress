"""
网表 → 简化 ANF 转换器

输入: AND/NOT 网表或直接 ANF
输出: (g, M, b) 满足 f(x) = g(Mx⊕b), 其中 g 是最简 ANF

核心策略:
1. 每个节点存 (g, M, b) 三元组: f(x) = g(Mx⊕b)
2. AND 时若 M 相同则同空间运算, 否则组合变换空间
3. T 超阈值时 simplify 压缩
"""
import numpy as np
from anf_factor import SparseANF, simplify
import time, re


# ===== ANF dict 运算 =====

def anf_not(d):
    if not d:
        return {0: 1}
    r = {0: 1}
    for m, c in d.items():
        r[m] = r.get(m, 0) ^ c
    return {k: v for k, v in r.items() if v}


def anf_xor(d1, d2):
    """d1 ⊕ d2"""
    r = dict(d1)
    for m, c in d2.items():
        r[m] = r.get(m, 0) ^ c
    return {k: v for k, v in r.items() if v}


def anf_and(d1, d2):
    if not d1 or not d2:
        return {}
    return SparseANF._anf_multiply(d1, d2)


def anf_or(d1, d2):
    """d1 OR d2 = d1 ⊕ d2 ⊕ d1·d2 (Boolean ring)"""
    return anf_xor(anf_xor(d1, d2), anf_and(d1, d2))


# ===== 变换组合 =====

def compose(M_sub, b_sub, M_acc, b_acc):
    """组合变换: z = M_sub·(M_acc·x ⊕ b_acc) ⊕ b_sub = M_new·x ⊕ b_new

    M_sub: m₂×m₁, M_acc: m₁×n  →  M_new: m₂×n
    b_sub: m₂,     b_acc: m₁    →  b_new: m₂
    """
    M_new = (np.array(M_sub, dtype=np.int64) @ np.array(M_acc, dtype=np.int64)) % 2
    b_sub_a = np.array(b_sub, dtype=np.int64).reshape(-1, 1)
    b_acc_a = np.array(b_acc, dtype=np.int64)
    b_new = (b_sub_a + np.array(M_sub, dtype=np.int64) @ b_acc_a.reshape(-1, 1)) % 2
    return M_new.tolist(), b_new.flatten().tolist()


def identity_M(n):
    return np.eye(n, dtype=np.int64).tolist()


def zero_b(n):
    return [0] * n


def same_matrix(M1, M2):
    """检查两个矩阵是否相同"""
    if len(M1) != len(M2):
        return False
    return all(all(a == b for a, b in zip(r1, r2)) for r1, r2 in zip(M1, M2))


# ===== 变换节点 =====

class TNode:
    """f(x) = g(Mx⊕b)

    g: ANF dict 在 m 变量空间中
    M: m×n 矩阵
    b: m×1 向量
    """
    __slots__ = ('g', 'M', 'b', 'm', 'n')

    def __init__(self, g, M, b, n):
        self.g = g
        self.M = M
        self.b = b
        self.m = len(M)
        self.n = n

    def copy(self):
        return TNode(dict(self.g), list(self.M), list(self.b), self.n)

    def __repr__(self):
        return f"TNode(m={self.m}, n={self.n}, T={len(self.g)})"


# ===== 网表解析器 =====

class CircuitSimplify:
    """解析 AND/NOT/OR 网表, 增量简化

    用法:
        cs = CircuitSimplify(threshold=4096)
        cs.parse(netlist_text)
        g, M, b = cs.result('om_0')
    """

    def __init__(self, threshold=4096, verbose=False):
        self.nodes = {}
        self.n = 0
        self.inputs = []
        self.outputs = []
        self.threshold = threshold
        self.verbose = verbose
        self.stats = {'simplify_calls': 0, 'simplify_time': 0.0}

    def add_input(self, name, n):
        """添加输入变量 (用于直接 ANF 输入而非网表)"""
        self.nodes[name] = TNode({1 << len(self.inputs): 1},
                                 identity_M(n), zero_b(n), n)

    # ---- 内部 AND 实现 ----

    def _compact_and_simplify(self, terms, m, M, b):
        """压缩变量空间后 simplify，返回 (g_s, M_s, b_s)"""
        used = set()
        for mask in terms:
            if mask == 0:
                continue
            t = mask
            while t:
                lsb = t & -t
                used.add((lsb.bit_length() - 1))
                t ^= lsb
        if 0 < len(used) < m:
            sorted_used = sorted(used)
            old_to_new = {old: new for new, old in enumerate(sorted_used)}
            compact_g = {}
            for mask, v in terms.items():
                new_mask = 0
                for old_k in sorted_used:
                    if (mask >> old_k) & 1:
                        new_mask |= (1 << old_to_new[old_k])
                compact_g[new_mask] = compact_g.get(new_mask, 0) ^ v
            compact_g = {k: v for k, v in compact_g.items() if v}
            compact_M = np.array([M[idx] for idx in sorted_used], dtype=np.int64)
            compact_b = [b[idx] for idx in sorted_used]
            g_s, M_s, b_s = simplify(SparseANF(compact_g, len(sorted_used)), verbose=False)
            M_new, b_new = compose(M_s.tolist(), b_s.tolist(),
                                    compact_M.tolist(), compact_b)
            return g_s, M_new, b_new
        else:
            g_s, M_s, b_s = simplify(SparseANF(terms, m), verbose=False)
            return g_s, M_s.tolist(), b_s.tolist()

    def _ensure_simplified(self, name):
        """简化并更新存储的节点（原地修改 nodes dict）"""
        node = self.nodes[name]
        if len(node.g) < self.threshold:
            return
        g_s, M_new, b_new = self._compact_and_simplify(node.g, node.m, node.M, node.b)
        self._stats_simplify(g_s.T(), len(node.g), node.m, g_s.n)
        self.nodes[name] = TNode(g_s.terms, M_new, b_new, node.n)

    def _and(self, name, name_A, name_B):
        # 预简化: 如果操作数太大, 先压缩再 AND（更新存储节点）
        self._ensure_simplified(name_A)
        self._ensure_simplified(name_B)
        A = self.nodes[name_A]
        B = self.nodes[name_B]

        # 情况 1: M 相同 → 同空间 AND
        if same_matrix(A.M, B.M):
            h = anf_and(A.g, B.g)
            h_T = len(h)
            if h_T > self.threshold:
                g_s, M_new, b_new = self._compact_and_simplify(h, A.m, A.M, A.b)
                self._stats_simplify(g_s.T(), h_T, A.m, g_s.n)
                self.nodes[name] = TNode(g_s.terms, M_new, b_new, A.n)
            else:
                self.nodes[name] = TNode(h, A.M, A.b, A.n)
            return

        # 情况 2: M 不同 → 组合变换空间
        mA, mB = A.m, B.m
        M_stacked = A.M + B.M
        b_stacked = A.b + B.b

        gB_shifted = {mk << mA: cv for mk, cv in B.g.items()}
        h = anf_and(A.g, gB_shifted)
        h_T = len(h)

        if h_T > self.threshold:
            g_s, M_new, b_new = self._compact_and_simplify(h, mA + mB, M_stacked, b_stacked)
            self._stats_simplify(g_s.T(), h_T, mA + mB, g_s.n)
            self.nodes[name] = TNode(g_s.terms, M_new, b_new, A.n)
        else:
            self.nodes[name] = TNode(h, M_stacked, b_stacked, A.n)

    def _stats_simplify(self, T_after, T_before, m_before, m_after):
        self.stats['simplify_calls'] += 1
        if self.verbose:
            ratio = (T_before - T_after) / T_before * 100 if T_before else 0
            print(f"  ⚡ simplify: T={T_before}→{T_after} ({ratio:.0f}%), "
                  f"m={m_before}→{m_after}")

    # ---- 网表解析 ----

    def parse(self, text):
        raw = text.replace('\n', ' ').replace('\r', '')
        raw = re.sub(r'\s+', ' ', raw).strip()

        minp = re.search(r'INORDER\s*=\s*([^;]+);', raw)
        if minp:
            self.inputs = minp.group(1).strip().split()
            self.n = len(self.inputs)
            I = identity_M(self.n)
            b0 = zero_b(self.n)
            for j, nm in enumerate(self.inputs):
                self.nodes[nm] = TNode({1 << j: 1}, I, b0, self.n)

        mout = re.search(r'OUTORDER\s*=\s*([^;]+);', raw)
        if mout:
            self.outputs = mout.group(1).strip().split()

        body_start = max(
            raw.find(';', raw.find('INORDER')),
            raw.find(';', raw.find('OUTORDER')) + 1
        )
        for stmt in raw[body_start:].split(';'):
            stmt = stmt.strip()
            if not stmt:
                continue
            m = re.match(r'(\w+)\s*=\s*(.+)$', stmt)
            if not m:
                continue
            name, expr = m.group(1), m.group(2).strip()
            self._parse_stmt(name, expr)
        return self

    def _parse_stmt(self, name, expr):
        # 直通: om_0 = i0
        if expr in self.nodes:
            self.nodes[name] = self.nodes[expr].copy()
            return

        # NOT: n65 = !i0
        mn = re.match(r'!(\w+)$', expr)
        if mn:
            src = self.nodes[mn.group(1)]
            self.nodes[name] = TNode(anf_not(src.g), src.M, src.b, src.n)
            return

        # AND: om_1 = i1 * n65
        ma = re.match(r'(\w+)\s*\*\s*(\w+)$', expr)
        if ma:
            self._and(name, ma.group(1), ma.group(2))
            return

        # OR: out = a + b
        mo = re.match(r'(\w+)\s*\+\s*(\w+)$', expr)
        if mo:
            name_A, name_B = mo.group(1), mo.group(2)
            self._ensure_simplified(name_A)
            self._ensure_simplified(name_B)
            A, B = self.nodes[name_A], self.nodes[name_B]
            if same_matrix(A.M, B.M):
                h = anf_or(A.g, B.g)
                if len(h) > self.threshold:
                    g_s, M_new, b_new = self._compact_and_simplify(h, A.m, A.M, A.b)
                    self._stats_simplify(g_s.T(), len(h), A.m, g_s.n)
                    self.nodes[name] = TNode(g_s.terms, M_new, b_new, A.n)
                else:
                    self.nodes[name] = TNode(h, A.M, A.b, A.n)
            else:
                raise NotImplementedError("OR with different M not implemented")
            return

        raise ValueError(f"语法错误: {name} = {expr}")

    # ---- 结果获取 ----

    def result(self, name):
        """获取简化结果 (g, M, b)"""
        node = self.nodes[name]
        g_s, M_f, b_f = self._compact_and_simplify(node.g, node.m, node.M, node.b)
        return g_s, M_f, b_f

    def result_raw(self, name):
        """获取未进一步简化的结果"""
        node = self.nodes[name]
        return SparseANF(node.g, node.m), node.M, node.b

    # ---- 验证 ----

    def verify(self, name, n_tests=200):
        """验证 f(x) = g(Mx⊕b) 正确性（直接求值，避免 expand）"""
        g, M, b = self.result(name)
        M_arr = np.array(M, dtype=np.int64)
        b_arr = np.array(b, dtype=np.int64)

        node = self.nodes[name]
        g_node = SparseANF(node.g, node.m)
        M_node = np.array(node.M, dtype=np.int64)
        b_node = np.array(node.b, dtype=np.int64)

        errors = 0
        for _ in range(n_tests):
            x = np.random.randint(0, 2, self.n)

            # f(x) = g_node(M_node · x ⊕ b_node)
            z_node = (M_node @ x + b_node) % 2
            z_node_mask = sum(int(z_node[j]) << j for j in range(node.m))
            f_val = g_node.eval_mask(z_node_mask)

            # g_result(M_result · x ⊕ b_result)
            z = (M_arr @ x + b_arr) % 2
            z_mask = sum(int(z[j]) << j for j in range(len(M)))

            g_val = g.eval_mask(z_mask)
            if f_val != g_val:
                errors += 1

        return errors


# ===== 输出格式化 =====

def fmt_result(g, M, b, x_names=None, output_name="f", orig_T=None):
    """格式化为可读文本"""
    n = len(M[0]) if M and M[0] else 0
    m = g.n
    x_names = x_names or [f"x{i}" for i in range(n)]
    z_names = [f"z{i}" for i in range(m)]

    lines = []
    lines.append("=" * 70)
    lines.append(f"  简化结果: {output_name}")
    lines.append("=" * 70)
    lines.append("")

    lines.append(f"  n={n}, m={m}")
    if orig_T is not None:
        lines.append(f"  T₀={orig_T} → T(g)={g.T()}")
        if orig_T > 0:
            lines.append(f"  压缩率: {(1-g.T()/orig_T)*100:.1f}%")
    lines.append("")

    lines.append("z = Mx ⊕ b:")
    for i in range(m):
        terms = [x_names[j] for j in range(n) if int(M[i][j]) % 2]
        if int(b[i]) % 2:
            terms.append("1")
        lines.append(f"  {z_names[i]} = {' ⊕ '.join(terms) if terms else '0'}")
    lines.append("")

    # g(z) 表达式
    parts = []
    for mask, coeff in sorted(g.terms.items(), key=lambda x: bin(x[0])):
        if coeff == 0:
            continue
        if mask == 0:
            parts.append("1")
        else:
            monom = [z_names[j] for j in range(m) if mask & (1 << j)]
            parts.append("·".join(monom) if monom else "1")
    g_str = " ⊕ ".join(parts) if parts else "0"
    lines.append(f"g(z) = {g_str}")
    lines.append("")

    return "\n".join(lines)


# ===== 直接 ANF 输入 =====

def simplify_anf(f, threshold=4096):
    """直接对 SparseANF 进行简化"""
    if f.T() > threshold:
        raise ValueError(f"ANF 太大 ({f.T()} > {threshold})，请使用电路输入方式")
    return simplify(f)


# ===== 测试 =====

def demo_priority_encoder():
    """测试 32 位 Priority Encoder（使用线性链，避免指数爆炸）"""
    n = 32
    inputs = [f"i{j}" for j in range(n)]
    outputs = [f"om_{j}" for j in range(n)]

    print("=" * 60)
    print(f"Priority Encoder ({n}输入, {n}输出)")
    print("=" * 60)

    cs = CircuitSimplify(threshold=2048, verbose=True)

    # 使用线性链构建：nₖ = NOT(iₖ) · nₖ₋₁, omₖ = iₖ · nₖ₋₁
    stmts = [f"INORDER = {' '.join(inputs)};",
             f"OUTORDER = {' '.join(outputs)};"]
    stmts.append(f"n0 = !{inputs[0]};")
    stmts.append(f"{outputs[0]} = {inputs[0]};")
    for k in range(1, n):
        stmts.append(f"t{k} = !{inputs[k]}; n{k} = t{k} * n{k-1};")
        stmts.append(f"{outputs[k]} = {inputs[k]} * n{k-1};")

    cs.parse(" ".join(stmts))
    x_names = inputs

    print(f"\n结果:")
    print(f"{'输出':>8} {'T(g)':>6} {'m':>4}  {'z 变换'}")
    for nm in cs.outputs:
        g, M, b = cs.result(nm)
        n_vars = g.n
        non_trivial = sum(1 for i in range(n_vars)
                          if any(int(M[i][j]) % 2 for j in range(n))
                          or int(b[i]) % 2)
        print(f"{nm:>8} {g.T():>6} {n_vars:>4}  "
              f"{non_trivial}/{n_vars} 行非平凡")

    # 详细输出 om_31
    g, M, b = cs.result('om_31')
    txt = fmt_result(g, M, b, x_names, "om_31")
    with open("/home/wangfz/bool/priority_result.txt", "w") as f:
        f.write(txt)
    print(f"\nom_31 结果已输出到 priority_result.txt")
    print(f"stats: {cs.stats['simplify_calls']} simplify 调用")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "demo":
        demo_priority_encoder()
    else:
        print("用法: python circuit_simplify.py demo")
        print("     或导入使用: from circuit_simplify import CircuitSimplify")
