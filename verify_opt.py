"""
验证优化电路等价性 — 直接求值，不经过 mask 解析。

验证方法：
  orig(x) = 原始电路输出（从 .txt 直接模拟）
  z = Mx⊕b（变换）
  opt(z) = 优化电路输出（从 _expr.poly 直接求值）

验证：orig(x) == opt(Mx⊕b) 对所有 x 成立
"""

import sys, os, re, random

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ========== 原始电路求值 ==========

def make_circuit_func(txt_path):
    """返回 (input_names, output_names, f(x_dict)→out_list)."""
    with open(txt_path) as f:
        text = f.read()
    raw = text.replace('\n', ' ').replace('\r', '')
    raw = re.sub(r'\s+', ' ', raw).strip()

    inp_m = re.search(r'INORDER\s*=\s*([^;]+);', raw)
    out_m = re.search(r'OUTORDER\s*=\s*([^;]+);', raw)
    inputs = inp_m.group(1).split()
    outputs = out_m.group(1).split()

    stmts = []
    body_start = max(raw.find(';', raw.find('INORDER')),
                     raw.find(';', raw.find('OUTORDER')) + 1)
    for stmt in raw[body_start:].split(';'):
        stmt = stmt.strip()
        if not stmt:
            continue
        m = re.match(r'(\w+)\s*=\s*(.+)$', stmt)
        if m:
            stmts.append((m.group(1), m.group(2)))

    def f(x_dict):
        v = dict(x_dict)
        for name, expr in stmts:
            if name in v:
                continue
            if expr == '0':
                v[name] = 0
            elif expr == '1':
                v[name] = 1
            elif expr[0] == '!':
                v[name] = 1 ^ v.get(expr[1:], 0)
            elif '*' in expr:
                ops = re.findall(r'\b(\w+)\b', expr)
                r = 1
                for op in ops:
                    r &= v.get(op, 0)
                v[name] = r
            elif '+' in expr:
                ops = re.findall(r'\b(\w+)\b', expr)
                r = 0
                for op in ops:
                    r ^= v.get(op, 0)
                v[name] = r
            else:
                v[name] = v.get(expr, 0)
        return [v[o] for o in outputs]

    return inputs, outputs, f


# ========== 变换求值 z = Mx⊕b ==========

def make_transform_func(trans_path, inputs):
    """返回 z_func(x_list)→z_list."""
    if not os.path.exists(trans_path) or os.path.getsize(trans_path) == 0:
        return None
    name_to_idx = {n: i for i, n in enumerate(inputs)}
    with open(trans_path) as f:
        lines = f.read().strip().splitlines()
    rows = []  # 每行: [(var_idx, ...), const]
    for line in lines:
        line = line.strip()
        if not line or line[0] == '#':
            continue
        m = re.match(r'z_\d+\s*=\s*(.+)', line)
        if not m:
            continue
        expr = m.group(1).strip()
        if expr == '0':
            rows.append(([], 0))
            continue
        terms = expr.split('+')
        idxs = []
        c = 0
        for t in terms:
            t = t.strip()
            if t == '1':
                c ^= 1
            elif t in name_to_idx:
                idxs.append(name_to_idx[t])
            elif t[0] == '!' and t[1:] in name_to_idx:
                c ^= 1
                idxs.append(name_to_idx[t[1:]])
        rows.append((idxs, c))

    def transform(x):
        return [((sum(x[i] for i in idxs) & 1) ^ c) for idxs, c in rows]
    return transform


# ========== 优化电路求值（预编译 ANF 表达式） ==========

def compile_anf(expr_str):
    """预编译 ANF 表达式为 (term_list, const) 结构。

    expr: 'z_7 + z_0*z_7 + z_1*z_7' 或 '0'
    返回: ([(i1,i2,...), ...], const) 或 (None, 0) 若 expr='0'
    """
    if expr_str == '0':
        return None, 0
    terms = []
    const = 0
    for mon in expr_str.split('+'):
        mon = mon.strip()
        if not mon:
            continue
        if mon == '1':
            const ^= 1
            continue
        idxs = tuple(int(x[2:]) for x in mon.split('*'))
        terms.append(idxs)
    return terms, const


def eval_anf_compiled(compiled, z):
    """求值预编译 ANF。"""
    terms, const = compiled
    if terms is None:
        return const  # always 0 for '0'
    val = const
    for idxs in terms:
        term = 1
        for i in idxs:
            term &= z[i]
        val ^= term
    return val


def make_opt_func(expr_path, outputs):
    """返回 g(z_list)→out_list（预编译 ANF 表达式）。"""
    if not os.path.exists(expr_path):
        return None
    with open(expr_path) as f:
        lines = f.read().strip().splitlines()
    compiled_of = {}
    for line in lines:
        line = line.strip()
        if not line or line[0] == '#':
            continue
        m = re.match(r'(\w+)\s*=\s*(.+)', line)
        if m:
            compiled_of[m.group(1)] = compile_anf(m.group(2).strip())

    def g(z):
        return [eval_anf_compiled(compiled_of.get(o, (None, 0)), z) for o in outputs]
    return g


# ========== 验证 ==========

def verify_instance(inst_name, base_path):
    txt_path = os.path.join(base_path, inst_name, f"{inst_name}.txt")
    if not os.path.exists(txt_path):
        return True, 0, {}

    inputs, outputs, circuit_f = make_circuit_func(txt_path)
    n = len(inputs)
    exhaustive = n <= 16

    ok = True
    total = 0
    results = {}

    for strategy in ['opt1', 'opt2']:
        trans_path = os.path.join(base_path, inst_name, f"{inst_name}_{strategy}_trans.poly")
        expr_path = os.path.join(base_path, inst_name, f"{inst_name}_{strategy}_expr.poly")

        if not os.path.exists(expr_path):
            continue

        transform_f = make_transform_func(trans_path, inputs)
        opt_f = make_opt_func(expr_path, outputs)
        if transform_f is None or opt_f is None:
            print(f"  {inst_name} ({strategy}): 文件缺失，跳过")
            continue

        if exhaustive:
            test_range = range(1 << n)
        else:
            rng = random.Random(f"v_{inst_name}_{strategy}")
            N = min(50000, 1 << n)
            test_range = [rng.randint(0, (1 << n) - 1) for _ in range(N)]

        mismatches = 0
        tested = 0
        for x_bits in test_range:
            x_dict = {inputs[i]: (x_bits >> i) & 1 for i in range(n)}
            x_list = [(x_bits >> i) & 1 for i in range(n)]

            orig = circuit_f(x_dict)
            z = transform_f(x_list)
            opt = opt_f(z)

            for j, o in enumerate(outputs):
                if orig[j] != opt[j]:
                    mismatches += 1
                    if mismatches <= 3:
                        print(f"    MISMATCH {strategy}: x={x_bits}, "
                              f"{o}: orig={orig[j]}, opt={opt[j]}, z={z}")
                    break
            tested += 1

        if mismatches:
            print(f"  ❌ {inst_name} ({strategy}): {mismatches}/{tested} mismatches")
            ok = False
            results[strategy] = False
        else:
            ts = f"全部{1<<n}" if exhaustive else f"{tested}次随机"
            print(f"  ✅ {inst_name} ({strategy}): {ts} 测试通过")
            total += tested
            results[strategy] = True

    return ok, total, results


def main():
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'examples')
    instances = ['hd08', 'hd07', 'hd03', 'hd04', 'ctrl', 'dec',
                 'int2float', 'hd01', 'hd02', 'cavlc']

    print("=" * 60)
    print("  验证优化电路等价性（直接求值）")
    print("=" * 60)

    all_ok = True
    grand_total = 0
    summary = []

    for inst in instances:
        txt_path = os.path.join(base, inst, f"{inst}.txt")
        inputs, outputs, _ = make_circuit_func(txt_path)
        n = len(inputs)
        print(f"\n--- {inst} (n={n}, 输出={len(outputs)}) ---")
        ok, nt, results = verify_instance(inst, base)
        if not ok:
            all_ok = False
        grand_total += nt

        r = []
        for s in ['opt1', 'opt2']:
            if s in results:
                r.append(f"{s}={'✅' if results[s] else '❌'}")
        summary.append(f"  {inst}: {', '.join(r)}")

    print(f"\n{'='*60}")
    if all_ok:
        print(f"  ✅ 全部通过 (共 {grand_total} 个测试)")
    else:
        print(f"  ❌ 有失败的实例")
    print(f"{'='*60}")
    print("\n汇总:")
    for s in summary:
        print(s)

    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())
