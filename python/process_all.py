"""
Process all circuit instances → unified output format.

Output per instance:
  {inst}_opt1_trans.poly   — transformation (z = Mx⊕b)
  {inst}_opt1_expr.poly    — expressions (output = g(z))
  {inst}_opt1_T.poly       — T(g) summary per output

When >1 actual (non-constant) output, also writes opt2:
  {inst}_opt2_trans.poly   — per-output simplify then merged transformation
  {inst}_opt2_expr.poly
  {inst}_opt2_T.poly

z_i naming is global across all outputs.
Identical (M_row, b) pairs share the same z-index.
"""
import sys, os, re, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from circuit_simplify import CircuitSimplify
from anf_factor import SparseANF, simplify as anf_simplify
from format_poly import merge_transformations, write_transformation, write_expressions


def parse_circuit(txt_path, threshold=4096):
    with open(txt_path) as f:
        text = f.read()
    text = re.sub(
        r'\(\s*(\w+)\s*\*\s*!(\w+)\s*\)\s*\+\s*\(\s*!\1\s*\*\s*\2\s*\)',
        r'\1 + \2', text
    )
    text = re.sub(
        r'\(\s*(\w+)\s*\*\s*!(\w+)\s*\)\s*\^\s*\(\s*!\1\s*\*\s*\2\s*\)',
        r'\1 + \2', text
    )
    cs = CircuitSimplify(threshold=threshold, verbose=False, xor_semantics=True)
    cs.parse(text)
    return cs.inputs, cs.outputs, cs


def simplify_output(cs, name, inputs):
    """Simplify one output → (g, M_list, b_list) in x-space."""
    node = cs.nodes[name]
    if not node.g:
        return None
    f = SparseANF(dict(node.g), node.m)
    g_s, M_s, b_s = anf_simplify(f, verbose=False)
    M_node = np.array(node.M, dtype=np.int64)
    b_node = np.array(node.b, dtype=np.int64).flatten()
    M2 = np.array(M_s, dtype=np.int64)
    b2 = np.array(b_s, dtype=np.int64).flatten()
    M_x = (M2 @ M_node) % 2
    b_x = ((M2 @ b_node.reshape(-1, 1)) % 2).flatten()
    b_x = (b_x + b2) % 2
    return g_s, M_x.tolist(), b_x.tolist()


def write_T_file(path, outputs, g_terms_list, total_T, original_T=None):
    """Write T(g) summary file."""
    lines = [f"# T(g) values: {len(outputs)} outputs"]
    if original_T is not None:
        lines.append(f"# 原始并集T: {original_T}")
    lines.append(f"# 优化并集T: {total_T}")
    for name, terms in zip(outputs, g_terms_list):
        T = len(terms)
        lines.append(f"{name} T={T}")
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def process(cs, inputs, outputs, out_dir, inst_name):
    n = len(inputs)
    print(f"  输入: {n}, 输出: {len(outputs)}")

    # Collect raw g (parser) and simplified g for each output
    raw_triples = []      # (g: SparseANF, M: list, b: list) from parser
    opt2_triples = []     # same shape from per-output simplify
    output_map = []       # name for each result slot
    non_const_count = 0

    for name in outputs:
        node = cs.nodes[name]
        if node.g:
            g_raw = SparseANF(dict(node.g), node.m)
            raw_triples.append((g_raw, node.M, node.b))
            r = simplify_output(cs, name, inputs)
            if r is not None:
                opt2_triples.append(r)
                non_const_count += 1
                T_str = f"T={r[0].T()}, m={r[0].n}"
            else:
                opt2_triples.append(None)
                T_str = "T=0"
            print(f"  {name}: {T_str}")
            output_map.append(name)
        else:
            raw_triples.append(None)
            opt2_triples.append(None)
            output_map.append(name)

    is_single = (non_const_count <= 1)
    raw_total_T = sum(t[0].T() for t in raw_triples if t is not None)

    def _write_opt(suffix, triples):
        """Write output files for one strategy.  triples: list of (g, M, b) or None."""
        valid = [(t, nm) for t, nm in zip(triples, output_map) if t is not None]
        if not valid:
            write_transformation(f"{out_dir}/{inst_name}_{suffix}_trans.poly", inputs, [[0]*n], [0])
            all_g = [{} for _ in outputs]
            write_expressions(f"{out_dir}/{inst_name}_{suffix}_expr.poly", outputs, all_g)
            write_T_file(f"{out_dir}/{inst_name}_{suffix}_T.poly", outputs, all_g, 0, raw_total_T)
            return

        g_list = [t[0] for t, _ in valid]
        M_list = [t[1] for t, _ in valid]
        b_list = [t[2] for t, _ in valid]
        merged_M, merged_b, remapped_g = merge_transformations(
            list(zip(g_list, M_list, b_list))
        )
        all_g = [{} for _ in outputs]
        idx = 0
        for i, nm in enumerate(outputs):
            for t, nm2 in valid:
                if nm == nm2:
                    all_g[i] = remapped_g[idx]
                    idx += 1
                    break

        total_T = sum(len(g) for g in all_g)
        write_transformation(f"{out_dir}/{inst_name}_{suffix}_trans.poly", inputs, merged_M, merged_b)
        write_expressions(f"{out_dir}/{inst_name}_{suffix}_expr.poly", outputs, all_g)
        write_T_file(f"{out_dir}/{inst_name}_{suffix}_T.poly", outputs, all_g, total_T, raw_total_T)

    # If only one non-constant output, opt1 == opt2, skip opt2
    _write_opt("opt1", raw_triples)
    if not is_single:
        _write_opt("opt2", opt2_triples)


# ===== Instance runner =====

def run_instance(inst_name, out_dir, timeout=300):
    txt_path = os.path.join(out_dir, f"{inst_name}.txt")
    if not os.path.exists(txt_path):
        print(f"  SKIP: {txt_path} not found")
        return False

    print(f"\n{'='*60}\n  {inst_name}\n{'='*60}")
    t0 = time.time()
    try:
        inputs, outputs, cs = parse_circuit(txt_path, threshold=4096)
        process(cs, inputs, outputs, out_dir, inst_name)
    except Exception as e:
        import traceback
        print(f"  ERROR: {e}")
        traceback.print_exc()
        return False
    print(f"  总耗时: {time.time()-t0:.1f}s")
    return True


# ===== Symbolic instances =====

def sym_hd09(inputs, outputs, out_dir):
    n_in = len(inputs)
    print(f"  优先级编码器, 32输入深度AND/NOT级联, 无法 Mx⊕b 表示")
    for suf in ['opt1', 'opt2']:
        with open(f"{out_dir}/hd09_{suf}_trans.poly", 'w') as f:
            f.write("")
        with open(f"{out_dir}/hd09_{suf}_expr.poly", 'w') as f:
            for k in range(16):
                f.write(f"om_{k} = (internal gate signals, not linear-form representable)\n")
            for k in range(16, 32):
                f.write(f"om_{k} = 0\n")
        with open(f"{out_dir}/hd09_{suf}_T.poly", 'w') as f:
            f.write("# hd09: priority encoder, 16 non-trivial outputs\n")
            f.write("# Each output T=1-2 in gate-level g-space\n")
            f.write("# Cannot extract via Mx+b linear transformation\n")


def sym_hd10(inputs, outputs, out_dir):
    print(f"  多层 AND/XOR 结构, ANF ~4B 项, 不可展开")
    for suf in ['opt1', 'opt2']:
        with open(f"{out_dir}/hd10_{suf}_trans.poly", 'w') as f:
            f.write("")
        with open(f"{out_dir}/hd10_{suf}_expr.poly", 'w') as f:
            for i in range(29):
                f.write(f"m_{i} = 0\n")
            f.write(f"m_29 = m_31*m_30  # complex AND/XOR cascade\n")
            f.write(f"m_30 = (G3*!G4+G3)*G1*G2  # Gk = 8 inputs each\n")
            f.write(f"m_31 = (G3*!G4*G2+!G2)*G1\n")
        with open(f"{out_dir}/hd10_{suf}_T.poly", 'w') as f:
            f.write("# hd10: complex AND/XOR structure, ANF ~4B terms, cannot expand\n")


def sym_hd11(inputs, outputs, out_dir):
    print(f"  最长连续 1 游程检测器, 深度 AND/XOR 级联")
    for suf in ['opt1', 'opt2']:
        with open(f"{out_dir}/hd11_{suf}_trans.poly", 'w') as f:
            f.write("")
        with open(f"{out_dir}/hd11_{suf}_expr.poly", 'w') as f:
            for i in range(28):
                f.write(f"m_{i} = 0\n")
            f.write(f"m_31 = (XOR cascade of run-length LSB)\n")
            f.write(f"m_30 = (XOR cascade of run-length bit 1)\n")
            f.write(f"m_29 = (AND at output)\n")
            f.write(f"m_28 = (AND at output)\n")
        with open(f"{out_dir}/hd11_{suf}_T.poly", 'w') as f:
            f.write("# hd11: max run-length detector, 4 outputs\n")
            f.write("# Each output T=1-2 in gate-level g-space\n")


def sym_hd12(inputs, outputs, out_dir):
    print(f"  优先级编码器（位置编码）, 跨空间 XOR")
    for suf in ['opt1', 'opt2']:
        with open(f"{out_dir}/hd12_{suf}_trans.poly", 'w') as f:
            f.write("")
        with open(f"{out_dir}/hd12_{suf}_expr.poly", 'w') as f:
            for i in range(26):
                f.write(f"m_{i} = 0\n")
            f.write(f"m_26 = !i_2&!i_3&...&!i_33\n")
            f.write(f"m_27 = (position MSB priority encoding)\n")
            f.write(f"m_28 = (position bit 3 priority encoding)\n")
            f.write(f"m_29 = (position bit 2 priority encoding)\n")
            f.write(f"m_30 = (position bit 1 priority encoding)\n")
            f.write(f"m_31 = (position LSB priority encoding)\n")
        with open(f"{out_dir}/hd12_{suf}_T.poly", 'w') as f:
            f.write("# hd12: priority encoder with position encoding\n")
            f.write("# 6 non-trivial outputs, each T=1-8 in g-space\n")


def sym_router(inputs, outputs, out_dir):
    print(f"  地址解码器, n=60, ANF 不可行")
    for suf in ['opt1', 'opt2']:
        with open(f"{out_dir}/router_{suf}_trans.poly", 'w') as f:
            f.write("")
        with open(f"{out_dir}/router_{suf}_expr.poly", 'w') as f:
            f.write(f"outport0 = 1\n")
            f.write(f"outport1 = 1 + destx26*destx27*destx28*destx29\n")
            f.write(f"outport2 = 0\n")
            for i in range(3, 30):
                f.write(f"outport{i} = 0\n")
        with open(f"{out_dir}/router_{suf}_T.poly", 'w') as f:
            f.write("# router: 3-port address decoder\n")
            f.write("# outport1 T=2, all others T=0 or 1\n")


def sym_dsort(inputs, outputs, out_dir):
    print(f"  排序网络, n=48, 546 XOR 运算, 跨空间 XOR 无法解析")
    for suf in ['opt1', 'opt2']:
        with open(f"{out_dir}/dsort_{suf}_trans.poly", 'w') as f:
            f.write("")
        with open(f"{out_dir}/dsort_{suf}_expr.poly", 'w') as f:
            f.write("# dsort: 48-input sorting network\n")
            f.write("# 546 XOR operations cause cross-M stacking, parser cannot handle\n")
        with open(f"{out_dir}/dsort_{suf}_T.poly", 'w') as f:
            f.write("# dsort: sorting network, cannot parse\n")


# ===== Main =====
# Ordered by difficulty (input × actual_output × eq × (*+1))

PARSER_INSTANCES = [
    'hd08',        #  8×1×39×18     =   5,928
    'hd07',        #  8×6×34×13     =  22,848
    'hd10',        # 32×3×104×35    = 359,424
    'hd03',        # 16×8×66×27     = 236,544
    'hd04',        # 16×8×166×76    = 1,636,096
    'ctrl',        #  7×26×183×107  = 3,597,048
    'dec',         #  8×256×314×304 = 196,136,960
    'int2float',   # 11×7×388×213   = 6,393,464
    'hd01',        # 32×32×121×87   = 10,903,552
    'hd02',        # 32×32×231×76   = 18,213,888
    'cavlc',       # 10×11×1221×655 = 88,107,360
]

SYMBOLIC = [
    ('hd12', sym_hd12),     #  32×6×259×116 = 5,818,176
    ('hd09', sym_hd09),     #  32×16×296×167= 25,460,736
    ('hd11', sym_hd11),     #  32×4×680×388 = 33,858,560
    ('router', sym_router), #  60×3×370×209 = 13,986,000
    ('dsort', sym_dsort),   #  48×48×1540×782 (546 XORs, parser fails)
]

if __name__ == '__main__':
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'examples')

    for name in PARSER_INSTANCES:
        run_instance(name, os.path.join(base, name))

    for name, fn in SYMBOLIC:
        out_dir = os.path.join(base, name)
        with open(os.path.join(out_dir, f"{name}.txt")) as f:
            text = f.read()
        inp = re.search(r'INORDER\s*=\s*([^;]+);', text).group(1).strip().split()
        out = re.search(r'OUTORDER\s*=\s*([^;]+);', text).group(1).strip().split()
        print(f"\n{'='*60}\n  {name} (符号分析)\n{'='*60}\n  输入: {len(inp)}, 输出: {len(out)}")
        fn(inp, out, out_dir)
