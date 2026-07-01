"""
dec.txt 网表 → ANF → 优化
8 输入 256 输出译码器（纯 AND/NOT 网表）
"""
import numpy as np
import time
from circuit_simplify import CircuitSimplify
from anf_factor import SparseANF
from anf_factor import simplify as anf_simplify
from simplify_poly import compute_union_T, _apply_merge_to_all


def extract_anf(node, n):
    M = np.array(node.M, dtype=np.int64)
    b = np.array(node.b, dtype=np.int64)
    g = SparseANF(node.g, node.m)
    tt_len = 1 << n
    tt = np.zeros(tt_len, dtype=np.int64)
    for x_mask in range(tt_len):
        x_arr = np.array([(x_mask >> i) & 1 for i in range(n)], dtype=np.int64)
        z = (M @ x_arr + b) % 2
        z_mask = sum(int(z[j]) << j for j in range(node.m))
        tt[x_mask] = g.eval_mask(z_mask)
    anf_mask = tt.copy()
    for i in range(n):
        step = 1 << i
        for j in range(0, tt_len, step << 1):
            for k in range(step):
                anf_mask[j + k + step] ^= anf_mask[j + k]
    return {mask: 1 for mask, val in enumerate(anf_mask) if val}


def main():
    out_dir = '/home/wangfz/bool'

    print("=" * 60)
    print("dec.txt 网表 → 优化")
    print("=" * 60)

    # ---- Parse netlist ----
    print("\n阶段 1: 解析网表...")
    t0 = time.time()
    cs = CircuitSimplify(threshold=4096, verbose=False)
    cs.parse(open('/home/wangfz/bool/dec.txt').read())
    n, inputs, outputs = cs.n, cs.inputs, cs.outputs
    print(f"  输入: {len(inputs)}, 输出: {len(outputs)}, 用时: {time.time()-t0:.2f}s")

    # ---- Extract ANF ----
    print("\n阶段 2: 提取 ANF...")
    t0 = time.time()
    funcs = {}
    total_T0 = 0
    union_masks = set()
    for name in outputs:
        anf = extract_anf(cs.nodes[name], n)
        funcs[name] = anf
        total_T0 += len(anf)
        union_masks.update(anf.keys())
    print(f"  ΣT₀ = {total_T0}, union T₀ = {len(union_masks)}, 用时: {time.time()-t0:.2f}s")

    # ---- Strategy 2: per-output ----
    print("\n策略 2（各自变换）...")
    t0 = time.time()
    total_T1 = 0
    simplified = {}  # cache
    for name in outputs:
        anf = funcs[name]
        if not anf:
            simplified[name] = None
            continue
        f = SparseANF(dict(anf), n)
        g, M_arr, b_arr = anf_simplify(f, verbose=False)
        simplified[name] = (g, M_arr, b_arr)
        total_T1 += g.T()
    t2 = time.time()
    pct2 = (total_T0 - total_T1) / total_T0 * 100
    print(f"  ΣT: {total_T0} → {total_T1} ({pct2:.1f}%↓), 用时: {t2-t0:.2f}s")

    # Verify strategy 2 (sample 5 outputs)
    errors = 0
    rng = np.random.default_rng(42)
    for name in outputs[:5]:
        item = simplified.get(name)
        if item is None:
            continue
        g, M_arr, b_arr = item
        f = SparseANF(dict(funcs[name]), n)
        m = g.n
        for _ in range(50):
            x = rng.integers(0, 2, n)
            xm = sum(int(x[j]) << j for j in range(n))
            z = (M_arr @ x + b_arr.flatten()) % 2
            zm = sum(int(z[j]) << j for j in range(m))
            if f.eval_mask(xm) != g.eval_mask(zm):
                errors += 1
    print(f"  验证: {errors}/250 错误")

    # ---- Strategy 1: shared transform ----
    print("\n策略 1（共享变换）...")
    func_list = [SparseANF(dict(funcs[name]), n) for name in outputs]
    orig_union = compute_union_T(func_list)

    # Complement search (n=8 → 256)
    t0 = time.time()
    best_comp, best_union = 0, orig_union
    for comp in range(1 << n):
        cand = []
        for f in func_list:
            M = np.eye(n, dtype=np.int64)
            b = np.array([(comp >> j) & 1 for j in range(n)], dtype=np.int64).reshape(-1, 1)
            cand.append(f.substitute_affine(M, b, verify=False))
        u = compute_union_T(cand)
        if u < best_union:
            best_union, best_comp = u, comp
    t1 = time.time()

    # Apply best complement
    M_s = np.eye(n, dtype=np.int64)
    b_s = np.array([(best_comp >> j) & 1 for j in range(n)], dtype=np.int64)
    cur = []
    for f in func_list:
        b_arr = np.array([(best_comp >> j) & 1 for j in range(n)], dtype=np.int64).reshape(-1, 1)
        cur.append(f.substitute_affine(np.eye(n, dtype=np.int64), b_arr, verify=False))

    # Greedy merge
    improved = True
    while improved:
        improved = False
        cur_u = compute_union_T(cur)
        best_i = best_j = -1
        best_mu = cur_u
        active = sorted(set(
            i for f in cur for mask in f.terms
            for i in range(n) if mask & (1 << i)
        ))
        for i in active:
            for j in active:
                if i == j:
                    continue
                cand = _apply_merge_to_all(cur, i, j, n)
                u = compute_union_T(cand)
                if u < best_mu:
                    best_mu, best_i, best_j = u, i, j
        if best_mu < cur_u:
            cur = _apply_merge_to_all(cur, best_i, best_j, n)
            M_s[best_i] = (M_s[best_i] + M_s[best_j]) % 2
            b_s[best_i] = (b_s[best_i] + b_s[best_j]) % 2
            improved = True
    t2 = time.time()
    final_union = compute_union_T(cur)
    pct1 = (orig_union - final_union) / orig_union * 100

    print(f"  union T: {orig_union} → {final_union} ({pct1:.1f}%↓)")
    print(f"  complement: {time.time()-t0:.2f}s + merge: {t2-t1:.2f}s")

    # Verify strategy 1
    errors = 0
    for idx in range(min(10, len(outputs))):
        f_orig = func_list[idx]
        f_new = cur[idx]
        for _ in range(20):
            x = rng.integers(0, 2, n)
            xm = sum(int(x[j]) << j for j in range(n))
            z = (M_s @ x + b_s) % 2
            zm = sum(int(z[j]) << j for j in range(n))
            if f_orig.eval_mask(xm) != f_new.eval_mask(zm):
                errors += 1
    print(f"  验证: {errors}/200 错误")

    # ---- Write output files ----
    with open(f"{out_dir}/dec_opt2.poly", 'w') as f:
        f.write(f"策略 2（各自变换）: dec.poly\n")
        f.write(f"  8 输入, 256 输出\n")
        f.write(f"  ΣT: {total_T0} → {total_T1} ({pct2:.1f}%↓)\n")

    with open(f"{out_dir}/dec_opt1.poly", 'w') as f:
        f.write(f"策略 1（共享变换）: dec.poly\n")
        f.write(f"  8 输入, 256 输出\n")
        f.write(f"  union T: {orig_union} → {final_union} ({pct1:.1f}%↓)\n")

    print(f"\n输出: dec_opt1.poly, dec_opt2.poly")
    print("=" * 60)


if __name__ == '__main__':
    main()
