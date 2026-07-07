#!/usr/bin/env python3
"""
poly2circuit.py — 从 .poly (ANF) 逆向生成 .txt 和 .eqn 电路文件。

.poly 格式：
  m               # 变量数
  k               # 输出数
  t0 t1 ... tk-1  # 每输出单项式数
  [b0, b1, ..., bm-1, coeff]  # 每输出所有单项式

输出：
  {inst}.txt — ISCAS 格式（AND=*, XOR=+, NOT=!）
  {inst}.eqn — 相同格式，输入为 x_0..x_{m-1}，输出为 y_0..y_{k-1}
  {inst}_stats.txt — 5 行数值（n k sum_T union_T max_deg）

用法：
  python3 scripts/poly2circuit.py preprocessed/aes_bool/aes_bool.poly
"""

import sys
import os


def parse_poly(path: str) -> tuple:
    """返回 (m, k, terms_per_output), terms[oi] = [mask0, mask1, ...]"""
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]

    m = int(lines[0])
    k = int(lines[1])
    term_counts = [int(x) for x in lines[2].split()]

    terms = [[] for _ in range(k)]
    idx = 3
    for oi in range(k):
        for _ in range(term_counts[oi]):
            raw = lines[idx].strip("[]")
            idx += 1
            bits = [int(x.strip()) for x in raw.split(",")]
            coeff = bits[-1]
            if coeff == 0:
                continue
            mask = sum(bits[j] << j for j in range(m))
            terms[oi].append(mask)
    return m, k, terms


def write_circuit(m: int, k: int, terms: list, path: str) -> None:
    """写 .txt/.eqn（同格式：*=AND, +=XOR, !=NOT）"""
    with open(path, "w") as f:
        f.write(f"INORDER = {' '.join(f'x_{j}' for j in range(m))};\n")
        f.write(f"OUTORDER = {' '.join(f'y_{j}' for j in range(k))};\n")

        gate_idx = 0

        def new_gate(expr: str) -> str:
            nonlocal gate_idx
            name = f"t_{gate_idx}"
            gate_idx += 1
            f.write(f"  {name} = {expr};\n")
            return name

        for oi in range(k):
            sigs = []
            has_const = False

            for mask in terms[oi]:
                if mask == 0:
                    has_const = True
                    continue

                vars_in = [f"x_{j}" for j in range(m) if (mask >> j) & 1]

                if len(vars_in) == 0:
                    continue
                elif len(vars_in) == 1:
                    sig = vars_in[0]
                else:
                    # 拆成 2-input AND 链（解析器不支持多输入 AND）
                    sig = vars_in[0]
                    for v in vars_in[1:]:
                        sig = new_gate(f"{sig} * {v}")
                sigs.append(sig)

            # XOR chain
            if not sigs:
                out = "0"
            elif len(sigs) == 1:
                out = sigs[0]
            else:
                acc = sigs[0]
                for s in sigs[1:]:
                    acc = new_gate(f"{acc} + {s}")
                out = acc

            # 常数 1 补码
            if has_const:
                if out == "0":
                    out = "1"
                else:
                    out = new_gate(f"!{out}")

            f.write(f"  y_{oi} = {out};\n")


def write_stats(m: int, k: int, terms: list, path: str) -> None:
    """写 _stats.txt：n k sum_T union_T max_deg"""
    sum_T = 0
    all_masks = set()
    max_deg = 0
    for oi in range(k):
        for mask in terms[oi]:
            deg = mask.bit_count()  # Python 3.8+: int.bit_count()
            if deg >= 2:
                sum_T += 1
                all_masks.add(mask)
                if deg > max_deg:
                    max_deg = deg
    union_T = len(all_masks)
    with open(path, "w") as f:
        f.write(f"{m}\n{k}\n{sum_T}\n{union_T}\n{max_deg}\n")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.poly> [--out-dir DIR]")
        sys.exit(1)

    poly_path = sys.argv[1]
    out_dir = None
    for i in range(2, len(sys.argv)):
        if sys.argv[i] == "--out-dir" and i + 1 < len(sys.argv):
            out_dir = sys.argv[i + 1]

    m, k, terms = parse_poly(poly_path)
    inst = os.path.splitext(os.path.basename(poly_path))[0]

    print(f"Parsed: {poly_path}")
    print(f"  n={m}, k={k}")
    for oi in range(k):
        const = " [1]" if 0 in terms[oi] else ""
        print(f"  Output {oi}: {len(terms[oi])} terms{const}")

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        txt_path = os.path.join(out_dir, f"{inst}.txt")
        eqn_path = os.path.join(out_dir, f"{inst}.eqn")
        stats_path = os.path.join(out_dir, f"{inst}_stats.txt")
    else:
        txt_path = f"{inst}.txt"
        eqn_path = f"{inst}.eqn"
        stats_path = f"{inst}_stats.txt"

    write_circuit(m, k, terms, txt_path)
    print(f"  Written: {txt_path}")
    write_circuit(m, k, terms, eqn_path)
    print(f"  Written: {eqn_path}")
    write_stats(m, k, terms, stats_path)
    print(f"  Written: {stats_path}")

    # 总计
    total_terms = sum(len(t) for t in terms)
    total_gates = sum(max(0, len(t) - (1 if 0 in t else 0) - 1) for t in terms)
    print(f"\n  Total ANF terms: {total_terms}")
    print(f"  Circuit gates (AND+XOR+NOT): ~{total_terms * 2}")


if __name__ == "__main__":
    main()
