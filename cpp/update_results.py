#!/usr/bin/env python3
"""Update results.md from _stats.txt, .affine, and _verify.txt files.

Usage:
  python3 update_results.py                # read from cpp/results/
  python3 update_results.py --dir <path>   # custom results directory
"""

import os, sys, glob, re

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
PREPROC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "preprocessed")
OUTPUT_FILE = os.path.join(RESULTS_DIR, "results.md")

# Instance metadata (n, k) for instances that may have no data files
INSTANCE_META = {
    "hd08": (8, 1), "hd07": (8, 8), "hd03": (16, 8), "hd04": (16, 8),
    "ctrl": (7, 26), "dec": (8, 256), "int2float": (11, 7), "cavlc": (10, 11),
    "hd10": (32, 32), "hd01": (32, 32), "hd02": (32, 32),
    "hd09": (32, 32), "hd11": (32, 32), "hd12": (32, 32),
}

# Strategy metadata: (direction, variant, display_name, applies_to)
# applies_to: "single" (k=1), "multi" (k>1, n≤16), "n32" (n=32)
STRATEGIES = [
    ("d1a", "opt",  "d1a_opt",  ["single"]),
    ("d3",  "opt",  "d3_opt",   ["single"]),
    ("d1b", "opt",  "d1b_opt",  ["single"]),
    ("d2",  "opt",  "d2_opt",   ["single"]),
    ("d1c", "opt",  "d1c_opt",  ["single"]),
    ("d1a", "opt1", "d1a_opt1", ["multi"]),
    ("d3",  "opt1", "d3_opt1",  ["multi"]),
    ("d1a", "opt2", "d1a_opt2", ["multi", "n32"]),
    ("d1b", "opt2", "d1b_opt2", ["multi", "n32"]),
    ("d2",  "opt2", "d2_opt2",  ["multi"]),
    ("d3",  "opt2", "d3_opt2",  ["multi", "n32"]),
    ("d1c", "opt2", "d1c_opt2", ["multi"]),
]

def read_stats(path):
    """Read a _stats.txt file, return (n, k, sum_T, union_T, max_deg) or None."""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]
    if len(lines) < 5:
        return None
    try:
        return (int(lines[0]), int(lines[1]), int(lines[2]),
                int(lines[3]), int(lines[4]))
    except (ValueError, IndexError):
        return None

def read_affine_m(path):
    """Read first line of .affine to get m (z variables count)."""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        line = f.readline().strip()
    try:
        return int(line)
    except ValueError:
        return None

def read_verify(path):
    """Read verify file, return '✅' if all pass, '❌' if any fail, '—' if missing."""
    if not os.path.isfile(path):
        return "—"
    with open(path) as f:
        content = f.read()
    if "All outputs PASS" in content or "Verified" in content:
        return "✅"
    if "FAIL" in content:
        return "❌"
    return "⏱"  # timeout or other

def get_strat_files(inst_dir, direction, variant):
    """Find the _stats.txt, .affine, _verify.txt for a strategy in an instance directory."""
    pat = f"{direction}_{variant}"
    stats = None
    affine = None
    verify = None
    for f in os.listdir(inst_dir):
        if f.endswith("_stats.txt") and pat in f and "_raw_" not in f:
            # Match exact direction_variant pattern
            # f is like "hd07_d3_opt2_stats.txt" — check the middle part
            base = f.replace("_stats.txt", "")
            parts = base.split("_", 1)  # split off instance name
            if len(parts) == 2 and parts[1] == pat:
                stats = os.path.join(inst_dir, f)
        if f.endswith(".affine") and pat in f:
            base = f.replace(".affine", "")
            parts = base.split("_", 1)
            if len(parts) == 2 and parts[1] == pat:
                affine = os.path.join(inst_dir, f)
        if f.endswith("_verify.txt") and pat in f:
            base = f.replace("_verify.txt", "")
            parts = base.split("_", 1)
            if len(parts) == 2 and parts[1] == pat:
                verify = os.path.join(inst_dir, f)
    return stats, affine, verify

def collect_instances():
    """Collect all instances, their raw stats, and strategy results."""
    instances = []
    inst_dirs = [d for d in INSTANCE_META if d in os.listdir(RESULTS_DIR)]

    for inst in inst_dirs:
        inst_path = os.path.join(RESULTS_DIR, inst)
        preproc_path = os.path.join(PREPROC_DIR, inst, f"{inst}_stats.txt")
        raw = read_stats(preproc_path)
        # For n=32, preprocessed stats are -1 (too large to compute)
        if raw and raw[0] >= 32 and raw[2] < 0:
            raw = None

        # Detect n/k: try raw stats, then any strategy stats, then metadata
        n, k = 0, 0
        if raw:
            n, k = raw[0], raw[1]
        if n == 0:
            for f in os.listdir(inst_path):
                if f.endswith("_stats.txt") and "_raw_" not in f:
                    s = read_stats(os.path.join(inst_path, f))
                    if s and s[0] > 0:
                        n, k = s[0], s[1]
                        break
        if n == 0 and inst in INSTANCE_META:
            n, k = INSTANCE_META[inst]

        # Classify instance
        if k == 1:
            strat_list = [s for s in STRATEGIES if "single" in s[3]]
        elif n >= 32:
            strat_list = [s for s in STRATEGIES if "n32" in s[3]]
        else:
            strat_list = [s for s in STRATEGIES if "multi" in s[3]]

        # Collect strategy data
        strategies = []
        for direction, variant, display, _ in strat_list:
            sp, ap, vp = get_strat_files(inst_path, direction, variant)
            st = read_stats(sp) if sp else None
            m_val = read_affine_m(ap) if ap else None
            v_stat = read_verify(vp) if vp else "—"

            strategies.append({
                "name": display,
                "m": m_val,
                "sum_T": st[2] if st else None,
                "union_T": st[3] if st else None,
                "verify": v_stat,
            })

        # Find best (min union_T, with ✅ > ❌ > —)
        verified_ok = [s for s in strategies if s["verify"] == "✅" and s["union_T"] is not None]
        verified_all = [s for s in strategies if s["union_T"] is not None]
        best = min(verified_ok, key=lambda s: s["union_T"]) if verified_ok else \
               (min(verified_all, key=lambda s: s["union_T"]) if verified_all else None)
        # If best is not verified, still show as "(unverified)"
        best_verified = best and best["verify"] == "✅"

        instances.append({
            "name": inst,
            "n": n,
            "k": k,
            "raw_sum_T": raw[2] if raw else None,
            "raw_union_T": raw[3] if raw else None,
            "strategies": strategies,
            "best": best,
            "best_verified": best_verified,
        })
    return instances

def fmt_raw(val, n):
    """Format value in raw stats table."""
    if val is None:
        return "—" if n >= 32 else "?"
    return str(val)

def fmt_opt(val):
    """Format a metric in optimization table."""
    if val is None:
        return "—"
    return str(val)

def verify_icon(v):
    return v

def compression_ratio(orig, opt):
    if orig is None or opt is None or opt == 0:
        return "—"
    ratio = orig / opt
    if ratio >= 10:
        return f"**{ratio:.1f}×**"
    elif ratio >= 1.01:
        return f"{ratio:.1f}×"
    else:
        return f"{ratio:.1f}×"

def generate_markdown(instances):
    lines = []
    lines.append("# ANF Optimization Results")
    lines.append("")
    lines.append("## Strategy Matrix")
    lines.append("")
    lines.append("| 方向 | opt (单输出) | opt1 (共享变换) | opt2 (各自变换) |")
    lines.append("|------|-------------|----------------|----------------|")
    lines.append("| d1a (真值表搜索) | ✅ 单输出 | ✅ 多输出 | ✅ 多输出 |")
    lines.append("| d1b (门级构建) | ✅ 单输出 | ❌ 不适用 | ✅ 多输出 |")
    lines.append("| d2 (稀疏g代数) | ✅ 单输出 | ❌ 不适用 | ✅ 多输出 |")
    lines.append("| d3 (Walsh+随机+爬山) | ✅ 单输出 | ✅ 多输出 | ✅ 多输出 |")
    lines.append("| d1c (complement-only) | ✅ 单输出 | ❌ 不适用 | ✅ 多输出 |")
    lines.append("")
    lines.append("单输出实例（hd08）：5 种策略 = d1a_opt, d1b_opt, d2_opt, d3_opt, d1c_opt")
    lines.append("多输出实例（其他）：7 种策略 = d1a_opt1, d1a_opt2, d1b_opt2, d2_opt2, d3_opt1, d3_opt2, d1c_opt2")
    lines.append("n=32 实例：仅 d1a_opt2, d1b_opt2, d3_opt2")
    lines.append("")
    lines.append("## Raw ANF Statistics")
    lines.append("")
    lines.append("| 实例 | n | k | 常量数 | 仿射数 | 非线性数 | 原始 sum_T | 原始 union_T |")
    lines.append("|------|---|----|--------|--------|----------|-----------|-------------|")
    for inst in instances:
        if inst["n"] >= 32:
            raw_sum = "—"
            raw_union = "—"
        else:
            raw_sum = fmt_raw(inst["raw_sum_T"], inst["n"])
            raw_union = fmt_raw(inst["raw_union_T"], inst["n"])
        # For n=32, use "—" for raw stats
        lines.append(f"| {inst['name']} | {inst['n']} | {inst['k']} | — | — | — | {raw_sum} | {raw_union} |")
    lines.append("")
    lines.append("注：n=32 实例的 ANF 项数过大，无法直接计算原始 sum_T/union_T（单项式空间 2^32 无法枚举）。")
    lines.append("")
    lines.append("## Optimization Results")
    lines.append("")
    lines.append("✅ = 验证通过  ❌ = 验证失败  — = 未运行/无结果  ⏱ = 超时")
    lines.append("")
    lines.append("sum_T/union_T 仅统计次数 ≥ 2 的单项式（剔除常数项和一次项）。")
    lines.append("")
    lines.append("---")
    lines.append("")

    for inst in instances:
        n = inst["n"]
        raw_ut = inst["raw_union_T"]
        raw_ut_str = fmt_raw(raw_ut, n)

        lines.append(f"### {inst['name']} (n={n}, k={inst['k']}, 原始 union_T={raw_ut_str})")
        lines.append("")
        lines.append("| 策略 | m | sum_T | union_T | 验证 |")
        lines.append("|------|---|-------|---------|------|")
        for s in inst["strategies"]:
            m_str = fmt_opt(s["m"])
            sum_str = fmt_opt(s["sum_T"])
            union_str = fmt_opt(s["union_T"])
            lines.append(f"| {s['name']} | {m_str} | {sum_str} | {union_str} | {s['verify']} |")

        # Best strategy line
        if inst["best"]:
            b = inst["best"]
            best_names = [s["name"] for s in inst["strategies"]
                         if s["verify"] == "✅" and s["union_T"] is not None
                         and s["union_T"] == b["union_T"]]
            if not best_names:
                # include unverified
                best_names = [s["name"] for s in inst["strategies"]
                             if s["union_T"] is not None and s["union_T"] == b["union_T"]]
            best_str = " / ".join(best_names)
            if inst["best_verified"]:
                lines.append(f"| **最优** | **{b['m']}** | **{b['sum_T']}** | **{b['union_T']}** | {best_str} |")
            else:
                lines.append(f"| **最优（未验证）** | **{b['m']}** | **{b['sum_T']}** | **{b['union_T']}** | {best_str} |")
        else:
            lines.append("| **最优** | — | — | — | 无有效结果 |")

        # Compression ratio
        if raw_ut is not None and inst["best"] and inst["best"]["union_T"] is not None:
            cr = compression_ratio(raw_ut, inst["best"]["union_T"])
            if cr != "—":
                lines.append("")
                lines.append(f"压缩比: {raw_ut}→{inst['best']['union_T']} = {cr}")

        lines.append("")
        lines.append("---")
        lines.append("")

    # Summary table
    lines.append("## Summary")
    lines.append("")
    lines.append("| 实例 | n | k | 原始 union_T | 最优 union_T | 最优 m | 压缩比 | 最佳策略 |")
    lines.append("|------|---|----|-------------|-------------|--------|--------|---------|")
    for inst in instances:
        raw_ut_str = fmt_raw(inst["raw_union_T"], inst["n"])
        if inst["best"]:
            best_ut = str(inst["best"]["union_T"])
            best_m = str(inst["best"]["m"]) if inst["best"]["m"] is not None else "—"
            if inst["best_verified"]:
                best_str = " / ".join(s["name"] for s in inst["strategies"]
                            if s["verify"] == "✅" and s["union_T"] is not None
                            and s["union_T"] == inst["best"]["union_T"])
            else:
                best_str_list = [s["name"] for s in inst["strategies"]
                               if s["union_T"] is not None
                               and s["union_T"] == inst["best"]["union_T"]]
                best_str = " / ".join(best_str_list) + " (未验证)"
            if inst["raw_union_T"] is not None and inst["best"]["union_T"] is not None:
                cr = compression_ratio(inst["raw_union_T"], inst["best"]["union_T"])
            else:
                cr = "—"
        else:
            best_ut = "—"
            best_m = "—"
            cr = "—"
            best_str = "无结果"
        lines.append(f"| {inst['name']} | {inst['n']} | {inst['k']} | {raw_ut_str} | **{best_ut}** | {best_m} | {cr} | {best_str} |")

    lines.append("")
    lines.append("## Key Observations")
    lines.append("")
    lines.append("1. **共享变换 (opt1) 最有效**: d1a_opt1 / d3_opt1 在 hd03(14×)、int2float(11.4×)、cavlc(4.3×) 取得最佳压缩比，且 m=n（z 空间维度与 x 空间相同）。")
    lines.append("2. **opt2 的 union_T = sum_T**: 各输出使用独立 z 变量时，不同输出的单项式在共享 z 空间中不重叠，union_T 自然等于 sum_T。")
    lines.append("3. **d1a_opt2 验证失败**: ctrl(61)、dec(57)、cavlc(864) 的 d1a_opt2 结果验证未通过，说明 main_large.cpp 的 z-space 评估可能存在近似误差。")
    lines.append("4. **n=32 难题**: 所有 n=32 实例均未找到通过验证的有效变换——d1a_opt2 验证失败、d1b_opt2 有 gate builder bug、d3_opt2 超时。")
    lines.append("5. **hd08**: d1a_opt/d3_opt/d1c_opt 都从 120→8（15×），d1b_opt/d2_opt 保持 120（无压缩）。hd08 只有一个非线性输出且次数很高（deg=7），因此压缩有限。")

    return "\n".join(lines) + "\n"

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Update results.md from stats files")
    parser.add_argument("--dir", help="Results directory (default: cpp/results/)")
    args = parser.parse_args()

    global RESULTS_DIR, PREPROC_DIR, OUTPUT_FILE
    if args.dir:
        RESULTS_DIR = args.dir
        OUTPUT_FILE = os.path.join(RESULTS_DIR, "results.md")

    if not os.path.isdir(RESULTS_DIR):
        print(f"ERROR: results directory not found: {RESULTS_DIR}", file=sys.stderr)
        sys.exit(1)

    instances = collect_instances()
    if not instances:
        print("No instance data found.", file=sys.stderr)
        sys.exit(1)

    md = generate_markdown(instances)
    with open(OUTPUT_FILE, "w") as f:
        f.write(md)
    print(f"Updated: {OUTPUT_FILE}")
    print(f"  Instances: {len(instances)}")
    total_strats = sum(len(inst["strategies"]) for inst in instances)
    print(f"  Strategies: {total_strats}")

if __name__ == "__main__":
    main()
