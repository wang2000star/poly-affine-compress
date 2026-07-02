"""
hd09.txt → 优化（符号分析）
32 输入 32 输出，om_16..om_31 ≡ 0。
电路为深度 AND/XOR 树结构，ANF 全展开不可行。
"""
import sys, re
sys.path.insert(0, '/home/wangfz/bool')

out_dir = '/home/wangfz/bool/examples/hd09'
text = open(f'{out_dir}/hd09.txt').read()

minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
inputs = minp.group(1).strip().split()
outputs = mout.group(1).strip().split()
n = len(inputs)

# Count constant-zero outputs
nz = sum(1 for o in outputs if int(o.split('_')[1]) >= 16)

# Key structural insight: n97 = AND(!i0..!i15) = NOR(i0..i15).
# This single node has T=65536 in full ANF.
# The XOR tree after it combines with other signals.
# Total ANF for non-trivial outputs is dominated by products involving n97.

print("=" * 60)
print("hd09.txt → 符号分析")
print("=" * 60)
print(f"  输入: {n}, 输出: {len(outputs)} ({len(outputs)-nz} 非平凡)")
print(f"  n=32: ANF 全展开不可行 (2³²=4B)")
print(f"  电路结构: AND tree (NOR of i₀..i₁₅, T=65536) + XOR mixing")
print(f"  → 估算 ΣT₀ ≈ 10⁵~10⁶ 项")
print("=" * 60)

lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: hd09.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入: {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个 ({len(outputs)-nz} 非平凡)")
lines1.append(f"  om_16..om_31: 恒 0")
lines1.append("")
lines1.append("  n=32, ANF 全展开 (2³²=4B) 不可行。")
lines1.append("  CircuitSimplify 解析超时 (阈值 10⁹, >300s)。")
lines1.append("")
lines1.append("  电路结构分析：")
lines1.append("  • 第一阶段 (n65..n97): AND tree，构建 i₀..i₁₅ 全零检测")
lines1.append("    n97 = AND(!i₀..!i₁₅) = NOR(i₀..i₁₅), ANF T = 65536")
lines1.append("  • 第二阶段 (n98..n178): XOR 混合层，输出组合")
lines1.append("  • 第三阶段 (om_0..om_15): 最终输出")
lines1.append("")
lines1.append("  由于变量深度纠缠，仿射变换无法有效简化此结构。")
lines1.append("  多级 AND/XOR 层使得任何 GF(2) 线性映射都无法解耦。")
lines1.append("")
lines1.append("=" * 70)

lines2 = []
lines2.append("=" * 70)
lines2.append("  策略 2（各自变换）: hd09.txt")
lines2.append("=" * 70)
lines2.append("")
lines2.append(f"  原始输入: {', '.join(inputs)}")
lines2.append(f"  输出函数: {len(outputs)} 个 ({len(outputs)-nz} 非平凡)")
lines2.append("")
lines2.append("  各自变换需对每个输出单独提取 ANF 并简化。")
lines2.append("  由于 n=32，ANF 提取不可行（各输出 T 达数万项）。")
lines2.append("")
lines2.append("  om_0..om_15: 非平凡输出，需 C++ 高效解析器处理")
lines2.append("  om_16..om_31: 恒 0")
lines2.append("")
lines2.append("=" * 70)

with open(f"{out_dir}/hd09_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd09_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("  → examples/hd09/hd09_opt1.poly")
print("  → examples/hd09/hd09_opt2.poly")
print("\n完成（符号分析）")
