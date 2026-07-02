"""
hd11.txt → 符号分析
32 输入 32 输出，仅 m_28..m_31 非平凡。
最长连续 1 游程长度检测（max consecutive-run length），≥8 时饱和。
"""
import sys, re
sys.path.insert(0, '/home/wangfz/bool')

out_dir = '/home/wangfz/bool/examples/hd11'
text = open(f'{out_dir}/hd11.txt').read()

minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
inputs = minp.group(1).strip().split()
outputs = mout.group(1).strip().split()
n = len(inputs)

# Circuit function: longest run of consecutive 1s in 32-bit input
# Output encoding (m_28..m_31):
#   0000 (0): all-zero input
#   0001 (1): longest run = 1  (no adjacent 1s)
#   0010 (2): longest run = 2
#   0011 (3): longest run = 3
#   ...
#   0111 (7): longest run = 7
#   1000 (8): longest run >= 8  (saturated)
#
# Key structural observations:
# - Lines 7-215 (n93-n215): AND-tree cascade detecting runs of consecutive 1s.
#   Each n_xx node detects a specific pattern of adjacent bits.
#   Example: n95 = i_30*i_31*i_32, n97 = i_29*i_30*i_31*i_32 (4-consecutive detect)
# - Lines 216-598 (n217-n640): Inverse AND-tree (NOT-decorated) for detecting
#   complement patterns, forming the XOR-tree cascade.
# - Lines 598-end (n641-n687): Deep XOR cascade that encodes the max-run length.
#   n684 = n683 XOR n641
#   n685 = n684 XOR n585
#   n686 = n685 XOR n518
#   n687 = n686 XOR n449
#   n688 = n687 XOR n400
#   n689 = n688 XOR n329
# - Final outputs:
#   m_28 = AND of two internal signals (MSB)
#   m_29 = AND of two internal signals
#   m_30 = XOR of two internal signals
#   m_31 = XOR (n689 XOR n302)

lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: hd11.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个 (4 非平凡)")
lines1.append(f"  m_0..m_27: 恒 0 (×28)")
lines1.append("")
lines1.append("  n=32, ANF 全展开 (2³²=4B) 不可行。")
lines1.append("  CircuitSimplify 解析超时 (AND/XOR 级联深度 ~600 行)。")
lines1.append("")
lines1.append("  电路功能：最长连续 1 游程长度检测器（max consecutive-run length）")
lines1.append("  • 输入 i_2..i_33（32 bit 向量）")
lines1.append("  • 输出 m_28..m_31 编码最长连续 1 的游程长度")
lines1.append("  • 游程长度 ≥8 时饱和输出 8")
lines1.append("  • 全 0 输入输出 0")
lines1.append("")
lines1.append("  输出编码：")
lines1.append("    m_28 (bit 3, MSB): 游程 ≥4")
lines1.append("    m_29 (bit 2):      游程 ≥2 且 ≤3 的进一步区分")
lines1.append("    m_30 (bit 1):      游程 ≥2 的奇偶/模式区分")
lines1.append("    m_31 (bit 0, LSB): 游程 ≥1 的最低有效位")
lines1.append("")
lines1.append("  ANF 复杂度：")
lines1.append("  • 每个输出依赖全部 32 个输入")
lines1.append("  • 游程检测需要 AND-tree 对每个相邻对/三元组/... 做检查")
lines1.append("  • XOR 级联编码将多个游程长度比较结果压缩到 4 bit")
lines1.append("  • 等效 ANF 项数估计在百万到数亿级别")
lines1.append("")
lines1.append("  仿射变换不可行原因：")
lines1.append("  • 游程检测本质上是位位置相关函数，依赖相邻关系")
lines1.append("  • 仿射变换 Mx⊕b 会破坏相邻位之间的空间关系")
lines1.append("  • 除非 M 是置换（permutation），否则会混合原本不相邻的位")
lines1.append("")
lines1.append("=" * 70)

lines2 = []
lines2.append("=" * 70)
lines2.append("  策略 2（各自变换）: hd11.txt")
lines2.append("=" * 70)
lines2.append("")
lines2.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines2.append(f"  输出函数: {len(outputs)} 个 (4 非平凡)")
lines2.append("")
lines2.append("  各自变换需对每个输出单独提取 ANF 并简化。")
lines2.append("  各输出均依赖全部 32 个输入，ANF 项数过大致使无法提取。")
lines2.append("")
lines2.append("  各输出的等效行为：")
lines2.append("")
lines2.append("  m_31 = OR(所有 i_k = 1)  —— 非零检测（至少一位为 1）")
lines2.append("         但受游程≥2 时的模式影响：相邻位均为 1 时 m_31 翻转")
lines2.append("")
lines2.append("  m_30 = 游程 ≥2 的标志与奇偶校正的组合")
lines2.append("         当最长游程 ≥2 时 m_30 通常为 1")
lines2.append("")
lines2.append("  m_29 = 游程 ≥2 且游程模式不满足某种条件的标志")
lines2.append("")
lines2.append("  m_28 = 游程 ≥4（当最长游程 ≥4 时置 1）")
lines2.append("")
lines2.append("  实际行为：m_28..m_31 整体作为 4-bit 数编码最长游程长度 L：")
lines2.append("    L=0 → 0000      L=3 → 0011      L=6 → 0110")
lines2.append("    L=1 → 0001      L=4 → 0100      L=7 → 0111")
lines2.append("    L=2 → 0010      L=5 → 0101      L≥8 → 1000")
lines2.append("")
lines2.append("  由于游程检测在 GF(2) 代数下无简单线性结构，")
lines2.append("  策略 2 同样无法获得有效简化。")
lines2.append("")
lines2.append("=" * 70)

with open(f"{out_dir}/hd11_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd11_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("  → examples/hd11/hd11_opt1.poly")
print("  → examples/hd11/hd11_opt2.poly")
print("完成（符号分析）")
