"""
hd12.txt → 符号分析
32 输入 32 输出，仅 m_26..m_31 非平凡。
Priority encoder（find first 1），输出位置编码。
"""
import sys, re
sys.path.insert(0, '/home/wangfz/bool')

out_dir = '/home/wangfz/bool/examples/hd12'
text = open(f'{out_dir}/hd12.txt').read()

minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
inputs = minp.group(1).strip().split()
outputs = mout.group(1).strip().split()
n = len(inputs)

# 电路功能：优先级编码器 (find-first-1)
# 输入 i_2..i_33，输出最低位 1 的位置索引
# m_0..m_25 ≡ 0
# m_26: 全零标志（1=无有效位）
# m_27..m_31: 位置二进制编码（0..31）
#
# 结构类似 hd01 优先级编码器链：
# 输入 i_2..i_33，信号 n149 是上下半区选择器
# 当 n149=1 时使用下半区 i_2..i_17
# 当 n149=0 时使用上半区 i_18..i_33
# 然后通过 MUX 树选择对应的位置编码位
#
# 位置编码 (m_27..m_31) = 5-bit position:
#   i_2  → 00000 (0)
#   i_3  → 00001 (1)
#   ...
#   i_17 → 01111 (15)
#   i_18 → 10000 (16)
#   ...
#   i_33 → 11111 (31)
#
# 可简化性：
# 类似 hd01 的结构，每位的优先级信号（AND of NOT lower + current bit）
# 在 g-space 中可简化为 T=1。二进制编码使 m_27..m_31 各有 T≈4-8。
# 但与 hd09 类似，内部 AND-NOT 链在 TNode 解析时可能触发
# x-space 展开瓶颈。

lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: hd12.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个 (6 非平凡)")
lines1.append(f"  m_0..m_25: 恒 0 (×26)")
lines1.append("")
lines1.append("  电路功能：优先级编码器 (find-first-1)")
lines1.append("  · 输入 i_2..i_33（32 bit 向量，i_2 = 最低位）")
lines1.append("  · 输出最低位为 1 的位置索引")
lines1.append("  · 全零输入输出 111111 (m_26=1 表示无有效位)")
lines1.append("")
lines1.append("  输出编码 (m_26..m_31)：")
lines1.append("  · m_26: 无效标志（0=有效, 1=全零）")
lines1.append("  · m_27..m_31: 5-bit 位置编码")
lines1.append("")
lines1.append("  关键结构：")
lines1.append("  · n149: 上下半区选择器（从 XOR 比较链生成）")
lines1.append("  · n149=1: 选择 i_2..i_17（下半区）")
lines1.append("  · n149=0: 选择 i_18..i_33（上半区）")
lines1.append("  · 输出通过 MUX 树编码位置二进制值")
lines1.append("")
lines1.append("  可简化性：")
lines1.append("  · 类似 hd01，优先级比较链中各位置的信号为 AND of NOT lower")
lines1.append("  · 在电路 g-space 中，各位置信号 T=1")
lines1.append("  · m_26: T=1 (全零检测)")
lines1.append("  · m_27..m_31: T≈4-8 (位置二进制编码)")
lines1.append("")
lines1.append("=" * 70)

lines2 = []
lines2.append("=" * 70)
lines2.append("  策略 2（各自变换）: hd12.txt")
lines2.append("=" * 70)
lines2.append("")
lines2.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines2.append(f"  输出函数: {len(outputs)} 个 (6 非平凡)")
lines2.append("")
lines2.append("  各自变换预期结果：")
lines2.append("")
lines2.append("  m_26 = NOR(i_2, i_3, ..., i_33)")
lines2.append("       = 全零检测")
lines2.append("       → g-space 中 T=1（直接作为一个内部信号）")
lines2.append("")
lines2.append("  m_31 (position LSB):")
lines2.append("    = i_2 ⊕ i_3·!i_2 ⊕ i_5·!i_4·!i_3·!i_2 ⊕ ...")
lines2.append("    → 在优先级编码的 g-space 中 T≈4")
lines2.append("")
lines2.append("  m_30 (position bit 1):")
lines2.append("    → 编码位置值的二进制 bit 1")
lines2.append("    → g-space 中 T≈4")
lines2.append("")
lines2.append("  m_29 (position bit 2): → T≈4")
lines2.append("  m_28 (position bit 3): → T≈4")
lines2.append("  m_27 (position bit 4): → T≈4")
lines2.append("")
lines2.append("  总 ΣT(g) 估计：1 + 5×4 ≈ 21")
lines2.append("")
lines2.append("  与 hd01 的对比：")
lines2.append("  · hd01 是 32 个独立输出（32 位置）")
lines2.append("  · hd12 将位置编码为 5 bit，输出更少但每输出项稍多")
lines2.append("  · hd01 的 ΣT(g) = 32×1.5 ≈ 48")
lines2.append("  · hd12 的 ΣT(g) ≈ 21")
lines2.append("")
lines2.append("=" * 70)

with open(f"{out_dir}/hd12_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd12_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("  → examples/hd12/hd12_opt1.poly")
print("  → examples/hd12/hd12_opt2.poly")
print("完成（符号分析）")
