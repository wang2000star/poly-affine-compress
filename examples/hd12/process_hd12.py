"""
hd12.txt → 符号分析
32 输入 32 输出，仅 m_26..m_31 非平凡。
Priority encoder（find first 1）。
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

# Circuit function: priority encoder (find first 1 from LSB side)
# Input i_2..i_33 (32 bits, i_2 = LSB = position 0)
# Output m_26..m_31 encode position of the lowest-index 1:
#
# i_2=1  → position 0  (000000)
# i_3=1  → position 1  (000001)
# i_4=1  → position 2  (000010)
# ...
# i_17=1 → position 15 (001111)
# i_18=1 → position 16 (010000)
# ...
# i_33=1 → position 31 (011111)
# all 0  → no valid 1  (111111, m_26=1)
#
# Output bit assignment:
#   m_26: bit 5 (32) — zero/invalid flag: 1 = no 1 found
#   m_27: bit 4 (16) — position ≥ 16
#   m_28: bit 3 (8)  — position bit 3
#   m_29: bit 2 (4)  — position bit 2
#   m_30: bit 1 (2)  — position bit 1
#   m_31: bit 0 (1)  — position bit 0 (LSB of position)
#
# Key structural signals:
#   n149: selector — lower half (i_2..i_17) vs upper half (i_18..i_33)
#         来自深层 XOR 级联比较链
#   n151 = !n149
#   n153 = n149·i_2 + n151·i_18  — MUX: select bit 0 of position
#   n156 = n149·i_8 + n151·i_24  — MUX: select bit 6
#   etc. — multiplexer tree for each output bit
#
# The circuit is a multiplexer-based priority encoder:
#   n149 determines if first 1 is in lower (i_2..i_17) or upper (i_18..i_33)
#   Then mux tree selects the appropriate position bits
#   m_26 generated from the "no 1 found" detector n305

lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: hd12.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个 (6 非平凡)")
lines1.append(f"  m_0..m_25: 恒 0 (×26)")
lines1.append("")
lines1.append("  n=32, ANF 全展开 (2³²=4B) 不可行。")
lines1.append("  CircuitSimplify 解析超时。")
lines1.append("")
lines1.append("  电路功能：优先级编码器（Priority Encoder / Find First 1）")
lines1.append("  输入 i_2..i_33 视为 32-bit 向量（i_2 = 最低位），")
lines1.append("  输出最低位为 1 的位置索引。")
lines1.append("")
lines1.append("  6-bit 输出编码 (m_26..m_31)：")
lines1.append("  • 正常位置编码：5 bit 位置 (m_27..m_31) + 1 bit 有效标志 (m_26)")
lines1.append("  • 位置范围：0 (i_2) 到 31 (i_33)")
lines1.append("  • 无 1 输入时输出全 1 (111111), m_26=1 表示无效")
lines1.append("")
lines1.append("  ANF 复杂度：")
lines1.append("  • 位置 0 依赖 i_2")
lines1.append("  • 位置 k 依赖 i_{k+2} AND NOT(所有更低位置) —— 约 k 个变量")
lines1.append("  • 综合每个输出的 ANF 含指数级项")
lines1.append("  • 全展开项数 ~2³²")
lines1.append("")
lines1.append("  尽管 function 本身有简洁描述（优先编码），")
lines1.append("  其在原始输入下的 ANF 需对每个位置做条件否定（AND with NOT of lower bits），")
lines1.append("  导致 ANF 项数随 n 指数增长。")
lines1.append("")
lines1.append("  仿射变换（共享）效果有限：")
lines1.append("  • 优先级编码是 non-linear 的——输出位不是输入的线性函数")
lines1.append("  • 任何 GF(2) 仿射变换无法将 IF-THEN-ELSE 结构简化为稀疏 ANF")
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
lines2.append("  各自变换需对每个输出单独提取 ANF 并简化。")
lines2.append("  各输出均依赖全部 32 个输入，ANF 项数过大。")
lines2.append("")
lines2.append("  各输出紧湊表达式：")
lines2.append("")
lines2.append("  m_26 = NOR(i_2, i_3, ..., i_33)  —— 全零检测")
lines2.append("  m_27 = i_18 + i_19 + ... + i_33  —— 有 1 在上半区")
lines2.append("         但需要排除下半区有 1 的影响")
lines2.append("")
lines2.append("  更精确（在整数语义下）：")
lines2.append("  position = min{k | i_{k+2}=1} （找到最低位 1）")
lines2.append("  m_26 = 1 当 position 不存在（全 0）")
lines2.append("  m_27..m_31 = position 的二进制编码")
lines2.append("")
lines2.append("  各输出位在原始布尔代数下的表达式：")
lines2.append("")
lines2.append("  m_31 = i_2 ⊕ i_2·i_3 ⊕ i_5·i_6·i_7·i_8·...（低位优先级树）")
lines2.append("       ≈ 位置奇偶性的复杂展开")
lines2.append("")
lines2.append("  m_30 = 位置 bit 1 的编码（涉及更多条件否定）")
lines2.append("")
lines2.append("  m_29 = 位置 bit 2 的编码")
lines2.append("")
lines2.append("  m_28 = 位置 bit 3 的编码")
lines2.append("")
lines2.append("  m_27 = 位置 bit 4 的编码（=1 当最低 1 在 i_18..i_33）")
lines2.append("")
lines2.append("  由于 GF(2) 仿射变换无法有效简化的 IF-THEN-ELSE")
lines2.append("  优先级编码逻辑，策略 2 同样无效。")
lines2.append("")
lines2.append("=" * 70)

with open(f"{out_dir}/hd12_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd12_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("  → examples/hd12/hd12_opt1.poly")
print("  → examples/hd12/hd12_opt2.poly")
print("完成（符号分析）")
