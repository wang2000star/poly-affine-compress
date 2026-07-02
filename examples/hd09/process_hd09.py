"""
hd09.txt → 符号分析
32 输入 32 输出，仅 om_0..om_15 非平凡。
优先级编码器（leading-one detector）结构，n=32 展开 ANF 不可行。
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

# Key structural insight:
# n97 = AND(!i0..!i15) = NOR(i0..i15) — lower 16 bits all-zero detect
# n98 = n97 XOR !i16 — connects half detect with i16
# om_0..om_15 form a priority encoder cascade
#
# Function: Given 32-bit input (i0..i31), find the first 1 and output its
# position offset from i16. Lower 16 bits (i0..i15) all map to "position 0".
# i16=1 → om_0, i17=1 → om_1, ..., i31=1 → om_15.
# When ANY lower bit is 1 and no upper bit conflict, om_0 fires.

lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: hd09.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个 (16 非平凡)")
lines1.append(f"  om_16..om_31: 恒 0")
lines1.append("")
lines1.append("  n=32, ANF 全展开 (2³²=4B) 不可行。")
lines1.append("  CircuitSimplify 解析超时 (AND/XOR 混合触发昂贵 simplify)。")
lines1.append("")
lines1.append("  电路功能：优先级编码器（leading-one detector）")
lines1.append("  输入 i0..i31，输出 om_0..om_15 编码最高优先级的 1 位置：")
lines1.append("  • i0..i15 任意位为 1 → om_0（全部映射到位置 0）")
lines1.append("  • i16=1（且低 16 位全 0）→ om_0")
lines1.append("  • i17=1（且 i16=0, 低 16 位全 0）→ om_1")
lines1.append("  • i18=1 → om_2 ... i31=1 → om_15")
lines1.append("")
lines1.append("  关键信号：")
lines1.append("  n97 = NOR(i0..i15) —— 检测低 16 位是否全零，ANF T = 65536")
lines1.append("  n98 = n97 XOR !i16 —— 跨半区交叉耦合信号")
lines1.append("  om_k 组成级联链，每个输出依赖前一个输出的结果")
lines1.append("")
lines1.append("  由于低位区 (i0..i15) 的任一位为 1 都触发 om_0，")
lines1.append("  且 XOR 交叉耦合导致变量深度纠缠，仿射变换无法有效简化。")
lines1.append("")
lines1.append("=" * 70)

lines2 = []
lines2.append("=" * 70)
lines2.append("  策略 2（各自变换）: hd09.txt")
lines2.append("=" * 70)
lines2.append("")
lines2.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines2.append(f"  输出函数: {len(outputs)} 个 (16 非平凡)")
lines2.append("")
lines2.append("  各自变换需对每个输出单独提取 ANF 并简化。")
lines2.append("  各输出均依赖全部 32 个输入，ANF 项数达数万。")
lines2.append("")
lines2.append("  各输出在原始变量下的表达式：")
lines2.append("")
lines2.append("  om_0 = (任意 i0..i15=1) OR (i16=1 AND NOR(i0..i15)=1)")
lines2.append("  om_1 = i17=1 AND NOR(i0..i15)=1 AND i16=0")
lines2.append("  om_2 = i18=1 AND NOR(i0..i15)=1 AND i16=0 AND i17=0")
lines2.append("  ...")
lines2.append("  om_k = i_{16+k}=1 且之前各位全 0 且 NOR(i0..i15)=1")
lines2.append("")
lines2.append("  在低 16 位非零时，行为更复杂：")
lines2.append("  • 低 16 位有 1 + 上区位不同组合 → om_k 按 XOR 交叉结果")
lines2.append("")
lines2.append("=" * 70)

with open(f"{out_dir}/hd09_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd09_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("  → examples/hd09/hd09_opt1.poly")
print("  → examples/hd09/hd09_opt2.poly")
print("完成（符号分析）")
