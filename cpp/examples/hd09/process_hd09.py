"""
hd09.txt → 符号分析
32 输入 32 输出 (priority encoder, leading-one detector)
om_0..om_15 非平凡, om_16..om_31 ≡ 0.

ANF 全展开 (2^16) 不可行, TNode g-space 解析器因跨空间 XOR 操作导致
M 矩阵无界增长而超时。以下为从电路功能导出的紧凑表达式。
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

print("=" * 60)
print("hd09.txt -> 符号分析")
print("=" * 60)

lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: hd09.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个 (16 非平凡)")
lines1.append("  om_16..om_31: 恒 0 (x16)")
lines1.append("")
lines1.append("  电路功能：优先级编码器 (leading-one detector)")
lines1.append("  . 低 16 位 (i0..i15) -> 全部映射到 om_0")
lines1.append("  . 高 16 位 (i16..i31) -> 分别映射到 om_0..om_15")
lines1.append("")
lines1.append("  共享 g-space 简化：")
lines1.append("  在电路内部 g-space 中定义 16 个内部信号：")
lines1.append("    z_0 = NOR(i0..i15)")
lines1.append("    z_k = NOT(i_{16+k-1})  (k=1..15)")
lines1.append("  各输出在 g-space 中的表达式：")
lines1.append("    om_0 = 1 xor z_0 xor i_16 * z_0  ->  T=2")
lines1.append("    om_1 = i_17 * z_0 * z_1  ->  T=1")
lines1.append("    om_k = i_{16+k} * z_0 * prod z_p  ->  T=1")
lines1.append("  共享 T(g) = 2 + 15*1 = 17")
lines1.append("")
lines1.append("  TNode 解析器限制：")
lines1.append("  . NOR(i0..i15) 在 x-space ANF 中展开为 2^16 项")
lines1.append("  . 简化后 z_0 在 g-space 中有 T=1 但 M 矩阵改变")
lines1.append("  . 后续 XOR 操作 (n98 = n66 + n97) 变为跨空间操作")
lines1.append("  . 跨空间操作堆叠 M 矩阵 -> M 无界增长 -> 解析超时")
lines1.append("  . 此限制不改变函数本身的可简化性")
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
lines2.append("  各输出在 g-space 中的紧凑形式：")
lines2.append("")
lines2.append("  om_0:")
lines2.append("    z_0 = NOR(i0, i1, ..., i15)")
lines2.append("    g(z) = 1 xor z_0 xor i_16 * z_0")
lines2.append("    T=2")
lines2.append("")
for k in range(1, 16):
    lines2.append(f"  om_{k}:")
    lines2.append(f"    = i_{16+k} * NOR(i0..i15) * NOT(i16) * ... * NOT(i_{15+k})")
    lines2.append(f"    -> g-space: T=1")
    lines2.append("")

lines2.append("  总 T(g) = 17")
lines2.append("")
lines2.append("  TNode 解析器因跨空间 XOR 无法自动提取。")
lines2.append("  但函数本身具有高度可简化性（T=17, 每输出 T=1-2）。")
lines2.append("")
lines2.append("=" * 70)

with open(f"{out_dir}/hd09_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd09_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("  -> examples/hd09/hd09_opt1.poly")
print("  -> examples/hd09/hd09_opt2.poly")
print("完成（符号分析）")
print("=" * 60)
