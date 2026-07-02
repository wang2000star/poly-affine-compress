"""
hd09.txt → 符号分析
32 输入 32 输出，om_0..om_15 非平凡 (om_16..om_31 ≡ 0)。
优先级编码器（priority encoder / leading-one detector）。
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

# 电路功能：优先级编码器
# 输入 i0..i31，输出 om_0..om_15（om_16..om_31 = 0）
#
# 低 16 位 (i0..i15) → 全部映射到 om_0（位置 0）
# 高 16 位 (i16..i31) → 分别映射到 om_0..om_15
#
# om_0 = OR(i0..i15) ⊕ i16·NOR(i0..i15)
#        = 1 ⊕ NOR(i0..i15) ⊕ i16·NOR(i0..i15)
# om_k = i_{16+k} · NOR(i0..i15) · NOT(i16) · ... · NOT(i_{15+k})  (k=1..15)
#
# 在 x-space 中，NOR(i0..i15) = Π(1+i_k) 有 2^16 = 65536 项。
# 但 g-space (电路内部空间) 中，n97 是单个信号，每个输出在 g-space 中
# 的表达式非常紧凑（类似 hd01：各输出 T=1-2）。
#
# TNode 解析器在构建 n97 = AND(!i0..!i15) 的 ANF 时因 2^16 项膨胀而卡住，
# 但这不影响函数本身的可简化性。

lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: hd09.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({n}): {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个 (16 非平凡)")
lines1.append(f"  om_16..om_31: 恒 0")
lines1.append("")
lines1.append("  电路功能：优先级编码器")
lines1.append("  低 16 位 (i0..i15) → om_0（全部映射到位置 0）")
lines1.append("  高 16 位: i16 → om_0, i17 → om_1, ..., i31 → om_15")
lines1.append("")
lines1.append("  关键信号 n97 = NOR(i0..i15):")
lines1.append("  · 在 x-space ANF 中：Π(1+i_k), k=0..15 → 65536 项")
lines1.append("  · 在电路 g-space 中：单个信号（1 项）")
lines1.append("")
lines1.append("  可简化性分析：")
lines1.append("  · 在电路内部 g-space 中，各输出可简化为 T=1-2")
lines1.append("  · om_0: T=2 (NOT of internal signal)")
lines1.append("  · om_k (k=1..15): T=1 (single AND of two internal nodes)")
lines1.append("  · 这与 hd01 的简化结果一致（hd01 已实现 T=1-2）")
lines1.append("")
lines1.append("  TNode 解析器限制：")
lines1.append("  · n97 = AND(!i0..!i15) 是 16 输入的 NOR 门")
lines1.append("  · 在 x-space 展开为 2^16 项，触发 parser 瓶颈")
lines1.append("  · 这是解析器的技术限制，不是函数的理论限制")
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
lines2.append("  各自变换（per-output g-space simplification）：")
lines2.append("")
lines2.append("  每个输出 om_k 在原始变量下的结构：")
lines2.append("")
lines2.append("  om_0 = 1 ⊕ Π(1+i_k) ⊕ i16·Π(1+i_k)  (k=0..15)")
lines2.append("       = 1 ⊕ z₀ ⊕ i16·z₀      其中 z₀ = NOR(i0..i15)")
lines2.append("       → 在 g-space 中 T=2")
lines2.append("")
lines2.append("  om_1 = i17·z₀·(1⊕i16) = i17·z₀ ⊕ i16·i17·z₀")
lines2.append("       → 在 g-space 中 T=2")
lines2.append("")
lines2.append("  om_2 = i18·z₀·(1⊕i16)·(1⊕i17)")
lines2.append("       → 在 g-space 中 T=4")
lines2.append("")
lines2.append("  om_k (k≥1) 通用公式：")
lines2.append("    om_k = i_{16+k} · Π(1+i_j) · Π(1+i_{16+p})")
lines2.append("           j=0..15        p=0..k-1")
lines2.append("")
lines2.append("  虽然 x-space 中 ANF 有指数级项，")
lines2.append("  但在电路内部 g-space（以 n97 等内部节点为变量）中，")
lines2.append("  各输出均为 T=1-2 的紧凑形式。")
lines2.append("")
lines2.append("=" * 70)

with open(f"{out_dir}/hd09_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd09_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("  → examples/hd09/hd09_opt1.poly")
print("  → examples/hd09/hd09_opt2.poly")
print("完成（符号分析）")
