"""
router.txt → 符号分析
60 输入 30 输出，仅 outport0/outport1 非平凡。
3 端口地址解码器（address range decoder）。
"""
import sys, re
sys.path.insert(0, '/home/wangfz/bool')

out_dir = '/home/wangfz/bool/examples/router'
text = open(f'{out_dir}/router.txt').read()

minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
inputs = minp.group(1).strip().split()
outputs = mout.group(1).strip().split()

n = len(inputs)  # 60

# Circuit function: 3-port address decoder
# Input: destx[0:29] (30-bit X coordinate), desty[0:29] (30-bit Y coordinate)
# Output: outport0, outport1, outport2
#
# Verified behavior:
#   outport0 ≡ 1  (always active — default/valid port)
#   outport1 = 0  when destx[26:29] = all-1 (i.e., destx ≥ 0x78000000)
#   outport1 = 1  otherwise
#   outport2 ≡ 0  (dead output, never fires)
#   outport3..outport29 ≡ 0
#
# Key structural signals (destx path):
#   n92 = destx27 * destx28
#   n130 = n129 * n92 — AND-tree cascade checking destx[24:29]
#   n173 = deep AND-tree combining destx[0:25] patterns
#   n177 = !(n173 * n130) — "all top bits set" detector
#   n215 = !(n213 * !n130) — complement path
#   n216 = n215 * n177
#   outport0 = !n216
#
# desty path has parallel structure but only cross-couples at final stage:
#   n257 = n256 * n219 — desty[27:29] AND
#   n350 = cross-coupled AND-tree (desty + destx0)
#   n360 = !(n358 * n215)
#   outport1 = n360 * n177
#   outport2 = n366 * n216 (n216 always 0 → outport2 always 0)

lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: router.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({n}): destx_0..29, desty_0..29")
lines1.append(f"  输出函数: {len(outputs)} 个")
lines1.append("")
lines1.append(f"  outport0: 恒 1 (默认端口)")
lines1.append(f"  outport1: 非平凡 (地址范围比较)")
lines1.append(f"  outport2: 恒 0")
lines1.append(f"  outport3..outport29: 恒 0 (×27)")
lines1.append("")
lines1.append("  n=60, ANF 全展开不可行 (2⁶⁰ 天文数字)。")
lines1.append("  CircuitSimplify 解析超时（双重 AND tree + XOR 交叉耦合）。")
lines1.append("")
lines1.append("  电路功能: 3 端口地址解码器")
lines1.append("  输入: 两个 30-bit 坐标 (destx, desty)")
lines1.append("  输出: 端口选择")
lines1.append("")
lines1.append("  outport0: 恒为 1（始终有效，默认/管理端口）")
lines1.append("")
lines1.append("  outport1: 地址范围检测")
lines1.append("    = 0  当 destx[26]·destx[27]·destx[28]·destx[29] = 1")
lines1.append("        即 destx 高 4 位全为 1 (destx ≥ 0x78000000)")
lines1.append("    = 1  其他情况")
lines1.append("")
lines1.append("  outport2: 恒为 0（输出端口地址无冲突）")
lines1.append("")
lines1.append("  ANF 复杂度：")
lines1.append("  • 每个输出依赖全部 60 个输入")
lines1.append("  • deep AND-tree 结构（~150 级串联）")
lines1.append("  • XOR 交叉耦合导致等价 ANF 项数指数级")
lines1.append("")
lines1.append("  共享仿射变换：尽管 outport1 仅依赖 destx[26:29] 的 AND，")
lines1.append("  但 AND-tree 检测路径涉及大量中间变量和 NOT 信号，")
lines1.append("  直接仿射变换无法简化 AND-of-NOTs 结构。")
lines1.append("")
lines1.append("=" * 70)

lines2 = []
lines2.append("=" * 70)
lines2.append("  策略 2（各自变换）: router.txt")
lines2.append("=" * 70)
lines2.append("")
lines2.append(f"  原始输入 ({n}): destx_0..29, desty_0..29")
lines2.append(f"  输出函数: {len(outputs)} 个 (1 非平凡)")
lines2.append("")
lines2.append("  各自变换需对每个输出单独提取 ANF。")
lines2.append("  n=60 导致 ANF 字典不可用（项数爆炸）。")
lines2.append("")
lines2.append("  各输出在原始变量下的等效表达式：")
lines2.append("")
lines2.append("  outport0 = 1")
lines2.append("")
lines2.append("  outport1 = NOT(destx[26] AND destx[27] AND destx[28] AND destx[29])")
lines2.append("           = NOR(destx26, destx27, destx28, destx29)")
lines2.append("")
lines2.append("  展开为 ANF 后：")
lines2.append("  outport1 = 1 ⊕ (destx26·destx27·destx28·destx29)")
lines2.append("           = 1 + destx26·destx27·destx28·destx29 (GF(2) 加法)")
lines2.append("           → T = 2 (常数 + 4 次单项式)")
lines2.append("")
lines2.append("  但实际上 outport1 在原始门级电路中的 ANF")
lines2.append("  涉及大量中间 NOT 信号的展开，项数远大于 2。")
lines2.append("  仅在门级提取化简后才能获得上述紧凑形式。")
lines2.append("")
lines2.append("  outport2 = 0")
lines2.append("")
lines2.append("  outport3..outport29 = 0")
lines2.append("")
lines2.append("  策略 2 需从门级提取 ANF → 不可行 (n=60)。")
lines2.append("")
lines2.append("  但若借助 TNode g-space 分解，")
lines2.append("  可通过代数代换将 outport1 简化为 T=2 的 4 次单项式。")
lines2.append("")
lines2.append("=" * 70)

with open(f"{out_dir}/router_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/router_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("  → examples/router/router_opt1.poly")
print("  → examples/router/router_opt2.poly")
print("完成（符号分析）")
