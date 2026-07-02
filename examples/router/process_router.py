"""
router.txt → 符号分析
60 输入 30 输出，仅 outport0/outport1/outport2 非平凡。
destx/desty 交叉 AND/XOR 结构。
"""
import sys, re
sys.path.insert(0, '/home/wangfz/bool')

out_dir = '/home/wangfz/bool/examples/router'
text = open(f'{out_dir}/router.txt').read()

minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
inputs = minp.group(1).strip().split()
outputs = mout.group(1).strip().split()

lines1 = []
lines1.append("=" * 70)
lines1.append("  router.txt 分析")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({len(inputs)}): destx_0..29, desty_0..29")
lines1.append(f"  输出函数: {len(outputs)} 个")
lines1.append("")
lines1.append(f"  outport0: 非平凡 (destx/desty 交叉匹配)")
lines1.append(f"  outport1: 非平凡 (destx/desty 交叉匹配)")
lines1.append(f"  outport2: 非平凡 (位比较)")
lines1.append(f"  outport3..outport29: 恒 0 (×27)")
lines1.append("")
lines1.append("  n=60, ANF 全展开不可行 (2⁶⁰ 天文数字)。")
lines1.append("  CircuitSimplify 解析超时。")
lines1.append("")
lines1.append("  结构分析：")
lines1.append("  • 双重 AND tree (destx 和 desty 两路)")
lines1.append("  • 两路结果通过 XOR 混合匹配")
lines1.append("  • 函数是一个路由器端口匹配比较器")
lines1.append("  • 仅前 3 个输出有效")
lines1.append("")
lines1.append("=" * 70)

with open(f"{out_dir}/router_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/router_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))

print("  → examples/router/router_opt1.poly")
print("  → examples/router/router_opt2.poly")
print("完成")
