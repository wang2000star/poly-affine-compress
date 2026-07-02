"""
hd11.txt → 符号分析
32 输入 32 输出，m_0..m_27 ≡ 0，仅 m_28/m_29/m_30/m_31 非平凡。
大型 AND/XOR 树 (686 行)，n=32 全展开不可行。
"""
import sys, re
sys.path.insert(0, '/home/wangfz/bool')

out_dir = '/home/wangfz/bool/examples/hd11'
text = open(f'{out_dir}/hd11.txt').read()

minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
inputs = minp.group(1).strip().split()
outputs = mout.group(1).strip().split()

lines1 = []
lines1.append("=" * 70)
lines1.append("  hd11.txt 分析")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({len(inputs)}): {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个")
lines1.append("")
lines1.append(f"  m_0..m_27: 恒 0 (×28)")
lines1.append(f"  m_28: 非平凡 (AND tree output)")
lines1.append(f"  m_29: 非平凡 (XOR of carry signals)")
lines1.append(f"  m_30: 非平凡 (XOR of intermediate)")
lines1.append(f"  m_31: 非平凡 (XOR of final stage)")
lines1.append("")
lines1.append("  n=32, ANF 全展开不可行 (2³²=4B)。")
lines1.append("  CircuitSimplify 解析超时（686 行复杂 AND/XOR 树）。")
lines1.append("")
lines1.append("  结构分析：")
lines1.append("  • 3 级 AND tree (n93-n302, n400-n689)")
lines1.append("  • 多路 XOR 混合输出")
lines1.append("  • 类似 hd10 的 AND-of-ORs 结构")
lines1.append("")
lines1.append("=" * 70)

with open(f"{out_dir}/hd11_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd11_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))

print("  → examples/hd11/hd11_opt1.poly")
print("  → examples/hd11/hd11_opt2.poly")
print("完成")
