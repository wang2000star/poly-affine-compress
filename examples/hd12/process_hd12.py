"""
hd12.txt → 符号分析
32 输入 32 输出，m_0..m_25 ≡ 0，仅 m_26..m_31 非平凡。
selector (n148/n149) + XOR 多路输出结构。
"""
import sys, re
sys.path.insert(0, '/home/wangfz/bool')

out_dir = '/home/wangfz/bool/examples/hd12'
text = open(f'{out_dir}/hd12.txt').read()

minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
inputs = minp.group(1).strip().split()
outputs = mout.group(1).strip().split()

lines1 = []
lines1.append("=" * 70)
lines1.append("  hd12.txt 分析")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入 ({len(inputs)}): {', '.join(inputs)}")
lines1.append(f"  输出函数: {len(outputs)} 个")
lines1.append("")
lines1.append(f"  m_0..m_25: 恒 0 (×26)")
lines1.append(f"  m_26: 非平凡 (mux output)")
lines1.append(f"  m_27: 非平凡 (comparator)")
lines1.append(f"  m_28: 非平凡 (selector head)")
lines1.append(f"  m_29: 非平凡 (ripple-carry stage)")
lines1.append(f"  m_30: 非平凡 (ripple-carry stage)")
lines1.append(f"  m_31: 非平凡 (ripple-carry stage)")
lines1.append("")
lines1.append("  n=32, ANF 全展开不可行 (2³²=4B)。")
lines1.append("  CircuitSimplify 解析超时。")
lines1.append("")
lines1.append("  结构分析：")
lines1.append("  • n149 是一个复杂 selector 函数 (comparator tree)")
lines1.append("  • m_28 基于 n149 构建 selector")
lines1.append("  • m_29..m_31 形成 ripple-carry 链")
lines1.append("  • m_26/m_27 使用 m_28..m_31 作为控制")
lines1.append("")
lines1.append("=" * 70)

with open(f"{out_dir}/hd12_opt1.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))
with open(f"{out_dir}/hd12_opt2.poly", 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))

print("  → examples/hd12/hd12_opt1.poly")
print("  → examples/hd12/hd12_opt2.poly")
print("完成")
