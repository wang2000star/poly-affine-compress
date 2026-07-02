"""
hd10.txt → 优化（符号分析版）
32 输入 32 输出，但 m_0..m_28 ≡ 0，仅 m_29/m_30/m_31 非平凡。
该电路为 AND-of-ORs 结构，ANF 全展开有 ~4B 项不可行。
本脚本直接做符号分析并写出紧凑表达式。
"""
import sys
sys.path.insert(0, '/home/wangfz/bool')
import time

inputs_32 = [f'i_{j}' for j in range(2, 34)]
n = 32


def group_or_anf(vars):
    """OR(vars) 的 ANF 字典：所有非空单项式"""
    terms = {}
    n_v = len(vars)
    # 生成所有 2^n_v - 1 个非空单项式
    for mask in range(1, 1 << n_v):
        monom = 0
        for k in range(n_v):
            if mask & (1 << k):
                bit = int(vars[k][2:]) - 2  # i_j → bit index
                monom |= 1 << bit
        terms[monom] = 1
    return terms


def multiply_anf_disjoint(A, B):
    """两个不相交变量组上的 ANF 相乘"""
    result = {}
    for ma, va in A.items():
        for mb, vb in B.items():
            mc = ma | mb
            result[mc] = result.get(mc, 0) ^ (va & vb)
    return {k: v for k, v in result.items() if v}


def anf_add(A, B):
    """ANF 相加（XOR）"""
    result = dict(A)
    for k, v in B.items():
        result[k] = result.get(k, 0) ^ v
    return {k: v for k, v in result.items() if v}


def anf_not(A):
    """ANF 取非"""
    result = dict(A)
    result[0] = result.get(0, 0) ^ 1
    return {k: v for k, v in result.items() if v}


def T(terms):
    return len(terms)


# ---- 构造各组的 OR ----
groups = {
    1: [f'i_{j}' for j in range(2, 10)],
    2: [f'i_{j}' for j in range(10, 18)],
    3: [f'i_{j}' for j in range(18, 26)],
    4: [f'i_{j}' for j in range(26, 34)],
}

print("=" * 60)
print("hd10.txt → 符号分析")
print("=" * 60)

t0 = time.time()

# 各组的 OR
OR_g = {}
for gid, vars in groups.items():
    OR_g[gid] = group_or_anf(vars)
    print(f"  OR(G{gid}) T = {T(OR_g[gid])}")

# m_29 = OR(G1) · OR(G2) · OR(G3) · OR(G4)
m29 = OR_g[1]
for gid in [2, 3, 4]:
    m29 = multiply_anf_disjoint(m29, OR_g[gid])

print(f"\n  m_29 ANF T = {T(m29)}")

# m_30 = n162 · n126
# n126 = OR(G1) · OR(G2)
n126 = multiply_anf_disjoint(OR_g[1], OR_g[2])
# n125 = OR(G1) (alias)
n125 = OR_g[1]
# n108 = NOR of G2 = NOT OR(G2)
n108 = anf_not(OR_g[2])
# n109 = NOT n108 = OR(G2)
n109 = OR_g[2]
# n157 = NOR of G3 = NOT OR(G3)
n157 = anf_not(OR_g[3])
# n158 = NOT n157 = OR(G3)
n158 = OR_g[3]
# n141 = NOR of G4 = NOT OR(G4)
n141 = anf_not(OR_g[4])
# n142 = NOT n141 = OR(G4)
n142 = OR_g[4]
# n159 = n158 · n142 = OR(G3) · OR(G4)
n159 = multiply_anf_disjoint(OR_g[3], OR_g[4])

# n161 = n158 · n141 = OR(G3) · NOT OR(G4)
n161 = multiply_anf_disjoint(OR_g[3], n141)
# n162 = n161 ⊕ n157 = (OR(G3)·NOR(G4)) ⊕ NOR(G3)
n162 = anf_add(n161, n157)

# m_30 = n162 · n126
m30 = multiply_anf_disjoint(n162, n126)
print(f"  m_30 ANF T = {T(m30)}")

# m_31 = n165 · n125
# n165 = n164 ⊕ n108
# n164 = n161 · n109 = (OR(G3)·NOR(G4)) · OR(G2)
n164 = multiply_anf_disjoint(n161, OR_g[2])
# n165 = n164 ⊕ n108
n165 = anf_add(n164, n108)
# m_31 = n165 · n125
m31 = multiply_anf_disjoint(n165, OR_g[1])
print(f"  m_31 ANF T = {T(m31)}")

dt = time.time() - t0
print(f"\n  符号分析耗时: {dt:.1f}s")

# ---- 输出 ----
# 策略 1 说明
lines1 = []
lines1.append("=" * 70)
lines1.append("  策略 1（共享变换）: hd10.txt")
lines1.append("=" * 70)
lines1.append("")
lines1.append(f"  原始输入: {', '.join(inputs_32)}")
lines1.append(f"  输出函数: 32 个")
lines1.append(f"  m_0..m_28 = 0 (常数)")
lines1.append("")
lines1.append("  ANF 全展开不可行（函数为 AND-of-ORs 结构，")
lines1.append(f"  非平凡输出 ANF 达 ~{T(m29)//10**9}B 项）。")
lines1.append("  以下为符号表示的 ANF 项数：")
lines1.append("")
lines1.append(f"  ΣT₀ (ANF 项数):")
lines1.append(f"    m_0..m_28: T = 0 (×29)")
lines1.append(f"    m_29:       T = {T(m29)}")
lines1.append(f"    m_30:       T = {T(m30)}")
lines1.append(f"    m_31:       T = {T(m31)}")
print(f"\n  union T₀ = {T(m29) + T(m30) + T(m31)}")
lines1.append("")
lines1.append(f"  union T₀ ≈ {T(m29) + T(m30) + T(m31)}")
lines1.append("")
lines1.append("  由于 ANF 过大，共享变换的搜索不可行。")
lines1.append("  函数本身已是乘积形式：")
lines1.append("")
lines1.append("  m_29 = OR(G₁) · OR(G₂) · OR(G₃) · OR(G₄)")
lines1.append("      其中 G₁={i_2..i_9}, G₂={i_10..i_17},")
lines1.append("           G₃={i_18..i_25}, G₄={i_26..i_33}")
lines1.append("")
lines1.append("  m_30 = [(OR(G₃)·NOR(G₄)) ⊕ NOR(G₃)] · OR(G₁)·OR(G₂)")
lines1.append("")
lines1.append("  m_31 = [(OR(G₃)·NOR(G₄)·OR(G₂)) ⊕ NOR(G₂)] · OR(G₁)")
lines1.append("")
lines1.append("  任何 GF(2) 仿射变换无法本质简化此结构。")
lines1.append("")
lines1.append("  各输出 ANF 展开:")
lines1.append("")

for name, terms in [("m_29", m29), ("m_30", m30), ("m_31", m31)]:
    lines1.append(f"  {name}: T={T(terms)}")
    # 输出前 20 项
    parts = []
    for mask in sorted(terms, key=lambda m: (bin(m).count('1'), m))[:20]:
        if mask == 0:
            parts.append("1")
        else:
            monom = [inputs_32[j] for j in range(n) if mask & (1 << j)]
            parts.append("·".join(monom))
    lines1.append(f"    前 20 项: {' ⊕ '.join(parts)}")
    lines1.append(f"    ... (共 {T(terms)} 项)")
    lines1.append("")

lines1.append("=" * 70)

with open('/home/wangfz/bool/hd10/hd10_opt1.poly', 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines1))

# 策略 2 类似
lines2 = []
lines2.append("=" * 70)
lines2.append("  策略 2（各自变换）: hd10.txt")
lines2.append("=" * 70)
lines2.append("")
lines2.append(f"  原始输入: {', '.join(inputs_32)}")
lines2.append(f"  输出函数: 32 个 (29 个恒 0)")
lines2.append("")
lines2.append("  ANF 全展开不可行（同上）。")
lines2.append("  各输出在原始变量下的 ANF 项数：")
lines2.append("")
for name in [f"m_{i}" for i in range(29)]:
    lines2.append(f"  {name}: 0")
lines2.append(f"  m_29: T = {T(m29)} (AND of 4 ORs)")
lines2.append(f"  m_30: T = {T(m30)} (XOR combination)")
lines2.append(f"  m_31: T = {T(m31)} (XOR combination)")
lines2.append("")
lines2.append("  策略 2 需要对每个输出单独找仿射变换。")
lines2.append("  由于 ANF 过大（~B 级），无法应用 simplify 管线。")
lines2.append("")
lines2.append("  各输出表达式（紧凑形式）:")
lines2.append("")
lines2.append("  m_29 = OR(i_2..i_9) · OR(i_10..i_17) · OR(i_18..i_25) · OR(i_26..i_33)")
lines2.append("")
lines2.append("  m_30 = [(OR(i_18..i_25) · NOR(i_26..i_33)) ⊕ NOR(i_18..i_25)]")
lines2.append("         · OR(i_2..i_9) · OR(i_10..i_17)")
lines2.append("")
lines2.append("  m_31 = [(OR(i_18..i_25) · NOR(i_26..i_33) · OR(i_10..i_17)) ⊕ NOR(i_10..i_17)]")
lines2.append("         · OR(i_2..i_9)")
lines2.append("")
lines2.append("=" * 70)

with open('/home/wangfz/bool/hd10/hd10_opt2.poly', 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines2))

print("\n  → hd10/hd10_opt1.poly")
print("  → hd10/hd10_opt2.poly")
print(f"\n完成")
print("=" * 60)
