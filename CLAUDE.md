# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

布尔函数 ANF（Algebraic Normal Form）简化研究。

核心问题：对于 n 元布尔函数 f(x) ∈ F₂ⁿ → F₂，寻找仿射变换 z = Mx ⊕ b（M ∈ F₂^{m×n}, b ∈ F₂^m，m 可小于 n，M 无需可逆或满秩）使得 g(z) = g(Mx ⊕ b) = f(x) 是纯代数多项式恒等式，且 g 的 ANF 单项式数 T(g) 相对于 T(f) 尽可能小。

扩展：向量布尔函数 F: F₂ⁿ → F₂ᵏ 的 ANF 简化（各分量联合优化）。

## 核心原理

- 纯代数代换：将 z = Mx⊕b 代入 g 得到 f，是多项式环 F₂[x]/(x_i²−x_i) 中的恒等式
- 无信息损失、无纤维条件、无秩限制：M 可以是任意 m×n 矩阵
- Boole 环乘法：x_i² = x_i，单项式相乘用 OR (|) 而非 XOR
- deg(g) = deg(f) 对任意 M,b 成立

## 文件结构

- `bool_anf.py` — 核心 library：真值表/ANF 互转、GF(2) 线性代数、Walsh 变换、向量布尔函数类
- `anf_factor.py` — 稀疏 ANF 实现（SparseANF 类）：代数代换 substitute_affine、线性系统求解 solve_sparsest_g、Walsh 方向搜索
- `benchmark_anf.py` — 基准测试：n=8 到 n=64 的性能与压缩率测试（导入 anf_factor 的 SparseANF）
- `solution-anf-simplify.md` — 理论解答文档

## SparseANF 类

位于 anf_factor.py。dict {mask: coeff} 存储 ANF，n 可达 64。

关键方法：
- `substitute_affine(M, b)` — 对任意 M,b 进行纯代数代换得到 g(z) = f(x)。m<n 时自动验证正确性
- `solve_sparsest_g(M, b, max_m=10)` — 线性代数方法求最稀疏的 g（适合 m≤10）
- `eval_mask(x)` — 在输入 x 处求值

## 常用命令

```bash
# 运行演示
python anf_factor.py

# 运行全部基准测试（含 n=8 到 n=64）
python benchmark_anf.py

# 仅运行速度基准
python -c "from benchmark_anf import bench_transform_time; bench_transform_time()"

# 仅运行 n=8 正确性验证
python -c "from benchmark_anf import verify_n8; verify_n8()"
```

## 测试函数生成

- `SparseANF.random_cubic(n, density, seed)` — 随机 ≤3 次函数
- `SparseANF.from_linear_structure(n, inner_n, seed)` — F(x)=G(Ax)，G 有 inner_n 个变量
- `benchmark_anf.make_power_function(n, e)` — F(x)=LSB(x^e) over F_{2^n}
- `benchmark_anf.make_hidden_structure(n, inner_n, seed)` — G 通过随机 H 隐藏

## 搜索策略

`search_affine_simplification(f, max_m, top_k, n_random)`:
1. Walsh 谱方向（单行候选）
2. 前 k 个 Walsh 方向组合（多行）
3. 随机 M,b（降维搜索）

返回按 T(g) 排序的候选列表。
