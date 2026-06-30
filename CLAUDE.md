# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

布尔函数 ANF（Algebraic Normal Form）简化研究。

核心问题：对于 n 元布尔函数 f(x) ∈ F₂ⁿ → F₂，寻找仿射变换 z = Mx ⊕ b（M ∈ F₂^{m×n}, b ∈ F₂^m）使得 g(z) = g(Mx ⊕ b) = f(x)，且 g 的 ANF 单项式数 T(g) 相对于 T(f) 尽可能小。

## 研究方向

- 代数次数在仿射变换下的不变性
- 仿射等价分类与标准型
- ANF 稀疏化的高效算法
- 特殊函数类（二次型、对称函数等）的最优压缩率

## Git 工作流

- 每次完成有意义的研究进展后，使用 `git commit` 保存
- 提交信息格式：中文简要描述做了什么

## 数学符号约定

- F₂: 二元域 {0, 1}，加法为 XOR
- ANF: Algebraic Normal Form，布尔函数的多项式表示
- T(f): f 的 ANF 中单项式个数
- deg(f): f 的代数次数
- 变量列向量：x = (x₀, ..., x_{n-1})ᵀ ∈ F₂ⁿ
