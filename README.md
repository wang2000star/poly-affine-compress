# 布尔函数 ANF 简化 / Boolean Function ANF Simplification

寻找仿射变换 $z = Mx \oplus b$（$M \in \mathbb{F}_2^{m \times n}, b \in \mathbb{F}_2^m$）使得 $g(z) = g(Mx \oplus b) = f(x)$ 且 $g$ 的 ANF 单项式数 $T(g)$ 尽可能小。

## 文件说明

### 核心库

| 文件 | 说明 |
|------|------|
| `bool_anf.py` | 基础库：真值表/ANF 互转、Walsh-Hadamard 变换、GF(2) 线性代数 |
| `anf_factor.py` | **单布尔函数 ANF 简化**：SparseANF 类、贪心 XOR 合并、Walsh 方向搜索、互补搜索、简化 pipeline |
| `int_poly_factor.py` | **单整系数多项式简化**：IntPoly 类、梯度引导合并、贪心合并、简化 pipeline |
| `vector_anf.py` | **向量布尔函数联合简化**：共享变换（策略 1）和各自变换（策略 2），支持 `--seed` 和 `--walsh-trials` |
| `vector_int_poly.py` | **向量整系数多项式联合简化**：向量版贪心合并、各自变换 |

### 基准测试

| 文件 | 说明 |
|------|------|
| `benchmark_anf.py` | 布尔函数 ANF 基准测试 |
| `benchmark_int.py` | 整系数多项式基准测试 |

### 文档

| 文件 | 说明 |
|------|------|
| `anf-simple.md` | 问题定义文档 |
| `solution-anf-simplify.md` | 理论解答 |
| `CLAUDE.md` | Claude Code 项目配置 |

### 旧文件

| 文件 | 说明 |
|------|------|
| `_misc/` | 历史参考文件 |

## 快速开始

```bash
# 单布尔函数简化演示
python anf_factor.py

# 向量布尔函数联合简化
python vector_anf.py --seed 0 --walsh-trials 10

# 向量整系数多项式联合简化
python vector_int_poly.py --seed 0

# 基准测试
python benchmark_anf.py
python benchmark_int.py
```

## 两种策略

- **策略 1（共享变换）**：所有分量使用同一个 $M, b$，计算并集单项式数
- **策略 2（各自变换）**：每个分量独立寻找自己的 $M_i, b_i$，可进一步合并公共线性形式计算并集
