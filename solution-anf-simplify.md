# 布尔函数 ANF 仿射简化问题解答

## 问题回顾

设 $f(x)$ 是 $n$ 元布尔函数，$x\in\mathbb{F}_2^n$，$T(f)$ 是 $f$ 的 ANF 形式单项式个数。

存在 $M\in\mathbb{F}_2^{m\times n}, b\in\mathbb{F}_2^m$ 使得 $z=Mx\oplus b\in\mathbb{F}_2^m$ 且
$g(z)=g(Mx\oplus b)=f(x)$。

- **问题 1**：$\deg(g)$ 与 $\deg(f)$ 是否相同？
- **问题 2**：如何构造 $M,b$ 使得 $\dfrac{T(g)}{T(f)}$ 尽可能小？是否存在高效算法？

---

## 问题 1：仿射变换下的代数次数不变性

**结论**：对于可逆仿射变换 $z = Mx \oplus b$（即 $M\in\operatorname{GL}(n,\mathbb{F}_2)$），
有 $\deg(g) = \deg(f)$。

### 证明

设 $f(x) = g(Mx\oplus b)$ 其中 $M$ 可逆。

1. **$\deg(g) \ge \deg(f)$**: 由于 $M$ 可逆，有 $x = M^{-1}(z\oplus b)$。每个 $x_i$ 是 $z$ 的
   一次仿射函数。因此 $g(z) = f(M^{-1}(z\oplus b))$ 是 $f$ 与一次仿射函数之复合。
   多项式复合不提高代数次数，故 $\deg(g) \le \deg(f)$。

2. **$\deg(g) \ge \deg(f)$**: 每个 $z_j = (Mx)_j \oplus b_j$ 是 $x$ 的一次仿射函数。
   $f(x) = g(Mx\oplus b)$ 是 $g$ 与一次函数之复合，故 $\deg(f) \le \deg(g)$。

综上 $\deg(g) = \deg(f)$。

$\square$

**注**：若 $M$ 不可逆（即 $m < n$ 或 $M$ 秩不足），则 $g$ 可视为 $f$ 在子空间上的限制，
可能出现 $\deg(g) < \deg(f)$。但对于 $g(z) = f(x)$ 对任意 $x$ 成立这一要求，
$M$ 必须有足够的秩来承载 $f$ 的所有信息。

---

## 问题 2：最小化 $T(g)/T(f)$ 的理论分析与算法

这是问题的核心。$T(f)$ 在仿射变换下**不是**不变量——这正是我们得以简化的原因。
一个单项式在 $x$ 变量下可能展开为多个 $z$ 单项式，反之多个 $x$ 单项式可能
"合并" 为更少的 $z$ 单项式。

### 2.1 理论基础

#### 2.1.1 ANF 的向量表示

将每个单项式 $x^s = \prod_{i:s_i=1} x_i$（其中 $s \in \{0,1\}^n$ 是指示向量）
视为 $\mathbb{F}_2$ 上的一个基向量。则 $f$ 可以写为：

$$f(x) = \bigoplus_{s \in \{0,1\}^n} a_s \cdot x^s,\quad a_s \in \mathbb{F}_2$$

$T(f)$ 就是非零系数 $a_s$ 的个数。

在仿射变换 $z = Mx \oplus b$ 下，每个 $x_i$ 变为 $z$ 的仿射函数：

$$x_i = (M^{-1}z)_i \oplus (M^{-1}b)_i = \bigoplus_{j=1}^n (M^{-1})_{ij} z_j \oplus c_i$$

#### 2.1.2 与仿射等价的联系

两个布尔函数 $f,g$ 称为**仿射等价**的，如果存在
$M \in \operatorname{GL}(n,\mathbb{F}_2), b \in \mathbb{F}_2^n$ 和一个仿射函数
$L: \mathbb{F}_2^n \to \mathbb{F}_2$ 使得：

$$f(x) = g(Mx \oplus b) \oplus L(x)$$

我们的情况是 $g$ 是 $f$ 在去掉 $L(x)$ 项后的像（允许 $L(x) = 0$），且允许
$m$ 与 $n$ 不同。但关键是：**$T(g)$ 的最小可能值取决于 $f$ 的仿射等价类**。

### 2.2 特殊函数类的精确解

#### 2.2.1 仿射函数 ($\deg(f) \le 1$)

若 $f$ 是仿射函数，则存在可逆仿射变换使其变为常数或单个变量：

$$T_{\min} = \begin{cases}
0 & \text{若 $f$ 为常数} \\
1 & \text{若 $f$ 为非平凡仿射函数}
\end{cases}$$

即 $T(g)/T(f) \to 0$（当 $T(f) > 1$ 时）。

#### 2.2.2 二次型函数 ($\deg(f) = 2$)

**这是唯一有完备分类的非平凡情形。**

二次布尔函数 $f(x) = x^\mathrm{T} A x \oplus L(x) \oplus c$，其中 $A$ 是严格上三角矩阵。
令 $B = A + A^\mathrm{T}$ 为 $\mathbb{F}_2$ 上的交错双线性型。

$B$ 的秩 $r$ 在可逆仿射变换下不变。根据秩 $r$ 的奇偶性，二次型有以下标准型：

$$
\begin{aligned}
\text{秩 $r$ 为偶数}: &\quad x_1x_2 \oplus x_3x_4 \oplus \cdots \oplus x_{r-1}x_r \oplus c \\
\text{秩 $r$ 为奇数}: &\quad x_1x_2 \oplus x_3x_4 \oplus \cdots \oplus x_{r-2}x_{r-1} \oplus x_r^2 \oplus c
\end{aligned}
$$

其中 $x_r^2 = x_r$（在 $\mathbb{F}_2$ 上），所以：

- 非退化二次型（$r = n$ 为偶）：$T_{\min} = n/2 + \delta$（$\delta=0$ 或 $1$ 取决于有无常数项）
- 更一般地，$T_{\min} \le \lfloor r/2 \rfloor + 2$

**算法**：由于 $B$ 的秩 $r$ 可以在 $O(n^3)$ 时间内通过高斯消元计算出，
同时对 $B$ 进行对称消元可直接得到变换矩阵。详见代码实现。

**示例**：$f(x) = x_0x_1 + x_0x_2 + x_1x_2 + x_0$
- $B = \begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix}$，秩 $r = 2$
- 标准型：$z_0z_1 + \cdots$（仅需 1-2 个单项式）

#### 2.2.3 具有线性结构的函数

若 $f(x \oplus a) = f(x)$ 对所有 $x$ 成立（$a \neq 0$），则 $a$ 是 $f$ 的一个线性结构。
存在仿射变换将 $a$ 对齐到某个坐标轴，从而减少一个变量。

**算法**：计算 $f$ 的自相关函数 $\mathcal{AC}_f(a) = \sum_x f(x) \oplus f(x\oplus a)$。
对每个 $a \neq 0$ 检查是否为线性结构。若存在 $k$ 个独立线性结构，
则可将 $n$ 元函数降为 $n-k$ 元。

### 2.3 一般情况的计算复杂性

对于 $\deg(f) \ge 3$ 的一般布尔函数，最小化 $T(g)/T(f)$ 问题被认为是计算困难的。

**相关已知结论**：

1. **仿射等价判定问题**：判断两个 $\deg \ge 3$ 的布尔函数是否仿射等价是
   $\mathrm{NP}$-难问题（已知归约自 Graph Isomorphism 问题的变种）。

2. **ANF 最小化**：找 $M,b$ 最小化 $T(g)$ 包含了仿射等价判定作为子问题
   （如果两个函数仿射等价，它们的 $T_{\min}$ 相同），因此至少是 $\mathrm{NP}$-难的。

3. **与算术电路复杂性的联系**：$T(f)$ 的最小值相关于 $f$ 的「乘法复杂度」
   （multiplicative complexity），即计算 $f$ 所需的最小 $\mathbb{F}_2$ 上乘积门数。

尽管一般问题是难的，以下启发式算法在实践中往往有效。

### 2.4 算法详述

#### 算法 A：暴力搜索（$n \le 5$）

对于 $n \le 5$，可穷举所有可逆仿射变换：

$$|\operatorname{GL}(n,\mathbb{F}_2)| = \prod_{i=0}^{n-1} (2^n - 2^i)$$

| $n$ | $|\operatorname{GL}(n,\mathbb{F}_2)|$ | 可行性 |
|-----|--------------------------------------|--------|
| 3   | 168                                  | 瞬间   |
| 4   | 20,160                               | 瞬间   |
| 5   | 9,999,360                            | 数秒   |
| 6   | 20,158,709,760                       | 不可行 |

对于每个 $(M,b)$ 计算 $g$ 的 ANF 并计数，取全局最小值。

#### 算法 B：二次型对角化（$\deg(f) \le 2$）

**输入**：$f(x) = x^\mathrm{T} A x \oplus Lx \oplus c$

**步骤**：
1. 计算 $B = A + A^\mathrm{T}$ 及其秩 $r$ 和左零空间 $N(B)$
2. 通过高斯消元求 $P$ 使 $P^\mathrm{T} B P = \begin{pmatrix} 0 & I_{r/2} \\ I_{r/2} & 0 \end{pmatrix} \oplus 0_{n-r}$（若 $r$ 偶）
3. 设在 $P$ 变换后的坐标 $(y_1,\ldots,y_n)$ 下，$f$ 变为：
   $$y_1y_2 + y_3y_4 + \cdots + y_{r-1}y_r + L'(y) + c$$
4. 通过仿射变换消去线性部分 $L'(y)$ 中与 $y_{r+1},\ldots,y_n$ 相关的项
5. 在零空间方向上，这些变量只以线性形式出现，可进一步化简

**复杂度**：$O(n^3)$。

#### 算法 C：贪心坐标合并（通用启发式）

**核心思想**：寻找 $x_i, x_j$ 使得合并它们（令 $z_k = x_i \oplus x_j$ 取代
$x_i$ 或 $x_j$ 作为新变量）能最大程度减少 $T$。

**步骤**：

```
输入: f 的 ANF 系数表 A: 2^n 个系数
输出: 变换矩阵 M, b

1. 初始化 M = I_n, b = 0
2. 对每个候选合并方向 d ∈ {0,1}^n \ {0}:
   a. 检查使用 d 作为某个新变量方向是否能减少 T
   b. 评分: ΔT(d) = T(g_after) - T(before)
3. 选择 ΔT 最小的 d 执行变换
4. 更新 M, b
5. 重复 2-4 直到再无改进或达到最大迭代次数
```

**评分函数的高效计算**：对每个候选方向 $d$，
新变量 $z_k = \langle d, x \rangle \oplus b_k$ 的引入将原 ANF 中的单项式
根据是否包含 $d$ 中的变量进行分组。利用 Mobius 变换的性质，
可在 $O(n \cdot 2^n)$ 时间内评估所有单方向变化。

#### 算法 D：基于 Walsh-Hadamard 谱的方法

**观察**：仿射变换对应 Walsh 谱的置换和符号变化。

- 计算 $f$ 的 Walsh-Hadamard 变换 $\hat{f}(\omega)$
- 寻找谱系数集中的方向——这暗示该方向上有低次数的依赖
- 谱集中的方向往往对应着可简化的变量组合

**步骤**：
1. 计算 Walsh 谱 $W_f(\omega) = \sum_x (-1)^{f(x) \oplus \langle \omega, x \rangle}$
2. 找出 $W_f(\omega)$ 绝对值较大的 $\omega$（排除 $\omega=0$）
3. 将这些 $\omega$ 作为新变量方向的候选
4. 测试每个候选对 ANF 复杂度的影响
5. 选择最优方向执行变换

**复杂度**：$O(n \cdot 2^n)$（Walsh 变换通过快速变换可达 $O(n \cdot 2^n)$）。

### 2.5 算法比较

| 算法 | 适用范围 | 复杂度 | 最优性保证 |
|------|---------|--------|-----------|
| A 暴力搜索 | $n \le 5$ | 指数 (穷举) | 全局最优 |
| B 二次型对角化 | $\deg(f) \le 2$ | $O(n^3)$ | 全局最优 |
| C 贪心合并 | 任意 $n$ | $O(n^3 \cdot 2^n)$ 启发式 | 无保证 |
| D Walsh 谱方法 | 任意 $n$ | $O(n \cdot 2^n)$ | 无保证 |

### 2.6 开问题与展望

1. **$\deg(f) \ge 3$ 的最小 $T$ 的刻画**：是否存在类似二次秩的不变量
   来刻画 $T_{\min}$？可能是某种高阶导数的秩或广义的多元相关系数。

2. **随机函数的 $T_{\min}$**：随机 $n$ 元布尔函数的 ANF 平均有 $2^{n-1}$ 个单项式。
   仿射变换能将这个数降低多少？是否存在普适的下界？

3. **与编码理论的联系**：将 ANF 系数向量视为 $\mathbb{F}_2^{2^n}$ 中的向量，
   仿射变换的作用对应某个线性群的表示。$T(f)$ 就是该向量的 Hamming 重量。
   问题转化为：求群作用下轨道中重量的最小值。这类似于寻找码的最小重量。

---

## 3. 扩展：向量布尔函数的 ANF 简化

### 3.1 问题定义

考虑向量布尔函数（S-盒）$F: \mathbb{F}_2^n \to \mathbb{F}_2^k$：

$$F(x) = (f^{(0)}(x), f^{(1)}(x), \ldots, f^{(k-1)}(x))$$

其中每个 $f^{(j)}$ 是 $n$ 元布尔函数。目标是找到仿射变换
$z = Mx \oplus b$（输入变换），使得变换后的函数
$G(z) = F(M^{-1}(z \oplus b))$（或等价地 $F(x) = G(Mx \oplus b)$）满足：

1. **次数约束**：$\deg(G^{(j)}) \le \deg(F^{(j)})$ 对所有 $j$ 成立
2. **项数极小化**：总项数 $T_{\text{total}}(G) = \sum_{j=0}^{k-1} T(G^{(j)})$ 显著小于 $T_{\text{total}}(F)$

### 3.2 关键差异分析与标量情况的对比

| 方面 | 标量布尔函数 | 向量布尔函数 |
|------|-------------|-------------|
| 优化目标 | 单个 $T(g)$ | 加权和 $\sum w_j \cdot T(g^{(j)})$ |
| 变换自由度 | 输入侧仿射 | 输入侧仿射 + 输出侧仿射（可选） |
| 次数约束 | 自动保持 | 需对各分量分别检查 |
| 最优解一致性 | 单一目标 | 各分量可能冲突 |

### 3.3 扩展算法

#### 算法 E1：加权贪心法

将各分量的 ANF 系数合并为加权和，在单一搜索中优化总目标。

**步骤**：
1. 计算 $F$ 中各分量 $f^{(j)}$ 的 ANF 系数
2. 构造加权 ANF 计数 $T_w(F) = \sum_{j=0}^{k-1} w_j \cdot T(f^{(j)})$
3. 对候选仿射变换 $(M,b)$，计算 $T_w(G) = \sum_j w_j \cdot T(g^{(j)})$
4. 选择使 $T_w$ 最小的变换

**权重选择**：
- 等权 $w_j = 1$：优化总项数
- 按次数加权 $w_j = \binom{n}{\deg(f^{(j)})}$：惩罚高次分量
- 自适应权重：迭代中根据改进难度调整

#### 算法 E2：主分量主导法

选择 $T(f^{(j)})$ 最大的分量作为主导分量，优先优化它：

1. 设 $j^* = \arg\max_j T(f^{(j)})$ 为主导分量
2. 对 $f^{(j^*)}$ 运行标量 ANF 简化算法，得到最优变换 $(M,b)$
3. 验证该变换对其他分量的影响：
   $$T_{\text{total}}(G) \stackrel{?}{\ll} T_{\text{total}}(F)$$
4. 若副作用过大，回退到加权优化

#### 算法 E3：输出仿射变换联合优化

在输入变换 $z = Mx \oplus b$ 的同时，考虑输出仿射变换
$y = N \cdot G(z) \oplus c$（$N \in \operatorname{GL}(k,\mathbb{F}_2)$）。

$$F(x) = N \cdot G(Mx \oplus b) \oplus c$$

输出仿射不改变各分量内部的 ANF 结构，但它：
- 将分量线性组合，可能**重新分配**各分量的复杂度
- 不同 $N$ 下各分量 $T(g^{(j)})$ 的和可能变化

**优化策略**：
1. 固定输入变换 $(M,b)$
2. 在输出侧枚举 $N \in \operatorname{GL}(k,\mathbb{F}_2)$ 和 $c \in \mathbb{F}_2^k$
3. 选取使 $\max_j T(g^{(j)})$ 最小或 $\sum_j T(g^{(j)})$ 最小的组合
4. 交替优化输入和输出变换

#### 算法 E4：矩阵分解法

当 $k$ 较大时，将 $F$ 视为矩阵值函数，
每个单项式的系数是一个 $k$ 维向量 $\vec{a}_s \in \mathbb{F}_2^k$：

$$F(x) = \bigoplus_{s \in \{0,1\}^n} \vec{a}_s \cdot x^s$$

仿射变换 $z = Mx \oplus b$ 重新排列和重组这些系数向量。
目标是找到 $M,b$ 使得系数矩阵中**非零行的数量**最小化。

### 3.4 实用性考量

对于密码学中常见的 S-盒（如 $n=8, k=8$）：

- 暴力搜索完全不可行（$|\operatorname{GL}(8,\mathbb{F}_2)| \approx 10^{17}$）
- 推荐策略：Walsh 谱分析 + 加权贪心迭代
- 可将 8 位 S-盒拆分为 $2 \times 4$ 位子盒分别优化
- 利用代数次数约束（通常 S-盒 $\deg \le 6$）来剪枝搜索空间

### 3.5 扩展版本代码实现

见 `bool_anf.py` 中以下函数：
- `VectorBooleanFunction` 类：向量布尔函数 ANF 表示
- `vector_weighted_minimize()`：加权贪心简化
- `vector_input_output_minimize()`：输入输出联合优化

---



## 参考文献

1. Carlet, C. (2010). *Boolean Functions for Cryptography and Error Correcting Codes*.
   In: Crama, Y., Hammer, P. (eds) Boolean Models and Methods in Mathematics,
   Computer Science, and Engineering.
2. Ryabov, A. (2018). *On Classification of Cubic Boolean Functions under
   Affine Transformations*. Journal of Applied and Industrial Mathematics.
3. MacWilliams, F. J., Sloane, N. J. A. (1977). *The Theory of Error-Correcting Codes*.
4. Bard, G. V. (2009). *Algorithms for Solving Polynomial Equations over Finite Fields*.
   (关于 Möbius 变换和 ANF 的基础算法)
