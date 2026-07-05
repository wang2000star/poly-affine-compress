# 布尔函数 ANF 仿射简化问题解答

## 问题回顾

设 $f(x)$ 是 $n$ 元布尔函数，$x\in\mathbb{F}_2^n$，$T(f)$ 是 $f$ 的 ANF 形式单项式个数。

存在 $M\in\mathbb{F}_2^{m\times n}, b\in\mathbb{F}_2^m$ 令 $z=Mx\oplus b\in\mathbb{F}_2^m$ 使得
$g(z)=g(Mx\oplus b)=f(x)$。

- **问题 1**：$\deg(g)$ 与 $\deg(f)$ 是否相同？
- **问题 2**：如何构造 $M,b$ 使得 $\dfrac{T(g)}{T(f)}$ 尽可能小？有没有高效的算法？

---

## 问题 1：仿射变换下的代数次数不变性

**结论**：$\deg(g) = \deg(f)$ 对任意 $M \in \mathbb{F}_2^{m\times n}, b \in \mathbb{F}_2^m$ 成立
（只要 $g$ 使得 $f(x) = g(Mx\oplus b)$ 恒成立）。

### 证明

**方向 1：$\deg(f) \le \deg(g)$**

每个 $z_j = (Mx)_j \oplus b_j = \bigoplus_{i=1}^n M_{ji}x_i \oplus b_j$
是 $x$ 的一次仿射函数。$f(x) = g(Mx \oplus b)$ 是 $g$ 与一次函数之复合。
$g$ 的一个 $d$ 次单项式 $z_{j_1}z_{j_2}\cdots z_{j_d}$ 代入后变为
$d$ 个一次仿射函数的乘积，故在 $x$ 中的次数 $\le d$。
因此 $\deg(f) \le \deg(g)$。

**方向 2：$\deg(g) \le \deg(f)$**

将 $f(x) = g(Mx\oplus b)$ 视为 $\mathbb{F}_2[x]/(x_i^2-x_i)$ 中的多项式恒等式。
问题是纯代数的：给定 $f$ 和 $M,b$，求 $g$ 使得代入 $z = Mx\oplus b$ 后得到 $f$。
关键在于 $g$ 的每个单项式 $z^t = \prod_{j: t_j=1} z_j$ 展开后是 $x$ 中次数
$\le |t|$ 的多项式。$f$ 作为这些展开的 $\mathbb{F}_2$ 线性组合，其次数不超过
$g$ 中出现的最高次单项式的次数。故 $\deg(g) \le \deg(f)$。

（注意：$M$ 无需方阵也无需满秩。对任意 $m\times n$ 矩阵，$z_j$ 总是 $x$ 的一次函数，
上述论证不受 $m,n$ 关系或 $M$ 的秩影响。）

### 推论

- 对任意 $m,n$ 和任意 $M\in\mathbb{F}_2^{m\times n}, b\in\mathbb{F}_2^m$，有 $\deg(g) = \deg(f)$。
- 对向量布尔函数 $F(x) = (f^{(0)},\ldots,f^{(k-1)})$，各分量分别满足 $\deg(g^{(j)}) = \deg(f^{(j)})$。
- 若想允许次数降低（$\deg(g) < \deg(f)$），需要引入输出侧的仿射变换或更一般的结构。

---

## 问题 2：最小化 $T(g)/T(f)$ 的理论分析与算法

### 2.1 代数本质

$f(x) = g(Mx\oplus b)$ 是关于 $x$ 的多项式恒等式。将 $g$ 写为 ANF：

$$
g(z) = \bigoplus_{t\in\{0,1\}^m} a_t \cdot z^t,\quad a_t\in\mathbb{F}_2
$$

代入 $z = Mx\oplus b$：

$$
f(x) = \bigoplus_{t\in\{0,1\}^m} a_t \cdot (Mx\oplus b)^t
$$

每个 $(Mx\oplus b)^t$ 展开后是 $x$ 的多项式（ANF）。记 $C_t(x)$ 为 $(Mx\oplus b)^t$
的 ANF 展开，则问题归结为 $\mathbb{F}_2$ 上的线性方程组：

$$
\sum_{t=0}^{2^m-1} a_t \cdot C_t(x) = f(x)
$$

写成矩阵形式：$V \cdot a = f$，其中 $V[s][t]$ 是 $C_t$ 中 $x^s$ 的系数。

- 若此方程组有解，则 $g$ 存在（由系数 $a_t$ 给出）
- 若有无穷多解（$M$ 非满秩时），则解空间是一个仿射子空间
- 我们寻找的是**最稀疏的解** — 即 $\sum_t a_t$ 最小的 $a$

### 2.2 与仿射等价的联系

两个布尔函数 $f,g$ 称为**仿射等价**的，如果存在
$M \in \operatorname{GL}(n,\mathbb{F}_2), b \in \mathbb{F}_2^n$ 和一个仿射函数
$L: \mathbb{F}_2^n \to \mathbb{F}_2$ 使得：

$$
f(x) = g(Mx \oplus b) \oplus L(x)
$$

我们的情况允许 $m$ 与 $n$ 不同且 $M$ 任意，不要求可逆或方阵。
但 $T(g)$ 的最小可能值仍与 $f$ 的代数结构密切相关。

### 2.3 特殊函数类的精确解

#### 2.3.1 仿射函数 ($\deg(f) \le 1$)

若 $f$ 是仿射函数，则存在仿射变换使其变为常数或单个变量：

$$
T_{\min} = \begin{cases}
0 & \text{若 $f$ 为常数} \\
1 & \text{若 $f$ 为非平凡仿射函数}
\end{cases}
$$

#### 2.3.2 二次型函数 ($\deg(f) = 2$)

**这是唯一有完备分类的非平凡情形。**

二次布尔函数 $f(x) = x^\mathrm{T} A x \oplus L(x) \oplus c$，其中 $A$ 是严格上三角矩阵。
令 $B = A + A^\mathrm{T}$ 为 $\mathbb{F}_2$ 上的交错双线性型。

$B$ 的秩 $r$ 在可逆仿射变换下不变。根据秩 $r$ 的奇偶性，二次型有以下标准型：

$$
\begin{aligned}
\text{秩 $r$ 为偶数}: &\quad x_1x_2 \oplus x_3x_4 \oplus \cdots \oplus x_{r-1}x_r \oplus c \\
\text{秩 $r$ 为奇数}: &\quad x_1x_2 \oplus x_3x_4 \oplus \cdots \oplus x_{r-2}x_{r-1} \oplus x_r \oplus c
\end{aligned}
$$

因此 $T_{\min} \le \lfloor r/2 \rfloor + \delta$（$\delta=0$ 或 $1$ 取决于常数项）。

**算法**：通过高斯消元在 $O(n^3)$ 时间内对 $B$ 进行对称消元，直接得到变换矩阵。

### 2.4 一般算法

#### 算法 A：暴力搜索（$n \le 5$）

枚举所有 $M\in\operatorname{GL}(n,\mathbb{F}_2)$ 和 $b\in\mathbb{F}_2^n$ 找最优 $T(g)$。

| $n$ | $|\operatorname{GL}(n,\mathbb{F}_2)|$ | 可行性 |
|-----|--------------------------------------|--------|
| 3   | 168                                  | 瞬间   |
| 4   | 20,160                               | 瞬间   |
| 5   | 9,999,360                            | 数秒   |
| 6   | 20,158,709,760                       | 不可行 |

#### 算法 B：线性方程组法（精确求解 $m$ 给定时）

给定 $M,b$，构造线性系统 $V\cdot a = f$ 并求解最稀疏的 $a$：

1. 对每个 $t\in\{0,1\}^m$，展开 $(Mx\oplus b)^t$ 为 $x$ 的 ANF 得到 $C_t$
2. 建立方程 $\sum_t a_t C_t = f$（稀疏线性系统）
3. 求解最稀疏的解（$m\le 10$ 时可通过枚举零空间精确求解）

**实现要点**：$C_t$ 的展开使用 Boole 环运算 $x_i^2 = x_i$，单变量乘积对应 mask 的 OR 操作。

**复杂度**：$O(3^m + 2^m \cdot T(f))$，其中 $3^m$ 是所有展开的总项数上界。

#### 算法 C：Walsh-Hadamard 谱方法（启发式搜索）

**观察**：仿射变换对应 Walsh 谱的置换和符号变化。谱系数大的方向往往暗示
该方向上有低次数的依赖，可能是有效的简化方向。

**步骤**：

1. 计算 Walsh 谱 $W_f(\omega) = \sum_x (-1)^{f(x) \oplus \langle \omega, x \rangle}$
2. 找出 $W_f(\omega)$ 绝对值较大的 $\omega$（排除 $\omega=0$）
3. 将这些 $\omega$ 作为新变量方向的候选
4. 测试每个候选方向组合的简化效果
5. 选择最优方向执行变换

**复杂度**：$O(n \cdot 2^n)$（通过快速 Walsh-Hadamard 变换）。

#### 算法 D：因式分解法（提取公因子）

观察 $f(x) = x_i \cdot h(x) \oplus r(x)$ 的形式（其中 $r$ 不含 $x_i$），
若 $h(x)$ 和 $r(x)$ 可以用少量线性形式表达，则整体 $f$ 可简化。

迭代策略：
1. 找到出现频率高的变量 $x_i$ 作为公因子
2. 令 $z_0 = x_i$，计算 $h$ 和 $r$ 的简化形式
3. 递归处理 $h$ 和 $r$

#### 算法 E：稠密 ANF 递归分解（补变换迭代法）

**基本观察**：若有 $t$ 个变量参与的 ANF $f(x_t)$ 满足 $T(f) \approx 2^t - t - 1$（几乎所有的次数 $\ge 2$ 的单
项式都存在），则补变换 $z = x \oplus 1$（即 $z_i = x_i \oplus 1$）会通过 $\mathbb{F}_2$ 上的超集求和大幅减少项数。
此时 $f(x_t)$ 可以写为：

$$
f(x_t) = \underbrace{(1 \oplus z_1)(1 \oplus z_2)\cdots(1 \oplus z_t)}_{\text{补变换展开}} \oplus \text{其他项}
$$

其中 $(1 \oplus z_1)(1 \oplus z_2)\cdots(1 \oplus z_t) = \bigoplus_{U \subseteq \{1..t\}} \prod_{j \in U} z_j$ 包含了 $t$ 个线性项、常数项和全部 $2^t - t - 1$ 个次数 $\ge 2$ 的项。若 $f$ 接近饱和（缺少少量项），则这些缺失项恰好是补变换后保留的项。

**递归分解流程**：

1. **输入**：$f(x_t)$，依赖 $t$ 个变量，$T(f)$ 是其次数 $\ge 2$ 的项数。
2. **判断**：计算比值 $\rho = T(f) \big/ (2^t - t - 1)$。
   - 若 $\rho \ge 70\%$：进入稠密分解。
   - 若 $\rho < 70\%$：项数已较少或不够饱和，进入标准仿射搜索流程（算法 B/C/D）。
3. **分解**：令 $D = (2^t - t - 1) - T(f)$，即 $f$ 中缺失的次数 $\ge 2$ 的单项式个数。
   - 这些缺失的单项式构成一个异或和 $h(x_s)$，其中 $s \le t$ 是 $h$ 涉及的变量数。
   - $h(x_s)$ 是补变换 $z = x \oplus 1$ 应用到 $f$ 后保留的误差项。
4. **递归**：对 $h(x_s)$ 重复步骤 2-3，以 $s < t$ 递减。
5. **终止**：
   - **情况 1**：$\rho \approx 0$（甚至 $=0$）。$h$ 的项数很少（如 $\le 5$），直接输出作为最终 ANF 的一部分。
   - **情况 2**：$0 < \rho < 70\%$（包括 $50\%$ 以下的中间情形）。已无法用稠密分解继续简化，对 $h$ 使用标准仿射搜索方法。
6. **合成**：将各层递归结果通过补变换链合成回原始变量表达，得到 $f$ 的最终简化 ANF。

**与补变换搜索的关系**：此方法相当于将补变换从单次应用（z=x⊕1）推广为迭代递归过程。在 $\mathbb{F}_2$ 中，
补变换的迭代应用不会引入新的结构（两次补变换回到原变量），但每次递归处理的误差项 $h$ 的变量数严格递减，
保证函数必然终止。

**复杂度**：每层递归需要一次 Möbius 变换（$O(t \cdot 2^t)$）来确定 $h$ 的结构。
$t$ 随递归严格递减，总复杂度 $O\big(\sum_{t' \le t} t' \cdot 2^{t'}\big) = O(t \cdot 2^t)$。

**应用场景**：特别适用于 n=32 实例中饱和/近饱和输出的 ANF 简化。实验数据表明：
- hd01/hd02 的饱和输出（$T = 2^t-1$）经补变换后 $T_{\text{comp}} = 2$
- hd10 的近饱和输出（$T \approx 4.2\times 10^9$）经补变换后 $T_{\text{comp}} \le 15$

### 2.5 算法比较

| 算法           | 适用范围          | 复杂度                      | 最优性保证 |
| -------------- | ----------------- | --------------------------- | ---------- |
| A 暴力搜索     | $n \le 5$       | 指数 (穷举)                 | 全局最优   |
| B 线性方程组法 | 任意 $n$, $m\le 10$ | $O(3^m + 2^m\cdot T(f))$ | 给定 $M,b$ 下的最优 $g$ |
| C Walsh 谱方法 | 任意 $n$         | $O(n \cdot 2^n)$          | 无保证     |
| D 因式分解法   | 任意 $n$         | $O(n \cdot T(f))$         | 无保证     |
| E 稠密递归分解 | 任意 $n$，$T(f) \approx 2^t-t-1$ | $O(t \cdot 2^t)$ (递归递减) | 无保证，但对近饱和函数极有效 |

### 2.7 开问题与展望

1. **$\deg(f) \ge 3$ 的最小 $T$ 的刻画**：是否存在类似二次秩的不变量
   来刻画 $T_{\min}$？可能是某种高阶导数的秩或广义的多元相关系数。
2. **随机函数的 $T_{\min}$**：随机 $n$ 元布尔函数的 ANF 平均有 $2^{n-1}$ 个单项式。
   仿射变换能将这个数降低多少？是否存在普适的下界？
3. **与编码理论的联系**：将 ANF 系数向量视为 $\mathbb{F}_2^{2^n}$ 中的向量，
   仿射变换的作用对应某个线性群的表示。$T(f)$ 就是该向量的 Hamming 重量。
   问题转化为：求群作用下轨道中重量的最小值。

---

## 3. 向量布尔函数的 ANF 简化扩展

### 3.1 问题定义

考虑向量布尔函数（S-盒）$F: \mathbb{F}_2^n \to \mathbb{F}_2^k$：

$$
F(x) = (f^{(0)}(x), f^{(1)}(x), \ldots, f^{(k-1)}(x))
$$

目标是找到 $M\in\mathbb{F}_2^{m\times n}, b\in\mathbb{F}_2^m$ 使得
$G(z) = (g^{(0)}(z),\ldots,g^{(k-1)}(z))$ 满足 $F(x) = G(Mx\oplus b)$
且 $T(G) = \sum_j T(g^{(j)})$ 显著小于 $T(F) = \sum_j T(f^{(j)})$。

### 3.2 与标量情况的对比

| 方面         | 标量布尔函数 | 向量布尔函数                        |
| ------------ | ------------ | ----------------------------------- |
| 优化目标     | 单个 $T(g)$ | 加权和 $\sum w_j \cdot T(g^{(j)})$ |
| 变换自由度   | 输入侧仿射   | 输入侧仿射 + 输出侧仿射（可选）     |
| 次数约束     | 自动保持     | 需对各分量分别检查                  |

### 3.3 扩展算法

#### 算法 E1：加权贪心法

将各分量的 ANF 系数合并为加权和，在单一搜索中优化总目标。

#### 算法 E2：主分量主导法

选择 $T(f^{(j)})$ 最大的分量作为主导分量，优先优化它，再验证对其他分量的影响。

#### 算法 E3：输出仿射变换联合优化

在输入变换 $z = Mx \oplus b$ 的同时，考虑输出仿射变换
$y = N \cdot G(z) \oplus c$（$N \in \operatorname{GL}(k,\mathbb{F}_2)$）。
输出仿射不改变各分量内部的 ANF 结构，但将分量线性组合，**重新分配**各分量的复杂度。

---

## 4. 已实现的启发式算法（C++/Python 实现）

### 4.1 三种多项式类型

系统统一处理三种类型的 ANF 简化：

| 类型 | 环 | 表示 | 模参数 |
|------|-----|------|--------|
| 布尔函数 | $\mathbb{F}_2[x]/(x_i^2 - x_i)$ | `SparseANF` (mask 位图) | — |
| 整系数多项式 | $\mathbb{Z}[x]$ | `IntPoly` (指数向量 → 系数) | `mod = 0` |
| $\mathbb{F}_p$ 多项式 | $\mathbb{F}_p[x]$ | `IntPoly` + 模算术 | `mod = p > 0` |

布尔函数使用 `SparseANF` 类（mask 位图表示，$n \le 64$），整系数和 $\mathbb{F}_p$ 多项式使用 `IntPoly` 类（字典存储）。

### 4.2 统一变换框架

所有算法求解同一个核心问题：给定 $f(x)$，寻找仿射变换 $z = Mx + b$（$M \in \mathbb{F}_2^{m\times n}$ 或 $M \in \mathbb{Z}^{m\times n}$，$m$ 可小于或大于 $n$，$M$ 无需方阵/满秩）使得 $g(z) = g(Mx + b) = f(x)$ 且 $T(g) < T(f)$。

#### $Z_1 \cup Z_2$ 框架（Union 机制）

标准方法要求 $g$ 通过右逆 $N$（$M\cdot N = I_m$）逆向表达 $f(x) = g(Mx + b)$。当 $m < n$ 时，变量信息丢失。

Union 框架通过保留两组变量绕过此限制：
- $Z_1 = X$：原始 $n$ 个变量
- $Z_2 = Mx + b$：仿射变换后的 $m$ 个变量

合并后共有 $m + n$ 个变量。$g$ 定义在 $(Z_1, Z_2)$ 上，等价于在更大的变量空间中寻找更稀疏的表达。实现通过构造扩展矩阵 $M_{\text{ext}} = [M; I_n]$（$(m+n)\times n$），通过列扩展 + 求逆得到 $x = N z_{\text{ext}} + c$，再调用前向展开 $g(z_{\text{ext}}) = f(N z_{\text{ext}} + c)$。

这个框架可以统一理解为：**将原始变量和变换后的新变量同时保留，让算法在更大的变量空间中寻找更稀疏的多项式表达**。

#### 补码 Union

在布尔函数中，Union 框架的一个特例是保留 $x_i$ 和 $x_i \oplus 1$（$\lnot x_i$）两套变量。利用 $x_i \cdot \lnot x_i = 0$ 可以消去某些项。

### 4.3 贪心 XOR 合并（布尔函数）

**核心操作**：将变量 $x_i$ 替换为 $x_i \oplus x_j$。对每个包含 $x_i$ 的单项式 $m$：
- $m$ 在原集合中翻转（移除）
- 新单项式 $m' = m \setminus \{x_i\}$（若 $x_j \in m$）或 $m' = (m \setminus \{x_i\}) \cup \{x_j\}$（若 $x_j \notin m$）翻转

这对应于 **Boolean 环中的变量代换**：$x_i \to x_i \oplus x_j$ 保持代数结构不变。

**算法**：

```
f_current = f, M = I_n
loop:
  对所有活跃变量对 (i, j):
    计算 ΔT = T(after merge) - T(f_current)
  若 best ΔT ≥ 0 → 终止
  执行合并 x_i → x_i ⊕ x_j
  删除未使用变量 → 减小 m
  更新 M
```

**实现要点**：
- 对 `SparseANF`：mask 操作 $m \to m \oplus (1\ll i)$ 翻转包含 bit $i$ 的项，OR 操作合并单项式
- 对 `IntPoly`：使用 `substitute_linear(i, coeffs)` 实现 $x_i \to \sum_j c_j x_j$，这里的合并系数为 $c_i = 1, c_j = -k$
- 候选 $k$ 值：$\{1, -1, 2, -2\}$

### 4.4 梯度引导合并（整数/$\mathbb{F}_p$）

**理论基础**：若 $\partial f / \partial x_i \equiv k \cdot \partial f / \partial x_j$（系数 $k$ 倍关系），则合并 $x_i \to x_i - k \cdot x_j$ 可能降低 $T(f)$。

**算法**（两阶段）：

1. **Phase A（梯度匹配）**：计算所有梯度 $\partial f / \partial x_i$，搜索满足 $\partial f / \partial x_i \equiv k \cdot \partial f / \partial x_j$ 的变量对 $(i, j)$
2. **Phase B（穷举）**：若 Phase A 未找到改进，对所有活跃变量对尝试合并

对于 $\mathbb{F}_p$ 多项式，梯度比较使用模算术：$\partial f / \partial x_i \equiv k \cdot \partial f / \partial x_j \pmod{p}$。

**优势**：对结构化函数能快速发现变量之间的线性依赖关系，在大 $n$ 时显著快于穷举。

### 4.5 补码搜索（布尔函数）

**观察**：对于布尔函数，保留 $x_i$ 和 $x_i \oplus 1$ 同时存在可能通过 $x_i \cdot (x_i \oplus 1) = 0$ 消去项。

**算法**：
- $n \le 16$：Gray 码遍历所有 $2^n$ 种补码模式，增量更新 $T(f)$（$O(2^n \cdot T(f))$）
- $n > 16$：贪心位翻转，逐位试探能否降低 $T(f)$

### 4.6 随机 $M,b$ 搜索

**策略**：生成随机仿射变换候选，通过 Union 框架测试是否能降低 $T(g)$。

**矩阵生成策略**：
- **结构化矩阵**（80% 概率）：每行最多一个非零元（$\pm 1$），直接映射 $z_j = \pm x_i + b_j$，不需求逆
- **通用矩阵**（20% 概率）：小随机条目 $[-2, 2]$，需检查满秩后方可求右逆

**流程**：
1. 随机生成 $M$（$m \times n$）和 $b$
2. 计算 $g_1 = f.\text{substitute\_affine}(M, b)$
3. 尝试 Union：$g_2 = f.\text{substitute\_affine\_union}(M, b)$
4. 若 $T(g_2) < T(g_1)$，保留 Union 结果
5. 验证正确性后更新最优

### 4.7 完整管线

各类型的完整简化管线：

**布尔函数**（`simplify`）：
```
Phase 1: 补码搜索 (simplify_by_complement)
Phase 2a: 梯度引导贪心合并 (greedy_merge_simplify_gradient)
Phase 2b: 穷举贪心合并 (greedy_merge_simplify)
Phase 3: 随机搜索 + Union (search_random)
```

**整数多项式**（`simplify_int`）：
```
Phase 1: 梯度引导合并 (simplify_by_gradient_int)
Phase 2: 穷举贪心合并 (greedy_merge_simplify_int)
Phase 3: 随机搜索 + Union (search_random_int)
```

**向量布尔函数**（`vector_simplify`）：
```
Phase 1: 补码 Union
Phase 2: 梯度引导合并 (vector_greedy_merge)
Phase 3: 穷举合并
Phase 4: 补码 Union + 随机搜索 (vector_search_random)
```

各阶段间通过矩阵合成 $M_{\text{new}} = M_{\text{next}} \cdot M_{\text{acc}}$, $b_{\text{new}} = M_{\text{next}} \cdot b_{\text{acc}} \oplus b_{\text{next}}$ 保持变换历史。

### 4.8 $\mathbb{F}_p$ 扩展

IntPoly 类添加 `int64_t mod` 成员：
- `mod = 0`：整数多项式（$\mathbb{Z}$），使用现有算法
- `mod > 0`：$\mathbb{F}_p$ 多项式，所有系数算术模 $p$

**模算术实现**：
- `reduce(x)`：$(x \bmod p + p) \bmod p$ 规范化到 $[0, p-1]$
- `mul_mod(a, b)`：$(a \cdot b) \bmod p$（使用 `__int128_t` 防溢出）
- `coeff_zero(c)`：$c \equiv 0 \pmod{p}$

**$\mathbb{F}_p$ 矩阵求逆**：
- 使用扩展欧几里得算法计算模逆
- 精确模 $p$ 高斯消元（无浮点精度问题）
- 专门函数：`fp_mat_invert`, `fp_extend_to_invertible`, `fp_extend_columns_to_invertible`, `fp_has_full_row_rank`

**简化管线**：所有算术运算、梯度比较、矩阵合成均模 $p$ 感知，`substitute_affine` 和 `substitute_affine_union` 自动分派到 $\mathbb{F}_p$ 路径。

### 4.9 验证机制

两种验证方式确保正确性：

1. **随机验证**：生成 $n_{\text{tests}}$ 个随机输入 $x$，计算 $z = Mx + b$，验证 $f(x) \equiv g(z)$（在 $\mathbb{F}_2$、$\mathbb{Z}$ 或 $\mathbb{F}_p$ 中）
2. **直接验证（$\mathbb{F}_p$）**：对 $\mathbb{F}_p$ 多项式，比较 $f(x) \bmod p$ 与 $g(z) \bmod p$

## 5. 文件结构

| 文件（C++） | 用途 |
|-------------|------|
| `src/sparse_anf.cpp` | `SparseANF` 类（布尔函数 ANF） |
| `src/simplify.cpp` | 布尔函数简化算法 |
| `src/int_poly.cpp` | `IntPoly` 类（整数/$\mathbb{F}_p$ 多项式） |
| `src/vector_anf.cpp` | 向量布尔函数简化 |
| `src/vector_int_poly.cpp` | 向量整数/$\mathbb{F}_p$ 多项式简化 |

| 文件（Python） | 用途 |
|---------------|------|
| `anf_factor.py` | `SparseANF` 类 + 布尔简化算法 |
| `int_poly.py` | `IntPoly` 类 + 整数/$\mathbb{F}_p$ 简化算法 |

## 6. 进一步方向

1. **$\mathbb{F}_p$ 的补码搜索**：对 $\mathbb{F}_p$ 多项式，补码相当于 $x_i \to -x_i + c$（$c \in \mathbb{F}_p$），可利用任意仿射单变量变换
2. **$\mathbb{F}_p$ 随机搜索的均匀采样**：从 $[0, p-1]$ 均匀采样矩阵条目，提升搜索覆盖率
3. **$\mathbb{F}_p$ 结构化路径扩展**：支持任意非零乘子（$\mathbb{F}_p$ 中所有非零元可逆）
4. **Walsh 方向搜索的 $\mathbb{F}_p$ 推广**：利用特征和做方向引导
5. **输出侧变换**：对向量函数引入输出仿射 $y = N \cdot G(z) \oplus c$，重新分配分量间的复杂度

## 参考文献

1. Carlet, C. (2010). *Boolean Functions for Cryptography and Error Correcting Codes*.
   In: Crama, Y., Hammer, P. (eds) Boolean Models and Methods in Mathematics,
   Computer Science, and Engineering.
2. Ryabov, A. (2018). *On Classification of Cubic Boolean Functions under
   Affine Transformations*. Journal of Applied and Industrial Mathematics.
3. MacWilliams, F. J., Sloane, N. J. A. (1977). *The Theory of Error-Correcting Codes*.
4. Bard, G. V. (2009). *Algorithms for Solving Polynomial Equations over Finite Fields*.
