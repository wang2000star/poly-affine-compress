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

### 2.5 算法比较

| 算法           | 适用范围          | 复杂度                      | 最优性保证 |
| -------------- | ----------------- | --------------------------- | ---------- |
| A 暴力搜索     | $n \le 5$       | 指数 (穷举)                 | 全局最优   |
| B 线性方程组法 | 任意 $n$, $m\le 10$ | $O(3^m + 2^m\cdot T(f))$ | 给定 $M,b$ 下的最优 $g$ |
| C Walsh 谱方法 | 任意 $n$         | $O(n \cdot 2^n)$          | 无保证     |
| D 因式分解法   | 任意 $n$         | $O(n \cdot T(f))$         | 无保证     |

### 2.6 开问题与展望

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

## 参考文献

1. Carlet, C. (2010). *Boolean Functions for Cryptography and Error Correcting Codes*.
   In: Crama, Y., Hammer, P. (eds) Boolean Models and Methods in Mathematics,
   Computer Science, and Engineering.
2. Ryabov, A. (2018). *On Classification of Cubic Boolean Functions under
   Affine Transformations*. Journal of Applied and Industrial Mathematics.
3. MacWilliams, F. J., Sloane, N. J. A. (1977). *The Theory of Error-Correcting Codes*.
4. Bard, G. V. (2009). *Algorithms for Solving Polynomial Equations over Finite Fields*.
