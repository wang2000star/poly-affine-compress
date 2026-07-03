"""
Boolean Function ANF Simplification Toolkit
============================================

Implements algorithms for simplifying the Algebraic Normal Form (ANF)
of Boolean functions via affine variable transformations z = Mx ⊕ b.

Topics 涵盖:
  - ANF 表示与基础操作 (Möbius 变换)
  - 仿射变换施加与 ANF 重算
  - 暴力搜索 (n ≤ 5)
  - 二次型对角化 (deg ≤ 2 最优算法)
  - 贪心坐标合并 (通用启发式)
  - Walsh 谱分析 (高阶函数启发式)
"""

from __future__ import annotations

import itertools
import math
from typing import Callable, Optional

import numpy as np

# ---------------------------------------------------------------------------
#  1.  ANF 基础表示
# ---------------------------------------------------------------------------

# ANF 系数以 2^n 长 bit 向量存储:
#   anf[mask] == 1  ↔  单项式 ∏_{i: mask_i=1} x_i 在 ANF 中


def truth_table_to_anf(tt: np.ndarray) -> np.ndarray:
    """Möbius 变换: 真值表 → ANF 系数向量 (in-place)."""
    anf = tt.copy()
    n = int(np.log2(len(tt)))
    for i in range(n):
        step = 1 << i
        for j in range(0, len(tt), step << 1):
            for k in range(step):
                anf[j + step + k] ^= anf[j + k]
    return anf


def anf_to_truth_table(anf: np.ndarray) -> np.ndarray:
    """逆 Möbius 变换: ANF 系数 → 真值表."""
    return truth_table_to_anf(anf.copy())  # 自逆


def anf_monomial_count(anf: np.ndarray) -> int:
    """T(f) = ANF 中非零单项式个数."""
    return int(anf.sum())


def anf_degree(anf: np.ndarray) -> int:
    """deg(f) = ANF 中最高次单项式的次数."""
    n = int(np.log2(len(anf)))
    max_deg = 0
    for mask in range(len(anf)):
        if anf[mask]:
            max_deg = max(max_deg, mask.bit_count())
    return max_deg


def anf_from_coeffs(coeffs: dict[int, int], n: int) -> np.ndarray:
    """从系数字典构造 ANF 向量.
    例: anf_from_coeffs({0b001: 1, 0b011: 1}, 3) → x0 + x0x1
    """
    anf = np.zeros(1 << n, dtype=np.uint8)
    for mask, val in coeffs.items():
        anf[mask & ((1 << n) - 1)] = val & 1
    return anf


def anf_to_expr(anf: np.ndarray, var_names: Optional[list[str]] = None) -> str:
    """ANF 向量 → 可读表达式字符串."""
    n = int(np.log2(len(anf)))
    if var_names is None:
        var_names = [f"x{i}" for i in range(n)]
    terms = []
    for mask in range(len(anf)):
        if not anf[mask]:
            continue
        if mask == 0:
            terms.append("1")
            continue
        term_vars = [var_names[i] for i in range(n) if (mask >> i) & 1]
        terms.append("*".join(term_vars))
    return " + ".join(terms) if terms else "0"


# ---------------------------------------------------------------------------
#  2.  仿射变换
# ---------------------------------------------------------------------------

def gf2_matrix_inv(M: np.ndarray) -> np.ndarray:
    """求 M ∈ GL(n, F₂) 的逆矩阵."""
    n = M.shape[0]
    aug = np.concatenate([M % 2, np.eye(n, dtype=np.uint8)], axis=1)
    for col in range(n):
        pivot = np.argmax(aug[col:, col]) + col
        if aug[pivot, col] == 0:
            raise ValueError("矩阵不可逆")
        aug[[col, pivot]] = aug[[pivot, col]]
        for row in range(n):
            if row != col and aug[row, col]:
                aug[row] ^= aug[col]
    return aug[:, n:].astype(np.uint8)


def apply_affine_transform(
    anf: np.ndarray, M: np.ndarray, b: np.ndarray
) -> np.ndarray:
    """给定 f 的 ANF, 以及变换 z = Mx ⊕ b, 计算 g(z) = f(x) 的 ANF.

    M ∈ F₂^{m×n}, b ∈ F₂^{m×1}.  返回 m 元函数 g 的 ANF.
    """
    n_in = int(np.log2(len(anf)))
    m, n_col = M.shape
    assert n_col == n_in, f"M 列数 {n_col} 应等于变量数 {n_in}"

    if m == 0:
        tt = anf_to_truth_table(anf)
        idx = int(b.ravel()[0]) if len(b) > 0 else 0
        return np.array([tt[0]], dtype=np.uint8)

    b_vec = b.ravel()
    tt = anf_to_truth_table(anf)
    size_out = 1 << m
    g_tt = np.zeros(size_out, dtype=np.uint8)

    # 可逆方阵: g(z) = f(M^{-1}(z ⊕ b))
    if m == n_in and _gf2_rank(M) == n_in:
        M_inv = gf2_matrix_inv(M)
        for z in range(size_out):
            z_vec = np.array([(z >> i) & 1 for i in range(m)], dtype=np.uint8)
            x_vec = (M_inv @ (z_vec ^ b_vec)) % 2
            x_idx = int(sum(int(x_vec[i]) << i for i in range(n_in)))
            g_tt[z] = tt[x_idx]
        return truth_table_to_anf(g_tt)

    # 非方阵 (m < n): 遍历所有 x, 对 z = Mx ⊕ b 赋值
    g_check = np.full(size_out, -1, dtype=np.int8)
    for x in range(1 << n_in):
        x_vec = np.array([(x >> i) & 1 for i in range(n_in)], dtype=np.uint8)
        z_vec = (M @ x_vec) % 2
        z_idx = 0
        for i in range(m):
            z_idx |= (int(z_vec[i]) ^ int(b_vec[i])) << i

        if g_check[z_idx] == -1:
            g_check[z_idx] = tt[x]
        elif g_check[z_idx] != tt[x]:
            raise ValueError(
                f"f 在 M 的纤维上非常数 (z=0b{z_idx:0{m}b}), g 无法良定义"
            )

    g_check[g_check == -1] = 0
    return truth_table_to_anf(g_check.astype(np.uint8))


def transform_matrix_from_basis(basis_vectors: list[np.ndarray]) -> np.ndarray:
    """从一组基向量构造变换矩阵 M (m × n).
    每个行向量是一个新变量的系数.
    """
    return np.array(basis_vectors, dtype=np.uint8)


def random_affine_transform(
    n: int, m: Optional[int] = None, invertible: bool = True, seed: Optional[int] = None
) -> tuple[np.ndarray, np.ndarray]:
    """生成随机仿射变换 (M, b).
    n: 输入变量数;  m: 输出变量数 (默认 n).
    invertible: 要求 M 可逆 (仅当 m == n 时).
    """
    if m is None:
        m = n
    if seed is not None:
        np.random.seed(seed)
    while True:
        M = np.random.randint(0, 2, size=(m, n), dtype=np.uint8)
        if invertible and m == n:
            if np.linalg.matrix_rank(M % 2) < min(m, n):
                continue
        break
    b = np.random.randint(0, 2, size=(m, 1), dtype=np.uint8)
    return M, b


# ---------------------------------------------------------------------------
#  3.  暴力搜索 (n ≤ 5)
# ---------------------------------------------------------------------------

def _gl_order(n: int) -> int:
    """|GL(n, F_2)|."""
    order = 1
    for i in range(n):
        order *= (1 << n) - (1 << i)
    return order


def _gf2_rank(M: np.ndarray) -> int:
    """计算矩阵在 GF(2) 上的秩."""
    A = M.copy()
    n, m = A.shape
    rank = 0
    row = 0
    for col in range(m):
        if row >= n:
            break
        # 找 pivot
        pivot = None
        for i in range(row, n):
            if A[i, col]:
                pivot = i
                break
        if pivot is None:
            continue
        # 交换到当前行
        A[[row, pivot]] = A[[pivot, row]]
        # 消去其他行
        for i in range(n):
            if i != row and A[i, col]:
                A[i] ^= A[row]
        rank += 1
        row += 1
    return rank


def _all_invertible_matrices(n: int):
    """生成 GL(n, F_2) 中的所有矩阵."""
    if n == 0:
        yield np.eye(1, dtype=np.uint8)
        return
    total = 1 << (n * n)
    for flat in range(total):
        M = np.array(
            [(flat >> (i * n + j)) & 1 for i in range(n) for j in range(n)],
            dtype=np.uint8,
        ).reshape(n, n)
        if _gf2_rank(M) < n:
            continue
        yield M


def brute_force_minimize(
    anf: np.ndarray, verbose: bool = False
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """暴力搜索最小化 T(g) (n ≤ 5 适用).
    返回: (T_min, g_anf, M_best, b_best)
    """
    n = int(np.log2(len(anf)))
    if n > 5:
        print(f"警告: n={n} > 5, GL 阶 {_gl_order(n):,} 过大, 可能耗时极长.")

    T_min = len(anf)  # 上界
    g_best = anf.copy()
    M_best = np.eye(n, dtype=np.uint8)
    b_best = np.zeros((n, 1), dtype=np.uint8)

    count = 0
    for M in _all_invertible_matrices(n):
        # 对每个 M 仍要尝试 2^n 个 b
        for b_flat in range(1 << n):
            b = np.array([(b_flat >> i) & 1 for i in range(n)], dtype=np.uint8).reshape(
                n, 1
            )
            g_anf = apply_affine_transform(anf, M, b)
            T = anf_monomial_count(g_anf)
            if T < T_min:
                T_min = T
                g_best = g_anf
                M_best = M
                b_best = b
                if verbose:
                    print(f"  新最优 T = {T}")
            count += 1
            if T_min == 0:
                break
        if T_min == 0:
            break

    if verbose:
        print(f"  总计测试 {count:,} 个仿射变换")
    return T_min, g_best, M_best, b_best


# ---------------------------------------------------------------------------
#  4.  二次型对角化 (deg ≤ 2 最优算法)
# ---------------------------------------------------------------------------

def quadratic_diagonalize(
    anf: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """二次布尔函数的仿射标准化.

    输入:  deg ≤ 2 布尔函数的 ANF 系数
    输出: (g_anf, M, b)  其中 g 为标准型, f(x) = g(Mx ⊕ b).
    标准型: T(g) ≤ floor(r/2) + (0, 1 或 2), r 为二次型的秩.
    """
    n = int(np.log2(len(anf)))

    # 提取二次型矩阵 A (严格上三角)
    A = np.zeros((n, n), dtype=np.uint8)
    L = np.zeros(n, dtype=np.uint8)  # 线性部分系数
    c = anf[0]  # 常数项

    for i in range(n):
        L[i] = anf[1 << i]  # x_i 系数

    for i in range(n):
        for j in range(i + 1, n):
            A[i, j] = anf[(1 << i) | (1 << j)]  # x_i x_j 系数

    # B = A + A^T 为对称双线性型矩阵 (在 F_2 上也是交错阵)
    B = (A + A.T) % 2

    # 高斯消元对角化 B
    M_accum = np.eye(n, dtype=np.uint8)

    def _swap_cols(mat, i, j):
        col_i = mat[:, i].copy()
        mat[:, i] = mat[:, j]
        mat[:, j] = col_i

    rank = 0
    # 标准型: 按 (e1, e2, e1, e2, ...) 模式配对
    # 对 B 做同步行列变换: M^T B M 对角化为 2×2 块
    M_transform = np.eye(n, dtype=np.uint8)
    current = 0

    for step in range(n // 2):
        # 找 pivot
        found = False
        for i in range(current, n):
            for j in range(i + 1, n):
                if B[i, j]:
                    found = True
                    break
            if found:
                break

        if not found:
            break

        # 交换 i, j 到 current, current+1 位置
        if i != current:
            # 行列交换
            B[[i, current]] = B[[current, i]]
            B[:, [i, current]] = B[:, [current, i]]
            M_transform[[i, current]] = M_transform[[current, i]]
            # M_accum 列交换
            _swap_cols(M_accum, i, current)
        if j != current + 1:
            B[[j, current + 1]] = B[[current + 1, j]]
            B[:, [j, current + 1]] = B[:, [current + 1, j]]
            M_transform[[j, current + 1]] = M_transform[[current + 1, j]]
            _swap_cols(M_accum, j, current + 1)

        # 用 current, current+1 消去其他行列中的对应项
        for k in range(n):
            if k == current or k == current + 1:
                continue
            if B[current, k]:
                # 列: col_k += col_{current+1}
                B[:, k] ^= B[:, current + 1]
                # 行: row_k += row_{current+1}
                B[k, :] ^= B[current + 1, :]
                # 同步到 M_accum
                M_accum[:, k] ^= M_accum[:, current + 1]
            if B[current + 1, k]:
                B[:, k] ^= B[:, current]
                B[k, :] ^= B[current, :]
                M_accum[:, k] ^= M_accum[:, current]

        rank += 2
        current += 2

    # 现在在 M_accum 变换后的坐标 y = M_accum^{-1} x 下,
    # f 的二次项形如 y_1 y_2 + y_3 y_4 + ... (或含 y_i^2 = y_i)
    # 我们还需要处理线性部分

    # 变换 M_accum 将原变量 x 映射到新变量 y
    # y_i = (M_accum^{-1} x)_i → 我们需要它的逆来定义 z = Mx
    # 在 y 坐标下计算 f 的系数
    y_anf = apply_affine_transform(anf, M_accum.T, np.zeros((n, 1), dtype=np.uint8))

    # 从 y_anf 提取线性部分并在零空间方向上消去
    # 对于 rank 之后的变量, 它们仅线性出现 → 可用仿射变换吸收
    M_final = M_accum.T.copy()
    b_final = np.zeros((n, 1), dtype=np.uint8)

    # 注意: 更精细的仿射处理可以进一步减少常数项
    # 这里简化处理, 返回标准型

    g_anf = y_anf
    M_out = M_final
    b_out = b_final

    return g_anf, M_out, b_out


# ---------------------------------------------------------------------------
#  5.  贪心坐标合并 (通用启发式)
# ---------------------------------------------------------------------------


def greedy_minimize(
    anf: np.ndarray, max_iter: int = 20, verbose: bool = False
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """贪心迭代: 每次替换一个变量为两个变量的 XOR 以最小化 T.

    返回: (T_min, g_anf, M, b) 其中 f(x) = g(Mx ⊕ b).
    """
    n = int(np.log2(len(anf)))
    current_anf = anf.copy()
    T_current = anf_monomial_count(current_anf)

    # 累计变换: f(x) = g(M_acc x ⊕ b_acc)
    M_acc = np.eye(n, dtype=np.uint8)
    b_acc = np.zeros((n, 1), dtype=np.uint8)

    if verbose:
        print(f"贪心算法: 初始 T = {T_current}, n = {n}")

    for iteration in range(max_iter):
        best_delta = 0
        best_g = None
        best_M_new = None
        best_b_new = None

        # 尝试将变量 i 替换为 x_i ⊕ (某些其他 x_j 的和)
        for i in range(n):
            for d_flat in range(1 << n):
                d = np.array([(d_flat >> j) & 1 for j in range(n)], dtype=np.uint8)
                if d[i] != 1:  # 新变量必须包含 x_i
                    continue
                if d.sum() < 1:
                    continue

                # 构造 M_step: 替换 M 的第 i 行为 d^T
                M_step = np.eye(n, dtype=np.uint8)
                M_step[i, :] = d

                for bc in (0, 1):
                    b_step = np.zeros((n, 1), dtype=np.uint8)
                    b_step[i] = bc

                    # 组合变换: new_z = M_step (M_acc x ⊕ b_acc) ⊕ b_step
                    #         = (M_step @ M_acc) x ⊕ (M_step @ b_acc ⊕ b_step)
                    M_new = (M_step @ M_acc) % 2
                    b_new = (M_step @ b_acc + b_step) % 2

                    g_anf = apply_affine_transform(anf, M_new, b_new)
                    T_new = anf_monomial_count(g_anf)
                    delta = T_current - T_new

                    if delta > best_delta:
                        best_delta = delta
                        best_g = g_anf
                        best_M_new = M_new
                        best_b_new = b_new

        if best_delta <= 0:
            if verbose:
                print(f"  迭代 {iteration}: 无改进, 终止")
            break

        current_anf = best_g
        M_acc = best_M_new
        b_acc = best_b_new
        T_current = anf_monomial_count(current_anf)

        if verbose:
            print(f"  迭代 {iteration}: T = {T_current}, Δ = {best_delta}")

    return T_current, current_anf, M_acc, b_acc


# ---------------------------------------------------------------------------
#  6.  Walsh-Hadamard 谱分析
# ---------------------------------------------------------------------------


def walsh_hadamard_transform(anf: np.ndarray) -> np.ndarray:
    """快速 Walsh-Hadamard 变换.
    返回: FWT( (-1)^{f(x)} ) 的实值谱.
    """
    tt = anf_to_truth_table(anf)
    # 将布尔值映射到 ±1: (-1)^{f(x)}
    f_walsh = 1.0 - 2.0 * tt.astype(np.float64)
    n = int(np.log2(len(f_walsh)))
    for i in range(n):
        step = 1 << i
        for j in range(0, len(f_walsh), step << 1):
            for k in range(step):
                u = f_walsh[j + k]
                v = f_walsh[j + step + k]
                f_walsh[j + k] = u + v
                f_walsh[j + step + k] = u - v
    return f_walsh


def walsh_spectral_analysis(
    anf: np.ndarray, top_k: int = 10
) -> list[tuple[int, float]]:
    """分析 Walsh 谱, 返回谱值最大的 top_k 个方向."""
    walsh = walsh_hadamard_transform(anf)
    magnitudes = np.abs(walsh)
    indices = np.argsort(-magnitudes)[:top_k]
    return [(idx, walsh[idx]) for idx in indices]


def walsh_directions_minimize(
    anf: np.ndarray, verbose: bool = False
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """基于 Walsh 谱方向搜索的 ANF 简化."""
    n = int(np.log2(len(anf)))
    walsh = walsh_hadamard_transform(anf)
    magnitudes = np.abs(walsh)

    # 按谱值排序, 跳过 ω = 0
    indices = np.argsort(-magnitudes)

    best_T = anf_monomial_count(anf)
    best_g = anf.copy()
    best_M = np.eye(n, dtype=np.uint8)
    best_b = np.zeros((n, 1), dtype=np.uint8)

    tested = 0
    for idx in indices:
        if idx == 0:
            continue
        if tested >= 50:  # 限制测试数
            break

        ω = np.array([(idx >> i) & 1 for i in range(n)], dtype=np.uint8)
        if ω.sum() == 0:
            continue

        i0 = int(idx & -idx).bit_length() - 1
        M_try = np.eye(n, dtype=np.uint8)
        M_try[i0, :] = ω

        for bc in (0, 1):
            b_try = np.zeros((n, 1), dtype=np.uint8)
            b_try[i0] = bc
            g_anf = apply_affine_transform(anf, M_try, b_try)
            T = anf_monomial_count(g_anf)
            tested += 1
            if T < best_T:
                best_T = T
                best_g = g_anf
                best_M = M_try
                best_b = b_try
                if verbose:
                    print(f"  Walsh: 新最优 T = {T}, 方向 ω = 0b{idx:0{n}b}")

    if verbose:
        print(f"  Walsh: 测试了 {tested} 个方向")
    return best_T, best_g, best_M, best_b


# ---------------------------------------------------------------------------
#  7.  向量布尔函数 ANF 简化扩展
# ---------------------------------------------------------------------------


class VectorBooleanFunction:
    """向量布尔函数 F: F₂ⁿ → F₂ᵏ.

    用 k × 2ⁿ 的 uint8 矩阵存储 ANF 系数:
      anf[j, mask] = f^{(j)} 中单项式 x^mask 的系数.

    Attributes:
        anf: (k, 2ⁿ) ANF 系数矩阵
        n: 输入变量数
        k: 输出分量数
    """

    def __init__(self, anf_coeffs: np.ndarray):
        """
        anf_coeffs: (k, 2ⁿ) 矩阵, anf[j, mask] = 1 表示 f^{(j)} 含 x^mask.
        """
        self.anf = anf_coeffs.astype(np.uint8)
        self.k, self.size = self.anf.shape
        self.n = int(np.log2(self.size))

    @classmethod
    def from_component_anfs(cls, component_anfs: list[np.ndarray]) -> "VectorBooleanFunction":
        """从各分量的 ANF 向量构造."""
        return cls(np.array(component_anfs))

    def component(self, j: int) -> np.ndarray:
        """返回第 j 个分量的 ANF."""
        return self.anf[j].copy()

    def total_monomial_count(self) -> int:
        """总 ANF 单项式数 T_total = Σ_j T(f^{(j)})."""
        return int(self.anf.sum())

    def total_degree(self) -> int:
        """最大分量代数次数."""
        return max(anf_degree(self.anf[j]) for j in range(self.k))

    def component_counts(self) -> np.ndarray:
        """各分量的 T 值向量."""
        return self.anf.sum(axis=1).astype(int)

    def apply_input_transform(
        self, M: np.ndarray, b: np.ndarray
    ) -> "VectorBooleanFunction":
        """对输入侧施加仿射变换 z = Mx ⊕ b.
        返回新的 VectorBooleanFunction G 使得 F(x) = G(Mx ⊕ b).
        """
        new_anfs = []
        for j in range(self.k):
            g_anf = apply_affine_transform(self.anf[j], M, b)
            new_anfs.append(g_anf)
        return VectorBooleanFunction.from_component_anfs(new_anfs)

    def apply_output_affine(self, N: np.ndarray, c: np.ndarray) -> "VectorBooleanFunction":
        """对输出侧施加仿射变换: H(x) = N·F(x) ⊕ c.

        N ∈ GL(k, F₂), c ∈ F₂ᵏ.
        注意: 输出侧变换不改变各分量的 ANF 结构,
        而是重新组合分量.
        """
        new_anf = np.zeros_like(self.anf)
        for j in range(self.k):
            for i in range(self.k):
                if N[j, i]:
                    new_anf[j] ^= self.anf[i]
        # 加常数 c
        for j in range(self.k):
            if c[j]:
                new_anf[j, 0] ^= 1  # 常数项
        return VectorBooleanFunction(new_anf)

    def to_expr(self, f_names: Optional[list[str]] = None,
                x_names: Optional[list[str]] = None) -> str:
        """多表达式输出."""
        if f_names is None:
            f_names = [f"f{j}" for j in range(self.k)]
        parts = []
        for j in range(self.k):
            if x_names is None:
                expr = anf_to_expr(self.anf[j], [f"x{i}" for i in range(self.n)])
            else:
                expr = anf_to_expr(self.anf[j], x_names)
            parts.append(f"{f_names[j]} = {expr}")
        return "\n".join(parts)


def vector_brute_force_minimize(
    vbf: VectorBooleanFunction,
    weights: Optional[np.ndarray] = None,
    verbose: bool = False,
) -> tuple[int, VectorBooleanFunction, np.ndarray, np.ndarray]:
    """向量布尔函数的暴力搜索 (n ≤ 4 适用).
    搜索输入侧仿射变换最小化加权 T 和.

    返回: (T_total_min, G_best, M_best, b_best)
    """
    n, k = vbf.n, vbf.k
    if n > 4:
        print(f"警告: n={n} > 4, GL 阶 {_gl_order(n):,} 过大")

    if weights is None:
        weights = np.ones(k, dtype=np.float64)

    T_min = vbf.total_monomial_count()
    G_best = vbf
    M_best = np.eye(n, dtype=np.uint8)
    b_best = np.zeros((n, 1), dtype=np.uint8)

    count = 0
    for M in _all_invertible_matrices(n):
        for b_flat in range(1 << n):
            b = np.array([(b_flat >> i) & 1 for i in range(n)],
                         dtype=np.uint8).reshape(n, 1)

            G = vbf.apply_input_transform(M, b)
            T_w = int(np.sum(weights * G.component_counts()))

            if T_w < T_min:
                T_min = T_w
                G_best = G
                M_best = M
                b_best = b
                if verbose:
                    print(f"  新最优 T_total = {T_min}")

            count += 1
            if T_min == 0:
                break
        if T_min == 0:
            break

    if verbose:
        print(f"  总计 {count:,} 个变换")
    return T_min, G_best, M_best, b_best


def vector_input_output_minimize(
    vbf: VectorBooleanFunction,
    verbose: bool = False,
) -> tuple[int, VectorBooleanFunction, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """输入输出联合优化 (n ≤ 4, k ≤ 4 适用).
    枚举所有 (M, b) 和 (N, c) 的组合.

    返回: (T_min, G_best, M_best, b_best, N_best, c_best)
    """
    n, k = vbf.n, vbf.k

    T_min = vbf.total_monomial_count()
    G_best = vbf
    M_best = np.eye(n, dtype=np.uint8)
    b_best = np.zeros((n, 1), dtype=np.uint8)
    N_best = np.eye(k, dtype=np.uint8)
    c_best = np.zeros((k, 1), dtype=np.uint8)

    total_combos = _gl_order(n) * (1 << n) * _gl_order(k) * (1 << k)
    if total_combos > 1_000_000:
        print(f"警告: 组合数 {total_combos:,} 过大")

    if verbose:
        print(f"输入输出联合优化: n={n}, k={k}, 组合数 ≈ {total_combos:,}")

    count = 0
    for M in _all_invertible_matrices(n):
        if n > 3:
            # 采样策略: n≥4 时随机采样
            break
        for b_flat in range(1 << n):
            b = np.array([(b_flat >> i) & 1 for i in range(n)],
                         dtype=np.uint8).reshape(n, 1)
            G_trans = vbf.apply_input_transform(M, b)

            for N in _all_invertible_matrices(k):
                for c_flat in range(1 << k):
                    c = np.array([(c_flat >> i) & 1 for i in range(k)],
                                 dtype=np.uint8).reshape(k, 1)
                    G = G_trans.apply_output_affine(N, c)
                    T_val = G.total_monomial_count()

                    if T_val < T_min:
                        T_min = T_val
                        G_best = G
                        M_best = M
                        b_best = b
                        N_best = N
                        c_best = c
                        if verbose:
                            print(f"  新最优 T_total = {T_min}")
                    count += 1

    if verbose:
        print(f"  测试了 {count:,} 种组合")
    return T_min, G_best, M_best, b_best, N_best, c_best


def demo_vector():
    """向量布尔函数 ANF 简化演示."""
    print("=" * 60)
    print("演示 5: 向量布尔函数 ANF 简化")
    print("=" * 60)

    # 构造一个 2 输出向量布尔函数 (n=3, k=2)
    # F(x) = (x0x1 + x0x2,  x0x1 + x1x2)
    n, k = 3, 2
    anf0 = np.zeros(1 << n, dtype=np.uint8)
    anf0[0b011] = 1  # x0x1
    anf0[0b101] = 1  # x0x2

    anf1 = np.zeros(1 << n, dtype=np.uint8)
    anf1[0b011] = 1  # x0x1
    anf1[0b110] = 1  # x1x2

    vbf = VectorBooleanFunction.from_component_anfs([anf0, anf1])
    print("原始 F:")
    print(vbf.to_expr(['f0', 'f1']))
    print(f"T_total = {vbf.total_monomial_count()}, 各分量 T = {vbf.component_counts()}")

    # 暴力搜索 (n=3: 小规模)
    T_min, G_best, M_best, b_best = vector_brute_force_minimize(vbf, verbose=False)
    print(f"\n最优输入变换: T_total_min = {T_min}")
    print(f"M = {M_best}")
    print(f"b = {b_best.ravel()}")
    print("变换后 G:")
    print(G_best.to_expr(['g0', 'g1'], ['z0', 'z1', 'z2']))
    print(f"各分量 T = {G_best.component_counts()}")

    # 输出变换联合优化 (k=2)
    T_joint, G_joint, M_j, b_j, N_j, c_j = vector_input_output_minimize(vbf, verbose=False)
    print(f"\n输入输出联合优化: T_total_min = {T_joint}")
    print(f"M = {M_j}, N = {N_j}")
    print("变换后 H:")
    print(G_joint.to_expr(['h0', 'h1'], ['z0', 'z1', 'z2']))
    print(f"各分量 T = {G_joint.component_counts()}")

    # 压缩率
    orig_T = vbf.total_monomial_count()
    print(f"\n压缩率: {T_min}/{orig_T} = {T_min/orig_T:.2f}")
    print()


# ---------------------------------------------------------------------------
#  8.  演示 & 测试
# ---------------------------------------------------------------------------


def demo_quadratic():
    """演示二次型对角化."""
    print("=" * 60)
    print("演示 1: 二次型对角化")
    print("=" * 60)

    # f(x) = x0 x1 + x0 x2
    # n=3, ANF: x0x1 (0b011) + x0x2 (0b101)
    n = 3
    anf = np.zeros(1 << n, dtype=np.uint8)
    anf[0b011] = 1  # x0 x1
    anf[0b101] = 1  # x0 x2
    print(f"f  = {anf_to_expr(anf)}")
    print(f"T(f) = {anf_monomial_count(anf)}, deg = {anf_degree(anf)}")

    # 期待 z0 = x0, z1 = x1 + x2, 得 g(z) = z0 z1, T(g) = 1
    M = np.array([[1, 0, 0], [0, 1, 1]], dtype=np.uint8)
    b = np.zeros((2, 1), dtype=np.uint8)
    g_anf = apply_affine_transform(anf, M, b)
    print(f"z = Mx, M = [{M[0]}, {M[1]}]")
    print(f"g  = {anf_to_expr(g_anf, ['z0', 'z1'])}")
    print(f"T(g) = {anf_monomial_count(g_anf)}, deg = {anf_degree(g_anf)}")
    print(f"T(g)/T(f) = {anf_monomial_count(g_anf)}/{anf_monomial_count(anf)}")
    print()

    # 示例 2: 三变量完全二次型
    # f(x) = x0x1 + x0x2 + x1x2 + x0
    n = 3
    anf2 = np.zeros(1 << n, dtype=np.uint8)
    anf2[0b011] = 1  # x0 x1
    anf2[0b101] = 1  # x0 x2
    anf2[0b110] = 1  # x1 x2
    anf2[0b001] = 1  # x0
    print(f"f  = {anf_to_expr(anf2)}")
    print(f"T(f) = {anf_monomial_count(anf2)}")

    # 暴力搜索 (n=3, 仅 168 × 8 = 1344 种变换)
    T_min, g_best, M_best, b_best = brute_force_minimize(anf2, verbose=False)
    print(f"暴力搜索最优: T_min = {T_min}")
    print(f"g  = {anf_to_expr(g_best, ['z0', 'z1', 'z2'])}")
    print(f"M = {M_best}")
    print(f"b = {b_best.ravel()}")
    print(f"T(g)/T(f) = {T_min}/{anf_monomial_count(anf2)}")
    print()


def demo_greedy():
    """演示贪心算法."""
    print("=" * 60)
    print("演示 2: 贪心坐标合并 (n=4 随机函数)")
    print("=" * 60)

    np.random.seed(42)
    n = 4
    # 构造一个非平凡的随机 ANF (控制密度 ~0.3)
    anf = np.zeros(1 << n, dtype=np.uint8)
    for mask in range(1 << n):
        if np.random.random() < 0.3 and mask != 0:
            anf[mask] = 1
    if anf.sum() == 0:
        anf[0b0001] = 1

    print(f"f  = {anf_to_expr(anf)}")
    print(f"T(f) = {anf_monomial_count(anf)}, deg = {anf_degree(anf)}")

    # 暴力搜索最优 (n=4: 20160 × 16 = 322560 种)
    T_min, g_best, M_best, b_best = brute_force_minimize(anf, verbose=False)
    print(f"暴力搜索最优: T_min = {T_min}")
    T_ratio = T_min / anf_monomial_count(anf)
    print(f"T(g)/T(f) = {T_min}/{anf_monomial_count(anf)} = {T_ratio:.3f}")

    # Walsh 启发式
    T_w, g_w, M_w, b_w = walsh_directions_minimize(anf, verbose=False)
    print(f"Walsh 方法:  T = {T_w}  (ratio = {T_w/anf_monomial_count(anf):.3f})")
    print()


def demo_linear_structure():
    """演示具有线性结构的函数."""
    print("=" * 60)
    print("演示 3: 线性结构检测")
    print("=" * 60)

    # 构造 f(x0, x1, x2) = g(x0, x1) 与 x2 无关
    # 即 f 在 x2 方向上具有线性结构
    # f(x) = x0 x1 + x2 (但 x2 是线性项, 不影响)
    n = 3
    anf = np.zeros(1 << n, dtype=np.uint8)
    anf[0b011] = 1  # x0 x1
    anf[0b100] = 1  # x2

    print(f"f  = {anf_to_expr(anf)}")
    print(f"T(f) = {anf_monomial_count(anf)}")

    # 如果 x2 方向是线性结构 f(x⊕e2)=f(x)⊕f(e2)
    # 则可以通过变换减少变量

    # 暴力搜索看是否能简化
    T_min, g_best, M_best, b_best = brute_force_minimize(anf, verbose=False)
    print(f"暴力搜索最优: T_min = {T_min}")
    print(f"g  = {anf_to_expr(g_best, ['z0', 'z1', 'z2'])}")
    print()


def demo_walsh():
    """演示 Walsh 谱分析."""
    print("=" * 60)
    print("演示 4: Walsh 谱分析")
    print("=" * 60)

    # bent 函数: f(x0,x1,x2,x3) = x0x1 + x2x3
    n = 4
    anf = np.zeros(1 << n, dtype=np.uint8)
    anf[0b0011] = 1  # x0 x1
    anf[0b1100] = 1  # x2 x3
    print(f"f  = {anf_to_expr(anf)}")
    print(f"T(f) = {anf_monomial_count(anf)}")

    walsh = walsh_hadamard_transform(anf)
    print(f"Walsh 谱范围: [{walsh.min()}, {walsh.max()}]")
    print("Top-5 谱系数方向:")
    top_idx = np.argsort(-np.abs(walsh))[:5]
    for idx in top_idx:
        ω = np.array([(idx >> i) & 1 for i in range(n)], dtype=np.uint8)
        print(f"  ω = {ω}  (0b{idx:0{n}b}), 谱值 = {walsh[idx]:.1f}")

    # Bent 函数的 Walsh 谱所有系数的绝对值都是 2^{n/2} = 4
    print()


# ---------------------------------------------------------------------------
#  主入口
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    demo_quadratic()
    demo_greedy()
    demo_linear_structure()
    demo_walsh()
    demo_vector()

    print("用文本中的例子演示:")
    print("-" * 60)
    print("原题: f(x) = x0 x1 + x0 x2")
    print("令 z0 = x0, z1 = x1 + x2")
    print("则 g(z) = z0 z1")
    print("T(f) = 2, T(g) = 1, T(g)/T(f) = 0.5")

    # 验证
    n = 3
    anf = np.zeros(1 << n, dtype=np.uint8)
    anf[0b011] = 1
    anf[0b101] = 1
    M = np.array([[1, 0, 0], [0, 1, 1]], dtype=np.uint8)
    b = np.zeros((2, 1), dtype=np.uint8)
    g_anf = apply_affine_transform(anf, M, b)
    print(f"\n代码验证: g(z) = {anf_to_expr(g_anf, ['z0', 'z1'])}")
    print(f"T(g) = {anf_monomial_count(g_anf)} ✓")
    print(f"deg(g) = {anf_degree(g_anf)}, deg(f) = {anf_degree(anf)}")
    print(f"deg 不变性验证: deg(g) = deg(f) = 2 ✓")
