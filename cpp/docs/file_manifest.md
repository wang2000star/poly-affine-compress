# C++ 文件清单

所有 C++ 优化程序生成的输出文件。原始电路输入文件在 `examples/`，结果输出到 `results/{inst}/`。

## 输入文件

```
examples/cavlc.txt      10 输入   11 输出
examples/ctrl.txt        7 输入   26 输出
examples/dec.txt         8 输入  256 输出
examples/dsort.txt       ? 输入    ? 输出
examples/hd01.txt       32 输入   32 输出
examples/hd02.txt       32 输入   32 输出
examples/hd03.txt       16 输入    8 输出
examples/hd04.txt       16 输入    8 输出
examples/hd07.txt        8 输入    8 输出
examples/hd08.txt        8 输入    1 输出
examples/hd09.txt        ? 输入    ? 输出
examples/hd10.txt       32 输入   32 输出
examples/hd11.txt        ? 输入    ? 输出
examples/hd12.txt        ? 输入    ? 输出
examples/int2float.txt  11 输入    7 输出
examples/router.txt      ? 输入    ? 输出
```

## 输出文件模板

对每个实例 `{inst}`，输出目录 `results/{inst}/` 下可能生成以下文件：

### 原始指标（优化前基线）

```
results/{inst}/
  {inst}_original_anf.poly       ← 原始电路 ANF（仅 degree ≥ 2 项）
  {inst}_original_T.txt          ← 原始 T 值（每输出 T + sum_T）
```

### d1a_opt2 — 方向1 变体1a + 各自变换（optimize_anf.cpp）

```
results/{inst}/
  {inst}_d1a_opt2_trans.poly     ← 变换矩阵 z = Mx⊕b
  {inst}_d1a_opt2_anf.poly       ← g(z) ANF（degree ≥ 2 仅）
  {inst}_d1a_opt2_T.txt          ← 每输出 T 值 + sum_T
  {inst}_d1a_opt2_verify.txt     ← 随机点验证结果
```

### d1a_opt1 — 方向1 变体1a + 共享变换（optimize_anf.cpp）

```
results/{inst}/
  {inst}_d1a_opt1_trans.poly
  {inst}_d1a_opt1_anf.poly
  {inst}_d1a_opt1_T.txt
  {inst}_d1a_opt1_verify.txt
```

### d1b_opt2 — 方向1 变体1b + 各自变换（main_gate_builder.cpp）

```
results/{inst}/
  {inst}_d1b_opt2_trans.poly     ← 变换矩阵 z = Mx⊕b + c/d 行
  {inst}_d1b_opt2_anf.poly       ← g(z) ANF（c-correction 后，degree ≥ 2 仅）
  {inst}_d1b_opt2_T.txt          ← 每输出 T 值 + sum_T
  {inst}_d1b_opt2_verify.txt     ← 随机点验证结果
```

### d1b_opt1 — 方向1 变体1b + 共享变换（main_gate_builder.cpp，待实现）

```
results/{inst}/
  {inst}_d1b_opt1_trans.poly
  {inst}_d1b_opt1_anf.poly
  {inst}_d1b_opt1_T.txt
  {inst}_d1b_opt1_verify.txt
```

### d2_opt2 — 方向2 + 各自变换（待实现）

```
results/{inst}/
  {inst}_d2_opt2_trans.poly
  {inst}_d2_opt2_anf.poly
  {inst}_d2_opt2_T.txt
  {inst}_d2_opt2_verify.txt
```

### d3_opt2 — 方向3 + 各自变换（待实现）

```
results/{inst}/
  {inst}_d3_opt2_trans.poly
  {inst}_d3_opt2_anf.poly
  {inst}_d3_opt2_T.txt
  {inst}_d3_opt2_verify.txt
```

### d3_opt1 — 方向3 + 共享变换（待实现）

```
results/{inst}/
  {inst}_d3_opt1_trans.poly
  {inst}_d3_opt1_anf.poly
  {inst}_d3_opt1_T.txt
  {inst}_d3_opt1_verify.txt
```

### best — 逐输出取最优合并

```
results/{inst}/
  {inst}_best_trans.poly
  {inst}_best_anf.poly
  {inst}_best_T.txt
  {inst}_best_verify.txt
  {inst}_best_summary.txt       ← 各组合对比表
```

## 汇总

每个实例最大可能文件数：2（原始）+ 6 组合 × 4 + 5（best）= **31 个文件**。
当前 C++ 已实现的组合：d1a_opt2、d1b_opt2。
