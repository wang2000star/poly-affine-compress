# C++ 文件清单

原始电路输入文件在 `examples/{inst}.txt`，结果输出分两类：
- `preprocessed/{inst}/` — 一次性预处理生成，固定不变
- `results/{inst}/` — 优化输出，每组合一套

## 输入文件

```
examples/cavlc.txt      10 输入   11 输出
examples/ctrl.txt        7 输入   26 输出
examples/dec.txt         8 输入  256 输出
examples/hd01.txt       32 输入   32 输出
examples/hd02.txt       32 输入   32 输出
examples/hd03.txt       16 输入    8 输出
examples/hd04.txt       16 输入    8 输出
examples/hd07.txt        8 输入    8 输出
examples/hd08.txt        8 输入    1 输出
examples/hd10.txt       32 输入   32 输出
examples/int2float.txt  11 输入    7 输出
```

## 预处理文件（preprocessed/{inst}/）

```
  {inst}.eqn             重命名电路（x_i, y_j, t_k）
  {inst}.bij             变量名映射（原始→x_i/y_j/t_k）
  {inst}_stats.txt       统计（5 行纯数字）
  {inst}_const.txt       常量输出
  {inst}_affine.txt      仿射输出
  {inst}_nonlinear.txt   非线性输出
  {inst}.poly            原始多项式（n≤16 仅）
```

## 优化输出（results/{inst}/）

每组合生成 4 文件，combo 取值：d1a_opt2, d1a_opt1, d1b_opt2, d1b_opt1, d2_opt2, d3_opt2, d3_opt1

```
  {inst}_{combo}.affine      变换矩阵（Z=[Mx+b; Cx+d]）
  {inst}_{combo}.poly        g(z) 多项式
  {inst}_{combo}_stats.txt   T 统计
  {inst}_{combo}_verify.txt  验证结果

  {inst}_best.affine         最优变换
  {inst}_best.poly           最优多项式
  {inst}_best_stats.txt      最优 T
  {inst}_best_verify.txt     最优验证
  {inst}_best_summary.txt    组合对比表
```
