# sum_T / union_T 统计算法说明

## .poly 文件中的 ANF 单项式矩阵

优化结果（`{策略}.poly`）包含所有输出的 ANF 单项式。将全部单项式排列为一个矩阵 W：

- 每行对应一个 ANF 单项式
- 前 m 列：z 变量选择位（mask），第 j 位为 1 表示该单项式包含 zⱼ
- 最后一列：系数（恒为 1，系数为 0 的项不写入文件）
- 行数 = 所有输出单项式数之和

## 统计方法

### sum_T

遍历 W 所有行，统计前 m 列中非 0 元素 ≥ 2 的行数（即 deg ≥ 2 的单项式总数）：

```
sum_T = 0
for each row i in W:
    if popcount(W[i][0..m-1]) >= 2:
        sum_T++
```

### union_T

先按全部 m+1 列去重（相同行只保留一条，不计入的文件修改），再统计 deg ≥ 2 行数：

```
W1 = deduplicate(W)    // 考虑全部 m+1 列，含系数
union_T = 0
for each row i in W1:
    if popcount(W1[i][0..m-1]) >= 2:
        union_T++
```

注意：
- sum_T ≥ union_T 恒成立（union 去重后只减不增）
- 去重时考虑全部 m+1 列（含系数列），但统计时只看前 m 列的 popcount
- 系数为 0 的项不写入文件，故实际所有行系数均为 1
