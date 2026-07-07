# 超时处理设计

## 问题

当前 `timeout` 发送 SIGTERM 直接杀进程，搜索已经有最优解但没来得及保存就丢失了。

## 四种结果状态

任何实例×策略运行结束后，应能明确区分以下四种情况：

| # | 搜索是否结束 | 是否超时 | 是否有解 | 行为 |
|---|------------|---------|---------|------|
| ① | 否 | 是 | 否 | 输出"搜索未结束、超时、无解"，不保存 |
| ② | 否 | 是 | 是 | **保存当前最优解**，输出"搜索未结束、超时、有解" |
| ③ | 是 | 否 | 是 | 保存结果，输出"搜索结束、未超时、有解" |
| ④ | 是 | 否 | 否 | 不保存，输出"搜索结束、未超时、无解" |

## 方案

### C++ 侧

新增 `--time-budget N` 参数（秒）：
- 在搜索循环（随机搜索、爬山等）的每批候选后检查 `elapsed ≥ budget`
- 如果超时且当前已有候选（`best.union_T < INT64_MAX`），break 退出走正常保存路径
- 如果超时且无候选，也 break 退出，标记无结果
- 日志输出 `"time budget exhausted"` 字符串，供脚本侧检测

### 脚本侧（run_strat）

每次运行前，先写入一行摘要头：

```
# [inst] [strat] args="--random N --hill-climb N" timeout=600s
```

运行结束后，追加结果行：

```
# result: [search_status] [timeout_status] [solution_status]
```

搜索状态：`search_ended` / `search_not_ended`
超时状态：`timeout` / `no_timeout`
解状态：`has_solution` / `no_solution`

四种情况：
- ① → `search_not_ended timeout no_solution`
- ② → `search_not_ended timeout has_solution`
- ③ → `search_ended no_timeout has_solution`
- ④ → `search_ended no_timeout no_solution`

判断依据：
1. **exit code**：`timeout` 杀 → 124-137；内部 budget 退出或正常结束 → 0
2. **结果文件**：`.affine` 是否存在且非空
3. **C++ 日志关键词**：`"time budget exhausted"` 是否出现
