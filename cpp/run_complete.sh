#!/bin/bash
# ====================================================================
# run_complete.sh — 非 n=32 实例完整策略组合运行脚本
#
# 单输出 (k_nl=1):  5 种策略
#   d1a_opt  d1b_opt  d2_opt  d3_opt  d1c_opt
#
# 多输出 (k_nl>1):  7 种策略
#   d1a_opt1 d1a_opt2 d1b_opt2 d2_opt2 d3_opt1 d3_opt2 d1c_opt2
#
# 无效组合（不存在）：
#   d2_opt1  d1b_opt1  d1c_opt1
#
# n=32 实例参见 run_n32.sh
# ====================================================================
set -e

# Signal handling — Ctrl+C kills all child processes
trap 'echo ""; echo "Interrupted. Terminating..."; kill 0 2>/dev/null; exit 130' INT TERM

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
EXAMPLES_DIR="$SCRIPT_DIR/examples"
RESULTS_DIR="$SCRIPT_DIR/results"
PREPROCESS_DIR="$SCRIPT_DIR/preprocessed"

MODE="${1:-test}"

# Strategy timeout per run (seconds). Override: TIMEOUT=1800 ./run_complete.sh --full
TIMEOUT="${TIMEOUT:-86400}"

# 按实例×策略的超时配置（覆盖 $TIMEOUT）
[ -f "$SCRIPT_DIR/time.cfg" ] && source "$SCRIPT_DIR/time.cfg"

# ---- 参数集 ----
if [ "$MODE" = "--full" ]; then
    # 共享/单输出搜索（跑一次）→ 最大力度
    P_COMMON="--random 100000 --hill-climb 100000"
    # m≤6 枚举极快
    P_D2="--random 100000 --hill-climb 100000 --max-m 6"
    # 补码搜索只做 Phase 1b，不需要随机/爬山
    P_D1C="--walsh-k 0 --random 0 --hill-climb 0"

    # ---- per-output 搜索参数（按实例规模定制） ----
    # n=16 评估最慢，k=8 → 15000/15000
    # n=10-11 中等，k=7-11 → 50000/50000
    # n=7-8 评估极快，k=8-256 → 50000~100000/50000~100000
    declare -A OPT2_P  # d1a_opt2
    OPT2_P[hd03]="--random 15000 --hill-climb 15000 --max-m 12"
    OPT2_P[hd04]="--random 15000 --hill-climb 15000 --max-m 12"
    OPT2_P[int2float]="--random 50000 --hill-climb 50000 --max-m 12"
    OPT2_P[cavlc]="--random 50000 --hill-climb 50000 --max-m 12"
    OPT2_P[hd07]="--random 100000 --hill-climb 100000 --max-m 12"
    OPT2_P[ctrl]="--random 100000 --hill-climb 500 --max-m 12"
    OPT2_P[dec]="--random 50000 --hill-climb 200 --max-m 12"

    declare -A D3_P  # d3_opt2
    D3_P[hd03]="--random 15000 --hill-climb 15000"
    D3_P[hd04]="--random 15000 --hill-climb 15000"
    D3_P[int2float]="--random 50000 --hill-climb 50000"
    D3_P[cavlc]="--random 50000 --hill-climb 50000"
    D3_P[hd07]="--random 100000 --hill-climb 100000"
    D3_P[ctrl]="--random 100000 --hill-climb 500"
    D3_P[dec]="--random 50000 --hill-climb 200"

    P_DEFAULT_OPT2="--random 30000 --hill-climb 30000 --max-m 12"
elif [ "$MODE" = "--list" ]; then
    echo "=== Instance × Strategy Matrix (non-n32) ==="
    echo ""
    echo "单输出 (k_nl=1, 5 strategies):"
    echo "  hd08 → d1a_opt  d1b_opt  d2_opt  d3_opt  d1c_opt"
    echo ""
    echo "多输出 (k_nl>1, 7 strategies):"
    echo "  hd07 hd03 hd04 ctrl dec int2float cavlc"
    echo "  → d1a_opt1  d1a_opt2  d1b_opt2  d2_opt2  d3_opt1  d3_opt2  d1c_opt2"
    echo ""
    echo "n=32 实例参见 run_n32.sh"
    echo ""
    echo "---"
    echo "其他用法:"
    echo "  ./run_complete.sh              快速测试"
    echo "  ./run_complete.sh --full       完整运行"
    echo "  ./run_complete.sh --clean      清理 results/ 下除 results.md 外的所有文件"
    echo "  ./run_complete.sh --update-md  从数据文件重新生成 results.md"
    exit 0

elif [ "$MODE" = "--clean" ]; then
    echo "=== Cleaning results/ ==="
    for inst in hd08 hd07 hd03 hd04 ctrl dec int2float cavlc; do
        [ -d "$RESULTS_DIR/$inst" ] && trash "$RESULTS_DIR/$inst"/* 2>/dev/null
    done
    echo "  Done"
    exit 0

elif [ "$MODE" = "--update-md" ]; then
    echo "=== Regenerating results.md ==="
    python3 "$SCRIPT_DIR/update_results.py"
    exit $?
else
    # 快速测试
    P_COMMON="--random 20 --walsh-k 10 --hill-climb 2"
    P_D2="--random 20 --hill-climb 2 --max-m 4"
    P_D1C="--walsh-k 0 --random 0 --hill-climb 0"
    declare -A OPT2_P D3_P
    for inst in hd07 hd03 hd04 int2float cavlc ctrl dec; do
        OPT2_P[$inst]="--random 20 --hill-climb 2 --max-m 8"
        D3_P[$inst]="--random 20 --hill-climb 2"
    done
fi

# ---- 配色 ----
RED='\033[0;31m'; BLUE='\033[0;34m'; NC='\033[0m'

check_exe() {
    if [ ! -f "$BUILD_DIR/$1" ]; then
        echo -e "  ${RED}SKIP: $1 not found${NC}"
        return 1
    fi
    return 0
}

# ---- 运行单一策略 ----
# 参数: <实例名> <策略标签> <电路文件> [额外参数]
run_strat() {
    local inst=$1 strat=$2 circuit=$3 extra_args=$4
    local exe="optimize_anf_${strat}"
    local out_dir="$RESULTS_DIR/$inst"
    mkdir -p "$out_dir"

    if ! check_exe "$exe"; then
        return 1
    fi

    local is_gb=0
    case "$strat" in
        d1b_opt|d1b_opt2) is_gb=1 ;;
    esac

    local t="${STRAT_TIMEOUT[${inst}_${strat}]:-$TIMEOUT}"
    local logfile="$out_dir/${inst}_${strat}_run.log"

    # 删除旧输出文件
    rm -f "$out_dir/${inst}_${strat}.affine" "$out_dir/${inst}_${strat}.poly" "$out_dir/${inst}_${strat}_stats.txt"

    echo -e "  ${BLUE}[${strat}]${NC} $inst (timeout=${t}s)"

    # 运行程序，输出同时到终端和日志（不做任何解析，程序自己输出判断结果）
    set +e
    if [ "$is_gb" -eq 1 ]; then
        timeout "$((t + 60))" "$BUILD_DIR/$exe" "$circuit" $extra_args --time-budget "$t" --out-dir "$out_dir" 2>&1 | tee "$logfile"
    else
        timeout "$((t + 60))" "$BUILD_DIR/$exe" "$circuit" $extra_args --time-budget "$t" --save-results "$out_dir" 2>&1 | tee "$logfile"
    fi
    set -e
}

# ---- 单输出实例：5 种策略 ----
run_single_output() {
    local inst=$1 circuit=$2
    echo ""
    echo ">> $inst (5 strategies)"
    run_strat "$inst" "d1a_opt" "$circuit" "$P_COMMON"
    run_strat "$inst" "d1b_opt" "$circuit" ""
    run_strat "$inst" "d2_opt"  "$circuit" "$P_D2"
    run_strat "$inst" "d3_opt"  "$circuit" "$P_COMMON"
    run_strat "$inst" "d1c_opt" "$circuit" "$P_D1C"
}

# ---- 多输出实例：7 种策略 ----
run_multi_output() {
    local inst=$1 circuit=$2
    echo ""
    echo ">> $inst (7 strategies)"
    run_strat "$inst" "d1a_opt1" "$circuit" "$P_COMMON"
    run_strat "$inst" "d1a_opt2" "$circuit" "${OPT2_P[$inst]:-$P_DEFAULT_OPT2}"
    run_strat "$inst" "d1b_opt2" "$circuit" ""
    run_strat "$inst" "d2_opt2"  "$circuit" "$P_D2"
    run_strat "$inst" "d3_opt1"  "$circuit" "$P_COMMON"
    run_strat "$inst" "d3_opt2"  "$circuit" "${D3_P[$inst]:-$P_DEFAULT_OPT2}"
    run_strat "$inst" "d1c_opt2" "$circuit" "$P_D1C"
}

# ====================================================================
# 主流程
# ====================================================================
echo "============================================"
echo " ANF Optimization — Complete Strategy Matrix"
echo " Mode: ${MODE}"
echo " Instances: non-n32 (hd08 hd07 hd03 hd04 ctrl dec int2float cavlc)"
echo "============================================"
echo ""

# ---- 预处理 ----
echo "--- Preprocessing ---"
for inst in hd08 hd07 hd03 hd04 ctrl dec int2float cavlc; do
    mkdir -p "$PREPROCESS_DIR/$inst"
    "$BUILD_DIR/preprocess" "$EXAMPLES_DIR/${inst}.txt" --out-dir "$PREPROCESS_DIR/$inst" &> /dev/null && echo "  $inst ✓" || echo "  $inst ✗"
done
echo ""

# ---- 单输出 ----
echo "=== 单输出实例 (5 strategies) ==="
run_single_output "hd08" "$EXAMPLES_DIR/hd08.txt"

# ---- 多输出 ----
echo ""
echo "=== 多输出实例 (7 strategies) ==="
for inst in hd07 hd03 hd04 ctrl dec int2float cavlc; do
    run_multi_output "$inst" "$EXAMPLES_DIR/${inst}.txt"
done

exit 0
