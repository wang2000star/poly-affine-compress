#!/bin/bash
# ====================================================================
# run_n32.sh — n=32 实例优化脚本
#
# 适用 3 种策略（per-output 搜索）：
#   d1a_opt2  d1b_opt2  d3_opt2
#
# 实例：hd10 hd01 hd02 hd09 hd11 hd12
#
# 其他实例参见 run_complete.sh
# ====================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
EXAMPLES_DIR="$SCRIPT_DIR/examples"
RESULTS_DIR="$SCRIPT_DIR/results"
PREPROCESS_DIR="$SCRIPT_DIR/preprocessed"

MODE="${1:-test}"

# Strategy timeout per run (seconds). Override: TIMEOUT=1800 ./run_n32.sh --full
TIMEOUT="${TIMEOUT:-86400}"

# 按实例×策略的超时配置（覆盖 $TIMEOUT）
[ -f "$SCRIPT_DIR/time.cfg" ] && source "$SCRIPT_DIR/time.cfg"

# ---- 参数集 ----
if [ "$MODE" = "--full" ]; then
    declare -A N32_P
    for inst in hd10 hd01 hd02 hd09 hd11 hd12; do
        N32_P[$inst]="--random 20000 --hill-climb 20000 --n32-random 1000"
    done
elif [ "$MODE" = "--list" ]; then
    echo "=== n=32 Instance × Strategy Matrix ==="
    echo ""
    echo "实例：hd10 hd01 hd02 hd09 hd11 hd12"
    echo "策略：d1a_opt2  d1b_opt2  d3_opt2"
    echo ""
    echo "---"
    echo "其他用法:"
    echo "  ./run_n32.sh              快速测试"
    echo "  ./run_n32.sh --full       完整运行"
    echo "  ./run_n32.sh --clean      清理 results/ 下 n=32 实例文件"
    exit 0

elif [ "$MODE" = "--clean" ]; then
    echo "=== Cleaning n=32 results/ ==="
    for inst in hd10 hd01 hd02 hd09 hd11 hd12; do
        [ -d "$RESULTS_DIR/$inst" ] && trash "$RESULTS_DIR/$inst"/* 2>/dev/null
    done
    echo "  Done"
    exit 0
else
    # 快速测试
    declare -A N32_P
    for inst in hd10 hd01 hd02 hd09 hd11 hd12; do
        N32_P[$inst]="--random 10 --n32-random 5 --hill-climb 1"
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

# ====================================================================
# 主流程
# ====================================================================
echo "============================================"
echo " ANF Optimization — n=32 Instances"
echo " Mode: ${MODE}"
echo "============================================"
echo ""

# ---- 预处理 ----
echo "--- Preprocessing ---"
for inst in hd10 hd01 hd02 hd09 hd11 hd12; do
    mkdir -p "$PREPROCESS_DIR/$inst"
    "$BUILD_DIR/preprocess" "$EXAMPLES_DIR/${inst}.txt" --out-dir "$PREPROCESS_DIR/$inst" &> /dev/null && echo "  $inst ✓" || echo "  $inst ✗"
done
echo ""

# ---- n=32 实例：各 3 种策略 ----
echo "=== n=32 实例 (3 strategies each) ==="
for inst in hd10 hd01 hd02 hd09 hd11 hd12; do
    echo ""
    echo ">> $inst"
    run_strat "$inst" "d1a_opt2" "$EXAMPLES_DIR/${inst}.txt" "${N32_P[$inst]}"
    run_strat "$inst" "d1b_opt2" "$EXAMPLES_DIR/${inst}.txt" ""
    run_strat "$inst" "d3_opt2"  "$EXAMPLES_DIR/${inst}.txt" "${N32_P[$inst]}"
done

exit 0
