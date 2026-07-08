#!/bin/bash
# ====================================================================
# run_aes.sh — AES S-box ANF 优化（7 种策略）
#
# aes_bool: n=8, k=8, 1983 gates, ~1013 ANF terms
# 多输出实例，适用 7 种策略。
# ====================================================================
set -e

# Signal handling — Ctrl+C kills all child processes
trap 'echo ""; echo "Interrupted. Terminating..."; kill 0 2>/dev/null; sleep 1; kill -9 0 2>/dev/null; exit 130' INT TERM

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
EXAMPLES_DIR="$SCRIPT_DIR/examples"
PREPROCESS_DIR="$SCRIPT_DIR/preprocessed"
RESULTS_DIR="$SCRIPT_DIR/results"
INST="aes_bool"

MODE="${1:-test}"

TIMEOUT="${TIMEOUT:-86400}"

# 按实例×策略的超时配置（覆盖 $TIMEOUT）
[ -f "$SCRIPT_DIR/time.cfg" ] && source "$SCRIPT_DIR/time.cfg"

# ---- 参数集 ----
if [ "$MODE" = "--full" ]; then
    # n=8 评估极快 (2^8=256 穷举)，搜索可以很激进
    P_COMMON="--random 100000 --hill-climb 100000"
    P_D2="--random 100000 --hill-climb 100000 --max-m 12"
    P_D1C="--walsh-k 100 --random 10000 --hill-climb 10000"
    P_OPT2="--random 100000 --hill-climb 2000 --max-m 12"
elif [ "$MODE" = "--list" ]; then
    echo "=== AES S-box Strategy Matrix ==="
    echo ""
    echo "aes_bool: n=8, k=8, 1983 gates"
    echo "7 strategies: d1a_opt1 d1a_opt2 d1b_opt2 d2_opt2 d3_opt1 d3_opt2 d1c_opt2"
    echo ""
    echo "---"
    echo "  ./run_aes.sh              快速测试"
    echo "  ./run_aes.sh --full       完整运行"
    echo "  ./run_aes.sh --clean      清理 results/aes_bool/ 下文件"
    exit 0
elif [ "$MODE" = "--clean" ]; then
    echo "=== Cleaning results/aes_bool/ ==="
    [ -d "$RESULTS_DIR/$INST" ] && trash "$RESULTS_DIR/$INST"/* 2>/dev/null
    echo "  Done"
    exit 0
else
    # 快速测试
    P_COMMON="--random 20 --walsh-k 10 --hill-climb 2"
    P_D2="--random 20 --hill-climb 2 --max-m 4"
    P_D1C="--walsh-k 10 --random 10 --hill-climb 2"
    P_OPT2="--random 20 --hill-climb 2 --max-m 8"
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

run_strat() {
    local strat=$1 circuit=$2 extra_args=$3
    local exe="optimize_anf_${strat}"
    local out_dir="$RESULTS_DIR/$INST"
    mkdir -p "$out_dir"

    if ! check_exe "$exe"; then
        return 1
    fi

    local is_gb=0
    case "$strat" in
        d1b_opt|d1b_opt2) is_gb=1 ;;
    esac

    local t="${STRAT_TIMEOUT[${INST}_${strat}]:-$TIMEOUT}"
    local logfile="$out_dir/${INST}_${strat}_run.log"

    # 删除旧输出文件
    rm -f "$out_dir/${INST}_${strat}.affine" "$out_dir/${INST}_${strat}.poly" "$out_dir/${INST}_${strat}_stats.txt"

    echo -e "  ${BLUE}[${strat}]${NC} $INST (timeout=${t}s)"

    # 运行程序，输出同时到终端和日志（不做任何解析，程序自己输出判断结果）
    set +e
    if [ "$is_gb" -eq 1 ]; then
        timeout --foreground "$((t + 60))" "$BUILD_DIR/$exe" "$circuit" $extra_args --time-budget "$t" --out-dir "$out_dir" 2>&1 | tee "$logfile"
    else
        timeout --foreground "$((t + 60))" "$BUILD_DIR/$exe" "$circuit" $extra_args --time-budget "$t" --save-results "$out_dir" 2>&1 | tee "$logfile"
    fi
    set -e
}

# ====================================================================
# 主流程
# ====================================================================
echo "============================================"
echo " AES S-box ANF Optimization (7 strategies)"
echo " Mode: ${MODE}"
echo "============================================"
echo ""

# ---- 预处理 ----
echo "--- Preprocessing ---"
mkdir -p "$PREPROCESS_DIR/$INST"
if "$BUILD_DIR/preprocess" "$EXAMPLES_DIR/${INST}.txt" --out-dir "$PREPROCESS_DIR/$INST" &> /dev/null; then
    echo "  $INST ✓"
else
    echo "  $INST ✗"
    exit 1
fi
echo ""

# ---- 7 种策略 ----
echo "=== Running 7 strategies ==="
echo ""

echo ">> $INST (7 strategies)"
run_strat "d1a_opt1" "$EXAMPLES_DIR/${INST}.txt" "$P_COMMON"
run_strat "d1a_opt2" "$EXAMPLES_DIR/${INST}.txt" "$P_OPT2"
run_strat "d1b_opt2" "$EXAMPLES_DIR/${INST}.txt" ""
run_strat "d2_opt2"  "$EXAMPLES_DIR/${INST}.txt" "$P_D2"
run_strat "d3_opt1"  "$EXAMPLES_DIR/${INST}.txt" "$P_COMMON"
run_strat "d3_opt2"  "$EXAMPLES_DIR/${INST}.txt" "$P_OPT2"
run_strat "d1c_opt2" "$EXAMPLES_DIR/${INST}.txt" "$P_D1C"

exit 0
