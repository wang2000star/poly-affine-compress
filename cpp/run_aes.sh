#!/bin/bash
# ====================================================================
# run_aes.sh — AES S-box ANF 优化（7 种策略）
#
# aes_bool: n=8, k=8, 1983 gates, ~1013 ANF terms
# 多输出实例，适用 7 种策略。
# ====================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
EXAMPLES_DIR="$SCRIPT_DIR/examples"
PREPROCESS_DIR="$SCRIPT_DIR/preprocessed"
RESULTS_DIR="$SCRIPT_DIR/results"
INST="aes_bool"

TIMEOUT="${TIMEOUT:-86400}"

# 按实例×策略的超时配置（覆盖 $TIMEOUT）
[ -f "$SCRIPT_DIR/time.cfg" ] && source "$SCRIPT_DIR/time.cfg"

# ---- 参数 ----
# n=8 评估极快 (2^8=256 穷举)，搜索可以很激进
P_COMMON="--random 100000 --hill-climb 100000"
P_D2="--random 100000 --hill-climb 100000 --max-m 12"
P_D1C="--walsh-k 100 --random 10000 --hill-climb 10000"
P_OPT2="--random 100000 --hill-climb 2000 --max-m 12"

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
echo " AES S-box ANF Optimization (7 strategies)"
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
