#!/bin/bash
# ====================================================================
# run_complete.sh — 完整策略组合运行脚本
#
# 单输出 (k_nl=1):  5 种策略
#   d1a_opt  d1b_opt  d2_opt  d3_opt  d1c_opt
#
# 多输出 (k_nl>1):  7 种策略
#   d1a_opt1 d1a_opt2 d1b_opt2 d2_opt2 d3_opt1 d3_opt2 d1c_opt2
#
# 无效组合（不存在）：
#   d2_opt1  d1b_opt1  d1c_opt1
# ====================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
EXAMPLES_DIR="$SCRIPT_DIR/examples"
RESULTS_DIR="$SCRIPT_DIR/results"
PREPROCESS_DIR="$SCRIPT_DIR/preprocessed"

MODE="${1:-test}"

# Strategy timeout per run (seconds). Override: TIMEOUT=1800 ./run_complete.sh --full
TIMEOUT="${TIMEOUT:-3600}"

# ---- 参数集 ----
if [ "$MODE" = "--full" ]; then
    P_COMMON="--random 10000 --walsh-k 10000 --hill-climb 10000"
    P_N32="--random 10000 --walsh-k 10000 --n32-random 500 --hill-climb 10000"
    P_D1A_OPT2="--random 5000 --walsh-k 5000 --hill-climb 5000 --max-m 12"
    P_D2="--random 10000 --hill-climb 10000 --max-m 6"
    P_D1C="--walsh-k 0 --random 0 --hill-climb 0"
elif [ "$MODE" = "--list" ]; then
    echo "=== Instance × Strategy Matrix ==="
    echo ""
    echo "单输出 (k_nl=1, 5 strategies):"
    echo "  hd08 → d1a_opt  d1b_opt  d2_opt  d3_opt  d1c_opt"
    echo ""
    echo "多输出 (k_nl>1, 7 strategies):"
    echo "  hd07 hd03 hd04 ctrl dec int2float cavlc"
    echo "  → d1a_opt1  d1a_opt2  d1b_opt2  d2_opt2  d3_opt1  d3_opt2  d1c_opt2"
    echo ""
    echo "n=32 实例 (仅 d1a_opt2 / d3_opt2 / d1b_opt2):"
    echo "  hd10 hd01 hd02 hd09 hd11 hd12"
    echo ""
    echo "---"
    echo "其他用法:"
    echo "  ./run_complete.sh              快速测试"
    echo "  ./run_complete.sh --full       完整运行"
    echo "  ./run_complete.sh --clean      清理 results/ 下除 results.md 外的所有文件"
    echo "  ./run_complete.sh --update-md  从数据文件重新生成 results.md"
    exit 0

elif [ "$MODE" = "--clean" ]; then
    echo "=== Cleaning results/ (keeping results.md) ==="
    for inst in hd08 hd07 hd03 hd04 ctrl dec int2float cavlc hd10 hd01 hd02 hd09 hd11 hd12; do
        if [ -d "$RESULTS_DIR/$inst" ]; then
            find "$RESULTS_DIR/$inst" -type f ! -name 'results.md' -delete 2>/dev/null
            rmdir "$RESULTS_DIR/$inst" 2>/dev/null || true
        fi
    done
    echo "  Done (kept results/results.md)"
    exit 0

elif [ "$MODE" = "--update-md" ]; then
    echo "=== Regenerating results.md ==="
    python3 "$SCRIPT_DIR/update_results.py"
    exit $?
else
    # 快速测试
    P_COMMON="--random 20 --walsh-k 10 --hill-climb 2"
    P_N32="--random 10 --walsh-k 5 --n32-random 5 --hill-climb 1"
    P_D1A_OPT2="--random 20 --walsh-k 10 --hill-climb 2 --max-m 8"
    P_D2="--random 20 --hill-climb 2 --max-m 4"
    P_D1C="--walsh-k 0 --random 0 --hill-climb 0"
fi

# ---- 配色 ----
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
pass_cnt=0; fail_cnt=0; skip_cnt=0

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
        skip_cnt=$((skip_cnt + 1))
        return 1
    fi

    local is_gb=0
    case "$strat" in
        d1b_opt|d1b_opt2) is_gb=1 ;;
    esac

    local logfile="$out_dir/${inst}_${strat}_run.log"
    echo -e "  ${BLUE}[${strat}]${NC} $inst"

    set +e
    if [ "$is_gb" -eq 1 ]; then
        timeout "$TIMEOUT" "$BUILD_DIR/$exe" "$circuit" $extra_args --out-dir "$out_dir" &> "$logfile"
    else
        timeout "$TIMEOUT" "$BUILD_DIR/$exe" "$circuit" $extra_args --save-results "$out_dir" &> "$logfile"
    fi
    local rc=$?
    set -e

    if [ $rc -eq 0 ] && grep -q "Verified\|All outputs verified\|Done\.\|Phase" "$logfile" 2>/dev/null; then
        # Check for "no valid transform" in log — means search failed to find anything
        if grep -q "no valid transform" "$logfile" 2>/dev/null; then
            local nv=$(grep -c "no valid transform" "$logfile")
            echo -e "    ${RED}✗ ${nv} outputs have no valid transform${NC}"
            tail -3 "$logfile" | sed 's/^/      /'
            fail_cnt=$((fail_cnt + 1))
        else
            echo -e "    ${GREEN}✓ PASS${NC}"
            pass_cnt=$((pass_cnt + 1))
        fi
    elif [ $rc -ne 0 ]; then
        echo -e "    ${RED}✗ FAIL (exit $rc)${NC}"
        tail -5 "$logfile" | sed 's/^/      /'
        fail_cnt=$((fail_cnt + 1))
    else
        echo -e "    ${YELLOW}? CHECK (exit=0)${NC}"
        tail -3 "$logfile" | sed 's/^/      /'
        pass_cnt=$((pass_cnt + 1))
    fi
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
    run_strat "$inst" "d1a_opt2" "$circuit" "$P_D1A_OPT2"
    run_strat "$inst" "d1b_opt2" "$circuit" ""
    run_strat "$inst" "d2_opt2"  "$circuit" "$P_D2"
    run_strat "$inst" "d3_opt1"  "$circuit" "$P_COMMON"
    run_strat "$inst" "d3_opt2"  "$circuit" "$P_COMMON"
    run_strat "$inst" "d1c_opt2" "$circuit" "$P_D1C"
}

# ---- 验证 ----
run_verify() {
    local inst=$1
    local eqn="$PREPROCESS_DIR/$inst/$inst.eqn"
    [ -f "$eqn" ] || eqn="$EXAMPLES_DIR/${inst}.eqn"
    [ -f "$eqn" ] || return

    VERIFY_EXE="$BUILD_DIR/verify_anf"
    [ -f "$VERIFY_EXE" ] || return

    for aff in "$RESULTS_DIR/$inst"/*.affine; do
        [ -f "$aff" ] || continue
        base="${aff%.affine}"
        poly="${base}.poly"
        [ -f "$poly" ] || continue
        verify_out="${base}_verify.txt"
        # Re-verify regardless of existing verify file
        echo -e "    verify: $(basename $base)"
        "$VERIFY_EXE" "$eqn" "$aff" "$poly" 2000 --output "$verify_out" &> /dev/null && \
            echo -e "      ${GREEN}✓ PASS${NC}" || \
            echo -e "      ${RED}✗ FAIL${NC}"
    done
}

# ====================================================================
# 主流程
# ====================================================================
echo "============================================"
echo " ANF Optimization — Complete Strategy Matrix"
echo " Mode: ${MODE}"
echo "============================================"
echo ""

# ---- 预处理 ----
echo "--- Preprocessing ---"
for inst in hd08 hd07 hd03 hd04 ctrl dec int2float cavlc hd10 hd01 hd02 hd09 hd11 hd12; do
    mkdir -p "$PREPROCESS_DIR/$inst"
    # hd10/hd01/hd02 use .txt, others .txt as well
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

# ---- n=32 实例 ----
echo ""
echo "=== n=32 实例 ==="
for inst in hd10 hd01 hd02 hd09 hd11 hd12; do
    echo ""
    echo ">> $inst"
    run_strat "$inst" "d1a_opt2" "$EXAMPLES_DIR/${inst}.txt" "$P_N32"
    run_strat "$inst" "d1b_opt2" "$EXAMPLES_DIR/${inst}.txt" ""
    run_strat "$inst" "d3_opt2"  "$EXAMPLES_DIR/${inst}.txt" "$P_N32"
done

# ---- 验证 ----
echo ""
echo "=== 验证 f(x) = g(Mx+b) ==="
for inst in hd08 hd07 hd03 hd04 ctrl dec int2float cavlc hd10 hd01 hd02 hd09 hd11 hd12; do
    run_verify "$inst"
done

# ====================================================================
# 汇总
# ====================================================================
echo ""
echo "============================================"
echo -e " ${GREEN}${pass_cnt} passed${NC}, ${RED}${fail_cnt} failed${NC}, ${YELLOW}${skip_cnt} skipped${NC}"
echo "============================================"

echo ""
echo "--- Best union_T per instance ---"
for inst in hd08 hd07 hd03 hd04 ctrl dec int2float cavlc hd10 hd01 hd02 hd09 hd11 hd12; do
    best_ut=999999999
    best_st=0
    best_name=""
    for sf in "$RESULTS_DIR/$inst"/*_stats.txt; do
        [ -f "$sf" ] || continue
        ut=$(sed -n '4p' "$sf" 2>/dev/null | tr -d ' ')
        st=$(sed -n '3p' "$sf" 2>/dev/null | tr -d ' ')
        [[ "$sf" == *_raw_stats.txt ]] && continue
        name=$(basename "$sf" _stats.txt)
        if [ -n "$ut" ] && [ "$ut" -lt "$best_ut" ] 2>/dev/null; then
            best_ut=$ut; best_st=$st; best_name=$name
        fi
    done
    if [ "$best_ut" -ne 999999999 ]; then
        printf "  %-12s sum_T=%-8s union_T=%-5s (%s)\n" "$inst" "$best_st" "$best_ut" "$best_name"
    fi
done

exit $fail_cnt
