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

# ---- 参数 ----
# n=8 评估极快 (2^8=256 穷举)，搜索可以很激进
P_COMMON="--random 100000 --hill-climb 100000"
P_D2="--random 100000 --hill-climb 100000 --max-m 6"
P_D1C="--walsh-k 0 --random 0 --hill-climb 0"
P_OPT2="--random 100000 --hill-climb 2000 --max-m 12"

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

run_strat() {
    local strat=$1 circuit=$2 extra_args=$3
    local exe="optimize_anf_${strat}"
    local out_dir="$RESULTS_DIR/$INST"
    mkdir -p "$out_dir"

    if ! check_exe "$exe"; then
        skip_cnt=$((skip_cnt + 1))
        return 1
    fi

    local is_gb=0
    case "$strat" in
        d1b_opt|d1b_opt2) is_gb=1 ;;
    esac

    local logfile="$out_dir/${INST}_${strat}_run.log"
    echo -e "  ${BLUE}[${strat}]${NC} $INST"

    set +e
    if [ "$is_gb" -eq 1 ]; then
        timeout "$TIMEOUT" "$BUILD_DIR/$exe" "$circuit" $extra_args --out-dir "$out_dir" &> "$logfile"
    else
        timeout "$TIMEOUT" "$BUILD_DIR/$exe" "$circuit" $extra_args --save-results "$out_dir" &> "$logfile"
    fi
    local rc=$?
    set -e

    if [ $rc -eq 0 ] && grep -q "Verified\|All outputs verified\|Done\.\|Phase" "$logfile" 2>/dev/null; then
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

run_verify() {
    local eqn="$PREPROCESS_DIR/$INST/$INST.eqn"
    [ -f "$eqn" ] || return

    local VERIFY_EXE="$BUILD_DIR/verify_anf"
    [ -f "$VERIFY_EXE" ] || return

    for aff in "$RESULTS_DIR/$INST"/*.affine; do
        [ -f "$aff" ] || continue
        base="${aff%.affine}"
        poly="${base}.poly"
        [ -f "$poly" ] || continue
        verify_out="${base}_verify.txt"
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

# ---- 验证 ----
echo ""
echo "=== 验证 f(x) = g(Mx+b) ==="
run_verify

# ====================================================================
# 汇总
# ====================================================================
echo ""
echo "============================================"
echo -e " ${GREEN}${pass_cnt} passed${NC}, ${RED}${fail_cnt} failed${NC}, ${YELLOW}${skip_cnt} skipped${NC}"
echo "============================================"

echo ""
echo "--- Best union_T ---"
best_ut=999999999
best_st=0
best_name=""
for sf in "$RESULTS_DIR/$INST"/*_stats.txt; do
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
    printf "  %-12s sum_T=%-8s union_T=%-5s (%s)\n" "$INST" "$best_st" "$best_ut" "$best_name"
fi

exit $fail_cnt
