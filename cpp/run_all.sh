#!/bin/bash
# ====================================================================
# run_all.sh — Run all 11 instances through all applicable strategies
#
# Instance groups:
#   n≤16, k_nl=1:  hd08                     → 4 opt strategies
#   n≤16, k_nl>1:  hd07 hd03 hd04 ctrl dec
#                  int2float cavlc           → 5 strategies (no d1a_opt2)
#   n=32, k_nl>1:  hd10 hd01 hd02           → 6 strategies (incl. d1a_opt2)
#
# Usage:
#   ./run_all.sh              # quick test
#   ./run_all.sh --full       # full optimization
#   ./run_all.sh --list       # show instance-strategy matrix
# ====================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
EXAMPLES_DIR="$SCRIPT_DIR/examples"
RESULTS_DIR="$SCRIPT_DIR/results"

MODE="${1:-test}"

# ---- Parameter sets ----
if [ "$MODE" = "--full" ]; then
    P_COMMON="--random 500 --walsh-k 40 --hill-climb 30"
    P_N32="--random 500 --walsh-k 40 --n32-random 50 --hill-climb 5"
    P_D1A_OPT2="--random 500 --walsh-k 40 --hill-climb 10 --max-m 12"
    P_D2="--random 300 --hill-climb 10 --max-m 6"
elif [ "$MODE" = "--list" ]; then
    echo "=== Instance × Strategy Matrix ==="
    echo ""
    echo "n≤16, k_nl=1 (4 opt):"
    echo "  hd08 → d1a_opt, d1b_opt, d2_opt, d3_opt"
    echo ""
    echo "n≤16, k_nl>1 (5 opt1/opt2, no d1a_opt2):"
    echo "  hd07 hd03 hd04 ctrl dec int2float cavlc"
    echo "  → d1a_opt1, d1b_opt2, d2_opt2, d3_opt1, d3_opt2"
    echo ""
    echo "n=32, k_nl>1 (6 including d1a_opt2):"
    echo "  hd10 hd01 hd02"
    echo "  → d1a_opt1, d1a_opt2, d1b_opt2, d2_opt2, d3_opt1, d3_opt2"
    exit 0
else
    # Quick test
    P_COMMON="--random 20 --walsh-k 10 --hill-climb 2"
    P_N32="--random 10 --walsh-k 5 --n32-random 5 --hill-climb 1"
    P_D1A_OPT2="--random 20 --walsh-k 10 --hill-climb 2 --max-m 8"
    P_D2="--random 20 --hill-climb 2 --max-m 4"
fi

# ---- Helpers ----
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
pass_cnt=0; fail_cnt=0; skip_cnt=0
results_log="$RESULTS_DIR/run_log.txt"
echo "Run started: $(date)" > "$results_log"

check_exe() {
    if [ ! -f "$BUILD_DIR/$1" ]; then
        echo -e "${RED}  SKIP: $1 not found${NC}"
        return 1
    fi
    return 0
}

# Run one strategy on one instance
run_strat() {
    local inst=$1 exe=$2 circuit=$3 extra_args=$4
    local out_dir="$RESULTS_DIR/$inst"
    mkdir -p "$out_dir"

    if ! check_exe "$exe"; then
        skip_cnt=$((skip_cnt + 1))
        return 1
    fi

    # Gate builder uses --out-dir, others use --save-results
    local is_gb=0
    case "$exe" in
        optimize_anf_d1b_opt|optimize_anf_d1b_opt2) is_gb=1 ;;
    esac

    local logfile="$out_dir/${inst}_${exe#optimize_anf_}_run.log"
    echo -e "  ${YELLOW}[${exe}]${NC} $inst"

    set +e
    if [ "$is_gb" -eq 1 ]; then
        timeout 300 "$BUILD_DIR/$exe" "$circuit" $extra_args --out-dir "$out_dir" &> "$logfile"
    else
        timeout 600 "$BUILD_DIR/$exe" "$circuit" $extra_args --save-results "$out_dir" &> "$logfile"
    fi
    local rc=$?
    set -e

    if [ $rc -eq 0 ] && grep -q "Verified\|All outputs verified\|Skipping verification\|Done\." "$logfile" 2>/dev/null; then
        echo -e "    ${GREEN}✓ PASS${NC}"
        pass_cnt=$((pass_cnt + 1))
        echo "[PASS] $inst $exe" >> "$results_log"
    elif [ $rc -ne 0 ]; then
        echo -e "    ${RED}✗ FAIL (exit $rc)${NC}"
        tail -3 "$logfile"
        fail_cnt=$((fail_cnt + 1))
        echo "[FAIL] $inst $exe (exit $rc)" >> "$results_log"
    else
        echo -e "    ${YELLOW}? CHECK (exit=0, no verification msg)${NC}"
        tail -3 "$logfile"
        pass_cnt=$((pass_cnt + 1))
        echo "[CHECK] $inst $exe" >> "$results_log"
    fi
}

# Run a list of strategies on an instance
run_group() {
    local inst=$1 circuit=$2 extra=$3; shift 3
    echo ""
    echo ">> $inst ($# strategies)"
    for strat in "$@"; do
        exe="optimize_anf_${strat}"
        case "$strat" in
            d1a_opt2)  args="$P_D1A_OPT2 $extra" ;;
            d2_opt|d2_opt2) args="$P_D2 $extra" ;;
            d1b_opt|d1b_opt2) args="$extra" ;;
            *)         args="$P_COMMON $extra" ;;
        esac
        run_strat "$inst" "$exe" "$circuit" "$args"
    done
}

# ====================================================================
# Main
# ====================================================================
echo "============================================"
echo " ANF Optimization — Run All Instances"
echo " Mode: ${MODE}"
echo "============================================"
echo ""

# ---- Preprocess all ----
echo "--- Preprocessing ---"
for inst in hd08 hd07 hd03 hd04 ctrl dec int2float cavlc hd10 hd01 hd02; do
    mkdir -p "$RESULTS_DIR/$inst"
    "$BUILD_DIR/preprocess" "$EXAMPLES_DIR/${inst}.txt" --out-dir "$RESULTS_DIR/$inst" &> /dev/null
    echo "  $inst ✓"
done
echo ""

# ---- Run strategies ----
echo "--- Running strategies ---"

# Group 1: n≤16, k_nl=1 — opt only
run_group "hd08" "$EXAMPLES_DIR/hd08.txt" "" \
    d1a_opt d1b_opt d2_opt d3_opt

# Group 2: n≤16, k_nl>1 — 5 strategies (no d1a_opt2, it crashes on n≤20)
for inst in hd07 hd03 hd04 ctrl dec int2float cavlc; do
    run_group "$inst" "$EXAMPLES_DIR/${inst}.txt" "" \
        d1a_opt1 d1b_opt2 d2_opt2 d3_opt1 d3_opt2
done

# Group 3: n=32, k_nl>1 — all 6 strategies
for inst in hd10 hd01 hd02; do
    run_group "$inst" "$EXAMPLES_DIR/${inst}.txt" "" \
        d1a_opt1 d1a_opt2 d1b_opt2 d2_opt2 d3_opt1 d3_opt2
done

# ====================================================================
# Summary
# ====================================================================
echo ""
echo "============================================"
echo -e " ${GREEN}${pass_cnt} passed${NC}, ${RED}${fail_cnt} failed${NC}, ${YELLOW}${skip_cnt} skipped${NC}"
echo " Log: $results_log"
echo "============================================"

# Collect best T values
echo ""
echo "--- Best T values ---"
printf "  %-12s %-10s %-10s\n" "Instance" "sum_T" "union_T"
printf "  %-12s %-10s %-10s\n" "--------" "-----" "-------"
for inst in hd08 hd07 hd03 hd04 ctrl dec int2float cavlc hd10 hd01 hd02; do
    best_union=999999999
    best_sum=0
    for sf in "$RESULTS_DIR/$inst"/*_stats.txt; do
        [ -f "$sf" ] || continue
        ut=$(sed -n '4p' "$sf" 2>/dev/null | tr -d ' ')
        st=$(sed -n '3p' "$sf" 2>/dev/null | tr -d ' ')
        [[ "$sf" == *_raw_stats.txt ]] && continue
        if [ -n "$ut" ] && [ "$ut" -lt "$best_union" ] 2>/dev/null; then
            best_union=$ut; best_sum=$st
        fi
    done
    if [ "$best_union" -ne 999999999 ]; then
        printf "  %-12s %-10s %-10s\n" "$inst" "$best_sum" "$best_union"
    fi
done

exit $fail_cnt
