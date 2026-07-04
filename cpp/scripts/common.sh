#!/bin/bash
# ====================================================================
# common.sh — Shared config for per-instance scripts
# Source this from per-instance scripts:  source "$(dirname "$0")/../common.sh"
# ====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$PROJ_DIR/build"
EXAMPLES_DIR="$PROJ_DIR/examples"
PREPROCESS_DIR="$PROJ_DIR/preprocessed"
RESULTS_DIR="$PROJ_DIR/results"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

# ---- Strategy tools for multi-output (k_nl>1) ----
STRATS_OPT1=(
    "d1a_opt1:optimize_anf_d1a_opt1"
    "d3_opt1:optimize_anf_d3_opt1"
)
STRATS_OPT2=(
    "d1a_opt2:optimize_anf_d1a_opt2"
    "d1b_opt2:optimize_anf_d1b_opt2"
    "d2_opt2:optimize_anf_d2_opt2"
    "d3_opt2:optimize_anf_d3_opt2"
)
ALL_STRATS_MULTI=( "${STRATS_OPT1[@]}" "${STRATS_OPT2[@]}" )

# ---- Strategy tools for single output (k_nl=1) ----
ALL_STRATS_SINGLE=(
    "d1a_opt:optimize_anf_d1a_opt"
    "d1b_opt:optimize_anf_d1b_opt"
    "d2_opt:optimize_anf_d2_opt"
    "d3_opt:optimize_anf_d3_opt"
)

# ---- Check that all exes exist ----
check_exe() {
    local exe="$1"
    if [ ! -f "$BUILD_DIR/$exe" ]; then
        echo -e "${RED}  SKIP: $exe not found${NC}"
        return 1
    fi
    return 0
}

# ---- Run one strategy ----
# Usage: run_strategy <instance> <strategy_label> <executable> <circuit> <extra_args>
run_strategy() {
    local inst=$1 label=$2 exe=$3 circuit=$4 extra_args=$5
    local out_dir="$RESULTS_DIR/$inst/${inst}_${label}"

    if ! check_exe "$exe"; then return 1; fi

    local logfile="$out_dir/run.log"
    echo -e "  ${YELLOW}[${label}]${NC} $inst"

    # All tools treat their output arg as a directory and append {inst}_{tag}
    case "$exe" in
        optimize_anf_d1b_opt|optimize_anf_d1b_opt2)
            timeout 600 "$BUILD_DIR/$exe" "$circuit" $extra_args --out-dir "$out_dir" &> "$logfile"
            ;;
        *)
            timeout 600 "$BUILD_DIR/$exe" "$circuit" $extra_args --save-results "$out_dir" &> "$logfile"
            ;;
    esac
    local rc=$?

    if [ $rc -eq 0 ]; then
        echo -e "    ${GREEN}✓${NC}"
    else
        echo -e "    ${RED}✗ (exit $rc)${NC}"
        tail -5 "$logfile"
    fi
    return $rc
}

# ---- Verify one result ----
# Usage: verify_result <instance> <strategy_label> <eqn_file>
verify_result() {
    local inst=$1 label=$2 eqn=$3
    local dir="$RESULTS_DIR/$inst/${inst}_${label}"
    local aff="$dir/${inst}_${label}.affine"
    local poly="$dir/${inst}_${label}.poly"
    local verify_out="$dir/${inst}_${label}_verify.txt"

    if [ ! -f "$aff" ] || [ ! -f "$poly" ]; then
        echo -e "    ${YELLOW}${label}: missing .affine or .poly${NC}"
        return 1
    fi

    # Skip if already verified OK
    if [ -f "$verify_out" ] && grep -q "All outputs PASS\|Verification.*PASS" "$verify_out" 2>/dev/null; then
        echo -e "    ${GREEN}${label}: already verified ✓${NC}"
        return 0
    fi

    "$BUILD_DIR/verify_anf" "$eqn" "$aff" "$poly" --output "$verify_out" 2>/dev/null
    if grep -q "All outputs PASS\|Verification.*PASS" "$verify_out" 2>/dev/null; then
        echo -e "    ${GREEN}${label}: ✅ PASS${NC}"
        return 0
    else
        echo -e "    ${RED}${label}: ❌ FAIL${NC}"
        head -20 "$verify_out" | grep -E "FAIL|PASS" | head -5
        return 1
    fi
}
