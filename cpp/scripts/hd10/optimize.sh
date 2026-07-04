#!/bin/bash
source "$(dirname "$0")/../common.sh"
INST="hd10"

echo "=== Optimize hd10 (6 multi-output strategies, n=32) ==="
CIRCUIT="$EXAMPLES_DIR/${INST}.txt"

P_COMMON="--random 200 --walsh-k 30 --hill-climb 5 --max-m 12"
P_D1A_OPT2="--random 200 --walsh-k 30 --hill-climb 3 --max-m 10"
P_D2="--random 100 --hill-climb 3 --max-m 5"

for strat_entry in "${STRATS_OPT1[@]}"; do
    label="${strat_entry%%:*}"
    exe="${strat_entry##*:}"
    run_strategy "$INST" "$label" "$exe" "$CIRCUIT" "$P_COMMON"
done

for strat_entry in "${STRATS_OPT2[@]}"; do
    label="${strat_entry%%:*}"
    exe="${strat_entry##*:}"
    case "$label" in
        d1a_opt2) args="$P_D1A_OPT2" ;;
        d2_opt2)  args="$P_D2" ;;
        d1b_opt2) args="" ;;
        *)        args="$P_COMMON" ;;
    esac
    run_strategy "$INST" "$label" "$exe" "$CIRCUIT" "$args"
done

echo "=== Optimization complete ==="
