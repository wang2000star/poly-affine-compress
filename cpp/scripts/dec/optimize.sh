#!/bin/bash
source "$(dirname "$0")/../common.sh"
INST="dec"

echo "=== Optimize dec (6 multi-output strategies) ==="
CIRCUIT="$PREPROCESS_DIR/$INST/$INST.eqn"

P_COMMON="--random 500 --walsh-k 40 --hill-climb 30"
P_D1A_OPT2="--random 500 --walsh-k 40 --hill-climb 10 --max-m 12"
P_D2="--random 300 --hill-climb 10 --max-m 6"

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
