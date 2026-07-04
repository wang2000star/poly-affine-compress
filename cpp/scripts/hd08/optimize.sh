#!/bin/bash
source "$(dirname "$0")/../common.sh"
INST="hd08"

echo "=== Optimize hd08 (4 single-output strategies) ==="
CIRCUIT="$PREPROCESS_DIR/$INST/$INST.eqn"

P_OPT="--random 5000 --walsh-k 40 --hill-climb 50"
P_D2="--random 1000 --hill-climb 20 --max-m 6"

for strat_entry in "${ALL_STRATS_SINGLE[@]}"; do
    label="${strat_entry%%:*}"
    exe="${strat_entry##*:}"
    case "$label" in
        d2_opt)  args="$P_D2" ;;
        d1b_opt) args="" ;;
        *)       args="$P_OPT" ;;
    esac
    run_strategy "$INST" "$label" "$exe" "$CIRCUIT" "$args"
done

echo "=== Optimization complete ==="
