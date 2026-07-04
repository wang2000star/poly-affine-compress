#!/bin/bash
source "$(dirname "$0")/../common.sh"
INST="ctrl"
STRATEGIES=(d1a_opt1 d1a_opt2 d1b_opt2 d2_opt2 d3_opt1 d3_opt2)

EQN="$PREPROCESS_DIR/$INST/$INST.eqn"
if [ ! -f "$EQN" ]; then
    echo -e "${RED}ERROR: $EQN not found — run preprocess.sh first${NC}"
    exit 1
fi

echo "=== Verify ctrl ==="
echo "  Circuit: $EQN"
echo ""

pass=0; fail=0
for label in "${STRATEGIES[@]}"; do
    if verify_result "$INST" "$label" "$EQN"; then
        pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
done

echo ""
echo "=== Result: $pass passed, $fail failed ==="
exit $fail
