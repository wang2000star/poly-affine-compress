#!/bin/bash
source "$(dirname "$0")/../common.sh"
INST="hd08"
STRATEGIES=(d1a_opt d1b_opt d2_opt d3_opt)

EQN="$PREPROCESS_DIR/$INST/$INST.eqn"
if [ ! -f "$EQN" ]; then
    echo -e "${RED}ERROR: $EQN not found — run preprocess.sh first${NC}"
    exit 1
fi

echo "=== Verify hd08 ==="
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
