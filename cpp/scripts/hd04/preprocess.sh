#!/bin/bash
source "$(dirname "$0")/../common.sh"
INST="hd04"

echo "--- Preprocessing hd04 ---"
mkdir -p "$PREPROCESS_DIR/$INST"
"$BUILD_DIR/preprocess" "$EXAMPLES_DIR/${INST}.txt" --out-dir "$PREPROCESS_DIR/$INST"
echo "  Done: $PREPROCESS_DIR/$INST/"
