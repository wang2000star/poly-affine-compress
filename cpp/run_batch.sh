#!/bin/bash
# Batch optimization for all small instances (n≤16)
# Each run with: max_m=n, walsh_k=40, random=500

EXAMPLES=/home/wangfz/bool/examples
OPT="/home/wangfz/bool/cpp/build/optimize_anf"
LOGS=/home/wangfz/bool/cpp/logs

mkdir -p "$LOGS"

instances=(
    "hd07:8:$EXAMPLES/hd07/hd07.txt"
    "hd03:16:$EXAMPLES/hd03/hd03.txt"
    "hd04:16:$EXAMPLES/hd04/hd04.txt"
    "ctrl:7:$EXAMPLES/ctrl/ctrl.txt"
    "dec:8:$EXAMPLES/dec/dec.txt"
    "int2float:11:$EXAMPLES/int2float/int2float.txt"
    "cavlc:10:$EXAMPLES/cavlc/cavlc.txt"
)

for inst in "${instances[@]}"; do
    IFS=':' read -r name n path <<< "$inst"
    echo "Starting $name (n=$n) at $(date)"
    
    # Full rank search
    nohup "$OPT" "$path" --max-m "$n" --walsh-k 40 --random 500 \
        > "${LOGS}/${name}_opt.log" 2>&1 &
    
    # Also try smaller m for dimensionality reduction
    small_m=$(( n > 8 ? 8 : n / 2 ))
    nohup "$OPT" "$path" --max-m "$small_m" --walsh-k 40 --random 500 \
        > "${LOGS}/${name}_opt_m${small_m}.log" 2>&1 &
done

echo "All batch jobs started at $(date)"
