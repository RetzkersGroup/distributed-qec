#!/bin/bash
set -e

INPUT="in/circuits/hyperbolic_floquet/8.8.8/k=10_d3a22e5e/semi_hyperbolic_floquet_l=[123]_*.stim"
P_VALUES="0.0005 0.001 0.002 0.005 0.01"
OUT_DIR="out/circuits/bias_comparison"

echo "=== Generating circuits with bell_bias=1 ==="
python tools/make_circuits.py \
  --input "$INPUT" \
  --p-values $P_VALUES \
  --bell-bias 1 \
  --output-dir $OUT_DIR

echo ""
echo "=== Generating circuits with bell_bias=10 ==="
python tools/make_circuits.py \
  --input "$INPUT" \
  --p-values $P_VALUES \
  --bell-bias 10 \
  --output-dir $OUT_DIR

echo ""
echo "=== Decoding all circuits ==="
python tools/decode.py \
  --circuit-files $OUT_DIR/*.stim \
  --max-shots 1000000 --max-errors 200 \
  --print-progress \
  --output-preamble bias_comparison
