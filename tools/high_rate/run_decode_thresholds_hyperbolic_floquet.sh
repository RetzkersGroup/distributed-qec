#!/bin/bash
set -e

python tools/decode.py \
  --circuit-files out/circuits/thresholds_hyperbolic_floquet/*v3*.stim \
  --max-shots 100000 --max-errors 1000 \
  --print-progress \
  --output-preamble threshold_hyperbolic_floquet
