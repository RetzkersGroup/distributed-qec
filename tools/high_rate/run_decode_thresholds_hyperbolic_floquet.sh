#!/bin/bash
set -e

python tools/decode.py \
  --circuit-files out/circuits/thresholds_hyperbolic_floquet/*k=10*v3*.stim \
  --max-shots 10000 --max-errors 100 \
  --print-progress \
  --output-preamble threshold_hyperbolic_floquet
