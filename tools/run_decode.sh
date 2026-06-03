#!/bin/bash
set -e

python tools/decode.py \
  --circuit-files out/circuits/hyperbolic_floquet/*.stim \
  --max-shots 100000 --max-errors 200 \
  --print-progress \
  --output-preamble hyperbolic_scan
