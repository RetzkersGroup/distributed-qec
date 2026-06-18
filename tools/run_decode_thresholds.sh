#!/bin/bash
set -e

python tools/decode.py \
  --circuit-files out/circuits/thresholds/*.stim \
  --max-shots 100000 --max-errors 1000 \
  --print-progress \
  --output-preamble threshold_3d
