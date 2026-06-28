#!/bin/bash
set -e

python tools/equivalent/decode.py \
  --circuit-files out/circuits/thresholds_equivalent_majority/*.stim \
  --max-shots 100000 --max-errors 1000 \
  --print-progress \
  --output-preamble threshold_equivalent_majority
