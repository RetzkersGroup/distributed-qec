#!/bin/bash
set -e

python tools/decode.py \
  --circuit-files out/circuits/thresholds_bell/*.stim \
  --max-shots 100000 --max-errors 1000 \
  --print-progress \
  --output-preamble threshold_bell
