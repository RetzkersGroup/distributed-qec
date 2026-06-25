#!/bin/bash
set -e

python tools/decode.py \
  --circuit-files out/circuits/thresholds_stairway/*v3*.stim \
  --max-shots 100000 --max-errors 1000 \
  --decoder tesseract_short \
  --print-progress \
  --output-preamble threshold_stairway
