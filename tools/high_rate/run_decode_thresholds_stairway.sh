#!/bin/bash
set -e

python tools/decode.py \
  --circuit-files out/circuits/thresholds_stairway/*v3*p=0.000*.stim \
  --max-shots 1000 --max-errors 10 \
  --decoder tesseract_short \
  --print-progress \
  --output-preamble threshold_stairway
