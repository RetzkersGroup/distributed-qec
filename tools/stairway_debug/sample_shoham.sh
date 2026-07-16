#!/bin/bash
set -e

python tools/decode.py \
  --circuit-files in/circuits/stairway/*em3*.stim \
  --max-shots 1000000 --max-errors 10 \
  --decoder tesseract_short \
  --print-progress \
  --output-preamble threshold_stairway