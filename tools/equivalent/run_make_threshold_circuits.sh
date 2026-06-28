#!/bin/bash
set -e

python tools/equivalent/make_threshold_circuits.py \
  --input "in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim" \
  --output-dir out/circuits/thresholds_equivalent
  # --skip-existing \
