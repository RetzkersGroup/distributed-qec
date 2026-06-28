#!/bin/bash
set -e

python tools/equivalent/make_threshold_circuits.py \
  --scheme equivalent_v3_majority \
  --input "in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim" \
  --output-dir out/circuits/thresholds_equivalent_majority
  # --skip-existing \
