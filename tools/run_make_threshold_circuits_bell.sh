#!/bin/bash
set -e

python tools/make_threshold_circuits_bell.py \
  --input "in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim" \
  --output-dir out/circuits/thresholds_bell \
  --noise-model biased_bell_v3
  # --skip-existing \
