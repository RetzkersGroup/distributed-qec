#!/bin/bash
set -e

python tools/make_threshold_circuits.py \
  --input "in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim" \
  --output-dir out/circuits/thresholds \
  --skip-existing
