#!/bin/bash
set -e

python tools/make_threshold_circuits_bell.py \
  --input "in/circuits/stairway/*.stim" \
  --output-dir out/circuits/thresholds_stairway \
  --noise-model biased_bell_v3 \
  --biases 10 \
  --p-bell-z-values
