#!/bin/bash
set -e

python tools/make_threshold_circuits_bell.py \
  --input "in/circuits/stairway/*no_noise*.stim" \
  --output-dir out/circuits/thresholds_stairway_bias1/ \
  --noise-model biased_bell_v3 \
  --biases 1 \
  --p-bell-z-values \
  --noiseless-start-end
