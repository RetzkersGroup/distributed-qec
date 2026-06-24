#!/bin/bash
set -e

python tools/make_threshold_circuits_bell.py \
  --input "in/circuits/hyperbolic_floquet/8.8.8/**/*.stim" \
  --output-dir out/circuits/thresholds_hyperbolic_floquet \
  --noise-model biased_bell_v3 \
  --biases 100 \
  --p-bell-z-values
