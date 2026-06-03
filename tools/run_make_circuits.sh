#!/bin/bash
set -e

python tools/make_circuits.py \
  --input "in/circuits/hyperbolic_floquet/8.8.8/k=10_d3a22e5e/*.stim" \
  --p-values 0.0005 0.001 0.002 0.005 0.01 \
  --bell-bias 1000.0 \
  --output-dir out/circuits/hyperbolic_floquet
