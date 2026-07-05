#!/bin/bash
set -e

# Plot the high-rate decode data collected by
# tools/high_rate/run_decode_thresholds_{hyperbolic_floquet,stairway}.sh.
# Writes figures to out/figures/. Pass --bias N to pick the Z-bias.
python tools/plot_high_rate.py "$@"
