#!/bin/bash
set -e

# Generate biased_bell_v3 honeycomb circuits on a direct (p, p_bell_z) grid
# (independent axes), as opposed to the bias-ray scan in
# run_make_threshold_circuits_bell.sh. Analysis of this grid is TBD.
python tools/make_grid_circuits_bell.py \
  --input "in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim" \
  --output-dir out/circuits/grid_bell \
  --noise-model biased_bell_v3
  # --p-values ... --p-bell-z-values ...   # override the default grids
  # --skip-existing \
