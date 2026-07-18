#!/bin/bash
set -e

# Plot the stairway LER-vs-p_z figure from the shared slice CSV.
#
# Stairway codes have NO threshold, so we only care about the LER-vs-p_z panel
# (out/figures/grid_stairway_ler_vs_pz.png): one curve per distance at the fixed
# p_res slice. plot_grid.py also emits a threshold-boundary panel and a paper
# figure, but with a single p_res those degenerate — the LER panel is the output.
#
# --p-example must match the slice p_res used by make_circuits_stairway.sh (5e-4);
# plot_grid.py snaps it to the nearest sampled p_res anyway.
#
#   bash tools/stairway_codes/plot_stairway.sh

cd "$(dirname "$0")/../.."

python tools/_grid_scans/plot_grid.py \
  --csv out/collection/grid_stairway.csv \
  --prefix grid_stairway \
  --code-label "Stairway 2BGA" \
  --p-example 5e-4 \
  --no-debug

echo
echo "Figure: out/figures/grid_stairway_ler_vs_pz.png"
