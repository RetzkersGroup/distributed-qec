#!/bin/bash
set -e

# Full (p_res, p_z) grid pipeline for BOTH codes in one shot: sample the
# honeycomb code (base bell + 3-repetition majority vote) and the hyperbolic
# Floquet code, then produce every per-code and combined figure.
#
#   bash tools/run_grid_all.sh          # hyperbolic test set (d=2,4,6)
#   bash tools/run_grid_all.sh full     # hyperbolic full set (d=2..12, heavy)
#
# Sampling uses `decode.py --resume`: each family's collection CSV
# (out/collection/grid_{bell,hyperbolic}.csv, sinter format) is resumed from, so
# a re-run only samples the shortfall -- nothing, if already at budget -- and
# then re-plots. Circuit generation uses --skip-existing where possible. No CSV
# is ever deleted; raise --max-shots/--max-errors in the sub-scripts to collect
# more (resume tops up).

cd "$(dirname "$0")/.."

HYP_MODE="${1:-test}"

echo "############################################################"
echo "# 1-2. Sample honeycomb (base bell + majority vote) [resume]"
echo "############################################################"
bash tools/_grid_scans/run_grid_bell.sh bell
bash tools/_grid_scans/run_grid_bell.sh majority
bash tools/_grid_scans/run_grid_bell.sh majority5

echo "############################################################"
echo "# 3. Sample hyperbolic Floquet (${HYP_MODE}, 3d rounds) [resume]"
echo "#    base + 3-rep + 5-rep majority vote"
echo "############################################################"
bash tools/_grid_scans/run_grid_hyperbolic.sh "$HYP_MODE" base
bash tools/_grid_scans/run_grid_hyperbolic.sh "$HYP_MODE" majority
bash tools/_grid_scans/run_grid_hyperbolic.sh "$HYP_MODE" majority5

echo "############################################################"
echo "# 4. Plot"
echo "############################################################"
# Per-code figures (boundary + LER + debug):
uv run python tools/_grid_scans/plot_grid.py \
  --csv "out/collection/grid_bell.csv" --prefix grid_bell
uv run python tools/_grid_scans/plot_grid.py \
  --csv "out/collection/grid_hyperbolic.csv" --prefix grid_hyperbolic \
  --code-label "Hyperbolic Floquet code" \
  --majority-label "Hyperbolic Floquet + 3-rep. majority vote"
# Combined figures (honeycomb base + majority vs hyperbolic):
uv run python tools/_grid_scans/plot_grid_summary.py
uv run python tools/_grid_scans/plot_grid_boundary.py \
  --source "Honeycomb code=out/collection/grid_bell.csv" \
  --source "Hyperbolic Floquet code=out/collection/grid_hyperbolic.csv"

echo
echo "Done. Figures in out/figures/ (grid_bell_*, grid_hyperbolic_*, grid_summary.pdf, grid_combined_threshold.png)."
