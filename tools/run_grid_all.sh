#!/bin/bash
set -e

# Full (p_res, p_z) grid pipeline for BOTH codes in one shot: sample the
# honeycomb code (base bell scheme + 3-repetition majority vote) and the
# hyperbolic Floquet code, then produce every per-code and combined figure.
#
#   bash tools/run_grid_all.sh          # hyperbolic test set (d=2,4,6)
#   bash tools/run_grid_all.sh full     # hyperbolic full set (d=2..12, heavy)
#
# Sampling is delegated to the per-family scripts (which use `uv run`). To make
# reruns fast, each family's sampling is SKIPPED when its collection CSV already
# exists (the data is reused and we just re-plot). Pass "force" as the 2nd arg to
# re-sample regardless. This script never deletes any CSV.
#
#   bash tools/run_grid_all.sh              # test hyperbolic; skip sampling if CSVs exist
#   bash tools/run_grid_all.sh full         # full hyperbolic; skip sampling if CSVs exist
#   bash tools/run_grid_all.sh test force    # re-sample everything from scratch

cd "$(dirname "$0")/.."

HYP_MODE="${1:-test}"
FORCE="${2:-}"

have_csv() { compgen -G "$1" >/dev/null 2>&1; }

echo "############################################################"
echo "# 1-2. Sample honeycomb (base bell + majority vote)"
echo "############################################################"
if [ "$FORCE" != "force" ] && have_csv "out/collection/grid_bell,*.csv"; then
  echo "  found out/collection/grid_bell,*.csv -> reusing it, skipping sampling"
  echo "  (pass 'force' as the 2nd arg to re-sample)"
else
  bash tools/_grid_scans/run_grid_bell.sh bell
  bash tools/_grid_scans/run_grid_bell.sh majority
fi

echo "############################################################"
echo "# 3. Sample hyperbolic Floquet (${HYP_MODE}, 3d rounds)"
echo "############################################################"
if [ "$FORCE" != "force" ] && have_csv "out/collection/grid_hyperbolic,*.csv"; then
  echo "  found out/collection/grid_hyperbolic,*.csv -> reusing it, skipping sampling"
else
  bash tools/_grid_scans/run_grid_hyperbolic.sh "$HYP_MODE"
fi

echo "############################################################"
echo "# 4. Plot"
echo "############################################################"
# Per-code figures (boundary + LER + debug):
uv run python tools/_grid_scans/plot_grid.py \
  --csv "out/collection/grid_bell,*.csv" --prefix grid_bell
uv run python tools/_grid_scans/plot_grid.py \
  --csv "out/collection/grid_hyperbolic,*.csv" --prefix grid_hyperbolic \
  --code-label "Hyperbolic Floquet code"
# Combined figures (honeycomb base + majority vs hyperbolic):
uv run python tools/_grid_scans/plot_grid_summary.py
uv run python tools/_grid_scans/plot_grid_boundary.py \
  --source "Honeycomb code=out/collection/grid_bell,*.csv" \
  --source "Hyperbolic Floquet code=out/collection/grid_hyperbolic,*.csv"

echo
echo "Done. Figures in out/figures/ (grid_bell_*, grid_hyperbolic_*, grid_summary.pdf, grid_combined_threshold.png)."
