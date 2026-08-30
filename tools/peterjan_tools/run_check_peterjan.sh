#!/bin/bash
set -e

# CHEAP CHECK (run this first): pure-axis thresholds for the gauge-honeycomb
# phenomenological repeated-measurement circuits, without the full 2-D grid.
#
#   - Scan A: pure p_m axis, p (data noise) pinned to 0 -> the p_m-axis intercept
#     of the boundary, one curve per gauge factor g. -> peterjan_check_pz.csv
#   - Scan B: pure p axis, p_m pinned to 0             -> the p-axis intercept.
#                                                       -> peterjan_check_pres.csv
#
# Use the fitted p_m thresholds here to set the per-g PZ_HI ceilings in
# run_grid_peterjan.sh. Ends with analyze_peterjan_check.py (per-g LER panels +
# threshold table). Distances 4, 8, 12; single memory observable.
#
#   bash tools/peterjan_tools/run_check_peterjan.sh          # g = 1 2 3
#   bash tools/peterjan_tools/run_check_peterjan.sh 2         # a single g
#
# Budgets via env: MAX_SHOTS (default 100000), MAX_ERRORS (default 500).

cd "$(dirname "$0")/../.."
export PYTHONUNBUFFERED=1
PY="${PYTHON:-.venv/bin/python}"

# g>=4 requires tools/peterjan_tools/extend_timelike_distance_data.py to have run.
GAUGE="${1:-1 2 3 4 5}"
DISTANCES="4 8 12"
OBS="${OBS:-memory_x}"
MAX_SHOTS="${MAX_SHOTS:-100000}"
MAX_ERRORS="${MAX_ERRORS:-500}"
NPZ=25
NP=20
# The pinned axis is set to a tiny floor, NOT exactly 0: at p=0 the data-noise
# edges vanish and the DEM no longer decomposes into graphlike components (these
# periodic circuits need every channel present, like the paper's renormalized
# noise). 1e-6 matches the p_res floor used for the full_vs_local pipeline
# (tools/_grid_scans/run_pres_threshold_check.sh) and keeps the graph well-formed
# without shifting the threshold.
FLOOR=1e-6

# p_m scan ranges bracket the expected thresholds (they rise with g). Widen per
# g once the first fit lands; these are conservative starting windows.
for G in $GAUGE; do
  case "$G" in
    1) PM_LO=1e-3; PM_HI=8e-2 ;;
    2) PM_LO=3e-3; PM_HI=1.5e-1 ;;
    3) PM_LO=5e-3; PM_HI=2.5e-1 ;;
    4) PM_LO=5e-3; PM_HI=3e-1 ;;
    5) PM_LO=5e-3; PM_HI=3.5e-1 ;;
    *) PM_LO=1e-3; PM_HI=3.5e-1 ;;
  esac
  PM=$($PY -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10($PM_LO), np.log10($PM_HI), $NPZ)))")

  echo "=== [g=$G] Scan A: pure p_m (p=0), ${NPZ} points ==="
  $PY tools/peterjan_tools/make_grid_circuits_peterjan.py \
    --distances $DISTANCES --gauge-factors "$G" --observable "$OBS" \
    --p-values $FLOOR --pm-values $PM \
    --output-dir "out/circuits/peterjan_circuits/check_pm_g${G}" --skip-existing
  $PY tools/decode.py \
    --circuit-files "out/circuits/peterjan_circuits/check_pm_g${G}"/*.stim \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" --print-progress \
    --output-preamble peterjan_check_pz --resume

  # Scan B: pure data-noise axis, p_m = 0. Threshold ~ a few 1e-3, so a tighter
  # log window regardless of g.
  P=$($PY -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-4, -1.3, $NP)))")
  echo "=== [g=$G] Scan B: pure p (p_m=0), ${NP} points ==="
  $PY tools/peterjan_tools/make_grid_circuits_peterjan.py \
    --distances $DISTANCES --gauge-factors "$G" --observable "$OBS" \
    --p-values $P --pm-values $FLOOR \
    --output-dir "out/circuits/peterjan_circuits/check_p_g${G}" --skip-existing
  $PY tools/decode.py \
    --circuit-files "out/circuits/peterjan_circuits/check_p_g${G}"/*.stim \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" --print-progress \
    --output-preamble peterjan_check_pres --resume
done

echo
echo "=== Analyze: fit pure-p_m thresholds per g ==="
$PY tools/peterjan_tools/analyze_peterjan_check.py

echo
echo "Data: out/collection/peterjan_check_pz.csv, peterjan_check_pres.csv"
echo "Figure: out/figures/peterjan_check_pz.pdf"
