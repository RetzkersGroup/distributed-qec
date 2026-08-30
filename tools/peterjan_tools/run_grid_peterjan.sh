#!/bin/bash
set -e

# Full (p, p_m) grid for the gauge-honeycomb phenomenological repeated-measurement
# circuits, one grid per gauge factor g. p is the data-noise axis
# (PAULI_CHANNEL_1(p/3,p/3,p/3), the p_res analog); p_m is the measurement-flip
# axis (the p_z analog). Decodes into a shared out/collection/peterjan_grid.csv
# (curves keyed by noise_model = peterjan_pheno_g{g}); then draws the boundary
# figure p_m* vs p* via plot_peterjan_boundary.py.
#
# Run run_check_peterjan.sh FIRST and set the per-g PM_HI ceilings below from the
# fitted pure-p_m thresholds (the boundary crossing must sit below the ceiling).
#
#   bash tools/peterjan_tools/run_grid_peterjan.sh          # g = 1 2 3
#   bash tools/peterjan_tools/run_grid_peterjan.sh 2        # a single g
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
PREAMBLE="peterjan_grid"
NP=25; NPM=25
P_HI=2e-3   # data-noise ceiling: repetition mainly moves the p_m threshold.
P=$($PY -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-4, np.log10($P_HI), $NP)))")

for G in $GAUGE; do
  # p_m ceiling grows with g; set from run_check_peterjan.sh so the crossing is
  # bracketed (and slightly exceeded) rather than pinned against the grid edge.
  case "$G" in
    1) PM_HI=8e-2 ;;
    2) PM_HI=1.5e-1 ;;
    3) PM_HI=2.5e-1 ;;
    4) PM_HI=3e-1 ;;
    5) PM_HI=3.5e-1 ;;
    *) PM_HI=3.5e-1 ;;
  esac
  PM=$($PY -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-3, np.log10($PM_HI), $NPM)))")
  OUTPUT_DIR="out/circuits/peterjan_circuits/grid_g${G}"

  echo "=== [g=$G] generate (p, p_m) grid (${NP} p x ${NPM} p_m) ==="
  $PY tools/peterjan_tools/make_grid_circuits_peterjan.py \
    --distances $DISTANCES --gauge-factors "$G" --observable "$OBS" \
    --p-values $P --pm-values $PM \
    --output-dir "$OUTPUT_DIR" --skip-existing

  echo "=== [g=$G] decode (resume) -> out/collection/${PREAMBLE}.csv ==="
  $PY tools/decode.py \
    --circuit-files "$OUTPUT_DIR"/*.stim \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" --print-progress \
    --output-preamble "$PREAMBLE" --resume
done

echo
echo "Done. Data: out/collection/${PREAMBLE}.csv"
echo "Figure: $PY tools/peterjan_tools/plot_peterjan_boundary.py"
