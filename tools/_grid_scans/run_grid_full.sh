#!/bin/bash
set -e

# Full-decoder repeated-measurement (p, p_z) grid runs, honeycomb only.
#
# For each repetition count r, applies repeat_measurements (r copies of every
# mid-circuit measurement layer + comparison detectors) BEFORE the 1-rep
# equivalent_v3 translation, so pymatching decodes the full repetition
# structure (no inner majority vote). Compare against the majority-vote grids
# in out/collection/grid_bell.csv via plot_full_vs_local.py.
#
#   bash tools/_grid_scans/run_grid_full.sh          # all r: 2 3 4 5
#   bash tools/_grid_scans/run_grid_full.sh 3        # a single r
#
# Writes circuits to out/circuits/grid_full_rep{r}/ and decodes into a shared
# out/collection/grid_full.csv (curves distinguished by the noise_model column:
# biased_bell_v3_equivalent_fullrep{r}). Decoding uses --resume, so re-invoking
# only tops up toward the shot/error budget.

cd "$(dirname "$0")/../.."

export PYTHONUNBUFFERED=1

REPS="${1:-2 3 4 5}"
INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="grid_full"
MAX_SHOTS=100000
MAX_ERRORS=500

# p axis matches the majority grids (P_HI=2e-3): repetition mainly moves the
# p_z threshold, and the p_res threshold only drops with r.
NP=25; NPZ=25
P_HI=2e-3
P=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-4, np.log10($P_HI), $NP)))")

for R in $REPS; do
  # p_z ceiling grows with r, bracketing (and slightly exceeding) the paired
  # majority-vote thresholds since the full decoder should sit above them:
  # 3-rep majority ~7%, 5-rep ~11%.
  case "$R" in
    2) PZ_HI=8e-2 ;;
    3) PZ_HI=1.2e-1 ;;
    4) PZ_HI=1.6e-1 ;;
    5) PZ_HI=2e-1 ;;
    *) PZ_HI=2e-1 ;;
  esac
  PZ=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-3, np.log10($PZ_HI), $NPZ)))")
  OUTPUT_DIR="out/circuits/grid_full_rep${R}"

  echo "=== [r=$R] Stage 1: generate full-decoder grid circuits (${NP} p x ${NPZ} p_z) ==="
  python tools/_grid_scans/make_grid_circuits_bell.py \
    --input "$INPUT_GLOB" \
    --output-dir "$OUTPUT_DIR" \
    --scheme equivalent_v3 \
    --repeat "$R" \
    --p-values $P \
    --p-bell-z-values $PZ \
    --skip-existing

  echo "=== [r=$R] Stage 2: decode (resume) -> out/collection/${PREAMBLE}.csv ==="
  python tools/decode.py \
    --circuit-files "$OUTPUT_DIR"/*.stim \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
    --print-progress \
    --output-preamble "$PREAMBLE" --resume
done

echo
echo "Done. Data: out/collection/${PREAMBLE}.csv"
echo "Figure:   python tools/_grid_scans/plot_full_vs_local.py"
