#!/bin/bash
set -e

# Correlated-matching twin of the full-decoder grid runs: the same repeated-
# measurement honeycomb circuits (repeat_measurements applied structurally, then
# the 1-rep equivalent_v3 translation, so the decoder sees the full repetition
# structure with no inner majority vote) decoded with pymatching's native
# correlated MWPM into a SEPARATE CSV out/collection/grid_full_correlated.csv.
#
# Correlated matching is exactly where it should help here: the repeated-
# measurement DEM is non-graphlike (comparison detectors), which plain
# pymatching decomposes lossily.
#
# A separate CSV is mandatory (see run_grid_bell_correlated.sh for why).
#
#   bash tools/_grid_scans/run_grid_full_correlated.sh        # all r: 2 3 4 5
#   bash tools/_grid_scans/run_grid_full_correlated.sh 3      # a single r
#
# Reuses the existing out/circuits/grid_full_rep{r}/ circuits (--skip-existing);
# only the decoder differs from the run that produced grid_full.csv. Plot:
# plot_full_vs_local.py --correlated.

cd "$(dirname "$0")/../.."

export PYTHONUNBUFFERED=1

REPS="${1:-2 3 4 5}"
INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="grid_full_correlated"
DECODER="pymatching-correlated"
MAX_SHOTS=100000
MAX_ERRORS=500

NP=25; NPZ=25
P_HI=2e-3
P=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-4, np.log10($P_HI), $NP)))")

for R in $REPS; do
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

  echo "=== [r=$R] Stage 2: decode ($DECODER, resume) -> out/collection/${PREAMBLE}.csv ==="
  python tools/decode.py \
    --circuit-files "$OUTPUT_DIR"/*.stim \
    --decoder "$DECODER" \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
    --print-progress \
    --output-preamble "$PREAMBLE" --resume
done

echo
echo "Done. Data: out/collection/${PREAMBLE}.csv"
echo "Figure:   python tools/_grid_scans/plot_full_vs_local.py --correlated"
