#!/bin/bash
set -e

# CHEAP CHECK: pure p_z threshold for the full-decoder repeated circuits,
# without running the full (p_res, p_z) grid. Pin p_res to a negligible floor
# and scan p_z for each repetition count r — the p_z-axis intercept of the
# boundary, directly comparable to the majority-vote intercepts from
# run_pres_threshold_check.sh (pz_check.csv).
#
#   bash tools/_grid_scans/run_pres_threshold_check_full.sh          # r = 2 3 4 5
#   bash tools/_grid_scans/run_pres_threshold_check_full.sh 3        # a single r
#
# Writes circuits to out/circuits/full_check_pz_rep{r}/ and decodes into
# out/collection/full_check_pz.csv (schemes distinguished by noise_model =
# biased_bell_v3_equivalent_fullrep{r}). Ends with analyze_full_check.py,
# which prints the fitted thresholds next to the r=1 / majority references
# from pz_check.csv. plot_full_vs_local.py also picks this CSV up as the
# p_res=0 intercept stars.

cd "$(dirname "$0")/../.."

export PYTHONUNBUFFERED=1

REPS="${1:-2 3 4 5}"
INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="full_check_pz"
MAX_SHOTS=100000
MAX_ERRORS=500

# As in run_pres_threshold_check.sh Scan B: p_res is pinned to a tiny floor
# rather than exactly 0, so the matching graph keeps its boundary connections
# without shifting the threshold.
P_FLOOR=1e-6
NPZ=25

for R in $REPS; do
  # p_z scan ranges bracket the expected thresholds: the full decoder should
  # sit at or above the paired majority vote (3-rep ~7%, 5-rep ~11%), and
  # r=2 somewhere between the r=1 baseline (~1.4%) and the 3-rep value.
  case "$R" in
    2) PZ_LO=3e-3; PZ_HI=1e-1 ;;
    3) PZ_LO=5e-3; PZ_HI=1.5e-1 ;;
    4) PZ_LO=5e-3; PZ_HI=1.8e-1 ;;
    5) PZ_LO=1e-2; PZ_HI=2.5e-1 ;;
    *) PZ_LO=3e-3; PZ_HI=2e-1 ;;
  esac
  PZ=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10($PZ_LO), np.log10($PZ_HI), $NPZ)))")
  OUTPUT_DIR="out/circuits/full_check_pz_rep${R}"

  echo "=== [r=$R] generate p_z-scan circuits (p_res=$P_FLOOR, ${NPZ} p_z) ==="
  python tools/_grid_scans/make_grid_circuits_bell.py \
    --input "$INPUT_GLOB" \
    --output-dir "$OUTPUT_DIR" \
    --scheme equivalent_v3 \
    --repeat "$R" \
    --p-values $P_FLOOR \
    --p-bell-z-values $PZ \
    --skip-existing

  echo "=== [r=$R] decode (resume) -> out/collection/${PREAMBLE}.csv ==="
  python tools/decode.py \
    --circuit-files "$OUTPUT_DIR"/*.stim \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
    --print-progress \
    --output-preamble "$PREAMBLE" --resume
done

echo
echo "=== Analyze: fit the p_z threshold per repetition ==="
python tools/_grid_scans/analyze_full_check.py

echo
echo "Done. Data: out/collection/${PREAMBLE}.csv | figure: out/figures/full_check_pz.pdf"
