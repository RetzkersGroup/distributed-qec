#!/bin/bash
set -e

# Generate stairway (p_res, p_z) SLICE circuits for the LER-vs-p_z figure.
#
# Stairway codes have NO threshold, so we don't scan a full grid — we only need a
# single fixed p_res slice with p_z swept. That is exactly what plot_grid.py's
# LER-vs-p_z panel draws (one curve per distance; no crossing needed).
#
# p_res is pinned to 5e-4 (the --p-values single value). p_z sweeps a log grid.
#
# Stairway 2BGA circuits have a noiseless head/tail marked by
# TICK[start_noise]/TICK[stop_noise]; only the *no_noise* inputs carry those
# tags, so we translate those with --noiseless-start-end (Bell-prep noise is
# injected only inside the window).
#
#   bash tools/stairway_codes/make_circuits_stairway.sh            # baseline (equivalent_v3)
#   bash tools/stairway_codes/make_circuits_stairway.sh majority   # 3-rep majority vote
#   bash tools/stairway_codes/make_circuits_stairway.sh majority5  # 5-rep majority vote

cd "$(dirname "$0")/../.."

SCHEME="${1:-bell}"

# Stairway inputs: use the no_noise circuits (they carry start/stop_noise tags).
INPUT_GLOB="in/circuits/stairway/*no_noise*.stim"

# Slice axes ------------------------------------------------------------------
P_RES=3e-4                    # fixed residual-noise slice
NPZ=12; PZ_HI=3e-2            # p_z log grid: 1e-3 .. PZ_HI
PZ_LO=$P_RES

case "$SCHEME" in
  bell)
    OUTPUT_DIR="out/circuits/grid_stairway"
    SCHEME_ARGS="--scheme equivalent_v3"
    ;;
  majority)
    OUTPUT_DIR="out/circuits/grid_stairway_majority"
    SCHEME_ARGS="--scheme equivalent_v3_majority"
    ;;
  majority5)
    OUTPUT_DIR="out/circuits/grid_stairway_majority5"
    SCHEME_ARGS="--scheme equivalent_v3_majority5"
    ;;
  *)
    echo "unknown scheme: $SCHEME (use bell|majority|majority5)" >&2; exit 1 ;;
esac

PZ=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10($PZ_LO), np.log10($PZ_HI), $NPZ)))")

echo "=== Generate '$SCHEME' stairway slice circuits (p_res=$P_RES x ${NPZ} p_z) ==="
echo "p_z grid: $PZ"
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$INPUT_GLOB" \
  --output-dir "$OUTPUT_DIR" \
  $SCHEME_ARGS \
  --p-values $P_RES \
  --p-bell-z-values $PZ \
  --noiseless-start-end \
  --skip-existing
