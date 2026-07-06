#!/bin/bash
set -e

# Grid scan for the hyperbolic Floquet code (biased_bell_v3, no repetition) on
# the same (p_res, p_z) grid as the honeycomb bell run. Reuses the generic
# generator (make_grid_circuits_bell.py) and decoder, writing to its own CSV
# (grid_hyperbolic) so the data stays separate from the honeycomb runs.
#
#   bash tools/_grid_scans/run_grid_hyperbolic.sh          # test: d=2,4,6  (l=1,2,3)
#   bash tools/_grid_scans/run_grid_hyperbolic.sh full     # large: d=2..12 (l=1..6)
#
# Uses `uv run` so it works standalone.

cd "$(dirname "$0")/../.."

# Distance set. l=1,2,3 -> d=2,4,6 (small, for validating the pipeline);
# l=1..6 -> d=2,4,6,8,10,12 (the full run, much heavier to decode).
if [ "${1:-test}" = "full" ]; then
  DIST_GLOB='semi_hyperbolic_floquet_l=[123456]_*.stim'
else
  DIST_GLOB='semi_hyperbolic_floquet_l=[123]_*.stim'
fi

INPUT_GLOB="in/circuits/hyperbolic_floquet/8.8.8/**/$DIST_GLOB"
SCALED_DIR="out/circuits/hyperbolic_3d_inputs"
OUTPUT_DIR="out/circuits/grid_hyperbolic"
PREAMBLE="grid_hyperbolic"
ROUNDS_MULTIPLE=3          # run each code for ROUNDS_MULTIPLE * d rounds
MAX_SHOTS=1000000
MAX_ERRORS=5000

# Same axes as the honeycomb bell grid: p_res 1e-4..1e-2, p_z 1e-3..3e-2, 25x25.
NP=25; NPZ=25; P_HI=1e-2; PZ_HI=3e-2
P=$(uv run python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-4, np.log10($P_HI), $NP)))")
PZ=$(uv run python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-3, np.log10($PZ_HI), $NPZ)))")

echo "=== Stage 0: rescale each code to ${ROUNDS_MULTIPLE}d rounds ==="
rm -rf "$SCALED_DIR"          # avoid stale distances from a previous run
uv run python tools/_grid_scans/scale_rounds.py \
  --input "$INPUT_GLOB" \
  --output-dir "$SCALED_DIR" \
  --rounds-multiple "$ROUNDS_MULTIPLE"

echo "=== Stage 1: generate hyperbolic grid circuits ==="
rm -rf "$OUTPUT_DIR"          # clean slate so the decode glob never mixes runs
uv run python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$SCALED_DIR/*.stim" \
  --output-dir "$OUTPUT_DIR" \
  --noise-model biased_bell_v3 \
  --p-values $P \
  --p-bell-z-values $PZ

echo "=== Stage 2: decode -> out/collection/${PREAMBLE},<date>.csv ==="
uv run python tools/decode.py \
  --circuit-files "$OUTPUT_DIR"/*.stim \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --print-progress \
  --output-preamble "$PREAMBLE"
