#!/bin/bash
set -e

# Reproduce the honeycomb (p, p_z) grid runs: generate circuits then decode into
# a SHARED CSV (out/collection/grid_bell,<date>.csv) so every scheme's grid data
# lives in one place (distinguished by the noise_model metadata column).
#
#   bash tools/run_grid_bell.sh            # default: two-ancilla bell scheme
#   bash tools/run_grid_bell.sh majority   # repeated-measurement (majority vote)
#
# Uses `uv run` so it works standalone. Each scheme writes its circuits to its
# own directory (so the decode glob never crosses schemes) but shares the CSV.

cd "$(dirname "$0")/.."

SCHEME="${1:-bell}"
INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="grid_bell"          # one CSV for every scheme's grid data
MAX_SHOTS=1000000
MAX_ERRORS=5000

# p axis: log grid 1e-4 .. P_HI.  p_z axis: log grid 1e-3 .. PZ_HI.
# NP/NPZ = points per axis. Ranges are per-scheme.
case "$SCHEME" in
  bell)
    OUTPUT_DIR="out/circuits/grid_bell"
    SCHEME_ARGS="--noise-model biased_bell_v3"
    NP=25; NPZ=25; P_HI=1e-2; PZ_HI=3e-2
    ;;
  majority)
    OUTPUT_DIR="out/circuits/grid_bell_majority"
    SCHEME_ARGS="--scheme equivalent_v3_majority"
    NP=25; NPZ=25; P_HI=2e-3; PZ_HI=1e-1
    ;;
  *)
    echo "unknown scheme: $SCHEME (use 'bell' or 'majority')" >&2; exit 1 ;;
esac

P=$(uv run python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-4, np.log10($P_HI), $NP)))")
PZ=$(uv run python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-3, np.log10($PZ_HI), $NPZ)))")

echo "=== Stage 1: generate '$SCHEME' grid circuits (${NP} p x ${NPZ} p_z) ==="
echo "p_z grid: $PZ"
uv run python tools/make_grid_circuits_bell.py \
  --input "$INPUT_GLOB" \
  --output-dir "$OUTPUT_DIR" \
  $SCHEME_ARGS \
  --p-values $P \
  --p-bell-z-values $PZ

echo "=== Stage 2: decode -> out/collection/${PREAMBLE},<date>.csv ==="
uv run python tools/decode.py \
  --circuit-files "$OUTPUT_DIR"/*.stim \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --print-progress \
  --output-preamble "$PREAMBLE"
