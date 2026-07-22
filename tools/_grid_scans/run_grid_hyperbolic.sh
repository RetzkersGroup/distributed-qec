#!/bin/bash
set -e

# Grid scan for the hyperbolic Floquet code (biased_bell_v3, no repetition) on
# the same (p_res, p_z) grid as the honeycomb bell run. Reuses the generic
# generator (make_grid_circuits_bell.py) and decoder, writing to its own CSV
# (grid_hyperbolic) so the data stays separate from the honeycomb runs.
#
#   bash tools/_grid_scans/run_grid_hyperbolic.sh                 # base, d=2,4,6
#   bash tools/_grid_scans/run_grid_hyperbolic.sh full            # base, d=2..12
#   bash tools/_grid_scans/run_grid_hyperbolic.sh test majority   # 3-rep majority vote
#
# base (biased_bell_v3) and majority (equivalent_v3_majority = 3 repeated syndrome
# measurements + majority vote) both write to grid_hyperbolic.csv, distinguished
# by noise_model -- like the honeycomb bell + majority share grid_bell.csv.
# Uses plain `python` — run from an environment with the project installed.

cd "$(dirname "$0")/../.."

# Distance set. l=1,2,3 -> d=2,4,6 (small, for validating the pipeline);
# l=1..6 -> d=2,4,6,8,10,12 (the full run, much heavier to decode).
if [ "${1:-test}" = "full" ]; then
  DIST_GLOB='semi_hyperbolic_floquet_l=[123456]_*.stim'
else
  DIST_GLOB='semi_hyperbolic_floquet_l=[123]_*.stim'
fi

# Scheme: base or the 3-repetition majority-vote equivalent. Each scheme's
# circuits go to their own dir; both decode into grid_hyperbolic.csv.
case "${2:-base}" in
  base)
    # non-repeated baseline = single-ancilla 1-rep equivalent (equivalent_v3),
    # which reproduces the two-ancilla biased_bell_v3 exactly.
    OUTPUT_DIR="out/circuits/grid_hyperbolic"
    SCHEME_ARGS="--scheme equivalent_v3"
    PZ_HI=3e-2 ;;
  majority)
    # Repeated measurement lifts the p_z threshold well above the base ceiling,
    # so the majority scan needs a higher p_z top (cf. honeycomb majority).
    OUTPUT_DIR="out/circuits/grid_hyperbolic_majority"
    SCHEME_ARGS="--scheme equivalent_v3_majority"
    PZ_HI=1e-1 ;;
  majority5)
    # 5-rep majority vote: even higher p_z threshold, so a higher p_z ceiling.
    OUTPUT_DIR="out/circuits/grid_hyperbolic_majority5"
    SCHEME_ARGS="--scheme equivalent_v3_majority5"
    PZ_HI=2e-1 ;;
  *) echo "unknown scheme: ${2} (use base|majority|majority5)" >&2; exit 1 ;;
esac

INPUT_GLOB="in/circuits/hyperbolic_floquet/8.8.8/**/$DIST_GLOB"
SCALED_DIR="out/circuits/hyperbolic_3d_inputs"
PREAMBLE="grid_hyperbolic"
ROUNDS_MULTIPLE=3          # run each code for ROUNDS_MULTIPLE * d rounds
MAX_SHOTS=1000000
MAX_ERRORS=5000

# p_res 1e-4..1e-2, 25x25; p_z 1e-3..PZ_HI (per-scheme, set above).
NP=25; NPZ=25; P_HI=1e-2
P=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-4, np.log10($P_HI), $NP)))")
PZ=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-3, np.log10($PZ_HI), $NPZ)))")

echo "=== Stage 0: rescale each code to ${ROUNDS_MULTIPLE}d rounds ==="
rm -rf "$SCALED_DIR"          # avoid stale distances from a previous run
python tools/_grid_scans/scale_rounds.py \
  --input "$INPUT_GLOB" \
  --output-dir "$SCALED_DIR" \
  --rounds-multiple "$ROUNDS_MULTIPLE"

echo "=== Stage 1: generate hyperbolic grid circuits ==="
rm -rf "$OUTPUT_DIR"          # clean slate so the decode glob never mixes runs
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$SCALED_DIR/*.stim" \
  --output-dir "$OUTPUT_DIR" \
  $SCHEME_ARGS \
  --p-values $P \
  --p-bell-z-values $PZ

echo "=== Stage 2: decode (resume) -> out/collection/${PREAMBLE}.csv ==="
python tools/decode.py \
  --circuit-files "$OUTPUT_DIR"/*.stim \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --print-progress \
  --output-preamble "$PREAMBLE" --resume
