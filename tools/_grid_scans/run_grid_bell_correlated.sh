#!/bin/bash
set -e

# Correlated-matching twin of run_grid_bell.sh: same honeycomb (p, p_z) grid
# circuits (local / majority-vote schemes), decoded with pymatching's native
# correlated MWPM (--decoder pymatching-correlated) into a SEPARATE CSV
# out/collection/grid_bell_correlated.csv.
#
# A separate CSV is mandatory: sinter's strong_id encodes the decoder and
# plot_grid.load_grid groups by strong_id while the threshold fit filters only
# on noise_model, so correlated rows in grid_bell.csv would double every grid
# point and corrupt the fits.
#
#   bash tools/_grid_scans/run_grid_bell_correlated.sh            # bell (r=1 base)
#   bash tools/_grid_scans/run_grid_bell_correlated.sh majority   # majority vote r=3
#   bash tools/_grid_scans/run_grid_bell_correlated.sh majority5  # majority vote r=5
#
# Circuits are decoder-independent, so this reuses the existing grid_bell*
# circuit directories (--skip-existing regenerates only what's missing, e.g.
# the r=3 majority set). Plot: plot_full_vs_local.py --correlated.

cd "$(dirname "$0")/../.."

SCHEME="${1:-bell}"
INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="grid_bell_correlated"   # one CSV for every scheme's correlated grid data
DECODER="pymatching-correlated"
MAX_SHOTS=1000000
MAX_ERRORS=5000

case "$SCHEME" in
  bell)
    OUTPUT_DIR="out/circuits/grid_bell"
    SCHEME_ARGS="--scheme equivalent_v3"
    NP=25; NPZ=25; P_HI=1e-2; PZ_HI=3e-2
    ;;
  majority)
    OUTPUT_DIR="out/circuits/grid_bell_majority"
    SCHEME_ARGS="--scheme equivalent_v3_majority"
    NP=25; NPZ=25; P_HI=2e-3; PZ_HI=1e-1
    ;;
  majority5)
    OUTPUT_DIR="out/circuits/grid_bell_majority5"
    SCHEME_ARGS="--scheme equivalent_v3_majority5"
    NP=25; NPZ=25; P_HI=2e-3; PZ_HI=2e-1
    ;;
  *)
    echo "unknown scheme: $SCHEME (use bell|majority|majority5)" >&2; exit 1 ;;
esac

P=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-4, np.log10($P_HI), $NP)))")
PZ=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(-3, np.log10($PZ_HI), $NPZ)))")

echo "=== Stage 1: generate '$SCHEME' grid circuits (${NP} p x ${NPZ} p_z) ==="
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$INPUT_GLOB" \
  --output-dir "$OUTPUT_DIR" \
  $SCHEME_ARGS \
  --p-values $P \
  --p-bell-z-values $PZ \
  --skip-existing

echo "=== Stage 2: decode ($DECODER, resume) -> out/collection/${PREAMBLE}.csv ==="
python tools/decode.py \
  --circuit-files "$OUTPUT_DIR"/*.stim \
  --decoder "$DECODER" \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --print-progress \
  --output-preamble "$PREAMBLE" --resume
