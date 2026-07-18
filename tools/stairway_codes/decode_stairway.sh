#!/bin/bash
set -e

# Decode the stairway slice circuits into a SHARED CSV (out/collection/grid_stairway.csv).
# Every scheme's slice data lives in one CSV (distinguished by the noise_model
# metadata column), so plot_grid.py can overlay them.
#
# Decoding uses --resume, so re-invoking only tops up toward the shot/error budget.
# Stairway decoding is slow, so the shot budget is kept modest here — bump
# MAX_SHOTS/MAX_ERRORS if you need tighter error bars.
#
#   bash tools/stairway_codes/decode_stairway.sh            # baseline (equivalent_v3)
#   bash tools/stairway_codes/decode_stairway.sh majority   # 3-rep majority vote
#   bash tools/stairway_codes/decode_stairway.sh majority5  # 5-rep majority vote

cd "$(dirname "$0")/../.."

# Stream child stdout/stderr live so sinter's --print-progress table shows in
# real time (e.g. under nohup) instead of block-buffering.
export PYTHONUNBUFFERED=1

SCHEME="${1:-bell}"
PREAMBLE="grid_stairway"      # one CSV for every scheme's slice data
MAX_SHOTS=1000
MAX_ERRORS=10

case "$SCHEME" in
  bell)      OUTPUT_DIR="out/circuits/grid_stairway" ;;
  majority)  OUTPUT_DIR="out/circuits/grid_stairway_majority" ;;
  majority5) OUTPUT_DIR="out/circuits/grid_stairway_majority5" ;;
  *) echo "unknown scheme: $SCHEME (use bell|majority|majority5)" >&2; exit 1 ;;
esac

echo "=== Decode '$SCHEME' (resume) -> out/collection/${PREAMBLE}.csv ==="
uv run python tools/decode.py \
  --circuit-files "$OUTPUT_DIR"/*.stim \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --decoder tesseract_short \
  --print-progress \
  --output-preamble "$PREAMBLE" --resume
