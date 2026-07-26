#!/bin/bash
set -e

# Relay-BP pure-p_z check for the full decoder, r=1 and r=3 only.
#
# Under pymatching the full-decoder p_z gain saturates (~1.8x vs majority's
# 5x) because pure-p_z logical errors are non-graphlike; tesseract (MLE
# search) is intractable at threshold-density syndromes. Relay-BP decodes the
# raw (undecomposed) DEM natively and is fast, so this scan measures the p_z
# threshold the full decoder actually achieves with a hyperedge-aware decoder.
# Spot check (d=3, r=3, pz=6.4%): relay-bp 0/500 errors vs pymatching 179/500.
# r=1 is included as the fair same-decoder baseline.
#
#   bash tools/_grid_scans/run_pz_check_relay_bp.sh
#   MAX_SHOTS=1000 MAX_ERRORS=50 bash tools/_grid_scans/run_pz_check_relay_bp.sh   # quick sanity pass
#
# Reuses the full_check_pz_rep{1,3} circuits (--skip-existing) and decodes into
# out/collection/full_check_pz_relay_bp.csv (one CSV per decoder, never mixed).
# Decoding uses --resume, so a later full-budget run tops up the sanity pass.

cd "$(dirname "$0")/../.."

export PYTHONUNBUFFERED=1

INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="full_check_pz_relay_bp"
MAX_SHOTS="${MAX_SHOTS:-100000}"
MAX_ERRORS="${MAX_ERRORS:-500}"

P_FLOOR=1e-6
NPZ=25

# r=1 brackets the baseline threshold (~1.4% under pymatching).
PZ_R1=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10(1e-3), np.log10(3e-2), $NPZ)))")
# r=3: wide top — the relay-bp threshold may sit at or above the 3-rep
# majority value (~7%).
PZ_R3=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10(5e-3), np.log10(1.5e-1), $NPZ)))")

echo "=== [r=1] generate p_z-scan circuits (skip-existing) ==="
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$INPUT_GLOB" \
  --output-dir out/circuits/full_check_pz_rep1 \
  --scheme equivalent_v3 \
  --p-values $P_FLOOR \
  --p-bell-z-values $PZ_R1 \
  --skip-existing

echo "=== [r=3] generate p_z-scan circuits (skip-existing) ==="
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$INPUT_GLOB" \
  --output-dir out/circuits/full_check_pz_rep3 \
  --scheme equivalent_v3 \
  --repeat 3 \
  --p-values $P_FLOOR \
  --p-bell-z-values $PZ_R3 \
  --skip-existing

for DIR in out/circuits/full_check_pz_rep1 out/circuits/full_check_pz_rep3; do
  echo "=== decode $DIR with relay-bp (resume) -> out/collection/${PREAMBLE}.csv ==="
  python tools/decode.py \
    --circuit-files "$DIR"/*.stim \
    --decoder relay-bp \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
    --print-progress \
    --output-preamble "$PREAMBLE" --resume
done

echo
echo "=== Analyze ==="
python tools/_grid_scans/analyze_full_check.py \
  --csv out/collection/${PREAMBLE}.csv \
  --out out/figures/full_check_pz_relay_bp.pdf

echo
echo "Done. Data: out/collection/${PREAMBLE}.csv | figure: out/figures/full_check_pz_relay_bp.pdf"
