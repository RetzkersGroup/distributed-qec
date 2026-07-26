#!/bin/bash
set -e

# Tesseract (short beam) pure-p_z check for the full decoder, r=1 and r=3 only.
#
# The pymatching full-decoder gains saturate (~1.8x) because pure-p_z logical
# errors are non-graphlike: MWPM can't exploit the repetition chain after
# hyperedge decomposition. Tesseract decodes the hyperedges natively, so this
# scan tests whether the full decoder's DEM actually carries the majority-level
# (or better) suppression. r=1 is included as the fair same-decoder baseline.
#
#   bash tools/_grid_scans/run_pz_check_tesseract.sh
#
# p_res is pinned to the usual 1e-6 floor; only p_z is scanned. Decodes into
# out/collection/full_check_pz_tesseract.csv (separate CSV so pymatching and
# tesseract results never mix in one file). Tesseract is slow -- the shot
# budget is deliberately small; re-run to top up (--resume).

cd "$(dirname "$0")/../.."

export PYTHONUNBUFFERED=1

INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="full_check_pz_tesseract"
MAX_SHOTS=20000
MAX_ERRORS=300

P_FLOOR=1e-6
NPZ=25

# r=1 scan range brackets the baseline threshold (~1.4% under pymatching).
PZ_R1=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10(1e-3), np.log10(3e-2), $NPZ)))")
# r=3 range: wide top, in case tesseract lifts the threshold toward the
# 3-rep majority value (~7%) or beyond.
PZ_R3=$(python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10(5e-3), np.log10(1.5e-1), $NPZ)))")

echo "=== [r=1] generate p_z-scan circuits ==="
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
  echo "=== decode $DIR with tesseract_short (resume) -> out/collection/${PREAMBLE}.csv ==="
  python tools/decode.py \
    --circuit-files "$DIR"/*.stim \
    --decoder tesseract_short \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
    --print-progress \
    --output-preamble "$PREAMBLE" --resume
done

echo
echo "=== Analyze ==="
python tools/_grid_scans/analyze_full_check.py \
  --csv out/collection/${PREAMBLE}.csv \
  --out out/figures/full_check_pz_tesseract.pdf

echo
echo "Done. Data: out/collection/${PREAMBLE}.csv | figure: out/figures/full_check_pz_tesseract.pdf"
