#!/bin/bash
set -e

# Tesseract (short beam) pure-p_z check for the full decoder, r=1 and r=3,
# restricted to p_z < 1e-2 where syndromes are sparse enough for tesseract's
# MLE search to run in reasonable time (it stalls at threshold-density
# syndromes, and relay-bp shows a convergence-failure error floor at low p_z
# on the larger distances — tesseract is the trustworthy decoder down here).
#
#   bash tools/_grid_scans/run_pz_check_tesseract_lowpz.sh                 # 1k shots
#   MAX_SHOTS=10000 bash tools/_grid_scans/run_pz_check_tesseract_lowpz.sh # top up
#
# Reuses the full_check_pz_rep{1,3} circuits where their grids fall below
# 1e-2 (--skip-existing) and decodes into
# out/collection/full_check_pz_tesseract.csv (resume -> raising MAX_SHOTS
# later only samples the shortfall).

cd "$(dirname "$0")/../.."

export PYTHONUNBUFFERED=1

INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="full_check_pz_tesseract"
MAX_SHOTS="${MAX_SHOTS:-1000}"
MAX_ERRORS="${MAX_ERRORS:-500}"

P_FLOOR=1e-6
NPZ=25

# The same log grids as the relay-bp/pymatching checks, truncated at 1e-2 so
# the points line up with the other CSVs' grids for direct comparison.
PZ_R1=$(python -c "import numpy as np; g=np.logspace(np.log10(1e-3), np.log10(3e-2), $NPZ); print(' '.join(f'{x:.10g}' for x in g[g < 1e-2]))")
PZ_R3=$(python -c "import numpy as np; g=np.logspace(np.log10(5e-3), np.log10(1.5e-1), $NPZ); print(' '.join(f'{x:.10g}' for x in g[g < 1e-2]))")

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

# Decode only the below-1e-2 circuits: pz=0.00... names match p_z < 1e-2
# (the grids contain no p_z < 1e-3, so the pz=0.00 prefix selects exactly
# the truncated range).
for DIR in out/circuits/full_check_pz_rep1 out/circuits/full_check_pz_rep3; do
  echo "=== decode $DIR (p_z < 1e-2) with tesseract_short (resume) -> out/collection/${PREAMBLE}.csv ==="
  python tools/decode.py \
    --circuit-files "$DIR"/*pz=0.00*.stim \
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
