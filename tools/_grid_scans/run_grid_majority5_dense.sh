#!/bin/bash
set -e

# Densify the 5-rep majority-vote (p_res, p_z) grids around the p_z crossing,
# for BOTH codes (honeycomb + hyperbolic Floquet).
#
# Why: the 5-rep LER curves are so steep that on the standard 25-point p_z grid
# each distance jumps from pL~0.01 to saturation within 1-2 grid steps, leaving
# the FSS fit almost no points in the transition region (d=7 has a single one
# at p_res=1e-4 -- see notebooks/fit_debug.ipynb). This script ADDS p_z points
# inside the crossing band only, at the SAME p_res values as the existing scan
# (so the new points land on the existing LER panels), then tops up the shared
# CSVs via decode --resume. Also densifies the p_res=0 intercept scan
# (pz_check.csv, cf. run_pres_threshold_check.sh Scan B) through the same band.
# Nothing is deleted; re-running only samples the shortfall.
#
#   bash tools/_grid_scans/run_grid_majority5_dense.sh          # hyperbolic test set (d=2,4,6)
#   bash tools/_grid_scans/run_grid_majority5_dense.sh full     # hyperbolic full set (d=2..12)
#
# Uses plain `python` — run from an environment with the project installed.

cd "$(dirname "$0")/../.."

MAX_SHOTS=1000000
MAX_ERRORS=5000

# Crossing band to densify (5-rep crossings sit at p_z ~ 0.09-0.18 for both
# codes) and subdivision: each existing grid gap inside the band is split into
# N_SUB equal log steps (N_SUB-1 new interior points per gap).
PZ_BAND_LO=0.05
PZ_BAND_HI=0.2
N_SUB=3

# New p_z values: interior subdivision points of the existing majority5 grid
# (logspace(-3, log10(0.2), 25)) inside the band. Emitting only interior points
# means no existing circuit is regenerated.
PZ_DENSE=$(python -c "
import numpy as np
grid = np.logspace(-3, np.log10(0.2), 25)
out = []
for lo, hi in zip(grid[:-1], grid[1:]):
    if lo < $PZ_BAND_LO or hi > $PZ_BAND_HI * 1.0001:
        continue
    out.extend(np.geomspace(lo, hi, $N_SUB + 1)[1:-1])
print(' '.join(f'{x:.10g}' for x in out))
")

# p_res values: the SAME grids as the original runs (so panels line up),
# restricted to where the crossing exists (protection is lost near the p_z=0
# intercept p_res* ~ 2.5e-3 honeycomb / 1.5e-3 hyperbolic; beyond that the fit
# has nothing to find).
P_HONEY=$(python -c "import numpy as np; p = np.logspace(-4, np.log10(2e-3), 25); print(' '.join(f'{x:.10g}' for x in p))")
P_HYP=$(python -c "import numpy as np; p = np.logspace(-4, -2, 25); print(' '.join(f'{x:.10g}' for x in p[p <= 1.6e-3]))")

echo "dense p_z values: $PZ_DENSE"
echo

# ── Honeycomb ────────────────────────────────────────────────────────────────
HONEY_INPUT="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
HONEY_DIR="out/circuits/grid_bell_majority5"       # same dir as run_grid_bell.sh majority5

echo "=== Honeycomb majority5: generate dense-band circuits ==="
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$HONEY_INPUT" \
  --output-dir "$HONEY_DIR" \
  --scheme equivalent_v3_majority5 \
  --p-values $P_HONEY \
  --p-bell-z-values $PZ_DENSE \
  --skip-existing

echo "=== Honeycomb majority5: decode (resume) -> out/collection/grid_bell.csv ==="
python tools/decode.py \
  --circuit-files "$HONEY_DIR"/*.stim \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --print-progress \
  --output-preamble grid_bell --resume

# ── Hyperbolic Floquet ───────────────────────────────────────────────────────
# Distance set mirrors run_grid_hyperbolic.sh: test = l=1,2,3 (d=2,4,6, matches
# the currently collected data), full = l=1..6 (d=2..12).
if [ "${1:-test}" = "full" ]; then
  DIST_GLOB='semi_hyperbolic_floquet_l=[123456]_*.stim'
else
  DIST_GLOB='semi_hyperbolic_floquet_l=[123]_*.stim'
fi
HYP_INPUT="in/circuits/hyperbolic_floquet/8.8.8/**/$DIST_GLOB"
SCALED_DIR="out/circuits/hyperbolic_3d_inputs"
HYP_DIR="out/circuits/grid_hyperbolic_majority5"   # same dir as run_grid_hyperbolic.sh majority5
ROUNDS_MULTIPLE=3

echo "=== Hyperbolic majority5: rescale inputs to ${ROUNDS_MULTIPLE}d rounds ==="
rm -rf "$SCALED_DIR"          # avoid stale distances from a previous run
python tools/_grid_scans/scale_rounds.py \
  --input "$HYP_INPUT" \
  --output-dir "$SCALED_DIR" \
  --rounds-multiple "$ROUNDS_MULTIPLE"

echo "=== Hyperbolic majority5: generate dense-band circuits ==="
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$SCALED_DIR/*.stim" \
  --output-dir "$HYP_DIR" \
  --scheme equivalent_v3_majority5 \
  --p-values $P_HYP \
  --p-bell-z-values $PZ_DENSE \
  --skip-existing

echo "=== Hyperbolic majority5: decode (resume) -> out/collection/grid_hyperbolic.csv ==="
python tools/decode.py \
  --circuit-files "$HYP_DIR"/*.stim \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --print-progress \
  --output-preamble grid_hyperbolic --resume

# ── p_res = 0 intercept (pz_check) ───────────────────────────────────────────
# The pure-p_z intercept scan (run_pres_threshold_check.sh Scan B) has the same
# sparsity problem through the majority5 crossing, so densify it too. Mirrors
# the check exactly: p_res pinned to its P_FLOOR (1e-6, not 0 -- exact 0 breaks
# the matching graph), dense points computed from ITS majority5 p_z grid
# (logspace 5e-3..0.2), decoded into pz_check.csv. New circuit dirs (the check
# script rm -rf's its own dirs on re-run; these must survive that).
P_FLOOR=1e-6
PZ_CHECK_DENSE=$(python -c "
import numpy as np
grid = np.logspace(np.log10(5e-3), np.log10(0.2), 25)
out = []
for lo, hi in zip(grid[:-1], grid[1:]):
    if lo < $PZ_BAND_LO or hi > $PZ_BAND_HI * 1.0001:
        continue
    out.extend(np.geomspace(lo, hi, $N_SUB + 1)[1:-1])
print(' '.join(f'{x:.10g}' for x in out))
")
echo
echo "dense pz_check values: $PZ_CHECK_DENSE"

echo "=== Honeycomb majority5 intercept: generate dense-band circuits ==="
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$HONEY_INPUT" \
  --output-dir "out/circuits/pz_check_honeycomb_majority5_dense" \
  --scheme equivalent_v3_majority5 \
  --p-values $P_FLOOR \
  --p-bell-z-values $PZ_CHECK_DENSE \
  --skip-existing

echo "=== Honeycomb majority5 intercept: decode (resume) -> out/collection/pz_check.csv ==="
python tools/decode.py \
  --circuit-files out/circuits/pz_check_honeycomb_majority5_dense/*.stim \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --print-progress \
  --output-preamble pz_check --resume

echo "=== Hyperbolic majority5 intercept: generate dense-band circuits ==="
# Reuses $SCALED_DIR (same 3d-rounds rescale as the check's own stage 0).
python tools/_grid_scans/make_grid_circuits_bell.py \
  --input "$SCALED_DIR/*.stim" \
  --output-dir "out/circuits/pz_check_hyperbolic_majority5_dense" \
  --scheme equivalent_v3_majority5 \
  --p-values $P_FLOOR \
  --p-bell-z-values $PZ_CHECK_DENSE \
  --skip-existing

echo "=== Hyperbolic majority5 intercept: decode (resume) -> out/collection/pz_check.csv ==="
python tools/decode.py \
  --circuit-files out/circuits/pz_check_hyperbolic_majority5_dense/*.stim \
  --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
  --print-progress \
  --output-preamble pz_check --resume

echo
echo "Done. Re-run the plot stage (tools/run_grid_all.sh) to refit with the dense band."
