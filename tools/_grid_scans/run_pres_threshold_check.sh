#!/bin/bash
set -e

# SANITY CHECK: probe the two axis-intercepts of the (p_res, p_z) boundary by
# pinning one axis to 0 and scanning the other, for all four code/scheme combos:
#
#   1. Honeycomb code                          (biased_bell_v3)
#   2. Honeycomb + 3-rep majority vote          (equivalent_v3_majority)
#   3. Hyperbolic Floquet code                 (biased_bell_v3)
#   4. Hyperbolic Floquet + 3-rep majority vote
#
# Scan A  (p_z = 0,   scan p_res) -> pure p_res threshold  -> out/collection/pres_check.csv
# Scan B  (p_res = 0, scan p_z)   -> pure p_z   threshold  -> out/collection/pz_check.csv
#
# These are the p_z->0 and p_res->0 intercepts of the boundary curves, meant to
# be reused in the final boundary figure. Expected p_z thresholds: ~1.5% (base),
# ~7% (3-rep majority vote).
#
#   bash tools/_grid_scans/run_pres_threshold_check.sh          # hyperbolic d=2,4,6
#   bash tools/_grid_scans/run_pres_threshold_check.sh full     # hyperbolic d=2..12
#
# Each scan keeps its own dedicated CSV (separate from every previous run) so it
# can be inspected directly; the four combos coexist in each CSV and stay
# distinguishable by input_circuit + noise_model. Decoding uses --resume, so
# re-invoking only tops up toward the budget.

cd "$(dirname "$0")/../.."

# Stream child stdout/stderr live (so sinter's --print-progress table shows in
# real time, e.g. under nohup) instead of block-buffering.
export PYTHONUNBUFFERED=1

MODE="${1:-test}"

# Circuit families ------------------------------------------------------------
HONEY_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"   # d=3,5,7

if [ "$MODE" = "full" ]; then
  HYP_DIST='semi_hyperbolic_floquet_l=[123456]_*.stim'   # d=2..12
else
  HYP_DIST='semi_hyperbolic_floquet_l=[123]_*.stim'      # d=2,4,6
fi
HYP_GLOB="in/circuits/hyperbolic_floquet/8.8.8/**/$HYP_DIST"

# Sampling --------------------------------------------------------------------
ROUNDS_MULTIPLE=3              # hyperbolic: run each code for 3*d rounds (as grid)
MAX_SHOTS=1000000
MAX_ERRORS=5000

# Scan B pins p_res to a negligible floor rather than exactly 0: at p_res=0 the
# majority gadget's error set leaves part of the matching graph with no boundary
# connection, so pymatching hits "no perfect matching". A tiny p_res (~4 orders
# below the p_z threshold) regularizes the graph without shifting the threshold.
P_FLOOR=1e-6

# Scan A axis: p_res log grid (brackets base ~2-3e-3, majority ~1-2e-3).
NP=25
P=$(uv run python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10(1e-4), np.log10(5e-3), $NP)))")
# Scan B axis: p_z log grids. Base thresholds ~1.5%, majority ~7%, so the two
# schemes get different tops (as in the (p_res, p_z) grid runs).
NPZ=25
PZ_BASE=$(uv run python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10(1e-3), np.log10(3e-2), $NPZ)))")
PZ_MAJ=$(uv run python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10(3e-3), np.log10(1.2e-1), $NPZ)))")
# 5-rep majority sits even higher (~10-13%), so a wider p_z top.
PZ_MAJ5=$(uv run python -c "import numpy as np; print(' '.join(f'{x:.10g}' for x in np.logspace(np.log10(5e-3), np.log10(2e-1), $NPZ)))")

# Scale the hyperbolic codes to 3d rounds once (honeycomb is used as-is, exactly
# like the grid runs, so the intercepts line up with the boundary curves).
HYP_SCALED="out/circuits/pres_check_hyperbolic_3d_inputs"
echo "=== Stage 0: rescale hyperbolic codes to ${ROUNDS_MULTIPLE}d rounds ==="
rm -rf "$HYP_SCALED"
uv run python tools/_grid_scans/scale_rounds.py \
  --input "$HYP_GLOB" \
  --output-dir "$HYP_SCALED" \
  --rounds-multiple "$ROUNDS_MULTIPLE"

# One combo: generate circuits over a single scan line, then decode/resume.
#   $1=label  $2=preamble(CSV)  $3=input-glob  $4=output-dir  $5..=make_grid args
COMBO_N=12          # total run_combo calls (6 per scan x 2 scans), for progress
COMBO_I=0
run_combo() {
  COMBO_I=$((COMBO_I + 1))
  local label="$1" preamble="$2" input="$3" outdir="$4"; shift 4
  local tag="[combo ${COMBO_I}/${COMBO_N}] [$label -> $preamble]"
  echo
  echo "=== $tag generate circuits ($(date +%H:%M:%S)) ==="
  rm -rf "$outdir"
  uv run python tools/_grid_scans/make_grid_circuits_bell.py \
    --input "$input" --output-dir "$outdir" "$@"
  echo "=== $tag decode (resume, live sinter progress) ($(date +%H:%M:%S)) ==="
  uv run python tools/decode.py \
    --circuit-files "$outdir"/*.stim \
    --max-shots "$MAX_SHOTS" --max-errors "$MAX_ERRORS" \
    --print-progress \
    --output-preamble "$preamble" --resume
}

echo
echo "############################################################"
echo "# Scan A: p_z = 0, scan p_res  -> pres_check.csv"
echo "############################################################"
run_combo "honeycomb 1-rep"      pres_check "$HONEY_GLOB"        out/circuits/pres_check_honeycomb_single     --scheme equivalent_v3           --p-values $P --p-bell-z-values 0.0
run_combo "honeycomb majority"   pres_check "$HONEY_GLOB"        out/circuits/pres_check_honeycomb_majority   --scheme equivalent_v3_majority  --p-values $P --p-bell-z-values 0.0
run_combo "honeycomb majority5"  pres_check "$HONEY_GLOB"        out/circuits/pres_check_honeycomb_majority5  --scheme equivalent_v3_majority5 --p-values $P --p-bell-z-values 0.0
run_combo "hyperbolic 1-rep"     pres_check "$HYP_SCALED/*.stim" out/circuits/pres_check_hyperbolic_single    --scheme equivalent_v3           --p-values $P --p-bell-z-values 0.0
run_combo "hyperbolic majority"  pres_check "$HYP_SCALED/*.stim" out/circuits/pres_check_hyperbolic_majority  --scheme equivalent_v3_majority  --p-values $P --p-bell-z-values 0.0
run_combo "hyperbolic majority5" pres_check "$HYP_SCALED/*.stim" out/circuits/pres_check_hyperbolic_majority5 --scheme equivalent_v3_majority5 --p-values $P --p-bell-z-values 0.0

echo
echo "############################################################"
echo "# Scan B: p_res = 0, scan p_z  -> pz_check.csv"
echo "############################################################"
run_combo "honeycomb 1-rep"      pz_check "$HONEY_GLOB"        out/circuits/pz_check_honeycomb_single     --scheme equivalent_v3           --p-values $P_FLOOR --p-bell-z-values $PZ_BASE
run_combo "honeycomb majority"   pz_check "$HONEY_GLOB"        out/circuits/pz_check_honeycomb_majority   --scheme equivalent_v3_majority  --p-values $P_FLOOR --p-bell-z-values $PZ_MAJ
run_combo "honeycomb majority5"  pz_check "$HONEY_GLOB"        out/circuits/pz_check_honeycomb_majority5  --scheme equivalent_v3_majority5 --p-values $P_FLOOR --p-bell-z-values $PZ_MAJ5
run_combo "hyperbolic 1-rep"     pz_check "$HYP_SCALED/*.stim" out/circuits/pz_check_hyperbolic_single    --scheme equivalent_v3           --p-values $P_FLOOR --p-bell-z-values $PZ_BASE
run_combo "hyperbolic majority"  pz_check "$HYP_SCALED/*.stim" out/circuits/pz_check_hyperbolic_majority  --scheme equivalent_v3_majority  --p-values $P_FLOOR --p-bell-z-values $PZ_MAJ
run_combo "hyperbolic majority5" pz_check "$HYP_SCALED/*.stim" out/circuits/pz_check_hyperbolic_majority5 --scheme equivalent_v3_majority5 --p-values $P_FLOOR --p-bell-z-values $PZ_MAJ5

echo
echo "############################################################"
echo "# Analyze: fit the threshold along each scan's axis"
echo "############################################################"
uv run python tools/_grid_scans/analyze_axis_threshold.py --csv out/collection/pres_check.csv
uv run python tools/_grid_scans/analyze_axis_threshold.py --csv out/collection/pz_check.csv

echo
echo "Done."
echo "  Scan A raw: out/collection/pres_check.csv | figure: out/figures/pres_check.pdf"
echo "  Scan B raw: out/collection/pz_check.csv   | figure: out/figures/pz_check.pdf"
