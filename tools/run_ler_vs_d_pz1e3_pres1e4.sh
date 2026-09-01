#!/bin/bash
set -e

# LER vs distance at a single fixed noise point (p_z = 1e-3, p_res = 1e-4),
# comparing schemes x decoders on the honeycomb-Gidney inputs (d = 3, 5, 7).
# Companion figure: notebooks/ler_vs_d_pz1e-3_pres1e-4.ipynb.
#
# All results accumulate in ONE sinter-format resume CSV
# out/collection/$PREAMBLE.csv (strong_id encodes the decoder, and the
# notebook groups curves by (noise_model, decoder), so sharing a CSV is safe
# here -- unlike the grid-threshold pipeline). Re-running only samples the
# shortfall, so raise the budgets below and re-run to extend statistics.
#
# To add a scheme or decoder later, append a line to SCHEMES / DECODERS.

cd "$(dirname "$0")/.."

# Noise point, overridable from the environment (use scientific notation so
# the derived CSV preamble stays readable):
#   P_Z=1e-2 bash tools/run_ler_vs_d_pz1e3_pres1e4.sh pymatching
P_Z="${P_Z:-1e-3}"     # --p-bell-z-values (Bell Z error rate)
P_RES="${P_RES:-1e-4}" # --p-values (residual noise p)
INPUT_GLOB="in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim"
PREAMBLE="ler_vs_d_pz${P_Z}_pres${P_RES}"
PY=.venv/bin/python

# name|circuit_dir|scheme flags for make_grid_circuits_bell.py
# "bell" follows the project convention: equivalent_v3 reproduces the
# two-ancilla biased_bell_v3 exactly (see run_grid_bell.sh).
SCHEMES=(
  "bell|out/circuits/grid_bell|--scheme equivalent_v3"
  "majority|out/circuits/grid_bell_majority|--scheme equivalent_v3_majority"
)

# decoder|max_shots|max_errors -- fast decoders first. belief_matching and
# tesseract_short are orders of magnitude slower per shot, so they get
# smaller first-pass budgets; bump and re-run to extend.
DECODERS=(
  "pymatching|5000000|2000"
  "pymatching-correlated|5000000|2000"
  "belief_matching|500000|500"
  "tesseract_short|500000|500"
)

# Optional: pass decoder names as args to run only those, e.g. on a bigger
# remote machine: `bash tools/run_ler_vs_d_pz1e3_pres1e4.sh tesseract_short`.
# A fresh remote checkout has no resume CSV, so its output contains only the
# decoders it ran; copy it back next to the local one as
# out/collection/ler_vs_d_pz1e-3_pres1e-4,<host>.csv -- the notebook globs
# ler_vs_d_pz1e-3_pres1e-4*.csv and sums batches per strong_id. (Don't copy
# the local CSV to the remote first, or its rows would be double counted.)
if [ "$#" -gt 0 ]; then
  filtered=()
  for dentry in "${DECODERS[@]}"; do
    for want in "$@"; do
      [[ "$dentry" == "$want|"* ]] && filtered+=("$dentry")
    done
  done
  [ "${#filtered[@]}" -gt 0 ] || { echo "no matching decoders: $*" >&2; exit 1; }
  DECODERS=("${filtered[@]}")
fi

# Circuit filenames embed the grid point via %.5g formatting.
POINT="pz=$(printf '%.5g' "$P_Z"),p=$(printf '%.5g' "$P_RES")"

echo "=== Stage 1: generate circuits at ($POINT) ==="
for entry in "${SCHEMES[@]}"; do
  IFS='|' read -r name dir scheme_flags <<< "$entry"
  echo "--- scheme: $name -> $dir"
  $PY tools/_grid_scans/make_grid_circuits_bell.py \
    --input "$INPUT_GLOB" \
    --output-dir "$dir" \
    $scheme_flags \
    --p-values "$P_RES" \
    --p-bell-z-values "$P_Z" \
    --skip-existing
done

echo "=== Stage 2: decode into out/collection/$PREAMBLE.csv ==="
for dentry in "${DECODERS[@]}"; do
  IFS='|' read -r decoder max_shots max_errors <<< "$dentry"
  for entry in "${SCHEMES[@]}"; do
    IFS='|' read -r name dir scheme_flags <<< "$entry"
    echo "--- decoder: $decoder, scheme: $name (max_shots=$max_shots, max_errors=$max_errors)"
    $PY tools/decode.py \
      --circuit-files "$dir"/*",$POINT.stim" \
      --decoder "$decoder" \
      --max-shots "$max_shots" \
      --max-errors "$max_errors" \
      --output-preamble "$PREAMBLE" \
      --print-progress
  done
done

echo "=== Done: out/collection/$PREAMBLE.csv ==="
