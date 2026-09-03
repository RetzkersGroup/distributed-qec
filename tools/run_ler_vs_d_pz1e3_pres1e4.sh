#!/bin/bash
set -e

# LER vs distance at a single fixed noise point (p_z = 1e-3, p_res = 1e-4),
# comparing schemes x decoders on the honeycomb-Gidney inputs (d = 3, 5, 7 by
# default; override with the DISTANCES env var).
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
PREAMBLE="ler_vs_d_pz${P_Z}_pres${P_RES}"
PY=.venv/bin/python

# Which code distances to include, overridable from the environment; each maps
# to a honeycomb-Gidney input circuit (rounds r). Default keeps the original
# d = 3, 5, 7. The larger inputs (d = 9, 11) were added later.
#   DISTANCES="3 5 7 9 11" bash tools/run_ler_vs_d_pz1e3_pres1e4.sh
DISTANCES="${DISTANCES:-3 5 7}"
INPUT_DIR="in/circuits/honeycomb_gidney"

# distance -> rounds (r) for the honeycomb-Gidney inputs. case keeps this
# portable to bash 3.2 (macOS), which has no associative arrays.
d_to_r() {
  case "$1" in
    3) echo 10 ;; 5) echo 16 ;; 7) echo 22 ;; 9) echo 28 ;; 11) echo 34 ;;
    *) return 1 ;;
  esac
}

INPUTS=()  # input .stim paths for the selected distances
RVALS=()   # matching rounds, aligned by position with DISTANCES
for d in $DISTANCES; do
  r="$(d_to_r "$d")" || { echo "unknown distance: $d (known: 3 5 7 9 11)" >&2; exit 1; }
  match=("$INPUT_DIR"/model=SDEM3,r="$r",*obs=H*.stim)
  [ -e "${match[0]}" ] || { echo "no input circuit for d=$d (r=$r)" >&2; exit 1; }
  INPUTS+=("${match[@]}")
  RVALS+=("$r")
done

# name|circuit_dir|scheme flags for make_grid_circuits_bell.py
# "bell" follows the project convention: equivalent_v3 reproduces the
# two-ancilla biased_bell_v3 exactly (see run_grid_bell.sh).
# "bell_rep3" is structural 3x measurement repetition with full decoding
# (repeat_measurements before translation; dir shared with run_grid_full.sh).
SCHEMES=(
  "bell|out/circuits/grid_bell|--scheme equivalent_v3"
  "majority|out/circuits/grid_bell_majority|--scheme equivalent_v3_majority"
  "bell_rep3|out/circuits/grid_full_rep3|--scheme equivalent_v3 --repeat 3"
)

# decoder|max_shots|max_errors -- fast decoders first. belief_matching,
# tesseract_short and the beam* configs are orders of magnitude slower per shot,
# so they get smaller first-pass budgets; bump and re-run to extend.
#
# The beam* configs are IonQ's BP + beam-search qLDPC decoder (Table I of
# arXiv:2512.07057; see BEAMSEARCH_CONFIGS in tools/decode.py for the params).
# It is NOT pip-installable: clone it once and build the extension into .venv
# (kept out of git):
#   git clone https://github.com/ionq-publications/BeamSearchDecoder.git BeamSearchDecoder
#   .venv/bin/python -m pip install "Cython>=3.0.0" "setuptools>=60.0.0"
#   (cd BeamSearchDecoder/decoder && ../../.venv/bin/python setup.py build_ext --inplace)
#   SP=$(.venv/bin/python -c "import site; print(site.getsitepackages()[0])")
#   echo "$PWD/BeamSearchDecoder" > "$SP/beamsearchdecoder.pth"
DECODERS=(
  "pymatching|5000000|2000"
  "pymatching-correlated|5000000|2000"
  "belief_matching|500000|500"
  "tesseract_short|500000|500"
  "tesseract_long|500|500"
  "beam64_32res_640iters|500|50"
)

# Beam-search configs (Table I, arXiv:2512.07057) selectable on the command line
# even when not pre-listed in DECODERS above; keep in sync with BEAMSEARCH_CONFIGS
# in tools/decode.py. When passed as an arg they run with this default budget.
BEAM_CONFIGS="beam8_230iters beam32_340iters beam64_640iters beam64_32res_640iters"
BEAM_SHOTS=500
BEAM_ERRORS=50

# Optional: pass scheme and/or decoder names as args to run only those, e.g.
#   bash tools/run_ler_vs_d_pz1e3_pres1e4.sh tesseract_short          # one decoder, all schemes
#   bash tools/run_ler_vs_d_pz1e3_pres1e4.sh bell_rep3                # one scheme, all decoders
#   bash tools/run_ler_vs_d_pz1e3_pres1e4.sh bell_rep3 pymatching     # one scheme x one decoder
# Copy CSVs from remote machines back next to the local one as
# out/collection/${PREAMBLE},<date>,<host>.csv -- the notebook globs
# ${PREAMBLE}*.csv and sums batches per strong_id.
if [ "$#" -gt 0 ]; then
  sel_schemes=()
  sel_decoders=()
  for want in "$@"; do
    hit=""
    for entry in "${SCHEMES[@]}"; do
      if [[ "$entry" == "$want|"* ]]; then sel_schemes+=("$entry"); hit=1; fi
    done
    for dentry in "${DECODERS[@]}"; do
      if [[ "$dentry" == "$want|"* ]]; then sel_decoders+=("$dentry"); hit=1; fi
    done
    # Beam-search config not pre-listed in DECODERS: synthesize a default budget.
    if [ -z "$hit" ]; then
      for cfg in $BEAM_CONFIGS; do
        if [ "$cfg" = "$want" ]; then
          sel_decoders+=("$want|$BEAM_SHOTS|$BEAM_ERRORS"); hit=1
        fi
      done
    fi
    if [ -z "$hit" ]; then echo "unknown scheme/decoder: $want" >&2; exit 1; fi
  done
  if [ "${#sel_schemes[@]}" -gt 0 ]; then SCHEMES=("${sel_schemes[@]}"); fi
  if [ "${#sel_decoders[@]}" -gt 0 ]; then DECODERS=("${sel_decoders[@]}"); fi
fi

# Circuit filenames embed the grid point via %.5g formatting.
POINT="pz=$(printf '%.5g' "$P_Z"),p=$(printf '%.5g' "$P_RES")"

echo "=== Stage 1: generate circuits at ($POINT) for d = $DISTANCES ==="
for entry in "${SCHEMES[@]}"; do
  IFS='|' read -r name dir scheme_flags <<< "$entry"
  echo "--- scheme: $name -> $dir"
  # --input takes one path at a time (no brace expansion), so generate per input.
  for in_stim in "${INPUTS[@]}"; do
    $PY tools/_grid_scans/make_grid_circuits_bell.py \
      --input "$in_stim" \
      --output-dir "$dir" \
      $scheme_flags \
      --p-values "$P_RES" \
      --p-bell-z-values "$P_Z" \
      --skip-existing
  done
done

echo "=== Stage 2: decode into out/collection/$PREAMBLE.csv ==="
for dentry in "${DECODERS[@]}"; do
  IFS='|' read -r decoder max_shots max_errors <<< "$dentry"
  for entry in "${SCHEMES[@]}"; do
    IFS='|' read -r name dir scheme_flags <<< "$entry"
    # Collect only the selected distances (output names carry the input stem,
    # e.g. "...,r=10,w=...__<scheme>,pz=...,p=....stim").
    files=()
    for r in "${RVALS[@]}"; do
      for f in "$dir"/*,r="$r",*",$POINT.stim"; do
        [ -e "$f" ] && files+=("$f")
      done
    done
    [ "${#files[@]}" -gt 0 ] || { echo "no circuits for scheme $name at ($POINT)" >&2; exit 1; }
    echo "--- decoder: $decoder, scheme: $name (${#files[@]} circuits, max_shots=$max_shots, max_errors=$max_errors)"
    $PY tools/decode.py \
      --circuit-files "${files[@]}" \
      --decoder "$decoder" \
      --max-shots "$max_shots" \
      --max-errors "$max_errors" \
      --output-preamble "$PREAMBLE" \
      --print-progress
  done
done

echo "=== Done: out/collection/$PREAMBLE.csv ==="
