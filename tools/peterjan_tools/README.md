# Gauge floquet phenomenological threshold pipeline

The grid-scan → threshold → boundary pipeline (the analog of
`out/figures/full_vs_local.pdf`) applied to the **native repeated-measurement**
circuits of Derks et al., *"Improving the spacetime overhead of dynamic
stabilizer codes by repeating measurements"* (vendored in
`floquet_colour_codes_numerics/`).

Unlike the `src/repeat_measurements` full-decoder route — which produces
non-graphlike circuits that MWPM decodes lossily (see
`docs/repeated_measurements_summary.md`) — here the **gauge factor `g` natively
repeats every check `g` times with comparison detectors built in**, so the
circuits stay graphlike and plain pymatching thresholds are meaningful. Larger
`g` suppresses measurement noise at the cost of more rounds.

## What varies

- **Code**: `GaugeHoneycombCode` (periodic honeycomb — *not* the planar honeycomb
  used elsewhere in this repo, so thresholds are not comparable to the planar
  numbers).
- **Axes**: `p` = data noise `PAULI_CHANNEL_1(p/3, p/3, p/3)` (x, the `p_res`
  analog); `p_m` = measurement-flip probability (y, the `p_z` analog — stored in
  metadata under `p_bell_z` so the existing fit stack works unchanged).
- **Gauge factor** `g ∈ {1, 2, 3, 4, 5}` — one boundary curve each (`g=1` = baseline).
- **Distances** 4, 8, 12; single memory observable (`memory_x` by default).

### Gauge factors ≥ 4 need the extended table (one-time)

`GaugeHoneycombCode.get_number_of_rounds_for_timelike_distance` reads a
precomputed timelike-distance table shipped with Kandel
(`.../timelike_distance_data/hcc_graphlike_td_data_phenomenological.json`) that
only covers gauge factors 1–3 (matching the paper's own scan); `g ≥ 4` raises
`KeyError`. Run once:

```bash
python tools/peterjan_tools/extend_timelike_distance_data.py --max-g 5   # ~5 min
```

This computes the `(4,4,4)`/`(5,5,5)` entries with Kandel's own logic and writes
a superset to `data/hcc_graphlike_td_data_phenomenological_extended.json`;
`make_grid_circuits_peterjan.py` loads it and monkeypatches `get_timelike_distance`
(the installed Kandel package is never mutated, `g ≤ 3` is unchanged). Requesting
`g ≥ 4` before this runs fails with a pointer to the extender.

## Setup

`main` is the OzYossarian/Kandel package (the paper's stim compiler); it is not
on PyPI and its pins are ancient. Install the optional group **with `--no-deps`**
against the current environment:

```bash
pip install -e ".[peterjan]" --no-deps          # kandel + stimcirq + networkx + alive-progress
# or, if the git dep is fussy: clone OzYossarian/Kandel and
#   pip install -e /path/to/Kandel --no-deps && pip install stimcirq networkx alive-progress
```

Smoke test: `python -c "from main.codes.tic_tac_toe.gauge.GaugeHoneycombCode import GaugeHoneycombCode"`.

## Run order

```bash
# 0. one-time, only for g >= 4:
python tools/peterjan_tools/extend_timelike_distance_data.py --max-g 5

# 1. Cheap 1-D check scans: locate the axis intercepts, calibrate grid ceilings.
bash tools/peterjan_tools/run_check_peterjan.sh        # -> peterjan_check_{pz,pres}.csv
                                                        #    peterjan_check_pz.pdf + table

# 2. Set the per-g PM_HI ceilings in run_grid_peterjan.sh from the fitted p_m*,
#    then run the full (p, p_m) grids and decode.
bash tools/peterjan_tools/run_grid_peterjan.sh         # -> peterjan_grid.csv

# 3. Boundary figure p_m* vs p*, one curve per g.
python tools/peterjan_tools/plot_peterjan_boundary.py  # -> out/figures/peterjan_boundary.pdf
```

Budgets are env-tunable: `MAX_SHOTS` (default 100000), `MAX_ERRORS` (default 500),
`OBS` (memory_x / memory_z), `PYTHON` (default `.venv/bin/python`). All decode
steps use `tools/decode.py --resume`, so re-running only tops up the budget.

## Files

| file | role |
|---|---|
| `extend_timelike_distance_data.py` | one-time: add g=4,5 entries to Kandel's timelike-distance table → `data/` |
| `make_grid_circuits_peterjan.py` | compile GaugeHoneycombCode per (d, g, p, p_m) via Kandel; write .stim + sidecar .json |
| `run_check_peterjan.sh` | pure-p_m and pure-p 1-D scans + `analyze_peterjan_check.py` |
| `run_grid_peterjan.sh` | full (p, p_m) grid per g + decode |
| `analyze_peterjan_check.py` | per-g LER-vs-p_m panels + fitted threshold table |
| `plot_peterjan_boundary.py` | boundary p_m* vs p* per g (reuses `tools/_grid_scans` fitters) |

Outputs: circuits under `out/circuits/peterjan_circuits/`, CSVs
`out/collection/peterjan_*.csv`, figures `out/figures/peterjan_*.pdf`.
