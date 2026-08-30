# Gauge-honeycomb phenomenological threshold pipeline

*(August 2026. Code: `tools/peterjan_tools/`. Reuses the fit/plot stack in
`tools/_grid_scans/` and the sampler in `tools/decode.py`.)*

## Goal

Reproduce the grid-scan → threshold → boundary workflow that produces
`out/figures/full_vs_local.pdf`, but for the **native repeated-measurement**
circuits of Derks et al., *"Improving the spacetime overhead of dynamic
stabilizer codes by repeating measurements"* (vendored in
`floquet_colour_codes_numerics/`). Output figure:
`out/figures/peterjan_boundary.pdf` — threshold boundary `p_m*` vs `p*`, one
curve per gauge factor.

## Why a separate route from `src/repeat_measurements`

`docs/repeated_measurements_summary.md` documents the full-decoder transform in
`src/repeat_measurements/`: it repeats measurement layers structurally and adds
comparison detectors, but the resulting circuits are **non-graphlike**, so MWPM
decodes a lossy decomposition (~1.8x cap) and the accurate decoders were either
intractable (tesseract at threshold density) or broken at large distance
(relay-bp floor).

The gauge floquet codes sidestep this: the **gauge factor `g` repeats each check
`g` times with the comparison detectors built into the code definition**, and the
compiler is constructed so the circuit stays **graphlike** (the paper's own
generator asserts `shortest_graphlike_error() >= d`). So plain pymatching gives
meaningful thresholds, and larger `g` cleanly trades more rounds for measurement-
noise suppression — the comparison this repo could not make with the
non-graphlike route.

## Circuit family and noise axes

- **Code**: `GaugeHoneycombCode` — the **periodic** honeycomb (no spatial
  boundary). Thresholds are therefore **not** comparable to the *planar*
  honeycomb numbers elsewhere in this repo.
- **Experiment**: single memory observable (`memory_x` by default; `memory_z`
  via `OBS=`), distances 4, 8, 12. Rounds are scaled per `(d, g)` so the
  timelike distance matches `d` (Kandel's phenomenological round-scaling).
- **x-axis `p`** — data-qubit noise `PAULI_CHANNEL_1(p/3, p/3, p/3)` (the
  `p_res` analog).
- **y-axis `p_m`** — ancilla measurement-flip probability (the `p_z` analog —
  what repetition suppresses). Stored in circuit metadata under the key
  `p_bell_z` so `plot_grid.load_grid` and the FSS fitter sweep it with **no
  changes to the shared code**; an honest `p_m` alias is stored alongside.

Noise is applied at **absolute** rates — no `per`-renormalization (unlike the
paper's `generate_circuits.py`, which fixes an average physical rate and
renormalizes `px:py:pz:pm` weights). Independent axes are exactly what the
2-D threshold fit needs.

### The `p ≈ 0` / `p_m ≈ 0` floor

At exactly `p = 0` the data-noise hyperedges vanish and the detector error model
no longer decomposes into graphlike components (these periodic circuits need
every channel present, which the paper's renormalized noise always had). The
1-D check scans therefore pin the "off" axis to a **1e-6 floor**, matching the
`p_res` floor in `tools/_grid_scans/run_pres_threshold_check.sh`. The floor does
not shift the threshold; it just keeps the matching graph well-formed. The full
2-D grid starts its `p` axis at `1e-4` (as `run_grid_full.sh` does), so it never
hits `p = 0` either.

## Gauge factor and the timelike-distance table

Gauge factors **1–5** are supported. Kandel computes the number of rounds for a
target timelike distance from a **precomputed graphlike-timelike-distance table**
shipped inside the package
(`.../timelike_distance_data/hcc_graphlike_td_data_phenomenological.json`), which
only covers uniform gauge factors 1–3 (the paper's scanned range). `g >= 4`
raises `KeyError('(4, 4, 4)')`.

`tools/peterjan_tools/extend_timelike_distance_data.py` computes the missing
`(4,4,4)` / `(5,5,5)` entries using Kandel's **own** generation logic
(`create_timelike_distance_data.get_td_bulk_and_boundary`) and writes a superset
to `tools/peterjan_tools/data/hcc_graphlike_td_data_phenomenological_extended.json`.
`make_grid_circuits_peterjan.py` loads that table at import time and
monkeypatches `GaugeHoneycombCode.get_timelike_distance` to read it — the
installed Kandel package is never mutated, and the `g <= 3` entries are identical
to the shipped ones, so their behaviour is unchanged. Requesting `g >= 4` before
running the extender fails with a pointer to it.

## Dependencies (`peterjan` optional group)

The paper's compiler is the `main` package = **OzYossarian/Kandel** (not on
PyPI; its pins are ancient). Install with `--no-deps` against the current env:

```bash
pip install -e ".[peterjan]" --no-deps      # kandel + stimcirq + networkx + alive-progress
```

Runtime deps beyond `main` itself: `stimcirq` (pulls `cirq-core`), `networkx`,
`alive-progress`. Verified working against the repo's current `stim` (1.16).

## Pipeline / run order

```bash
# 0. one-time (only needed for g >= 4)
python tools/peterjan_tools/extend_timelike_distance_data.py --max-g 5

# 1. cheap 1-D check scans -> axis intercepts, and PM_HI calibration
bash tools/peterjan_tools/run_check_peterjan.sh          # peterjan_check_{pz,pres}.csv
                                                          # + peterjan_check_pz.pdf, table

# 2. full (p, p_m) grid per g -> decode
bash tools/peterjan_tools/run_grid_peterjan.sh           # peterjan_grid.csv

# 3. boundary figure
python tools/peterjan_tools/plot_peterjan_boundary.py    # out/figures/peterjan_boundary.pdf
```

Env knobs: `MAX_SHOTS` (100000), `MAX_ERRORS` (500), `OBS` (memory_x/memory_z),
`PYTHON` (`.venv/bin/python`). All decode steps use `tools/decode.py --resume`,
so re-running only tops up toward the budget.

## Data flow

| stage | script | writes |
|---|---|---|
| generate | `make_grid_circuits_peterjan.py` (Kandel per point, parallel) | `out/circuits/peterjan_circuits/**/*.stim` + sidecar `.json` |
| decode | `tools/decode.py --resume` (pymatching) | `out/collection/peterjan_{grid,check_pz,check_pres}.csv` |
| fit + panels | `analyze_peterjan_check.py` | `out/figures/peterjan_check_pz.pdf` + per-`g` threshold table |
| boundary | `plot_peterjan_boundary.py` | `out/figures/peterjan_boundary.pdf` |

Metadata keys consumed by the shared fitter: `p`, `p_bell_z` (=`p_m`),
`expected_distance`, `k=1`, `rounds`, `noise_model = peterjan_pheno_g{g}`,
`input_circuit` (its filename carries `_n=`/`_d=`/`_k=` for the `[[n,k,d]]`
labels). Threshold fits reuse `plot_grid.fit_pz_threshold` (FSS across
distances) with a per-`g` `p_m` window cap (`_PM_MAX`), since the threshold rises
with `g` and would otherwise pin against the shared 0.09 default.

## Status / caveats

- Pipeline validated end-to-end on a mini-run (generate → decode → analyze →
  boundary). Rough mini-run pure-`p_m` thresholds rose monotonically with `g`
  (≈1.2% → 2.5% → 3.0% for g=1/2/3 at d=4,8 only, low shots), confirming the
  qualitative repetition benefit; production numbers need the full d=4,8,12 grid
  at the default budgets.
- The `_PM_MAX` window caps and per-`g` `PM_HI` scan ceilings in the run scripts
  are starting estimates — after the first check run, set `PM_HI` from the
  fitted `p_m*` so the boundary crossing sits below the grid ceiling.
- Cross-check against the paper's own thresholds via
  `floquet_colour_codes_numerics/phenomenological_noise/out/data/` +
  `heatmap_plot.ipynb` (periodic code, so not comparable to planar honeycomb).
