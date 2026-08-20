# Repeated measurements: full decoder vs majority vote — summary

*(July 2026. Code: `src/repeat_measurements/`, `tools/_grid_scans/`. Spec:
`src/repeat_measurements/add_measurements_to_circuit.md`.)*

## Goal

Reduce measurement errors in Floquet circuits by repeating each pairwise
measurement r times, and compare two ways of decoding the repetitions:

- **Local (majority vote)**: the existing `equivalent_v3_majority`/`majority5`
  schemes — the r-repetition vote is folded into the noise block, the outer
  decoder sees one effective measurement. Only odd r makes sense.
- **Full decoder**: repeat the measurements structurally in the circuit
  (r copies per mid-circuit layer + comparison detectors between consecutive
  copies) and decode everything together. Works for any r; should in principle
  beat the majority vote since the decoder sees strictly more information.

## What was built

1. **`repeat_measurements(circuit, r, split_detectors=True)`**
   (`src/repeat_measurements/repeat.py`): repeats each mid-circuit measurement
   layer r times, inserts comparison detectors, rewrites all
   DETECTOR/OBSERVABLE_INCLUDE rec offsets. To preserve matchability, each
   original annotation's rec targets are classified per target: *past-facing*
   (annotation in the same time step as the measurement) → first copy;
   *future-facing* (later time step) → last copy. Multiple future-facing uses
   are allowed (normal in Floquet codes); multiple past-facing uses raise
   (escape hatch: `split_detectors=False`). Respects
   `TICK[start_noise]`/`TICK[stop_noise]` windows; final readout layer is
   never repeated. 28 tests in `tests/test_repeat_measurements.py`.

2. **Pipeline integration**: `make_grid_circuits_bell.py --repeat r` applies
   the transform before `translate_circuit` and names the scheme
   `biased_bell_v3_equivalent_fullrep{r}` (metadata: `repetitions`,
   `decoder_kind`). Guard: refuses `--repeat` combined with majority schemes
   (would double-count repetition).

3. **Run scripts** (all resume-able, honeycomb d=3,5,7, p_res pinned to 1e-6
   for the p_z scans):
   - `run_grid_full.sh` — full (p_res, p_z) grids per r → `grid_full.csv`.
   - `run_pres_threshold_check_full.sh` — cheap pure-p_z scan per r
     → `full_check_pz.csv` (pymatching).
   - `run_pz_check_relay_bp.sh` — same scan, relay-bp → `full_check_pz_relay_bp.csv`.
   - `run_pz_check_tesseract_lowpz.sh` — same scan capped at p_z < 1e-2,
     tesseract → `full_check_pz_tesseract.csv`.
   - `run_pz_check_tesseract.sh` — uncapped tesseract scan (abandoned, see below).
   - Analyzers/figures: `analyze_full_check.py` (per-r LER panels + threshold
     table), `plot_full_vs_local.py` (boundary figure p_z* vs p_res*,
     majority curves from `grid_bell.csv` vs fullrep curves from `grid_full.csv`).

4. **Decoders wired into `tools/decode.py`**: `relay-bp`, `relay-bp-light`
   (10-set ensemble), `relay-bp-ultralight` (4-set), `mem-bp`, `msl-bp`
   (via the new `bp` optional-dependency group: relay-bp, beliefmatching,
   stimbposd), alongside the existing pymatching / tesseract_short / tesseract_long.

## Results

### 1. Pymatching (MWPM): full decoder gain saturates at ~1.8x

Pure-p_z thresholds (p_res ≈ 0, honeycomb, FSS fit):

| scheme                | p_z*    | vs r=1 |
|-----------------------|---------|--------|
| r=1 baseline          | 1.37 %  | 1.00x  |
| full decoder, r=2     | 2.16 %  | 1.57x  |
| full decoder, r=3     | 2.35 %  | 1.71x  |
| full decoder, r=4     | 2.38 %  | 1.73x  |
| full decoder, r=5     | 2.41 %  | 1.76x  |
| 3-rep majority vote   | 6.97 %  | 5.08x  |
| 5-rep majority vote   | 11.2 %  | 8.13x  |

The majority thresholds match the combinatorial prediction exactly
(3p*² → 6.8 %, 10p*³ → 11.1 %). The full decoder under MWPM captures only a
fraction of the repetition benefit and saturates by r=3.

**Diagnosis**: pure-p_z logical errors in these circuits are non-graphlike —
`shortest_graphlike_error()` finds nothing at any r. Honeycomb measurements
are referenced by up to 4 detectors (150 measurements have 3 future-facing +
1 past-facing use), so single measurement flips fire >2 detectors and MWPM
decodes a lossy decomposition. The DEM carries the information; matching
can't extract it.

### 2. Tesseract at threshold density: intractable

`tesseract_short` on the uncapped p_z scan: ~4 shots completed in 13 minutes
(minutes per shot). Not a CSS issue — MLE beam search degrades with detection
event *density*, and a threshold scan sits exactly in the dense regime.
Killed.

### 3. Relay-BP: fast, huge win in spot checks, but a convergence-failure floor

Spot check (d=3, r=3, 500 shots): at p_z=2.7 % relay-bp 0 errors vs
pymatching 36; at 6.4 % relay-bp 0 vs pymatching 179. Default ensemble too
slow at d=7 (~2.5 s/shot) → added `relay-bp-light` (10 sets, stop_nconv=1,
~5x faster, no accuracy loss in spot checks at d=3/d=7-few-shots).

But the 1k-shot scan showed LERs *worse than pymatching* at d=5/7, growing
with distance — BP convergence failures, not physics. Three-decoder spot
check at p_z=1e-3, r=1 confirmed:

| decoder          | d=3        | d=5     | d=7     |
|------------------|------------|---------|---------|
| pymatching (1M)  | 3.0e-4     | 1e-6    | 0       |
| relay-bp-light (1k) | 0       | 1e-3    | 2e-3    |
| tesseract (1k)   | 0          | 0       | 0       |

relay-bp-light has an error floor ≥ ~1e-3 at d≥5 that inverts the distance
ordering. Its fitted thresholds would be meaningless. (Untested: whether the
full 60-set default or heavier ensembles remove the floor.)

### 4. Tesseract below p_z = 1e-2: essentially perfect, and far above MWPM

10k shots per point, p_z < 1e-2 (~53 min): **zero errors in every cell**
(r=1 and r=3, all distances) except 3 isolated d=3 events around
p_z ≈ 0.6–0.8 %. At p_z=0.55 %, d=3, pymatching sits at 1.3e-2 — tesseract
at < 1e-4. So:

- The MWPM "r=1 baseline = 1.37 %" is itself a decoder artifact; the true
  r=1 threshold with a good decoder is **> 1 %** and likely much higher.
- Every threshold in table 1 is decoder-limited, not code-limited.
- The r=1 vs r=3 comparison is still open: neither curve lifts off below 1e-2.

## Where this stands

The full-decoder *transform* is correct and verified (deterministic
detectors, decomposable DEM, r=1 identity, exact combinatorial agreement for
the majority reference). The open problem is purely **decoder choice in the
crossing region** (p_z ~ 2–10 %+, dense syndromes):

- MWPM: fast but lossy (~1.8x cap).
- Tesseract: accurate but intractable above ~1e-2.
- relay-bp-light: fast but broken at d≥5.

## Next steps

1. Benchmark stronger relay-bp ensembles (default 60-set, heavier) at the
   sub-1e-2 points against tesseract's zeros; if some config removes the d=7
   floor at acceptable ms/shot, rerun the r=1 + r=3 p_z scans with it.
2. Otherwise: push tesseract point-by-point toward ~2e-2 (d=3/5 are cheap;
   d=7 dominates cost) to bracket at least the r=1 crossing.
3. The full (p_res, p_z) boundary comparison (`plot_full_vs_local.py`,
   `run_grid_full.sh`) is ready to re-run once a trustworthy fast decoder is
   chosen; r=2 grid data (pymatching) is complete, r=3 partial.
4. Intended production run on soca-purple-job:
   `pip install -e ".[sampling,bp]"` then the chosen check script
   (all scripts resume into their CSVs; budgets via MAX_SHOTS/MAX_ERRORS env vars).

## Artifact map

| artifact | content |
|---|---|
| `out/collection/full_check_pz.csv` | pymatching p_z scans, r=2..5 (complete) |
| `out/collection/full_check_pz_relay_bp.csv` | relay-bp-light r=1,3 (1k sanity; floor-contaminated at d≥5) |
| `out/collection/full_check_pz_tesseract.csv` | tesseract r=1,3, p_z<1e-2, 10k shots (clean) |
| `out/collection/grid_full.csv` | full-decoder (p_res,p_z) grid, r=2 complete (pymatching) |
| `out/collection/spot_check_pz1e3.csv` | 3-decoder spot check at p_z=1e-3 |
| `out/figures/full_check_pz.pdf` | pymatching per-r LER panels + fitted thresholds |
| `out/figures/full_check_pz_tesseract.pdf` | tesseract panels (flat ~0 below 1e-2) |
| `out/figures/full_vs_local.pdf` | boundary figure: majority vs fullrep (pymatching-fitted) |
