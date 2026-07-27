# distributed-ghz

Tools for comparing different ways to perform syndrome extraction on a distributed
quantum computer using stim.

See [context.md](./context.md) for the full project scope and design notes.

## Installation

Requires Python 3.10+.

```bash
git clone <repo-url> distributed-ghz
cd distributed-ghz
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Optional extras:

```bash
pip install -e ".[sampling]"    # sinter (for sampling/decoding)
pip install -e ".[notebooks]"   # jupyter, matplotlib, pandas
pip install -e ".[all]"         # everything
```

## Project layout

```
distributed-ghz/
├── src/                     # library code
│   ├── bell_translation.py          # Bell-pair measurement translation + noise
│   ├── equivalent_noise_calculator.py  # single-ancilla equivalent noise model
│   ├── threshold_fit.py             # quadratic finite-size-scaling threshold fit
│   ├── bell_threshold_analysis.py   # load/fit/plot backend for the Bell threshold scans
│   ├── high_rate_analysis.py        # high-rate (hyperbolic Floquet) analysis
│   └── noise_utils.py               # shared noise helpers
├── tools/                   # standalone scripts
│   ├── translate.py                 # CLI for bell_translation
│   ├── make_threshold_circuits_bell.py  # generate Bell threshold-scan circuits
│   ├── decode.py                    # sinter sampling/decoding into CSVs
│   ├── _grid_scans/                 # (p_res, p_z) grid pipeline (Fig. 4)
│   └── stairway_codes/              # stairway LER-vs-p_z pipeline (Fig. 5)
├── notebooks/               # analysis and figure notebooks
├── tests/                   # pytest suite
├── in/                      # input circuits + metadata
├── out/                     # output circuits, CSVs, and figures
└── pyproject.toml
```

## Tools

### `tools/translate.py` — Bell-pair translation

Replaces every 2-body `MPP` / `MXX` / `MZZ` / `MYY` measurement in a stim
circuit with the equivalent Bell-mediated circuit, optionally adding the
`biased_bell` noise model. Existing noise on the input circuit is stripped
first.

Each pairwise measurement becomes 2 ancilla measurements; `DETECTOR` and
`OBSERVABLE_INCLUDE` record references are rewritten automatically so the
detectors and observables of the translated circuit reproduce those of the
original.

#### biased_bell noise parameters

| Flag             | Channel                                | Where applied                                                  |
|------------------|----------------------------------------|----------------------------------------------------------------|
| `--p-init-x/y/z` | `PAULI_CHANNEL_1(px, py, pz)`          | One ancilla, right after the Bell-pair CX                      |
| `--p-internal`   | `DEPOLARIZE2(p)`                       | After every CX between an ancilla and a data qubit             |
| `--p-m`          | Measurement-flip probability           | Argument of `M(p)` / `MX(p)` (all measurements)                |
| `--p-1`          | `DEPOLARIZE1(p)`                       | After every single-qubit gate (S, S_DAG, H, ...)               |
| `--p-init`       | `X_ERROR(p)` / `Z_ERROR(p)`            | After every data-qubit `R` / `RX` (state preparation)          |

**Allowed input operations.** The translator has a strict allowlist: it
accepts `R` / `RX`, single-qubit Cliffords, single-qubit Z- or X-basis
measurements (`M`, `MR`, `MX`, `MRX`), pairwise measurements (`MPP` /
`MXX` / `MZZ` / `MYY`), and annotations (`TICK`, `DETECTOR`,
`OBSERVABLE_INCLUDE`, `QUBIT_COORDS`, `SHIFT_COORDS`, `REPEAT`). Two-qubit
gates between data qubits (CX, CZ, ...) are rejected with an error — in
the distributed model, all data-qubit interactions must go through the
Bell-mediated measurements.

`--input` accepts a single `.stim` file, a directory (every `.stim` file
inside is processed, non-recursively), or a glob pattern. **Quote glob
patterns** so the shell doesn't expand them — argparse only takes a single
value, and Python expands the glob itself.

#### Examples

Noisy translation, single file:

```bash
python tools/translate.py \
  --input "in/circuits/hyperbolic_floquet/circuit.stim" \
  --output-dir out/ \
  --p-init-x 0.001 --p-init-y 0.001 --p-init-z 0.001 \
  --p-internal 0.001 --p-m 0.001 --p-1 0.0001
```

Translate every circuit in a directory (or matching a glob):

```bash
python tools/translate.py \
  --input "in/circuits/test/*.stim" \
  --output-dir out/circuits/ \
  --noiseless
```

Noiseless run (useful for verifying that the translation preserves logical
structure — the resulting circuit should have a deterministic DEM):

```bash
python tools/translate.py \
  --input "in/circuits/hyperbolic_floquet/circuit.stim" \
  --output-dir out/ \
  --output-suffix biased_bell_noiseless \
  --noiseless
```

Specifying a metadata file explicitly (defaults to `<input-stem>.json`,
only valid for a single-file `--input`):

```bash
python tools/translate.py \
  --input "in/circuits/some.stim" \
  --metadata "in/circuits/other_metadata.json" \
  --output-dir out/ \
  --p-internal 0.001 --p-m 0.001
```

#### Outputs

For an input `<stem>.stim` and `--output-suffix biased_bell`, the tool writes:

- `out/<stem>_biased_bell.stim` — the translated circuit
- `out/<stem>_biased_bell.json` — input metadata merged with translation
  metadata (noise model name, noise params, paths, qubit / measurement /
  detector / observable counts)

## Library usage

```python
import stim
from bell_translation import translate_circuit

src = stim.Circuit.from_file("in/circuits/hyperbolic_floquet/circuit.stim")

# Noiseless
translated = translate_circuit(src)

# With biased_bell noise
noisy = translate_circuit(src, {
    "p_init_x": 0.001, "p_init_y": 0.001, "p_init_z": 0.001,
    "p_internal": 0.001, "p_m": 0.001, "p_1": 0.0001,
})
```

## Reproducing the figures

**Figure 4** — honeycomb and hyperbolic Floquet threshold grids:

```bash
bash tools/run_grid_all.sh
```

**Figure 5** — stairway-code logical error rate vs. `p_z`:

```bash
bash tools/stairway_codes/make_circuits_stairway.sh
bash tools/stairway_codes/decode_stairway.sh
bash tools/stairway_codes/plot_stairway.sh
```
