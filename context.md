The goal of this project is to compare different ways to perform SE on a distributed quantum computer.

Mostly using stim.

## First approach - Bell measurements:

This fits floquet codes where the operations are only pairwise measurements like MZZ or MXX (or mpp X0*X1).

In this case, you are using a bell state like this:

```
MXX 0 1
```

turns into
```
RX 2
R 3
CX 2 3
CX 2 0
CX 3 1
MX 2 3
```
where 2 and 3 are ancilla qubits (and they can have other indexes, like, the first free index in the circuit so they won't collide).

For MZZ in `biased_bell_v3`, the data qubits are H-conjugated so the same ancilla→data CX pattern as MXX can be reused:
```
MZZ 0 1
```
turns into
```
RX 2
R 3
CX 2 3
H 1 0
CX 2 0
CX 3 1
H 1 0
MX 2 3
```
(H gates on data conjugate Z↔X, so XX-style readout on the Bell pair measures ZZ on the original data. See Appendix B for the v2 form, which used reversed-direction CX and a Z-basis measurement.)

For MYY, S gates on the data qubits conjugate Y to X, allowing the Bell pair to measure in the Y basis:
```
MYY 0 1
```
turns into
```
R 3
RX 2
CX 2 3
S 0 1
CX 3 1 2 0
MX 3 2
S_DAG 0 1
```

## Second approach - shor-style measurements (GHZ-based)

This is for multi-qubit stabilizer codes, like Bivariate Bicycle codes.

You generate a GHZ state (somehow, to be determined), and use it to do shor-style measurement on the stabilizer. 

TBD - will be implemented after the Bell approach is working.


# Implementation details

The script gets a stim circuit from the in/ folder. The input circuit may be noisy — we always strip noise first using `circuit.without_noise()`.

The circuits should have either MXX, MZZ, MYY, or MPP instructions. MXX/MZZ/MYY are converted to MPP form before translation. 

There should be two translation methods, one for each approach. Translation method also adds noise.

## Bell measurements method:

For a pairwise measurement (MXX, MZZ or MPP with targets at the length of two) the script should replace the instruction with an equivalent circuit that uses Bell measurements (as described above).

## GHZ-based measurement:

Here we will use MPP with more than two-body operators. The order of the operators might be important, since the ghz creation method might not be symmetric for all qubits.

TBD.

## Circuit metadata

Each circuit (input and output) should include a json file with the same name that represents the metadata of the circuit.

**Input metadata:** Can contain anything relevant to the code (typically `code_name`, `distance`, `rounds`), but no strict schema required.

**Output metadata:** Should include everything from the input metadata plus translation-relevant fields: noise model name, noise parameters, path to the output stim circuit, and any other relevant translation details.

Look [here](https://github.com/shohamjac/stairway-codes/tree/main) for example in another project I did. 

## Noise

The noise model is called `biased_bell_v3`.

The general approach is that each data qubit `q` sits on a "site" which holds exactly one dedicated ancilla at qubit index `num_data + q`. A pairwise measurement on data qubits `(q1, q2)` always uses ancillas `(num_data + q1, num_data + q2)`. The same ancilla is reused across all measurement rounds (reset at the start of each Bell block). The Bell state is generated between the two ancilla sites and used to jointly-measure the data qubits.

The only difference between `biased_bell_v3` and the v2 model (Appendix B) is the **MZZ gadget**: v3 sandwiches the ancilla↔data CX block in `H` gates on the data qubits and reads out in the X basis (so the Bell measurement form is identical to the MXX gadget). v2 instead reversed the CX direction and read out in the Z basis. The noise channels are unchanged.

### Noise parameterization

The model is controlled by two parameters: **`p`** (physical error rate) and **`p_bell_z`** (Z-component of Bell-pair prep noise, set independently). The individual noise channels are derived as:

| Channel | Value | Description |
|---------|-------|-------------|
| `p_bell_x` | `p / 2` | X-component of Bell-pair prep noise |
| `p_bell_y` | `p / 2` | Y-component of Bell-pair prep noise |
| `p_bell_z` | `p_bell_z` | Z-component of Bell-pair prep noise (set independently) |
| `p_data_z` | `p` | Z_ERROR on the two participating data qubits during Bell-pair prep |
| `p_internal` | `0` | 2-qubit depolarizing after ancilla↔data CX (specified parameter, default 0) |
| `p_m` | `p` | Measurement flip probability |
| `p_1` | `p / 10` | 1-qubit depolarizing after single-qubit gates (including the MZZ-gadget H gates) |
| `p_init` | `p` | Data-qubit state-prep noise (X_ERROR after R, Z_ERROR after RX) |

After Bell-pair creation (and before the ancilla↔data CX gates), apply:

- `PAULI_CHANNEL_1(p_bell_x, p_bell_y, p_bell_z)` on one ancilla qubit, and
- `Z_ERROR(p)` on the two data qubits participating in the measurement.

The ancilla resets inside the Bell-state preparation are *not* affected by `p_init`; their noise is captured by `p_bell_{x,y,z}`. `p_internal` is retained as a specified parameter (defaulting to `0`) so the v1 channel can still be turned on if desired. The H gates introduced for the MZZ gadget pick up `DEPOLARIZE1(p_1)` like any other single-qubit gate.

### GHZ noise model

For the GHZ-based (Shor-style) measurement approach, the noise model is called `biased_ghz`. After the GHZ state on the `k` ancilla sites is prepared (one ancilla per data qubit participating in the stabilizer), a per-qubit Pauli channel is applied to each ancilla in the GHZ. The other channels match the Bell model.

| Channel | Value | Description |
|---------|-------|-------------|
| `p_ghz_x` | `p / 2` | X error on each qubit of the GHZ |
| `p_ghz_y` | `p / 2` | Y error on each qubit of the GHZ |
| `p_ghz_z` | `p_bell_z` | Z error on each qubit of the GHZ (set independently), almost equivalent to `number_of_qubits*p_bell_z` on one |
| `p_internal` | `p` | 2-qubit depolarizing after ancilla↔data CX |
| `p_m` | `p` | Measurement flip probability |
| `p_1` | `p / 10` | 1-qubit depolarizing after single-qubit gates |
| `p_init` | `p` | Data-qubit state-prep noise (X_ERROR after R, Z_ERROR after RX) |

Applied as `PAULI_CHANNEL_1(p_ghz_x, p_ghz_y, p_ghz_z)` on every ancilla qubit of the GHZ after preparation. As with the Bell model, ancilla resets inside the GHZ preparation are *not* affected by `p_init` — their noise is captured by `p_ghz_{x,y,z}`.

## Allowed input operations

The Bell translation has a strict allowlist — any instruction in the input circuit that doesn't fall into one of these categories is rejected with a clear error:

- State prep: `R`, `RX`
- Single-qubit Clifford gates (H, S, S_DAG, etc.)
- Single-qubit Z- or X-basis measurements (`M`, `MR`, `MX`, `MRX`)
- Pairwise measurements: 2-body `MPP`, `MXX`, `MZZ`, `MYY`
- Annotations / control: `TICK`, `DETECTOR`, `OBSERVABLE_INCLUDE`, `QUBIT_COORDS`, `SHIFT_COORDS`, `REPEAT` blocks

In particular, **two-qubit gates between data qubits (CX, CZ, etc.) are not allowed** — under the distributed model, all data-qubit interactions must go through the Bell-mediated measurements.

## Pipeline

The overall pipeline:

input circuit + metadata → translate method (adds noise, updates metadata) → output circuit + metadata → sampling script (using sinter) → data CSV (with metadata) → plots using notebooks.

### Milestones

**Milestone 1:** Bell translation with `biased_bell` noise (no sampling yet). ✓
- Translate MXX/MZZ/MPP(2-body) → Bell measurement circuits
- Apply `biased_bell` noise model
- Output translated circuit + metadata

**Milestone 2:** Sampling with sinter. ✓
- `tools/make_circuits.py` — generates translated circuits for a sweep of `p` values (with configurable `p_bell_z`)
- `tools/decode.py` — runs sinter.collect() with pymatching, outputs CSV with metadata
- End-to-end: input circuit → make_circuits (translate + noise) → decode (sinter) → CSV

**Milestone 3:** GHZ-based measurements (TBD).


# Threshold plot

We will create threshold plots for the `biased_bell_v3` noise model, following the same workflow used for `flexible_pairwise` (see [threshold_compare.md](./threshold_compare.md) and `tools/make_threshold_circuits.py`).

## Parameter space

Unlike `flexible_pairwise`, which has three independent noise channels and therefore lives on a 2-simplex of directions, `biased_bell_v3` has only **two** parameters:

- `p` — the overall physical error rate that scales most channels
- `p_bell_z` — the Z-component of the Bell-pair prep noise, set independently of `p`

So the parameter space is the 2D plane `(p, p_bell_z)`. To match the threshold-scan idiom, we pick a small set of "directions" through this plane and sweep `p` along each direction:

- We parameterize each direction by a **bias** `β = p_bell_z / p`. For a chosen `β`, we set `p_bell_z = β · p` and sweep `p`.
- A useful set of biases covers `β ∈ {0, 1, 10, 100, ...}`, including the unbiased case (`β = 0`, plain Z-prep noise off) and several Z-biased regimes that probe the noise asymmetry the model is designed to study.
- Optionally we also include a few **fixed `p_bell_z`** scans (e.g., `p_bell_z = 0.001` or `0.01`) where `p_bell_z` is held constant while `p` varies — these are easier to interpret experimentally because `p_bell_z` is often set by hardware independently of the rest of the circuit.

The other channels follow the v2 derivation in the noise section above (`p_bell_x = p_bell_y = p/2`, `p_data_z = p`, `p_internal = 0`, `p_m = p`, `p_1 = p/10`, `p_init = p`).

## Pipeline

```
input circuit (.stim, noiseless, MPP/MXX/MZZ/MYY)
    → translate_circuit() with biased_bell_v3 noise at (p, p_bell_z)
    → write .stim + .json (sidecar metadata) per (circuit, bias, p)
    → sinter.collect() with pymatching → CSV with metadata
    → plots: logical error rate vs p, grouped by (distance, bias)
              and threshold(p_bell_z) curves
```

This reuses `tools/make_circuits.py` (already takes `--p-values` and `--p-bell-z`) under the hood. The threshold-scan driver — `tools/make_threshold_circuits_bell.py` (analogous to `tools/make_threshold_circuits.py`) — iterates over a list of biases (or fixed `p_bell_z` values) and a `p` sweep, calling the Bell translator for every `(input circuit, bias, p)` triple. Filenames embed `bias` (or `pbz`) and `p` so downstream decoding can group runs.

## What the plots show

Two complementary views are produced from the same data:

1. **Per-bias threshold curves** — for each bias `β`, plot logical error rate vs `p` for several distances; the crossing point is the threshold `p*(β)`. This is the same shape as `notebooks/threshold_comparison.ipynb`, just with `β` (or `p_bell_z`) playing the role of the simplex direction.
2. **Threshold vs bias** — plot `p*(β)` (or `p*(p_bell_z)` at fixed `p_bell_z`) on a single axis. This is a 2D version of the 3D simplex surface plot in `notebooks/threshold_3d_surface.ipynb` and is the headline figure: it shows how robust the code is to Z-biased Bell-pair prep noise.

See `notebooks/threshold_comparison.ipynb` (per-bias curves) and a new threshold-vs-bias notebook for the full workflow.


# Appendix

## Appendix A: `biased_bell` (v1) noise model

The original Bell-mediated noise model, superseded by `biased_bell_v3` above. Kept here for reference and reproducibility of earlier results.

The general approach is that each data qubit `q` sits on a "site" which holds exactly one dedicated ancilla at qubit index `num_data + q`. A pairwise measurement on data qubits `(q1, q2)` always uses ancillas `(num_data + q1, num_data + q2)`. The same ancilla is reused across all measurement rounds (reset at the start of each Bell block). The Bell state is generated between the two ancilla sites and used to jointly-measure the data qubits.

### Noise parameterization

The model is controlled by two parameters: **`p`** (physical error rate) and **`p_bell_z`** (Z-component of Bell-pair prep noise, set independently). The individual noise channels are derived as:

| Channel | Value | Description |
|---------|-------|-------------|
| `p_bell_x` | `p / 2` | X-component of Bell-pair prep noise |
| `p_bell_y` | `p / 2` | Y-component of Bell-pair prep noise |
| `p_bell_z` | `p_bell_z` | Z-component of Bell-pair prep noise (set independently) |
| `p_internal` | `p` | 2-qubit depolarizing after ancilla↔data CX |
| `p_m` | `p` | Measurement flip probability |
| `p_1` | `p / 10` | 1-qubit depolarizing after single-qubit gates |
| `p_init` | `p` | Data-qubit state-prep noise (X_ERROR after R, Z_ERROR after RX) |

Applied as `PAULI_CHANNEL_1(p_bell_x, p_bell_y, p_bell_z)` on one ancilla qubit after Bell pair creation. The ancilla resets inside the Bell-state preparation are *not* affected by `p_init`; their noise is captured by `p_bell_{x,y,z}`.

## Appendix B: `biased_bell_v2` noise model

Predecessor of `biased_bell_v3`. Kept for reference and reproducibility of earlier scans. Only the **MZZ gadget** differs from v3 — the noise channels are identical.

### MZZ gadget (v2)

```
MZZ 0 1
```
turns into
```
RX 2
R 3
CX 2 3
CX 0 2
CX 1 3
M 2 3
```
(CNOT direction is reversed: data→ancilla, and measurement is in Z basis. v3 instead H-conjugates the data qubits and reuses the MXX-style ancilla→data CX pattern with an X-basis readout.)

The MXX and MYY gadgets are the same as in v3.

### Noise parameterization (v2)

Same as v3:

| Channel | Value | Description |
|---------|-------|-------------|
| `p_bell_x` | `p / 2` | X-component of Bell-pair prep noise |
| `p_bell_y` | `p / 2` | Y-component of Bell-pair prep noise |
| `p_bell_z` | `p_bell_z` | Z-component of Bell-pair prep noise (set independently) |
| `p_data_z` | `p` | Z_ERROR on the two participating data qubits during Bell-pair prep |
| `p_internal` | `0` | 2-qubit depolarizing after ancilla↔data CX (specified parameter, default 0) |
| `p_m` | `p` | Measurement flip probability |
| `p_1` | `p / 10` | 1-qubit depolarizing after single-qubit gates |
| `p_init` | `p` | Data-qubit state-prep noise (X_ERROR after R, Z_ERROR after RX) |

After Bell-pair creation (and before the ancilla↔data CX gates), apply `PAULI_CHANNEL_1(p_bell_x, p_bell_y, p_bell_z)` on one ancilla qubit and `Z_ERROR(p)` on the two participating data qubits. Ancilla resets inside the Bell-state preparation are *not* affected by `p_init`.
