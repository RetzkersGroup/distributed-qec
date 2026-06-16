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

For MZZ, the bell state is the same (|00>+|11> = |++>+|-->), but the CNOTs are flipped and measurements are in the Z basis:
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
(CNOT direction is reversed: data→ancilla, and measurement is in Z basis.)

For MYY, S gates on the data qubits conjugate Y to X, allowing the Bell pair to measure in the Y basis:
```
MYY 0 2
```
turns into (qubits 0,2 are data; qubits 1,3 are ancilla):
```
R 3
RX 1
CX 1 3
S 0 2
CX 3 2 1 0
MX 3 1
S_DAG 0 2
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

The noise model is called `biased_bell`.

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

### GHZ noise model

For the GHZ-based (Shor-style) measurement approach, the noise model is called `biased_ghz`. After the GHZ state on the `k` ancilla sites is prepared (one ancilla per data qubit participating in the stabilizer), a per-qubit Pauli channel is applied to each ancilla in the GHZ. The other channels match the Bell model.

| Channel | Value | Description |
|---------|-------|-------------|
| `p_ghz_x` | `p_ghz_x` | X error on each qubit of the GHZ |
| `p_ghz_y` | `p_ghz_y` | Y error on each qubit of the GHZ |
| `p_ghz_z` | `p_ghz_z` | Z error on each qubit of the GHZ |
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
