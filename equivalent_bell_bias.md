This file defines two new noise models
`biased_bell_v2_equivalent`
`biased_bell_v2_equivalent_majority`
They are is similar to `biased_bell_v2` but instead of bell based measurements, we are going to do:

```
MXX 0 1
```

turns into
```
R 2
XCX 0 2
XCX 1 2
<NOISE-OPERATIONS>
M 2
```
where 2 is an ancilla qubits (and can have other indexes, like, the first free index in the circuit so they won't collide).

`<NOISE-OPERATIONS>` are going to be defined later.

For MZZ, the bell state is the same (|00>+|11> = |++>+|-->), but the CNOTs are flipped and measurements are in the Z basis:
```
MZZ 0 1
```
turns into
```
R 2
H 0 1
XCX 0 2
XCX 1 2
<NOISE-OPERATIONS>
H 0 1
M 2
```

For MYY, S gates on the data qubits conjugate Y to X, allowing the Bell pair to measure in the Y basis:
```
MYY 0 1
```
turns into (qubits 0,2 are data; qubits 1,3 are ancilla):
```
R 2
S 0 1
YCX 0 2
YCX 1 2
<NOISE-OPERATIONS>
S_DAG 0 1
M 2
```


## `biased_bell_v2_equivalent` noise:

This is a noise which is equivalent to the `biased_bell_v2` noise. 

### Noise parameterization

The model is controlled by two parameters: **`p`** (physical error rate) and **`p_bell_z`** (Z-component of Bell-pair prep noise, set independently). The individual noise channels are derived as:

| Channel | Value | Description |
|---------|-------|-------------|
| `p_m` | `p` | Measurement flip probability |
| `p_1` | `p / 10` | 1-qubit depolarizing after single-qubit gates |
| `p_init` | `p` | Data-qubit state-prep noise (X_ERROR after R, Z_ERROR after RX) |

and, replace `<NOISE-OPERATIONS>` with:

```
E(p) X2
E(p_bell_z) X2
E(p/2) X0 X2
E(p/2) X0
E(p) Z0 X2
E(p) Z1 X2
```
where 2 is the ancilla index, 0 and 1 are the ....

The ancilla resets inside the Bell-state preparation are *not* affected by `p_init`; their noise is captured by `p_bell_{x,y,z}`. `p_internal` is retained as a specified parameter (defaulting to `0` in v2) so the v1 channel can still be turned on if desired.

## Allowed input operations

The Bell translation has a strict allowlist — any instruction in the input circuit that doesn't fall into one of these categories is rejected with a clear error:

- State prep: `R`, `RX`
- Single-qubit Clifford gates (H, S, S_DAG, etc.)
- Single-qubit Z- or X-basis measurements (`M`, `MR`, `MX`, `MRX`)
- Pairwise measurements: 2-body `MPP`, `MXX`, `MZZ`, `MYY`
- Annotations / control: `TICK`, `DETECTOR`, `OBSERVABLE_INCLUDE`, `QUBIT_COORDS`, `SHIFT_COORDS`, `REPEAT` blocks

In particular, **two-qubit gates between data qubits (CX, CZ, etc.) are not allowed** — under the distributed model, all data-qubit interactions must go through the Bell-mediated measurements.

