Here, we present the general concepts for creating a repeated measurement floquet circuit.

# Goal

In order to reduce measurement errors, one can repeat the pairwise measurements composing a Floquet circuit.

The results can be decoded using majority vote (inner-decoder) and outer decoder (same as the normal circuit), or a full decoder that decodes everything together. The full decoder will present superior error rates and also allow for other scenarios such as 2-repetitions (which does not make sense with majority decoder).

For the full-decoder, we build a floquet circuit that has the normal detectors, as well as small detectors between the repeated measurements.

# Input circuit

A stim circuit. Following the project pipeline convention, the input is pre-processed with `circuit.without_noise()` (noise is re-added afterward by the noise model) and `.flattened()` (REPEAT blocks expanded, so measurement-index bookkeeping is simple — same as `translate_circuit` in `src/bell_translation.py`).

Note: measurement layers may be `MXX`/`MZZ`/`MYY` or `MPP` (the honeycomb example uses one `MPP` per round mixing 2-body and 1-body products). The transform repeats measurement *instructions* generically.

Example: `in/circuits/honeycomb_gidney/model=SDEM3,r=10,w=6,h=9,p=0.001,obs=H,sheared=False.stim`

## Scope: which layers are repeated

Only mid-circuit rounds are repeated. The final readout layer (the last measurement layer of the circuit) stays single.

# Output

## Multiply the measurement

Each mid-circuit measurement time step (between two TICKs) is preceded by r-1 new copies of the same measurement instruction, each in its own time step. For example, r=2:

```
...
TICK
MXX 0 1 2 3
DETECTOR rec[-1]
TICK
DETECTOR rec[-8]
```

should become:

```
TICK
MXX 0 1 2 3
TICK
MXX 0 1 2 3
DETECTOR rec[-1]
TICK
DETECTOR rec[-8]
```

## Add comparison detectors

Consecutive copies measure the same operator with nothing in between, so they must agree. For each measurement in the layer, add r-1 detectors comparing copy k with copy k+1. With a layer of N measurements and r=2:

```
TICK
MXX 0 1 2 3
TICK
MXX 0 1 2 3
# new comparison detectors (copy 1 vs copy 2):
DETECTOR rec[-1] rec[-3]
DETECTOR rec[-2] rec[-4]
DETECTOR rec[-1]
TICK
DETECTOR rec[-8]
```

## Fix records (detectors and observables)

Since stim rec targets are relative, every annotation after an inserted copy must be re-offset. Implementation: maintain a `measurement_map: list[list[int]]` from original measurement index to new measurement indices, and rewrite each annotation's rec targets — the same pattern as `_rewrite_annotation_to_buf()` in `src/bell_translation.py`.

### Which copy does an original annotation point to?

To keep the decoding graph matchable, the two original detectors that use a repeated measurement must be split between the first and last copy. Otherwise a measurement error on one copy would flip 3 detectors (the comparison detector plus both originals), which is not matchable.

**Decision: past-facing → first copy, future-facing → last copy.** Together with the comparison detectors, the r copies form a chain in time anchored at both ends; a measurement error on any single copy flips exactly 2 detectors, and data errors between copies stay local in time (the repetition window behaves like extra rounds inserted in time).

**Classification rule (per rec target, not per detector):** a reference to measurement `m` is

- **past-facing** if the annotation appears in the same time step as `m`'s layer (between `m`'s measurement instruction and the next TICK) → remap to the **first** copy;
- **future-facing** if it appears in any later time step → remap to the **last** copy.

A single detector may mix both: e.g. `DETECTOR(1.5, 4, 0) rec[-29] rec[-26] rec[-9]` right after its round's MPP gets `rec[-9]` (current layer) → first copy, and `rec[-29]`, `rec[-26]` (earlier rounds) → last copies of their layers.

### Assertions

The classification rule assumes past-facing detectors sit in the same time step as their measurement. Input circuits might not follow this, so when splitting is enabled the transform must verify it and fail loudly rather than produce a subtly broken graph:

1. **Past-facing usage check:** every measurement in a repeated layer must have at most one past-facing detector reference. Multiple past-facing uses would concentrate too many detectors on the first copy (original past-facing + comparison detectors), breaking matchability that we could otherwise preserve.
2. **Multiple future-facing uses are allowed.** In Floquet codes (e.g. honeycomb), a single measurement commonly contributes to multiple detectors in later rounds. All future-facing references map to the last copy — this preserves the original circuit's structure (the original already had that same measurement used by multiple detectors). `OBSERVABLE_INCLUDE` references are exempt from the past-facing count (observables do not contribute to hyperedge weight in the DEM) but are still classified by the same rule.

Error messages should point the user at `split_detectors=False`.

### Optional splitting

Splitting is controlled by a flag (`split_detectors=True` by default). With `split_detectors=False` — for inputs that are not matchable anyway — all references to a repeated measurement are remapped to one copy (the last), no assertions applied.

### OBSERVABLE_INCLUDE

Handled by the same rule as detectors: classify each rec target as past- or future-facing and remap to the first or last copy accordingly. Any single consistent copy is a valid observable (all copies are deterministically equal in the noiseless circuit; including all copies would change the parity for even r). Note the honeycomb circuit accumulates the observable incrementally — `OBSERVABLE_INCLUDE(1)` appears right after most rounds' MPP with current-layer offsets (past-facing → first copy), not just at the end.

# Pipeline: make_full_decoder_circuits

The full-decoder pipeline composes this transform with the existing noise machinery, in this order:

1. `repeat_measurements(circuit, r, split_detectors=True)` — structural, noiseless.
2. `translate_circuit(noise_params, ...)` from `src/bell_translation.py` — adds Bell gadgets + noise.

This composes cleanly: `translate_circuit` rewrites the comparison detectors' recs like any other detector (each rec expands to the two ancilla recs; matchability is preserved since each ancilla measurement flip still touches at most 2 detectors). Each repeated copy becomes a full Bell gadget with its own prep/internal/measurement noise — including `Z_ERROR(p_data_z)` on the two participating data qubits per copy, so repetition trades measurement-flip suppression against r× data-qubit Z noise per round.

Caveats:

- **Do not combine with the `equivalent_v3_majority*` schemes** — those already fold measurement repetition into the noise block, so composing would double-count it. Use `scheme="bell"` or plain `equivalent_v3`. (Conversely, full-decoder r=3 vs `equivalent_v3_majority` at the same parameters is a good sanity check.)
- **TICK tags must survive:** the stairway 2BGA circuits rely on `TICK[start_noise]` / `TICK[stop_noise]` for `translate_circuit(noiseless_start_end=True)`. `repeat_measurements` must preserve tags on existing TICKs and emit the inserted between-copy TICKs untagged (never duplicating a tagged TICK). When the tags are present, repetition should apply only to layers inside the noise window (repeating noiseless head/tail rounds adds measurements for no benefit).
- Metadata gains the repeat parameter `r` at step 1.

# Implementation notes

- `repeat_measurements(circuit, r, split_detectors=True) -> stim.Circuit` in `src/repeat_measurements/`, following the `buf: list[str]` accumulate-then-parse style of `bell_translation.py`. r=1 must be the identity.
- A `tools/` script applies the transform to input circuits and writes sidecar JSON metadata (input metadata + repeat parameters), following `in/circuits/honeycomb_gidney/gen_metadata.py` conventions.
- Tests mirror `tests/test_bell_translation.py`: r=1 identity, detectors remain deterministic, `detector_error_model(decompose_errors=True)` succeeds when `split_detectors=True`, assertion triggers on a hand-built circuit with a late detector, plus a pymatching smoke decode.
