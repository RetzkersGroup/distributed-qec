"""Bell-pair-based translation of 2-body pairwise measurements.

Replaces MPP / MXX / MZZ / MYY (2-body) instructions with circuits that
implement the same measurement using a Bell pair shared between two
ancilla qubits, optionally adding the `biased_bell_v3` noise model inline
(or v2 — selectable via `NoiseParams.mzz_gadget`; v3 is the default and
only changes the MZZ gadget's shape, not the noise channels).

A 2-body MPP measurement produces 2 ancilla measurement records (instead
of 1 logical record). To preserve the meaning of every DETECTOR /
OBSERVABLE_INCLUDE record reference, each original measurement index is
mapped to its corresponding new measurement indices, and `rec[-k]` is
rewritten accordingly. XOR-ing the two ancilla outcomes yields the
original Pauli-product outcome, so a single-target rec reference becomes
a two-target rec reference.

The input circuit is flattened first (REPEAT blocks expanded) so the
measurement-index map can be computed straightforwardly. Existing noise
on the input is stripped via `circuit.without_noise()`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import stim

PAIRWISE_MEASUREMENT_NAMES = {"MPP", "MXX", "MZZ", "MYY"}

# Single-qubit measurements supported in the input. Z-basis variants get
# X_ERROR semantics (handled via the (p) argument); X-basis variants get
# Z_ERROR semantics. We do not support Y-basis here.
_Z_BASIS_MEAS = {"M", "MR"}
_X_BASIS_MEAS = {"MX", "MRX"}
_SINGLE_QUBIT_MEAS_NAMES = _Z_BASIS_MEAS | _X_BASIS_MEAS

# Reset gates supported in the input.
_Z_BASIS_RESET = {"R"}
_X_BASIS_RESET = {"RX"}
_RESET_NAMES = _Z_BASIS_RESET | _X_BASIS_RESET

# Pure single-qubit Clifford gates that get DEPOLARIZE1(p_1) appended.
_SINGLE_QUBIT_GATE_NAMES = {
    "H", "H_XY", "H_XZ", "H_YZ",
    "S", "S_DAG", "SQRT_X", "SQRT_X_DAG", "SQRT_Y", "SQRT_Y_DAG",
    "SQRT_Z", "SQRT_Z_DAG", "X", "Y", "Z",
    "C_XYZ", "C_ZYX", "I",
}

# Annotations / control flow that pass through unchanged with no noise.
_PASS_THROUGH_NAMES = {
    "TICK", "QUBIT_COORDS", "SHIFT_COORDS",
}


@dataclass
class NoiseParams:
    p_bell_x: float = 0.0
    p_bell_y: float = 0.0
    p_bell_z: float = 0.0
    # Z_ERROR applied to the two participating data qubits during Bell-pair
    # prep (v2 channel; v1 left this implicit via p_internal).
    p_data_z: float = 0.0
    p_internal: float = 0.0
    p_m: float = 0.0
    p_1: float = 0.0
    # Data-qubit state-preparation noise: X_ERROR after R, Z_ERROR after RX.
    p_init: float = 0.0
    # MZZ gadget shape: "v2" -> data->ancilla CX + Z-basis readout
    # (the original biased_bell_v2 form);
    # "v3" -> H-conjugated data + ancilla->data CX + X-basis readout
    # (default; matches biased_bell_v3 in context.md).
    mzz_gadget: str = "v3"

    @classmethod
    def from_dict(cls, d: dict | None) -> "NoiseParams | None":
        if d is None:
            return None
        return cls(**d)

    def to_dict(self) -> dict:
        return {
            "p_bell_x": self.p_bell_x,
            "p_bell_y": self.p_bell_y,
            "p_bell_z": self.p_bell_z,
            "p_data_z": self.p_data_z,
            "p_internal": self.p_internal,
            "p_m": self.p_m,
            "p_1": self.p_1,
            "p_init": self.p_init,
            "mzz_gadget": self.mzz_gadget,
        }


def translate_circuit(
    circuit: stim.Circuit,
    noise_params: dict | NoiseParams | None = None,
    *,
    noiseless_start_end: bool = False,
) -> stim.Circuit:
    """Translate a circuit's pairwise measurements into Bell-mediated circuits.

    Pre-processing: strips noise via `circuit.without_noise()` and flattens
    REPEAT blocks (so measurement-index bookkeeping is simple).

    DETECTOR and OBSERVABLE_INCLUDE record refs are rewritten so the
    translated circuit retains the same logical detectors and observables.

    If `noiseless_start_end=True`, the input must contain matching
    `TICK[start_noise]` and `TICK[stop_noise]` tags. Bell-prep noise is
    only injected for instructions in the window `[start_noise,
    stop_noise)`; ops outside the window are translated structurally with
    no noise channels (useful for noiseless head/tail rounds, e.g. the
    stairway-codes 2BGA circuits). When the flag is off the tags are
    ignored and the whole circuit is noisy.
    """
    if isinstance(noise_params, dict):
        noise = NoiseParams.from_dict(noise_params)
    else:
        noise = noise_params

    cleaned = circuit.without_noise().flattened()
    return _translate_flat(cleaned, noise, noiseless_start_end=noiseless_start_end)


def _translate_flat(
    circuit: stim.Circuit,
    noise: NoiseParams | None,
    *,
    noiseless_start_end: bool = False,
) -> stim.Circuit:
    out = stim.Circuit()
    num_data = circuit.num_qubits
    # Each data qubit q has a dedicated ancilla at index (num_data + q).
    # A measurement on data qubits (q1, q2) uses ancillas (num_data+q1, num_data+q2).
    measurement_map: list[list[int]] = []
    new_meas_count = 0

    # Accumulate all instructions as strings, parse once at the end.
    buf: list[str] = []

    # Optional noiseless head/tail: when `noiseless_start_end` is True we
    # require matching TICK[start_noise] / TICK[stop_noise] tags and only
    # apply noise to ops in the window. When False (default) the tags are
    # ignored and the whole circuit is noisy.
    if noiseless_start_end:
        has_start = any(op.name == "TICK" and op.tag == "start_noise" for op in circuit)
        has_stop = any(op.name == "TICK" and op.tag == "stop_noise" for op in circuit)
        if not has_start or not has_stop:
            raise ValueError(
                "noiseless_start_end=True requires both TICK[start_noise] "
                f"and TICK[stop_noise] in the input; got start={has_start}, "
                f"stop={has_stop}."
            )
        active_noise: NoiseParams | None = None
    else:
        active_noise = noise

    for op in circuit:
        name = op.name
        if name in PAIRWISE_MEASUREMENT_NAMES:
            products = _extract_products(op)
            for basis, qubits in products:
                if len(qubits) == 2:
                    a, b = qubits
                    anc1 = num_data + a
                    anc2 = num_data + b
                    _bell_block_to_buf(buf, basis, a, b, anc1, anc2, active_noise)
                    measurement_map.append(
                        [new_meas_count, new_meas_count + 1]
                    )
                    new_meas_count += 2
                else:  # len(qubits) == 1: 1-body MPP product
                    _single_pauli_meas_to_buf(buf, basis, qubits[0], active_noise)
                    measurement_map.append([new_meas_count])
                    new_meas_count += 1
        elif name in _SINGLE_QUBIT_MEAS_NAMES:
            qubits = [t.qubit_value for t in op.targets_copy()]
            p_m = active_noise.p_m if active_noise is not None else 0.0
            q_str = " ".join(str(q) for q in qubits)
            if p_m > 0:
                buf.append(f"{name}({p_m}) {q_str}")
            else:
                buf.append(f"{name} {q_str}")
            for _ in qubits:
                measurement_map.append([new_meas_count])
                new_meas_count += 1
        elif name in _RESET_NAMES:
            qubits = [t.qubit_value for t in op.targets_copy()]
            q_str = " ".join(str(q) for q in qubits)
            buf.append(f"{name} {q_str}")
            _reset_noise_to_buf(buf, name, qubits, active_noise)
        elif name == "DETECTOR" or name == "OBSERVABLE_INCLUDE":
            _rewrite_annotation_to_buf(
                buf, name, op, measurement_map, new_meas_count
            )
        elif name in _SINGLE_QUBIT_GATE_NAMES:
            qubits = [t.qubit_value for t in op.targets_copy()]
            q_str = " ".join(str(q) for q in qubits)
            buf.append(f"{name} {q_str}")
            if active_noise is not None and active_noise.p_1 > 0 and qubits:
                buf.append(f"DEPOLARIZE1({active_noise.p_1}) {q_str}")
        elif name in _PASS_THROUGH_NAMES:
            buf.append(str(op))
            if noiseless_start_end and name == "TICK":
                if op.tag == "start_noise":
                    active_noise = noise
                elif op.tag == "stop_noise":
                    active_noise = None
        else:
            raise ValueError(
                f"Bell translation has no defined noise model for instruction "
                f"{name!r}. Operations on data qubits must be one of: "
                f"R, RX (state prep); single-qubit Clifford gates; "
                f"single-qubit Z- or X-basis measurements (M, MR, MX, MRX); "
                f"or pairwise measurements (MPP / MXX / MZZ / MYY). "
                f"Two-qubit data-qubit gates (CX, CZ, etc.) are not allowed."
            )

    return stim.Circuit("\n".join(buf))


def _extract_products(
    op: stim.CircuitInstruction,
) -> list[tuple[str, list[int]]]:
    """Return list of (basis, [qubit, ...]) for each product in the op.

    Bell translation supports 1-body and 2-body products. 2-body products
    must be same-basis (XX/YY/ZZ); 1-body products are translated to a
    direct single-qubit measurement.
    """
    name = op.name
    targets = op.targets_copy()

    if name == "MPP":
        # Within an MPP, each "product" is a sequence of pauli targets joined
        # by combiner targets. Products are separated by *no combiner* between
        # consecutive paulis.
        products: list[list[stim.GateTarget]] = [[]]
        prev_was_pauli = False
        for t in targets:
            if t.is_combiner:
                prev_was_pauli = False
                continue
            if prev_was_pauli:
                products.append([])
            products[-1].append(t)
            prev_was_pauli = True

        out: list[tuple[str, list[int]]] = []
        for prod in products:
            if len(prod) == 1:
                out.append((_pauli_of(prod[0]), [prod[0].qubit_value]))
            elif len(prod) == 2:
                basis = _pauli_of(prod[0])
                basis2 = _pauli_of(prod[1])
                if basis != basis2:
                    raise ValueError(
                        f"Bell translation only supports same-basis 2-body MPP "
                        f"products (XX, YY, ZZ); got {basis}{basis2}."
                    )
                out.append(
                    (basis, [prod[0].qubit_value, prod[1].qubit_value])
                )
            else:
                raise ValueError(
                    f"Bell translation only supports 1- or 2-body MPP "
                    f"products; got product of length {len(prod)}."
                )
        return out

    # MXX / MZZ / MYY: targets are bare qubits, paired up.
    basis = name[1]  # 'X', 'Y', or 'Z'
    if len(targets) % 2 != 0:
        raise ValueError(f"{name} requires an even number of targets.")
    return [
        (basis, [targets[i].qubit_value, targets[i + 1].qubit_value])
        for i in range(0, len(targets), 2)
    ]


def _pauli_of(t: stim.GateTarget) -> str:
    if t.is_x_target:
        return "X"
    if t.is_y_target:
        return "Y"
    if t.is_z_target:
        return "Z"
    raise ValueError(f"Expected a Pauli target, got {t!r}")


def _bell_block_to_buf(
    buf: list[str],
    basis: str, a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    if basis == "X":
        _xx_to_buf(buf, a, b, anc1, anc2, noise)
    elif basis == "Z":
        _zz_to_buf(buf, a, b, anc1, anc2, noise)
    elif basis == "Y":
        _yy_to_buf(buf, a, b, anc1, anc2, noise)
    else:
        raise AssertionError(basis)


def _xx_to_buf(
    buf: list[str], a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    buf.append(f"RX {anc1}")
    buf.append(f"R {anc2}")
    buf.append(f"CX {anc1} {anc2}")
    if noise is not None:
        _bell_prep_noise_to_buf(buf, anc2, a, b, noise)
    buf.append(f"CX {anc1} {a}")
    if noise is not None and noise.p_internal > 0:
        buf.append(f"DEPOLARIZE2({noise.p_internal}) {anc1} {a}")
    buf.append(f"CX {anc2} {b}")
    if noise is not None and noise.p_internal > 0:
        buf.append(f"DEPOLARIZE2({noise.p_internal}) {anc2} {b}")
    p_m = noise.p_m if noise is not None else 0.0
    if p_m > 0:
        buf.append(f"MX({p_m}) {anc1} {anc2}")
    else:
        buf.append(f"MX {anc1} {anc2}")


def _zz_to_buf(
    buf: list[str], a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    gadget = noise.mzz_gadget if noise is not None else "v3"
    if gadget == "v2":
        _zz_to_buf_v2(buf, a, b, anc1, anc2, noise)
    elif gadget == "v3":
        _zz_to_buf_v3(buf, a, b, anc1, anc2, noise)
    else:
        raise ValueError(
            f"NoiseParams.mzz_gadget must be 'v2' or 'v3'; got {gadget!r}."
        )


def _zz_to_buf_v2(
    buf: list[str], a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    """biased_bell_v2 MZZ gadget: data->ancilla CX + Z-basis readout."""
    buf.append(f"RX {anc1}")
    buf.append(f"R {anc2}")
    buf.append(f"CX {anc1} {anc2}")
    if noise is not None:
        _bell_prep_noise_to_buf(buf, anc2, a, b, noise)
    buf.append(f"CX {a} {anc1}")
    if noise is not None and noise.p_internal > 0:
        buf.append(f"DEPOLARIZE2({noise.p_internal}) {a} {anc1}")
    buf.append(f"CX {b} {anc2}")
    if noise is not None and noise.p_internal > 0:
        buf.append(f"DEPOLARIZE2({noise.p_internal}) {b} {anc2}")
    p_m = noise.p_m if noise is not None else 0.0
    if p_m > 0:
        buf.append(f"M({p_m}) {anc1} {anc2}")
    else:
        buf.append(f"M {anc1} {anc2}")


def _zz_to_buf_v3(
    buf: list[str], a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    """biased_bell_v3 MZZ gadget: H-conjugate data, ancilla->data CX, X-basis readout."""
    buf.append(f"RX {anc1}")
    buf.append(f"R {anc2}")
    buf.append(f"CX {anc1} {anc2}")
    if noise is not None:
        _bell_prep_noise_to_buf(buf, anc2, a, b, noise)
    buf.append(f"H {a} {b}")
    if noise is not None and noise.p_1 > 0:
        buf.append(f"DEPOLARIZE1({noise.p_1}) {a} {b}")
    buf.append(f"CX {anc1} {a}")
    if noise is not None and noise.p_internal > 0:
        buf.append(f"DEPOLARIZE2({noise.p_internal}) {anc1} {a}")
    buf.append(f"CX {anc2} {b}")
    if noise is not None and noise.p_internal > 0:
        buf.append(f"DEPOLARIZE2({noise.p_internal}) {anc2} {b}")
    buf.append(f"H {a} {b}")
    if noise is not None and noise.p_1 > 0:
        buf.append(f"DEPOLARIZE1({noise.p_1}) {a} {b}")
    p_m = noise.p_m if noise is not None else 0.0
    if p_m > 0:
        buf.append(f"MX({p_m}) {anc1} {anc2}")
    else:
        buf.append(f"MX {anc1} {anc2}")


def _yy_to_buf(
    buf: list[str], a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    buf.append(f"R {anc2}")
    buf.append(f"RX {anc1}")
    buf.append(f"CX {anc1} {anc2}")
    if noise is not None:
        _bell_prep_noise_to_buf(buf, anc2, a, b, noise)
    buf.append(f"S {a} {b}")
    if noise is not None and noise.p_1 > 0:
        buf.append(f"DEPOLARIZE1({noise.p_1}) {a} {b}")
    buf.append(f"CX {anc2} {b} {anc1} {a}")
    if noise is not None and noise.p_internal > 0:
        buf.append(f"DEPOLARIZE2({noise.p_internal}) {anc2} {b}")
        buf.append(f"DEPOLARIZE2({noise.p_internal}) {anc1} {a}")
    p_m = noise.p_m if noise is not None else 0.0
    if p_m > 0:
        buf.append(f"MX({p_m}) {anc2} {anc1}")
    else:
        buf.append(f"MX {anc2} {anc1}")
    buf.append(f"S_DAG {a} {b}")
    if noise is not None and noise.p_1 > 0:
        buf.append(f"DEPOLARIZE1({noise.p_1}) {a} {b}")


def _single_pauli_meas_to_buf(
    buf: list[str], basis: str, q: int, noise: NoiseParams | None,
) -> None:
    """1-body MPP product: emit a direct single-qubit measurement on `q`.

    No Bell pair, no ancilla. Uses M/MX/MY for Z/X/Y respectively, with
    the measurement-flip channel `p_m`.
    """
    inst = {"X": "MX", "Y": "MY", "Z": "M"}[basis]
    p_m = noise.p_m if noise is not None else 0.0
    if p_m > 0:
        buf.append(f"{inst}({p_m}) {q}")
    else:
        buf.append(f"{inst} {q}")


def _bell_prep_noise_to_buf(
    buf: list[str], anc_qubit: int, data_a: int, data_b: int, noise: NoiseParams,
) -> None:
    """Bell-pair-prep noise: PAULI_CHANNEL_1 on one ancilla + Z_ERROR on data."""
    if noise.p_bell_x != 0 or noise.p_bell_y != 0 or noise.p_bell_z != 0:
        buf.append(
            f"PAULI_CHANNEL_1({noise.p_bell_x}, {noise.p_bell_y}, {noise.p_bell_z}) {anc_qubit}"
        )
    if noise.p_data_z != 0:
        buf.append(f"Z_ERROR({noise.p_data_z}) {data_a} {data_b}")


def _reset_noise_to_buf(
    buf: list[str], reset_name: str, qubits: list[int],
    noise: NoiseParams | None,
) -> None:
    if noise is None or noise.p_init == 0 or not qubits:
        return
    q_str = " ".join(str(q) for q in qubits)
    if reset_name in _Z_BASIS_RESET:
        buf.append(f"X_ERROR({noise.p_init}) {q_str}")
    elif reset_name in _X_BASIS_RESET:
        buf.append(f"Z_ERROR({noise.p_init}) {q_str}")


def _rewrite_annotation_to_buf(
    buf: list[str],
    name: str,
    op: stim.CircuitInstruction,
    measurement_map: list[list[int]],
    new_meas_count: int,
) -> None:
    """Rewrite DETECTOR/OBSERVABLE_INCLUDE rec targets and append as string."""
    targets = op.targets_copy()
    args = op.gate_args_copy()
    n_old = len(measurement_map)

    rec_parts: list[str] = []
    for t in targets:
        if t.is_measurement_record_target:
            k = -t.value
            old_idx = n_old - k
            if old_idx < 0:
                raise ValueError(f"rec[-{k}] precedes any measurement.")
            for new_idx in measurement_map[old_idx]:
                offset = new_idx - new_meas_count
                rec_parts.append(f"rec[{offset}]")
        else:
            rec_parts.append(str(t))

    if name == "DETECTOR":
        if args:
            coords = ", ".join(str(a) for a in args)
            buf.append(f"DETECTOR({coords}) {' '.join(rec_parts)}")
        else:
            buf.append(f"DETECTOR {' '.join(rec_parts)}")
    else:
        idx = int(args[0]) if args else 0
        buf.append(f"OBSERVABLE_INCLUDE({idx}) {' '.join(rec_parts)}")
