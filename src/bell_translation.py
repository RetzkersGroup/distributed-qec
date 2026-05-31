"""Bell-pair-based translation of 2-body pairwise measurements.

Replaces MPP / MXX / MZZ / MYY (2-body) instructions with circuits that
implement the same measurement using a Bell pair shared between two
ancilla qubits, optionally adding the `biased_bell` noise model inline.

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
    p_init_x: float = 0.0
    p_init_y: float = 0.0
    p_init_z: float = 0.0
    p_internal: float = 0.0
    p_m: float = 0.0
    p_1: float = 0.0
    # Data-qubit state-preparation noise: X_ERROR after R, Z_ERROR after RX.
    p_init: float = 0.0

    @classmethod
    def from_dict(cls, d: dict | None) -> "NoiseParams | None":
        if d is None:
            return None
        return cls(**d)

    def to_dict(self) -> dict:
        return {
            "p_init_x": self.p_init_x,
            "p_init_y": self.p_init_y,
            "p_init_z": self.p_init_z,
            "p_internal": self.p_internal,
            "p_m": self.p_m,
            "p_1": self.p_1,
            "p_init": self.p_init,
        }


def translate_circuit(
    circuit: stim.Circuit,
    noise_params: dict | NoiseParams | None = None,
) -> stim.Circuit:
    """Translate a circuit's pairwise measurements into Bell-mediated circuits.

    Pre-processing: strips noise via `circuit.without_noise()` and flattens
    REPEAT blocks (so measurement-index bookkeeping is simple).

    DETECTOR and OBSERVABLE_INCLUDE record refs are rewritten so the
    translated circuit retains the same logical detectors and observables.
    """
    if isinstance(noise_params, dict):
        noise = NoiseParams.from_dict(noise_params)
    else:
        noise = noise_params

    cleaned = circuit.without_noise().flattened()
    return _translate_flat(cleaned, noise)


def _translate_flat(
    circuit: stim.Circuit, noise: NoiseParams | None
) -> stim.Circuit:
    out = stim.Circuit()
    num_data = circuit.num_qubits
    # Each data qubit q has a dedicated ancilla at index (num_data + q).
    # A measurement on data qubits (q1, q2) uses ancillas (num_data+q1, num_data+q2).
    measurement_map: list[list[int]] = []
    new_meas_count = 0

    for op in circuit:
        name = op.name
        if name in PAIRWISE_MEASUREMENT_NAMES:
            pairs = _extract_pairs(op)
            for basis, a, b in pairs:
                anc1 = num_data + a
                anc2 = num_data + b
                _emit_bell_block(out, basis, a, b, anc1, anc2, noise)
                measurement_map.append([new_meas_count, new_meas_count + 1])
                new_meas_count += 2
        elif name in _SINGLE_QUBIT_MEAS_NAMES:
            qubits = [t.qubit_value for t in op.targets_copy()]
            p_m = noise.p_m if noise is not None else 0.0
            out.append(name, qubits, p_m)
            for _ in qubits:
                measurement_map.append([new_meas_count])
                new_meas_count += 1
        elif name in _RESET_NAMES:
            qubits = [t.qubit_value for t in op.targets_copy()]
            out.append(name, qubits)
            _add_reset_noise(out, name, qubits, noise)
        elif name == "DETECTOR":
            new_targets = _rewrite_rec_targets(
                op.targets_copy(), measurement_map, new_meas_count
            )
            out.append("DETECTOR", new_targets, op.gate_args_copy())
        elif name == "OBSERVABLE_INCLUDE":
            new_targets = _rewrite_rec_targets(
                op.targets_copy(), measurement_map, new_meas_count
            )
            out.append("OBSERVABLE_INCLUDE", new_targets, op.gate_args_copy())
        elif name in _SINGLE_QUBIT_GATE_NAMES:
            out.append(op)
            _add_1q_noise_after(out, op, noise)
        elif name in _PASS_THROUGH_NAMES:
            out.append(op)
        else:
            raise ValueError(
                f"Bell translation has no defined noise model for instruction "
                f"{name!r}. Operations on data qubits must be one of: "
                f"R, RX (state prep); single-qubit Clifford gates; "
                f"single-qubit Z- or X-basis measurements (M, MR, MX, MRX); "
                f"or pairwise measurements (MPP / MXX / MZZ / MYY). "
                f"Two-qubit data-qubit gates (CX, CZ, etc.) are not allowed."
            )
    return out


def _extract_pairs(op: stim.CircuitInstruction) -> list[tuple[str, int, int]]:
    """Return list of (basis, q1, q2) for each 2-body product in the op."""
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

        pairs: list[tuple[str, int, int]] = []
        for prod in products:
            if len(prod) != 2:
                raise ValueError(
                    f"Bell translation only supports 2-body MPP products; "
                    f"got product of length {len(prod)}."
                )
            basis = _pauli_of(prod[0])
            basis2 = _pauli_of(prod[1])
            if basis != basis2:
                raise ValueError(
                    f"Bell translation only supports same-basis 2-body MPP "
                    f"products (XX, YY, ZZ); got {basis}{basis2}."
                )
            pairs.append((basis, prod[0].qubit_value, prod[1].qubit_value))
        return pairs

    # MXX / MZZ / MYY: targets are bare qubits, paired up.
    basis = name[1]  # 'X', 'Y', or 'Z'
    if len(targets) % 2 != 0:
        raise ValueError(f"{name} requires an even number of targets.")
    return [
        (basis, targets[i].qubit_value, targets[i + 1].qubit_value)
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


def _emit_bell_block(
    out: stim.Circuit,
    basis: str, a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    if basis == "X":
        _emit_xx(out, a, b, anc1, anc2, noise)
    elif basis == "Z":
        _emit_zz(out, a, b, anc1, anc2, noise)
    elif basis == "Y":
        _emit_yy(out, a, b, anc1, anc2, noise)
    else:
        raise AssertionError(basis)


def _emit_xx(
    out: stim.Circuit, a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    # RX anc1; R anc2; CX anc1 anc2; CX anc1 a; CX anc2 b; MX anc1 anc2
    out.append("RX", [anc1])
    out.append("R", [anc2])
    out.append("CX", [anc1, anc2])
    if noise is not None:
        _add_init_noise(out, anc2, noise)
    out.append("CX", [anc1, a])
    if noise is not None and noise.p_internal > 0:
        out.append("DEPOLARIZE2", [anc1, a], noise.p_internal)
    out.append("CX", [anc2, b])
    if noise is not None and noise.p_internal > 0:
        out.append("DEPOLARIZE2", [anc2, b], noise.p_internal)
    p_m = noise.p_m if noise is not None else 0.0
    out.append("MX", [anc1, anc2], p_m)


def _emit_zz(
    out: stim.Circuit, a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    # RX anc1; R anc2; CX anc1 anc2; CX a anc1; CX b anc2; M anc1 anc2
    out.append("RX", [anc1])
    out.append("R", [anc2])
    out.append("CX", [anc1, anc2])
    if noise is not None:
        _add_init_noise(out, anc2, noise)
    out.append("CX", [a, anc1])
    if noise is not None and noise.p_internal > 0:
        out.append("DEPOLARIZE2", [a, anc1], noise.p_internal)
    out.append("CX", [b, anc2])
    if noise is not None and noise.p_internal > 0:
        out.append("DEPOLARIZE2", [b, anc2], noise.p_internal)
    p_m = noise.p_m if noise is not None else 0.0
    out.append("M", [anc1, anc2], p_m)


def _emit_yy(
    out: stim.Circuit, a: int, b: int, anc1: int, anc2: int,
    noise: NoiseParams | None,
) -> None:
    # R anc2; RX anc1; CX anc1 anc2; S a b; CX anc2 b anc1 a; MX anc2 anc1; S_DAG a b
    out.append("R", [anc2])
    out.append("RX", [anc1])
    out.append("CX", [anc1, anc2])
    if noise is not None:
        _add_init_noise(out, anc2, noise)
    out.append("S", [a, b])
    if noise is not None and noise.p_1 > 0:
        out.append("DEPOLARIZE1", [a, b], noise.p_1)
    out.append("CX", [anc2, b, anc1, a])
    if noise is not None and noise.p_internal > 0:
        out.append("DEPOLARIZE2", [anc2, b], noise.p_internal)
        out.append("DEPOLARIZE2", [anc1, a], noise.p_internal)
    p_m = noise.p_m if noise is not None else 0.0
    out.append("MX", [anc2, anc1], p_m)
    out.append("S_DAG", [a, b])
    if noise is not None and noise.p_1 > 0:
        out.append("DEPOLARIZE1", [a, b], noise.p_1)


def _add_init_noise(out: stim.Circuit, qubit: int, noise: NoiseParams) -> None:
    if noise.p_init_x == 0 and noise.p_init_y == 0 and noise.p_init_z == 0:
        return
    out.append(
        "PAULI_CHANNEL_1",
        [qubit],
        (noise.p_init_x, noise.p_init_y, noise.p_init_z),
    )


def _add_1q_noise_after(
    out: stim.Circuit, op: stim.CircuitInstruction, noise: NoiseParams | None
) -> None:
    if noise is None or noise.p_1 == 0:
        return
    qubits = [t.qubit_value for t in op.targets_copy() if t.qubit_value is not None]
    if qubits:
        out.append("DEPOLARIZE1", qubits, noise.p_1)


def _add_reset_noise(
    out: stim.Circuit, reset_name: str, qubits: list[int],
    noise: NoiseParams | None,
) -> None:
    """Add state-prep noise after a reset.

    R prepares |0>, an eigenstate of Z, so the dominant noise is a Pauli that
    anticommutes with Z — i.e. X_ERROR. Symmetrically, RX prepares |+> and
    gets Z_ERROR.
    """
    if noise is None or noise.p_init == 0 or not qubits:
        return
    if reset_name in _Z_BASIS_RESET:
        out.append("X_ERROR", qubits, noise.p_init)
    elif reset_name in _X_BASIS_RESET:
        out.append("Z_ERROR", qubits, noise.p_init)


def _rewrite_rec_targets(
    targets: Iterable[stim.GateTarget],
    measurement_map: list[list[int]],
    new_meas_count: int,
) -> list[stim.GateTarget]:
    """Map every rec[-k] target to the corresponding new rec targets.

    For a pairwise measurement that became two ancilla measurements, the
    XOR of those two ancilla outcomes equals the original outcome, so a
    single rec target becomes two rec targets.
    """
    n_old = len(measurement_map)
    out: list[stim.GateTarget] = []
    for t in targets:
        if t.is_measurement_record_target:
            k = -t.value  # k > 0
            old_idx = n_old - k
            if old_idx < 0:
                raise ValueError(f"rec[-{k}] precedes any measurement.")
            for new_idx in measurement_map[old_idx]:
                offset = new_idx - new_meas_count  # negative
                out.append(stim.target_rec(offset))
        else:
            out.append(t)
    return out
