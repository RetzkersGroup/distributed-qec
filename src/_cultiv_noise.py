import collections
from typing import Iterator, AbstractSet, DefaultDict, Any

import stim

"""
This module is a verbatim copy from Strilanc/magic-state-cultivation
"""


ANNOTATION_OPS = {
    "DETECTOR",
    "OBSERVABLE_INCLUDE",
    "QUBIT_COORDS",
    "SHIFT_COORDS",
    "TICK",
    "MPAD",
}
OP_TO_MEASURE_BASES = {
    "M": "Z",
    "MR": "Z",
    "MX": "X",
    "MY": "Y",
    "MZ": "Z",
    "MRX": "X",
    "MRY": "Y",
    "MRZ": "Z",
    "MXX": "XX",
    "MYY": "YY",
    "MZZ": "ZZ",
    "MPP": "*",
}


class NoiseRule:
    """Describes how to add noise to an operation."""

    def __init__(
        self,
        *,
        before: dict[str, float | tuple[float, ...]] | None = None,
        after: dict[str, float | tuple[float, ...]] | None = None,
        flip_result: float = 0,
    ):
        if after is None:
            after = {}
        if before is None:
            before = {}
        if not (0 <= flip_result <= 1):
            raise ValueError(f"not (0 <= {flip_result=} <= 1)")
        for k, p in [*after.items(), *before.items()]:
            gate_data = stim.gate_data(k)
            if gate_data.produces_measurements or not gate_data.is_noisy_gate:
                raise ValueError(f"not a pure noise channel: {k} from {after=}")
            if gate_data.num_parens_arguments_range == range(1, 2):
                if not isinstance(p, (int, float)) or not (0 <= p <= 1):
                    raise ValueError(f"not a probability: {p!r}")
            else:
                if len(p) not in gate_data.num_parens_arguments_range:
                    raise ValueError(f"Wrong number of arguments {p!r} for gate {k!r}")
                if not isinstance(p, (list, tuple)) or not (0 <= sum(p) <= 1):
                    raise ValueError(f"not a tuple of disjoint probabilities: {p!r}")
        self.before = before
        self.after = after
        self.flip_result = flip_result

    def append_noisy_version_of(
        self,
        *,
        split_op: stim.CircuitInstruction,
        out_during_moment: stim.Circuit,
        before_moments: DefaultDict[Any, stim.Circuit],
        after_moments: DefaultDict[Any, stim.Circuit],
        immune_qubit_indices: AbstractSet[int],
    ) -> None:
        targets = split_op.targets_copy()
        if immune_qubit_indices and any(
            (t.is_qubit_target or t.is_x_target or t.is_y_target or t.is_z_target)
            and t.value in immune_qubit_indices
            for t in targets
        ):
            out_during_moment.append(split_op)
            return

        args = split_op.gate_args_copy()
        if self.flip_result:
            gate_data = stim.gate_data(split_op.name)
            assert gate_data.produces_measurements
            assert gate_data.is_noisy_gate
            assert gate_data.num_parens_arguments_range == range(0, 2)
            assert len(args) == 0
            args = [self.flip_result]

        out_during_moment.append(split_op.name, targets, args)
        raw_targets = [t.value for t in targets if not t.is_combiner]
        for op_name, arg in self.before.items():
            before_moments[(op_name, arg)].append(op_name, raw_targets, arg)
        for op_name, arg in self.after.items():
            after_moments[(op_name, arg)].append(op_name, raw_targets, arg)


class NoiseModel:
    def __init__(
        self,
        idle_depolarization: float = 0,
        tick_noise: NoiseRule | None = None,
        additional_depolarization_waiting_for_m_or_r: float = 0,
        gate_rules: dict[str, NoiseRule] | None = None,
        measure_rules: dict[str, NoiseRule] | None = None,
        any_measurement_rule: NoiseRule | None = None,
        any_clifford_1q_rule: NoiseRule | float | None = None,
        any_clifford_2q_rule: NoiseRule | float | None = None,
        allow_multiple_uses_of_a_qubit_in_one_tick: bool = False,
    ):
        if isinstance(any_clifford_1q_rule, float):
            any_clifford_1q_rule = NoiseRule(
                after={"DEPOLARIZE1": any_clifford_1q_rule}
            )
        if isinstance(any_clifford_2q_rule, float):
            any_clifford_2q_rule = NoiseRule(
                after={"DEPOLARIZE2": any_clifford_2q_rule}
            )
        self.idle_depolarization = idle_depolarization
        self.tick_noise = tick_noise
        self.additional_depolarization_waiting_for_m_or_r = (
            additional_depolarization_waiting_for_m_or_r
        )
        self.gate_rules = {} if gate_rules is None else gate_rules
        self.measure_rules = measure_rules
        self.any_measurement_rule = any_measurement_rule
        self.any_clifford_1q_rule = any_clifford_1q_rule
        self.any_clifford_2q_rule = any_clifford_2q_rule
        self.allow_multiple_uses_of_a_qubit_in_one_tick = (
            allow_multiple_uses_of_a_qubit_in_one_tick
        )
        assert self.tick_noise is None or not self.tick_noise.flip_result

    @staticmethod
    def si1000(p: float) -> "NoiseModel":
        """Superconducting inspired noise.

        As defined in "A Fault-Tolerant Honeycomb Memory" https://arxiv.org/abs/2108.10457
        """
        return NoiseModel(
            idle_depolarization=p / 10,
            additional_depolarization_waiting_for_m_or_r=2 * p,
            any_clifford_1q_rule=NoiseRule(after={"DEPOLARIZE1": p / 10}),
            any_clifford_2q_rule=NoiseRule(after={"DEPOLARIZE2": p}),
            measure_rules={
                "Z": NoiseRule(after={"DEPOLARIZE1": p}, flip_result=p * 5),
                "ZZ": NoiseRule(after={"DEPOLARIZE2": p}, flip_result=p * 5),
            },
            gate_rules={
                "R": NoiseRule(after={"X_ERROR": p * 2}),
            },
        )

    @staticmethod
    def uniform_depolarizing(
        p: float, *, single_qubit_only: bool = False
    ) -> "NoiseModel":
        """Near-standard circuit depolarizing noise.

        Everything has the same parameter p.
        """
        dep2 = "DEPOLARIZE1" if single_qubit_only else "DEPOLARIZE2"
        return NoiseModel(
            idle_depolarization=p,
            any_clifford_1q_rule=NoiseRule(after={"DEPOLARIZE1": p}),
            any_clifford_2q_rule=NoiseRule(after={dep2: p}),
            measure_rules={
                "X": NoiseRule(after={"DEPOLARIZE1": p}, flip_result=p),
                "Y": NoiseRule(after={"DEPOLARIZE1": p}, flip_result=p),
                "Z": NoiseRule(after={"DEPOLARIZE1": p}, flip_result=p),
                "XX": NoiseRule(after={dep2: p}, flip_result=p),
                "XY": NoiseRule(after={dep2: p}, flip_result=p),
                "XZ": NoiseRule(after={dep2: p}, flip_result=p),
                "YX": NoiseRule(after={dep2: p}, flip_result=p),
                "YY": NoiseRule(after={dep2: p}, flip_result=p),
                "YZ": NoiseRule(after={dep2: p}, flip_result=p),
                "ZX": NoiseRule(after={dep2: p}, flip_result=p),
                "ZY": NoiseRule(after={dep2: p}, flip_result=p),
                "ZZ": NoiseRule(after={dep2: p}, flip_result=p),
            },
            gate_rules={
                "RX": NoiseRule(after={"Z_ERROR": p}),
                "RY": NoiseRule(after={"X_ERROR": p}),
                "R": NoiseRule(after={"X_ERROR": p}),
            },
        )

    @classmethod
    def biased_noise(cls, p: float, p_d: float, p_m: float) -> 'NoiseModel':
        return cls(
            idle_depolarization=p/10,
            any_clifford_1q_rule=NoiseRule(after={'DEPOLARIZE1': p/10}),
            any_clifford_2q_rule=NoiseRule(after={'DEPOLARIZE2': p}),
            measure_rules={
                'X': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
                'Y': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
                'Z': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
                'XX': NoiseRule(after={'DEPOLARIZE2': p_d}, flip_result=p_m),
                'YY': NoiseRule(after={'DEPOLARIZE2': p_d}, flip_result=p_m),
                'ZZ': NoiseRule(after={'DEPOLARIZE2': p_d}, flip_result=p_m),
            },
            gate_rules={
                'RX': NoiseRule(after={'Z_ERROR': p}),
                'RY': NoiseRule(after={'X_ERROR': p}),
                'R': NoiseRule(after={'X_ERROR': p}),
            }
        )

    @classmethod
    def normalized_pairwise_errors(cls, p: float) -> 'NoiseModel':
        """
        normalize the error of a pairwise error to be p per gate. Such that: p = p_d + p_m - p_d*p_m
        """
        p_m = 1 - (1-p)**0.5
        p_d = 1 - (1-p)**0.5
        return cls.biased_noise(p=p, p_d=p_d, p_m=p_m)

    @classmethod
    def flexible_pairwise(cls, p_m_bulk: float, p_spam: float, p_res: float) -> 'NoiseModel':
        """Pairwise-measurement noise with independent measurement, SPAM, and residual depolarizing rates.

        No idle noise, no single-qubit gate noise.
        """
        pairwise_after = {'DEPOLARIZE1': p_res} if p_res > 0 else {}
        spam_after_r = {'X_ERROR': p_spam} if p_spam > 0 else {}
        spam_after_rx = {'Z_ERROR': p_spam} if p_spam > 0 else {}
        return cls(
            idle_depolarization=0,
            any_clifford_1q_rule=None,
            any_clifford_2q_rule=None,
            allow_multiple_uses_of_a_qubit_in_one_tick=True,
            measure_rules={
                'X': NoiseRule(flip_result=p_spam),
                'Y': NoiseRule(flip_result=p_spam),
                'Z': NoiseRule(flip_result=p_spam),
                'XX': NoiseRule(after=pairwise_after, flip_result=p_m_bulk),
                'YY': NoiseRule(after=pairwise_after, flip_result=p_m_bulk),
                'ZZ': NoiseRule(after=pairwise_after, flip_result=p_m_bulk),
            },
            gate_rules={
                'R': NoiseRule(after=spam_after_r),
                'RX': NoiseRule(after=spam_after_rx),
            },
        )

    def _noise_rule_for_split_operation(
        self,
        *,
        split_op: stim.CircuitInstruction,
    ) -> NoiseRule | None:
        if occurs_in_classical_control_system(split_op):
            return None

        rule = self.gate_rules.get(split_op.name)
        if rule is not None:
            return rule

        gate_data = stim.gate_data(split_op.name)

        if (
            self.any_clifford_1q_rule is not None
            and gate_data.is_unitary
            and gate_data.is_single_qubit_gate
        ):
            return self.any_clifford_1q_rule
        if (
            self.any_clifford_2q_rule is not None
            and gate_data.is_unitary
            and gate_data.is_two_qubit_gate
        ):
            return self.any_clifford_2q_rule
        if self.measure_rules is not None:
            rule = self.measure_rules.get(_measure_basis(split_op=split_op))
            if rule is not None:
                return rule
        if self.any_measurement_rule is not None and gate_data.produces_measurements:
            return self.any_measurement_rule
        if gate_data.is_reset and gate_data.produces_measurements:
            m_name, r_name = {
                "MRX": ("MX", "RX"),
                "MRY": ("MY", "RY"),
                "MR": ("M", "R"),
            }[gate_data.name]
            r_noise = self._noise_rule_for_split_operation(
                split_op=stim.CircuitInstruction(r_name, split_op.targets_copy())
            )
            m_noise = self._noise_rule_for_split_operation(
                split_op=stim.CircuitInstruction(m_name, split_op.targets_copy())
            )
            return NoiseRule(
                before=r_noise.before if r_noise is not None else {},
                after=r_noise.after if r_noise is not None else {},
                flip_result=m_noise.flip_result if m_noise is not None else 0,
            )

    def _append_idle_error(
        self,
        *,
        moment_split_ops: list[stim.CircuitInstruction],
        out: stim.Circuit,
        system_qubit_indices: AbstractSet[int],
        immune_qubit_indices: AbstractSet[int],
    ) -> None:
        collapse_qubits = []
        clifford_qubits = []
        pauli_qubits = []
        for split_op in moment_split_ops:
            if occurs_in_classical_control_system(split_op):
                continue
            gate_data = stim.gate_data(split_op.name)
            if gate_data.is_reset or gate_data.produces_measurements:
                qubits_out = collapse_qubits
            elif split_op.name in "IXYZ":
                qubits_out = pauli_qubits
            elif gate_data.is_unitary:
                qubits_out = clifford_qubits
            for target in split_op.targets_copy():
                if not target.is_combiner:
                    qubits_out.append(target.value)

        usage_counts = collections.Counter(collapse_qubits + clifford_qubits)
        for p in pauli_qubits:
            if usage_counts[p] == 0:
                usage_counts[p] += 1
        qubits_used_multiple_times = {q for q, c in usage_counts.items() if c != 1}
        if (
            qubits_used_multiple_times
            and not self.allow_multiple_uses_of_a_qubit_in_one_tick
        ):
            moment = stim.Circuit()
            for op in moment_split_ops:
                moment.append(op)
            raise ValueError(
                f"Qubits were operated on multiple times without a TICK in between:\n"
                f"multiple uses: {sorted(qubits_used_multiple_times)}\n"
                f"moment:\n"
                f"{moment}"
            )

        collapse_qubits_set = set(collapse_qubits)
        clifford_qubits_set = set(clifford_qubits + pauli_qubits)
        idle = sorted(
            system_qubit_indices
            - collapse_qubits_set
            - clifford_qubits_set
            - immune_qubit_indices
        )
        if idle and self.idle_depolarization:
            out.append("DEPOLARIZE1", idle, self.idle_depolarization)

        waiting_for_mr = sorted(
            system_qubit_indices - collapse_qubits_set - immune_qubit_indices
        )
        if (
            collapse_qubits_set
            and waiting_for_mr
            and self.additional_depolarization_waiting_for_m_or_r
        ):
            out.append(
                "DEPOLARIZE1", idle, self.additional_depolarization_waiting_for_m_or_r
            )

        if self.tick_noise is not None:
            for k, p in self.tick_noise.before.items():
                out.append(k, system_qubit_indices - immune_qubit_indices, p)
            for k, p in self.tick_noise.after.items():
                out.append(k, system_qubit_indices - immune_qubit_indices, p)

    def _append_noisy_moment(
        self,
        *,
        moment_split_ops: list[stim.CircuitInstruction],
        out: stim.Circuit,
        system_qubits_indices: AbstractSet[int],
        immune_qubit_indices: AbstractSet[int],
    ) -> None:
        skip_pauli_targets = set()
        for split_op in moment_split_ops:
            gate_data = stim.gate_data(split_op.name)
            if (
                gate_data.is_unitary
                and gate_data.is_single_qubit_gate
                and not split_op.name in "IXYZ"
            ):
                for t in split_op.targets_copy():
                    skip_pauli_targets.add(t.qubit_value)

        before = collections.defaultdict(stim.Circuit)
        after = collections.defaultdict(stim.Circuit)
        grow = stim.Circuit()
        for split_op in moment_split_ops:
            rule = self._noise_rule_for_split_operation(split_op=split_op)
            if rule is None:
                grow.append(split_op)
            elif split_op.name in "IXYZ":
                new_targets = []
                skipped_targets = []
                for t in split_op.targets_copy():
                    if t.qubit_value in skip_pauli_targets:
                        skipped_targets.append(t)
                    else:
                        new_targets.append(t)
                        skip_pauli_targets.add(t.qubit_value)
                if skipped_targets:
                    grow.append(
                        stim.CircuitInstruction(
                            split_op.name, skipped_targets, split_op.gate_args_copy()
                        )
                    )
                if new_targets:
                    rule.append_noisy_version_of(
                        split_op=stim.CircuitInstruction(
                            split_op.name, new_targets, split_op.gate_args_copy()
                        ),
                        out_during_moment=grow,
                        before_moments=before,
                        after_moments=after,
                        immune_qubit_indices=immune_qubit_indices,
                    )
            else:
                rule.append_noisy_version_of(
                    split_op=split_op,
                    out_during_moment=grow,
                    before_moments=before,
                    after_moments=after,
                    immune_qubit_indices=immune_qubit_indices,
                )
        for k in sorted(before.keys()):
            out += before[k]
        out += grow
        for k in sorted(after.keys()):
            out += after[k]

        self._append_idle_error(
            moment_split_ops=moment_split_ops,
            out=out,
            system_qubit_indices=system_qubits_indices,
            immune_qubit_indices=immune_qubit_indices,
        )

    def noisy_circuit_skipping_mpp_boundaries(
        self,
        circuit: stim.Circuit,
        *,
        immune_qubit_indices: AbstractSet[int] | None = None,
    ) -> stim.Circuit:
        """Adds noise to the circuit except for MPP operations at the start/end."""
        allowed = {
            "TICK",
            "OBSERVABLE_INCLUDE",
            "DETECTOR",
            "MPP",
            "QUBIT_COORDS",
            "SHIFT_COORDS",
        }
        start = 0
        end = len(circuit)
        while start < len(circuit) and circuit[start].name in allowed:
            start += 1
        while end > 0 and circuit[end - 1].name in allowed:
            end -= 1
        while end < len(circuit) and circuit[end].name != "MPP":
            end += 1
        while end > 0 and circuit[end - 1].name == "TICK":
            end -= 1
        if end <= start:
            raise ValueError("end <= start")
        noisy = self.noisy_circuit(
            circuit[start:end], immune_qubit_indices=immune_qubit_indices
        )
        return circuit[:start] + noisy + circuit[end:]

    def noisy_circuit(
        self,
        circuit: stim.Circuit,
        *,
        system_qubit_indices: AbstractSet[int] | None = None,
        immune_qubit_indices: AbstractSet[int] | None = None,
    ) -> stim.Circuit:
        """Returns a noisy version of the given circuit, by applying the receiving noise model."""
        if system_qubit_indices is None:
            system_qubit_indices = set(range(circuit.num_qubits))
        if immune_qubit_indices is None:
            immune_qubit_indices = set()

        result = stim.Circuit()

        first = True
        for moment_split_ops in _iter_split_op_moments(
            circuit, immune_qubit_indices=immune_qubit_indices
        ):
            if first:
                first = False
            elif result and isinstance(result[-1], stim.CircuitRepeatBlock):
                pass
            else:
                result.append("TICK")
            if isinstance(moment_split_ops, stim.CircuitRepeatBlock):
                noisy_body = self.noisy_circuit(
                    moment_split_ops.body_copy(),
                    system_qubit_indices=system_qubit_indices,
                    immune_qubit_indices=immune_qubit_indices,
                )
                noisy_body.append("TICK")
                result.append(
                    stim.CircuitRepeatBlock(
                        repeat_count=moment_split_ops.repeat_count, body=noisy_body
                    )
                )
            else:
                self._append_noisy_moment(
                    moment_split_ops=moment_split_ops,
                    out=result,
                    system_qubits_indices=system_qubit_indices,
                    immune_qubit_indices=immune_qubit_indices,
                )

        return result


def occurs_in_classical_control_system(op: stim.CircuitInstruction) -> bool:
    if op.name in ANNOTATION_OPS:
        return True

    gate_data = stim.gate_data(op.name)
    if gate_data.is_unitary and gate_data.is_two_qubit_gate:
        targets = op.targets_copy()
        for k in range(0, len(targets), 2):
            a = targets[k]
            b = targets[k + 1]
            classical_0 = a.is_measurement_record_target or a.is_sweep_bit_target
            classical_1 = b.is_measurement_record_target or b.is_sweep_bit_target
            if not (classical_0 or classical_1):
                return False
        return True
    return False


def _split_targets_if_needed(
    op: stim.CircuitInstruction,
    immune_qubit_indices: AbstractSet[int],
) -> list[stim.CircuitInstruction]:
    gate_data = stim.gate_data(op.name)
    if gate_data.is_unitary and gate_data.is_two_qubit_gate:
        yield from _split_targets_if_needed_clifford_2q(op, immune_qubit_indices)
    elif op.name == "MPP":
        yield from _split_targets_if_needed_m_basis(op)
    elif op.name in ANNOTATION_OPS:
        yield op
    elif gate_data.is_noisy_gate and not gate_data.produces_measurements:
        yield op
    elif gate_data.is_single_qubit_gate:
        yield from _split_out_immune_targets_assuming_single_qubit_gate(
            op, immune_qubit_indices
        )
    elif gate_data.is_two_qubit_gate:
        yield from _split_out_immune_targets_assuming_two_qubit_gate(
            op, immune_qubit_indices
        )
    else:
        raise NotImplementedError(f"{op=}")


def _split_out_immune_targets_assuming_single_qubit_gate(
    op: stim.CircuitInstruction,
    immune_qubit_indices: AbstractSet[int],
) -> list[stim.CircuitInstruction]:
    if immune_qubit_indices:
        args = op.gate_args_copy()
        for t in op.targets_copy():
            yield stim.CircuitInstruction(op.name, [t], args)
    else:
        yield op


def _split_out_immune_targets_assuming_two_qubit_gate(
    op: stim.CircuitInstruction,
    immune_qubit_indices: AbstractSet[int],
) -> list[stim.CircuitInstruction]:
    if immune_qubit_indices:
        args = op.gate_args_copy()
        targets = op.targets_copy()
        for k in range(len(targets)):
            t1 = targets[k]
            t2 = targets[k + 1]
            yield stim.CircuitInstruction(op.name, [t1, t2], args)
    else:
        yield op


def _split_targets_if_needed_clifford_2q(
    op: stim.CircuitInstruction,
    immune_qubit_indices: AbstractSet[int],
) -> list[stim.CircuitInstruction]:
    gate_data = stim.gate_data(op.name)
    assert gate_data.is_unitary and gate_data.is_two_qubit_gate
    targets = op.targets_copy()
    if immune_qubit_indices or any(t.is_measurement_record_target for t in targets):
        args = op.gate_args_copy()
        for k in range(0, len(targets), 2):
            yield stim.CircuitInstruction(op.name, targets[k : k + 2], args)
    else:
        yield op


def _split_targets_if_needed_m_basis(
    op: stim.CircuitInstruction,
) -> list[stim.CircuitInstruction]:
    targets = op.targets_copy()
    args = op.gate_args_copy()
    k = 0
    start = k
    while k < len(targets):
        if k + 1 == len(targets) or not targets[k + 1].is_combiner:
            yield stim.CircuitInstruction(op.name, targets[start : k + 1], args)
            k += 1
            start = k
        else:
            k += 2
    assert k == len(targets)


def _iter_split_op_moments(
    circuit: stim.Circuit,
    *,
    immune_qubit_indices: AbstractSet[int],
) -> Iterator[stim.CircuitRepeatBlock | list[stim.CircuitInstruction]]:
    cur_moment = []

    for op in circuit:
        if isinstance(op, stim.CircuitRepeatBlock):
            if cur_moment:
                yield cur_moment
                cur_moment = []
            yield op
        elif isinstance(op, stim.CircuitInstruction):
            if op.name == "TICK":
                yield cur_moment
                cur_moment = []
            else:
                cur_moment.extend(
                    _split_targets_if_needed(
                        op, immune_qubit_indices=immune_qubit_indices
                    )
                )
    if cur_moment:
        yield cur_moment


def _measure_basis(*, split_op: stim.CircuitInstruction) -> str | None:
    result = OP_TO_MEASURE_BASES.get(split_op.name)
    if result == "*":
        result = ""
        targets = split_op.targets_copy()
        for k in range(0, len(targets), 2):
            t = targets[k]
            if t.is_x_target:
                result += "X"
            elif t.is_y_target:
                result += "Y"
            elif t.is_z_target:
                result += "Z"
            else:
                raise NotImplementedError(f"{targets=}")
    return result
