"""Core implementation of the repeat-measurements transform.

Given a stim circuit and a repetition count r, repeats each mid-circuit
measurement layer r times, inserts comparison detectors between consecutive
copies, and rewrites all DETECTOR / OBSERVABLE_INCLUDE rec references so the
circuit retains correct detector semantics and (when split_detectors=True)
matchability.

The input circuit is expected to be noiseless (pipeline convention: strip
noise first, re-add via the noise model afterward).
"""

from __future__ import annotations

import stim


# Measurement gate names we recognize and repeat.
_MEASUREMENT_NAMES = {
    "M", "MR", "MX", "MRX", "MY", "MRY",
    "MXX", "MYY", "MZZ", "MPP",
}

# Annotations that carry rec targets and must be rewritten.
_ANNOTATION_NAMES = {"DETECTOR", "OBSERVABLE_INCLUDE"}

# Instructions that are structural / pass-through.
_PASS_THROUGH_NAMES = {
    "TICK", "QUBIT_COORDS", "SHIFT_COORDS",
}


def repeat_measurements(
    circuit: stim.Circuit,
    r: int,
    *,
    split_detectors: bool = True,
) -> stim.Circuit:
    """Repeat each mid-circuit measurement layer r times.

    Pre-processing: strips noise via `circuit.without_noise()` and flattens
    REPEAT blocks (so measurement-index bookkeeping is simple).

    Parameters
    ----------
    circuit : stim.Circuit
        Input circuit (may contain noise and REPEAT blocks).
    r : int
        Repetition count (r=1 is identity).
    split_detectors : bool
        When True (default), remap past-facing rec targets to the first copy
        and future-facing to the last copy, preserving matchability. Raises
        ValueError if the circuit structure violates the required assumptions.
        When False, all references go to the last copy.

    Returns
    -------
    stim.Circuit
        The transformed circuit with repeated measurements and comparison
        detectors inserted.
    """
    if r < 1:
        raise ValueError(f"r must be >= 1; got {r}.")

    cleaned = circuit.without_noise().flattened()

    if r == 1:
        return cleaned

    # Parse the flat circuit into time steps (TICK-separated blocks).
    time_steps = _parse_time_steps(cleaned)

    # Identify which time steps contain measurements, and which is the final
    # measurement layer (excluded from repetition — "mid-circuit only").
    meas_step_indices = [
        i for i, ts in enumerate(time_steps)
        if any(_is_measurement(op) for op in ts)
    ]
    if not meas_step_indices:
        return cleaned

    last_meas_step = meas_step_indices[-1]

    # Determine which steps are inside the noise window (respect TICK tags).
    noise_window = _find_noise_window(time_steps)

    # Build the measurement map and output buffer.
    return _build_repeated_circuit(
        time_steps, r, split_detectors, last_meas_step, noise_window
    )


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_time_steps(circuit: stim.Circuit) -> list[list[stim.CircuitInstruction]]:
    """Split a flat circuit into time steps separated by TICK instructions.

    Each time step is a list of instructions between consecutive TICKs.
    TICK instructions themselves are stored as a single-element list forming
    their own "step" — this preserves TICK tags and ordering.
    """
    steps: list[list[stim.CircuitInstruction]] = [[]]
    for op in circuit:
        if op.name == "TICK":
            # Store the TICK as its own step so we can preserve its tag.
            steps.append([op])
            steps.append([])
        else:
            steps[-1].append(op)
    return steps


def _find_noise_window(
    time_steps: list[list[stim.CircuitInstruction]],
) -> tuple[int, int] | None:
    """Find the step indices of TICK[start_noise] and TICK[stop_noise].

    Returns (start_idx, stop_idx) or None if no tags are present.
    """
    start = stop = None
    for i, ts in enumerate(time_steps):
        if len(ts) == 1 and ts[0].name == "TICK":
            if ts[0].tag == "start_noise":
                start = i
            elif ts[0].tag == "stop_noise":
                stop = i
    if start is not None and stop is not None:
        return (start, stop)
    return None


def _is_measurement(op: stim.CircuitInstruction) -> bool:
    return op.name in _MEASUREMENT_NAMES


def _count_measurements_in_op(op: stim.CircuitInstruction) -> int:
    """Count how many measurement records a single instruction produces."""
    name = op.name
    if name == "MPP":
        # Each "product" in an MPP produces one record. Products are separated
        # by non-combiner boundaries.
        targets = op.targets_copy()
        count = 0
        prev_was_pauli = False
        for t in targets:
            if t.is_combiner:
                prev_was_pauli = False
                continue
            if prev_was_pauli:
                count += 1
            prev_was_pauli = True
        count += 1  # last product
        return count
    elif name in {"MXX", "MYY", "MZZ"}:
        return len(op.targets_copy()) // 2
    else:
        # M, MR, MX, MRX, MY, MRY — one record per qubit target.
        return len(op.targets_copy())


# ---------------------------------------------------------------------------
# Core transform
# ---------------------------------------------------------------------------

def _build_repeated_circuit(
    time_steps: list[list[stim.CircuitInstruction]],
    r: int,
    split_detectors: bool,
    last_meas_step: int,
    noise_window: tuple[int, int] | None,
) -> stim.Circuit:
    """Build the output circuit with repeated measurements."""
    # measurement_map[old_idx] = list of new measurement indices for that
    # original measurement. For non-repeated measurements: [new_idx].
    # For repeated measurements: [copy_0_idx, copy_1_idx, ..., copy_{r-1}_idx].
    measurement_map: list[list[int]] = []
    new_meas_count = 0

    # Track which measurements belong to the "current layer" for each step
    # (used for past/future classification).
    # layer_of[new_meas_idx] = step_index where it was emitted.
    layer_of: list[int] = []

    # For the usage-pattern assertion: count past-facing uses per original
    # measurement index. Only populated when split_detectors=True.
    past_uses: dict[int, int] = {}

    buf: list[str] = []
    most_recent_meas_step: int | None = None

    for step_idx, ts in enumerate(time_steps):
        # Is this step a TICK?
        if len(ts) == 1 and ts[0].name == "TICK":
            buf.append(str(ts[0]))
            continue

        # Should measurements in this step be repeated?
        should_repeat = (
            any(_is_measurement(op) for op in ts)
            and step_idx != last_meas_step
            and _in_repeat_window(step_idx, noise_window)
        )

        if should_repeat:
            # Separate measurement ops from annotations and other ops.
            meas_ops = [op for op in ts if _is_measurement(op)]
            annotations = [op for op in ts if op.name in _ANNOTATION_NAMES]
            other_ops = [
                op for op in ts
                if not _is_measurement(op) and op.name not in _ANNOTATION_NAMES
            ]

            # Count measurements in this layer.
            layer_meas_count = sum(
                _count_measurements_in_op(op) for op in meas_ops
            )

            # Emit other (non-measurement, non-annotation) ops once.
            for op in other_ops:
                buf.append(str(op))

            # Emit r copies of the measurement layer.
            first_copy_start = new_meas_count
            for copy_idx in range(r):
                if copy_idx > 0:
                    buf.append("TICK")
                for op in meas_ops:
                    buf.append(str(op))

            # Build measurement_map entries: each original measurement maps to
            # [first_copy_idx, ..., last_copy_idx].
            for m in range(layer_meas_count):
                copies = [
                    first_copy_start + copy * layer_meas_count + m
                    for copy in range(r)
                ]
                measurement_map.append(copies)
                for c in copies:
                    layer_of.append(step_idx)

            new_meas_count += r * layer_meas_count
            most_recent_meas_step = step_idx

            # Emit comparison detectors between consecutive copies.
            for copy_idx in range(r - 1):
                for m in range(layer_meas_count):
                    idx_a = first_copy_start + copy_idx * layer_meas_count + m
                    idx_b = first_copy_start + (copy_idx + 1) * layer_meas_count + m
                    offset_a = idx_a - new_meas_count
                    offset_b = idx_b - new_meas_count
                    buf.append(
                        f"DETECTOR rec[{offset_a}] rec[{offset_b}]"
                    )

            # Rewrite original annotations.
            for op in annotations:
                _rewrite_annotation(
                    buf, op, measurement_map, new_meas_count,
                    most_recent_meas_step, layer_of,
                    split_detectors, past_uses,
                )

        else:
            # Non-repeated step: emit everything, tracking measurements.
            for op in ts:
                if _is_measurement(op):
                    n = _count_measurements_in_op(op)
                    for _ in range(n):
                        measurement_map.append([new_meas_count])
                        layer_of.append(step_idx)
                        new_meas_count += 1
                    most_recent_meas_step = step_idx
                    buf.append(str(op))
                elif op.name in _ANNOTATION_NAMES:
                    _rewrite_annotation(
                        buf, op, measurement_map, new_meas_count,
                        most_recent_meas_step, layer_of,
                        split_detectors, past_uses,
                    )
                else:
                    buf.append(str(op))

    return stim.Circuit("\n".join(buf))


def _in_repeat_window(
    step_idx: int,
    noise_window: tuple[int, int] | None,
) -> bool:
    """Check whether a step is inside the repeat window.

    When noise window tags are present, only repeat measurements inside
    (start_noise, stop_noise). Otherwise repeat everything (except the final
    readout, handled by the caller).
    """
    if noise_window is None:
        return True
    start, stop = noise_window
    return start < step_idx < stop


def _rewrite_annotation(
    buf: list[str],
    op: stim.CircuitInstruction,
    measurement_map: list[list[int]],
    new_meas_count: int,
    most_recent_meas_step: int | None,
    layer_of: list[int],
    split_detectors: bool,
    past_uses: dict[int, int],
) -> None:
    """Rewrite a DETECTOR or OBSERVABLE_INCLUDE, choosing which copy per target."""
    name = op.name
    targets = op.targets_copy()
    args = op.gate_args_copy()
    n_old = len(measurement_map)

    rec_parts: list[str] = []
    for t in targets:
        if not t.is_measurement_record_target:
            rec_parts.append(str(t))
            continue

        k = -t.value
        old_idx = n_old - k
        if old_idx < 0:
            raise ValueError(f"rec[-{k}] precedes any measurement.")

        new_indices = measurement_map[old_idx]
        if len(new_indices) == 1:
            # Not a repeated measurement — straightforward remap.
            offset = new_indices[0] - new_meas_count
            rec_parts.append(f"rec[{offset}]")
        else:
            # Repeated measurement — pick which copy.
            meas_step = layer_of[new_indices[0]]

            if split_detectors:
                is_past = (meas_step == most_recent_meas_step)

                if is_past:
                    chosen = new_indices[0]  # first copy
                    # Assertion: at most one past-facing detector per
                    # repeated measurement. Multiple past-facing uses would
                    # concentrate too many detectors on the first copy
                    # (original past-facing + comparison detectors), breaking
                    # matchability that we could otherwise preserve.
                    if name == "DETECTOR":
                        past_uses[old_idx] = past_uses.get(old_idx, 0) + 1
                        if past_uses[old_idx] > 1:
                            raise ValueError(
                                f"Measurement at original index {old_idx} has "
                                f"multiple past-facing detector references. "
                                f"The input circuit may not be matchable. "
                                f"Try split_detectors=False."
                            )
                else:
                    chosen = new_indices[-1]  # last copy
            else:
                # No splitting: everything goes to the last copy.
                chosen = new_indices[-1]

            offset = chosen - new_meas_count
            rec_parts.append(f"rec[{offset}]")

    # Emit the rewritten annotation.
    if name == "DETECTOR":
        if args:
            coords = ", ".join(str(a) for a in args)
            buf.append(f"DETECTOR({coords}) {' '.join(rec_parts)}")
        else:
            buf.append(f"DETECTOR {' '.join(rec_parts)}")
    else:
        idx = int(args[0]) if args else 0
        buf.append(f"OBSERVABLE_INCLUDE({idx}) {' '.join(rec_parts)}")
