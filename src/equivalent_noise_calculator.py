"""Majority-vote equivalent noise block for `biased_bell_v3_equivalent_majority`.

The majority-vote model repeats the Bell-mediated syndrome measurement three
times and majority-votes the three outcomes. This module computes a *single-shot*
correlated `E(...)` noise block that is equivalent (to first order) to that
three-round gadget, so it can be dropped into the single-ancilla translation in
`bell_translation.py`.

Method (ported from notebooks/majority_vote.ipynb):

1. Build a small reference circuit: a data Bell pair (qubits 0,1 paired with
   10,11), then the three noisy Bell-mediated measurements, each tagged with a
   syndrome detector D0/D1/D2, followed by four logical-Pauli readout detectors
   D3..D6 (Z0, Z1, X0, X1). Take its detector error model.
2. Drop gauge mechanisms (p == 0.5 lines touching only the logical detectors).
3. Enumerate every single mechanism and every *pair* of mechanisms (first-order
   expansion: pair probability = p1 * p2), XOR-ing their detector sets.
4. Collapse each detector set to a flip pattern over {M, X0, X1, Z0, Z1}:
     - M (a measurement error) fires iff >= 2 of D0,D1,D2 flipped (majority vote),
     - X0/X1/Z0/Z1 are the logical detectors, carried through directly.
5. Aggregate probabilities per flip pattern and return the list, dropping the
   identity (no-flip) pattern.

The returned representation, `list[(probability, frozenset[label])]`, matches
`bell_translation._v3_equivalent_mechanisms`, so the same renderer emits the
`E(...)` lines.
"""

from __future__ import annotations

import collections
from itertools import combinations

import stim

# Detector layout in the reference circuit:
#   D0, D1, D2 — three repeated syndrome measurements (majority-voted into M)
#   D3 -> Z0,  D4 -> Z1,  D5 -> X0,  D6 -> X1
_SYNDROME_DETS = frozenset({0, 1, 2})
_LOGICAL_MAP = {3: "Z0", 4: "Z1", 5: "X0", 6: "X1"}

# How many times the syndrome measurement is repeated before majority voting.
_NUM_REPEATS = 3


def _reference_circuit(p_e: float, p_bell_z: float) -> stim.Circuit:
    """The 3-repeat majority-vote reference circuit (MXX gadget, v3 channels).

    Mirrors notebooks/majority_vote.ipynb cell 4. The three repeated blocks each
    carry the biased_bell_v3 noise channels (Z_ERROR(p_e) on data before and
    after the ancilla->data CX, PAULI_CHANNEL_1(p_e/2, p_e/2, p_bell_z) on the
    Bell ancilla, and MX(p_e) readout). The trailing MXX/MZZ readouts define the
    logical-Pauli detectors D3..D6.
    """
    blocks = []
    # Three repeated syndrome blocks on fresh ancilla pairs (2,3), (4,5), (6,7).
    for i, (anc_x, anc_z) in enumerate([(2, 3), (4, 5), (6, 7)]):
        blocks.append(f"""
Z_ERROR({p_e}) 0 1
RX {anc_x}
R {anc_z}
CX {anc_x} {anc_z}
PAULI_CHANNEL_1({p_e / 2}, {p_e / 2}, {p_bell_z}) {anc_z}
Z_ERROR({p_e}) 0 1
CX {anc_x} 0 {anc_z} 1
MX({p_e}) {anc_x} {anc_z}
DETECTOR({i}) rec[-1] rec[-2]""")

    return stim.Circuit(
        """
# initialize data qubits and pairs
RX 0 1
R 10 11
CX 0 10 1 11
MXX 0 1
"""
        + "".join(blocks)
        + """
# MEASURE PAULI ON THE DATA QUBITS
MXX 0 10
DETECTOR(1, 1) rec[-1]
MXX 1 11
DETECTOR(2, 1) rec[-1]
MZZ 0 10
DETECTOR(3, 1) rec[-1]
MZZ 1 11
DETECTOR(4, 1) rec[-1]
"""
    )


def _raw_mechanisms(circuit: stim.Circuit) -> list[tuple[float, frozenset]]:
    """DEM error mechanisms as (probability, detector-set), dropping gauges."""
    dem = circuit.detector_error_model(allow_gauge_detectors=True)
    raw: list[tuple[float, frozenset]] = []
    for instr in dem:
        if instr.type != "error":
            continue
        p = instr.args_copy()[0]
        if p == 0.5:  # gauge line
            continue
        dets = frozenset(
            t.val for t in instr.targets_copy() if t.is_relative_detector_id()
        )
        raw.append((p, dets))
    return raw


def _to_equivalent(dets: frozenset) -> frozenset:
    """Collapse a detector set to a flip pattern over {M, X0, X1, Z0, Z1}.

    M fires iff a majority (>= 2 of 3) of the syndrome detectors flipped; the
    logical detectors map through directly.
    """
    flips = set()
    if len(dets & _SYNDROME_DETS) >= (_NUM_REPEATS // 2 + 1):
        flips.add("M")
    for d in dets - _SYNDROME_DETS:
        flips.add(_LOGICAL_MAP[d])
    return frozenset(flips)


def get_equivalent_noise_operations(
    p_e: float, p_bell_z: float
) -> list[tuple[float, frozenset]]:
    """Majority-vote equivalent noise block for biased_bell_v3_equivalent_majority.

    Returns ``list[(probability, frozenset[label])]`` with labels drawn from
    {"M", "X0", "X1", "Z0", "Z1"}, aggregated per flip pattern and with the
    identity (no-flip) pattern dropped. Same representation as
    ``bell_translation._v3_equivalent_mechanisms``.
    """
    raw = _raw_mechanisms(_reference_circuit(p_e, p_bell_z))

    # First-order expansion: every single mechanism + every pair (p1 * p2).
    agg: dict[frozenset, float] = collections.defaultdict(float)
    for p, dets in raw:
        agg[_to_equivalent(dets)] += p
    for (p1, d1), (p2, d2) in combinations(raw, 2):
        agg[_to_equivalent(d1 ^ d2)] += p1 * p2

    # Drop the identity pattern; sort by descending probability for readability.
    return [
        (p, flips)
        for flips, p in sorted(agg.items(), key=lambda kv: -kv[1])
        if flips
    ]
