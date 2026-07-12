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

# Detector layout in the reference circuit (num_repeats = R):
#   D0 .. D(R-1)  — the R repeated syndrome measurements (voted into M)
#   D(R) -> Z0,  D(R+1) -> Z1,  D(R+2) -> X0,  D(R+3) -> X1
# (For R=3: syndrome dets {0,1,2}; logical D3->Z0, D4->Z1, D5->X0, D6->X1.)

# Default number of repeats before majority voting. R=1 is a single measurement
# (no vote); R=3 is the majority-vote model.
_NUM_REPEATS = 3


def _syndrome_dets(num_repeats: int) -> frozenset:
    """Detector ids of the repeated syndrome measurements (D0 .. D(R-1))."""
    return frozenset(range(num_repeats))


def _logical_map(num_repeats: int) -> dict[int, str]:
    """Detector id -> logical-Pauli label for the four trailing readouts."""
    r = num_repeats
    return {r: "Z0", r + 1: "Z1", r + 2: "X0", r + 3: "X1"}


def _reference_circuit(
    p_e: float, p_bell_z: float, num_repeats: int = _NUM_REPEATS
) -> stim.Circuit:
    """The ``num_repeats``-repeat vote reference circuit (MXX gadget, v3 channels).

    Mirrors notebooks/majority_vote.ipynb cell 4, generalized to R = num_repeats
    repeats. The R repeated blocks each carry the biased_bell_v3 noise channels,
    matching the deployed two-ancilla gadget (bell_translation._xx_to_buf) one for
    one: PAULI_CHANNEL_1(p_e/2, p_e/2, p_bell_z) on the Bell ancilla, a single
    Z_ERROR(p_e) on the data before the ancilla->data CX, and MX(p_e) readout. The
    trailing MXX/MZZ readouts define the logical-Pauli detectors D(R)..D(R+3).
    R=1 gives a single, un-voted measurement.

    Each repeat uses a fresh ancilla pair (2+2i, 3+2i); with data partners at
    10,11 this stays collision-free for num_repeats <= 4 (only 1 and 3 are used).
    """
    if not 1 <= num_repeats <= 4:
        raise ValueError(f"num_repeats must be in 1..4, got {num_repeats}")
    blocks = []
    # Repeated syndrome blocks on fresh ancilla pairs (2,3), (4,5), ...
    for i in range(num_repeats):
        anc_x, anc_z = 2 + 2 * i, 3 + 2 * i
        blocks.append(f"""
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
    # approximate_disjoint_errors: the gadget's PAULI_CHANNEL_1 has disjoint
    # Pauli mechanisms (e.g. at p_bell_z=0 only X/Y remain) that stim must be
    # told it may treat as independent. It's a no-op when not actually needed,
    # so it doesn't change any p_bell_z>0 result.
    dem = circuit.detector_error_model(
        allow_gauge_detectors=True, approximate_disjoint_errors=True)
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


def _to_equivalent(dets: frozenset, num_repeats: int = _NUM_REPEATS) -> frozenset:
    """Collapse a detector set to a flip pattern over {M, X0, X1, Z0, Z1}.

    M fires iff a majority (>= R//2 + 1 of R) of the syndrome detectors flipped;
    the logical detectors map through directly. For R=1 the "majority" is 1, so M
    fires iff the single syndrome detector flipped (an un-voted measurement).
    """
    syndrome = _syndrome_dets(num_repeats)
    logical_map = _logical_map(num_repeats)
    flips = set()
    if len(dets & syndrome) >= (num_repeats // 2 + 1):
        flips.add("M")
    for d in dets - syndrome:
        flips.add(logical_map[d])
    return frozenset(flips)


def get_equivalent_noise_operations(
    p_e: float, p_bell_z: float, num_repeats: int = _NUM_REPEATS
) -> list[tuple[float, frozenset]]:
    """Equivalent noise block for the repeated Bell-mediated measurement.

    With the default ``num_repeats=3`` this is the majority-vote model used by
    ``biased_bell_v3_equivalent_majority``; ``num_repeats=1`` gives the single,
    un-voted measurement (see ``get_equivalent_noise_operations_single``).

    Returns ``list[(probability, frozenset[label])]`` with labels drawn from
    {"M", "X0", "X1", "Z0", "Z1"}, aggregated per flip pattern and with the
    identity (no-flip) pattern dropped. Same representation as
    ``bell_translation._v3_equivalent_mechanisms``.
    """
    raw = _raw_mechanisms(_reference_circuit(p_e, p_bell_z, num_repeats))

    # First-order expansion: every single mechanism + every pair (p1 * p2).
    agg: dict[frozenset, float] = collections.defaultdict(float)
    for p, dets in raw:
        agg[_to_equivalent(dets, num_repeats)] += p
    for (p1, d1), (p2, d2) in combinations(raw, 2):
        agg[_to_equivalent(d1 ^ d2, num_repeats)] += p1 * p2

    if num_repeats == 1:
        # A single, unrepeated measurement error flips only the lone syndrome
        # detector, which is a gauge in the reference DEM (locally undetectable),
        # so stim drops it (reports p=0.5) and P(M) would lose its p_e/p_bell_z
        # dependence. Restore it with the physical single-measurement flip
        # probability -- two MX(p_e) readouts plus the Z-biased ancilla error --
        # matching the pure-{M} terms of bell_translation._v3_equivalent_mechanisms.
        agg[frozenset({"M"})] += 2 * p_e + p_bell_z

    # Drop the identity pattern; sort by descending probability for readability.
    return [
        (p, flips)
        for flips, p in sorted(agg.items(), key=lambda kv: -kv[1])
        if flips
    ]


def get_equivalent_noise_operations_single(
    p_e: float, p_bell_z: float
) -> list[tuple[float, frozenset]]:
    """Single-repetition (num_repeats=1) equivalent block: no majority vote.

    Same derivation machinery as the 3-repeat majority model, but with one
    measurement, so P(M) is *linear* in the physical error rate instead of the
    ~3 p^2 quadratic majority-vote suppression. The isolated measurement error is
    a gauge in the R=1 reference DEM, so it is injected by hand (see
    ``get_equivalent_noise_operations``) to keep the block faithful; the result is
    thus comparable to the static hand-written
    ``bell_translation._v3_equivalent_mechanisms``.
    """
    return get_equivalent_noise_operations(p_e, p_bell_z, num_repeats=1)


def _two_anc_reference_circuit(p: float, p_bell_z: float) -> stim.Circuit:
    """Reference circuit for the *deployed* two-ancilla Bell MXX gadget.

    Built with the real bell_translation._xx_to_buf so it faithfully matches the
    biased_bell_v3 ("bell") scheme's physical channels (PAULI_CHANNEL_1(p/2, p/2,
    p_bell_z) on the ancilla, Z_ERROR(p) on the data, MX(p) readout; p_internal=0,
    so no two-qubit depolarizing). Same scaffolding as _reference_circuit: a
    noiseless data Bell pair (0,1)<->(10,11), the single noisy measurement (its
    two-ancilla parity is the syndrome detector D0), then the four noiseless
    logical-Pauli readouts D1->Z0, D2->Z1, D3->X0, D4->X1.
    """
    from bell_translation import NoiseParams, _xx_to_buf  # lazy: avoid import cycle

    noise = NoiseParams(
        p_bell_x=p / 2, p_bell_y=p / 2, p_bell_z=p_bell_z,
        p_data_z=p, p_internal=0.0, p_m=p, p_1=p / 10, p_init=p,
        mzz_gadget="v3", scheme="bell",
    )
    buf = ["RX 0 1", "R 10 11", "CX 0 10 1 11", "MXX 0 1"]
    _xx_to_buf(buf, 0, 1, 2, 3, noise)            # deployed two-ancilla MXX gadget
    buf.append("DETECTOR(0) rec[-1] rec[-2]")     # syndrome = parity of the pair
    buf += ["MXX 0 10", "DETECTOR(1) rec[-1]",
            "MXX 1 11", "DETECTOR(2) rec[-1]",
            "MZZ 0 10", "DETECTOR(3) rec[-1]",
            "MZZ 1 11", "DETECTOR(4) rec[-1]"]
    return stim.Circuit("\n".join(buf))


def get_two_anc_noise_operations(
    p: float, p_bell_z: float
) -> list[tuple[float, frozenset]]:
    """Effective noise block of the deployed two-ancilla Bell (biased_bell_v3) gadget.

    Same {M, X0, X1, Z0, Z1} representation and derivation (DEM collapse + pair
    expansion + injected single-measurement flip) as the 1-rep equivalent, but
    from the physical two-ancilla gadget rather than the single-ancilla one, so
    the two can be compared mechanism-by-mechanism.
    """
    raw = _raw_mechanisms(_two_anc_reference_circuit(p, p_bell_z))
    agg: dict[frozenset, float] = collections.defaultdict(float)
    for pr, dets in raw:
        agg[_to_equivalent(dets, 1)] += pr
    for (p1, d1), (p2, d2) in combinations(raw, 2):
        agg[_to_equivalent(d1 ^ d2, 1)] += p1 * p2
    # Single measurement error is a gauge here too; inject it (see the R=1 note).
    agg[frozenset({"M"})] += 2 * p + p_bell_z
    return [
        (pr, flips)
        for flips, pr in sorted(agg.items(), key=lambda kv: -kv[1])
        if flips
    ]
