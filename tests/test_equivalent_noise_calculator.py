"""Tests for src/equivalent_noise_calculator.py (majority-vote noise block)."""

import pytest

from equivalent_noise_calculator import (
    get_equivalent_noise_operations,
    _reference_circuit,
    _raw_mechanisms,
    _to_equivalent,
)

_VALID_LABELS = {"M", "X0", "X1", "Z0", "Z1"}


class TestMechanismFormat:
    def test_returns_prob_label_pairs(self):
        mechs = get_equivalent_noise_operations(0.001, 0.01)
        assert mechs, "expected a non-empty mechanism list"
        for prob, labels in mechs:
            assert isinstance(prob, float) and prob > 0
            assert isinstance(labels, frozenset)
            assert labels <= _VALID_LABELS
            assert labels, "identity (no-flip) pattern must be dropped"

    def test_patterns_unique_and_sorted(self):
        mechs = get_equivalent_noise_operations(0.001, 0.01)
        patterns = [labels for _, labels in mechs]
        assert len(patterns) == len(set(patterns)), "patterns aggregated per flip set"
        probs = [p for p, _ in mechs]
        assert probs == sorted(probs, reverse=True)


class TestMajorityVote:
    def test_two_of_three_syndromes_give_M(self):
        assert "M" in _to_equivalent(frozenset({0, 1}))
        assert "M" in _to_equivalent(frozenset({0, 1, 2}))

    def test_one_syndrome_does_not_give_M(self):
        assert "M" not in _to_equivalent(frozenset({0}))

    def test_logical_detectors_map_through(self):
        assert _to_equivalent(frozenset({3})) == frozenset({"Z0"})
        assert _to_equivalent(frozenset({4})) == frozenset({"Z1"})
        assert _to_equivalent(frozenset({5})) == frozenset({"X0"})
        assert _to_equivalent(frozenset({6})) == frozenset({"X1"})

    def test_mixed_syndrome_and_logical(self):
        # Majority of syndromes + a logical flip.
        assert _to_equivalent(frozenset({0, 1, 5})) == frozenset({"M", "X0"})


class TestReferenceCircuit:
    def test_dem_builds_with_gauge(self):
        circ = _reference_circuit(0.001, 0.01)
        dem = circ.detector_error_model(allow_gauge_detectors=True)
        assert dem.num_detectors == 7  # D0..D6

    def test_raw_drops_gauges(self):
        circ = _reference_circuit(0.001, 0.01)
        raw = _raw_mechanisms(circ)
        assert raw, "expected non-gauge mechanisms"
        assert all(p != 0.5 for p, _ in raw)


class TestScaling:
    def test_zero_noise_gives_empty(self):
        assert get_equivalent_noise_operations(0.0, 0.0) == []

    def test_probabilities_increase_with_p(self):
        lo = dict((labels, p) for p, labels in get_equivalent_noise_operations(0.001, 0.01))
        hi = dict((labels, p) for p, labels in get_equivalent_noise_operations(0.01, 0.01))
        # The pure-measurement pattern should grow with p_e.
        assert hi[frozenset({"M"})] > lo[frozenset({"M"})]

    def test_pbell_z_feeds_measurement(self):
        # Raising p_bell_z (a Z error on the ancilla -> measurement flip) should
        # increase the M-flip weight.
        base = dict((labels, p) for p, labels in get_equivalent_noise_operations(0.001, 0.001))
        more = dict((labels, p) for p, labels in get_equivalent_noise_operations(0.001, 0.05))
        assert more[frozenset({"M"})] > base[frozenset({"M"})]
