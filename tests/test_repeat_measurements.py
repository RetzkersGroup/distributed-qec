"""Tests for src/repeat_measurements/."""

import pathlib

import pytest
import stim

from repeat_measurements import repeat_measurements


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_instruction(circuit: stim.Circuit, name: str) -> bool:
    return any(op.name == name for op in circuit.flattened())


def _instructions_with_name(circuit: stim.Circuit, name: str):
    return [op for op in circuit.flattened() if op.name == name]


def _count_measurements(circuit: stim.Circuit) -> int:
    return circuit.num_measurements


# ---------------------------------------------------------------------------
# 1. Identity: r=1 returns the cleaned (noiseless, flattened) circuit
# ---------------------------------------------------------------------------

class TestIdentity:
    def test_r1_is_identity(self):
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 1)
        assert result.num_measurements == src.num_measurements
        assert result.num_detectors == src.num_detectors
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_r1_strips_noise(self):
        src = stim.Circuit("""
            R 0 1
            DEPOLARIZE1(0.1) 0 1
            TICK
            MPP(0.01) Z0*Z1
            DETECTOR rec[-1]
        """)
        result = repeat_measurements(src, 1)
        assert not _has_instruction(result, "DEPOLARIZE1")

    def test_r1_flattens_repeat(self):
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            REPEAT 3 {
                TICK
                MPP Z0*Z1
                DETECTOR rec[-1] rec[-2]
            }
        """)
        result = repeat_measurements(src, 1)
        # No REPEAT blocks in output.
        assert not any(
            isinstance(op, stim.CircuitRepeatBlock) for op in result
        )
        assert result.num_detectors == 4

    def test_r0_raises(self):
        src = stim.Circuit("R 0 1\nTICK\nMPP Z0*Z1\nDETECTOR rec[-1]")
        with pytest.raises(ValueError, match="r must be >= 1"):
            repeat_measurements(src, 0)


# ---------------------------------------------------------------------------
# 2. Measurement count grows correctly
# ---------------------------------------------------------------------------

class TestMeasurementCount:
    def test_r2_doubles_mid_circuit_measurements(self):
        # 2 mid-circuit measurement layers (each with 1 meas) + 1 final.
        # Mid-circuit: 2 layers * 1 meas * r=2 = 4, final: 1. Total = 5.
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        # 2 mid-circuit layers * 2 copies + 1 final = 5.
        assert result.num_measurements == 5

    def test_r3_with_multi_measurement_layer(self):
        # Layer with 2 measurements (MXX on 4 qubits = 2 pairs).
        src = stim.Circuit("""
            RX 0 1 2 3
            TICK
            MXX 0 1 2 3
            DETECTOR rec[-1]
            DETECTOR rec[-2]
            TICK
            MXX 0 1 2 3
            DETECTOR rec[-1] rec[-3]
            DETECTOR rec[-2] rec[-4]
        """)
        result = repeat_measurements(src, 3)
        # First layer: 2 meas * 3 copies = 6. Final (not repeated): 2. Total = 8.
        assert result.num_measurements == 8


# ---------------------------------------------------------------------------
# 3. Comparison detectors are inserted
# ---------------------------------------------------------------------------

class TestComparisonDetectors:
    def test_r2_adds_comparison_detectors(self):
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        # Original: 2 detectors. r=2 adds 1 comparison per repeated layer.
        # 1 mid-circuit layer with 1 meas -> 1 comparison. Total = 2 + 1 = 3.
        assert result.num_detectors == 3

    def test_r3_adds_two_comparison_per_meas(self):
        src = stim.Circuit("""
            RX 0 1
            TICK
            MXX 0 1
            DETECTOR rec[-1]
            TICK
            MXX 0 1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 3)
        # 1 mid-circuit layer, 1 meas, r=3 -> 2 comparisons. Total = 2 + 2 = 4.
        assert result.num_detectors == 4


# ---------------------------------------------------------------------------
# 4. Deterministic detectors (noiseless -> DEM has 0 errors)
# ---------------------------------------------------------------------------

class TestDeterministic:
    def test_r2_mpp_zz(self):
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_r2_mxx(self):
        src = stim.Circuit("""
            RX 0 1
            TICK
            MXX 0 1
            DETECTOR rec[-1]
            TICK
            MXX 0 1
            DETECTOR rec[-1] rec[-2]
            TICK
            MXX 0 1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_r3_multi_round(self):
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 3)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_r2_mixed_one_and_two_body_mpp(self):
        # Like the honeycomb: mixed 1-body and 2-body products in one MPP.
        src = stim.Circuit("""
            R 0 1 2
            TICK
            MPP Z0 Z1*Z2
            DETECTOR rec[-1]
            DETECTOR rec[-2]
            TICK
            MPP Z0 Z1*Z2
            DETECTOR rec[-1] rec[-3]
            DETECTOR rec[-2] rec[-4]
            TICK
            MPP Z0 Z1*Z2
            DETECTOR rec[-1] rec[-3]
            DETECTOR rec[-2] rec[-4]
        """)
        result = repeat_measurements(src, 2)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_r2_with_repeat_block(self):
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            REPEAT 4 {
                TICK
                MPP Z0*Z1
                DETECTOR rec[-1] rec[-2]
            }
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0


# ---------------------------------------------------------------------------
# 5. Final readout is NOT repeated
# ---------------------------------------------------------------------------

class TestFinalReadout:
    def test_final_layer_single(self):
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 3)
        # 1 mid-circuit layer * 3 + 1 final = 4.
        assert result.num_measurements == 4

    def test_single_layer_circuit_not_repeated(self):
        # Only one measurement layer -> it's the final -> nothing repeated.
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
        """)
        result = repeat_measurements(src, 3)
        assert result.num_measurements == 1
        assert result.num_detectors == 1


# ---------------------------------------------------------------------------
# 6. OBSERVABLE_INCLUDE rewriting
# ---------------------------------------------------------------------------

class TestObservableInclude:
    def test_observable_rewritten_deterministic(self):
        src = stim.Circuit("""
            RX 0 1
            TICK
            MXX 0 1
            OBSERVABLE_INCLUDE(0) rec[-1]
            TICK
            MXX 0 1
            OBSERVABLE_INCLUDE(0) rec[-1]
            DETECTOR rec[-1] rec[-2]
            TICK
            MXX 0 1
            OBSERVABLE_INCLUDE(0) rec[-1]
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0
        obs_ops = _instructions_with_name(result, "OBSERVABLE_INCLUDE")
        assert len(obs_ops) == 3


# ---------------------------------------------------------------------------
# 7. split_detectors=False mode
# ---------------------------------------------------------------------------

class TestNoSplit:
    def test_no_split_still_deterministic(self):
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2, split_detectors=False)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_no_split_no_assertion_error(self):
        # A circuit where the same measurement is used by 3 detectors
        # (non-matchable) should not raise when split_detectors=False.
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        # Should not raise.
        result = repeat_measurements(src, 2, split_detectors=False)
        assert result.num_measurements > 0


# ---------------------------------------------------------------------------
# 8. Assertions when split_detectors=True
# ---------------------------------------------------------------------------

class TestAssertions:
    def test_multiple_past_facing_raises(self):
        # Two detectors in the same step referencing the same measurement.
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        with pytest.raises(ValueError, match="past-facing"):
            repeat_measurements(src, 2, split_detectors=True)

    def test_multiple_future_facing_allowed(self):
        # Multiple future-facing uses of the same measurement are allowed —
        # they all land on the last copy, consistent with the original circuit
        # already having that structure (common in Floquet codes).
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            DETECTOR rec[-2]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2, split_detectors=True)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0


# ---------------------------------------------------------------------------
# 9. TICK tags preserved
# ---------------------------------------------------------------------------

class TestTickTags:
    def test_tags_round_trip(self):
        src = stim.Circuit("""
            R 0 1
            TICK[start_noise]
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK[stop_noise]
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        tags = [op.tag for op in result if op.name == "TICK"]
        assert "start_noise" in tags
        assert "stop_noise" in tags

    def test_inserted_ticks_untagged(self):
        src = stim.Circuit("""
            R 0 1
            TICK[start_noise]
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK[stop_noise]
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 3)
        tick_tags = [op.tag for op in result if op.name == "TICK"]
        # Only the original start/stop tags; inserted TICKs have empty tag.
        tagged = [t for t in tick_tags if t]
        assert tagged == ["start_noise", "stop_noise"]

    def test_outside_noise_window_not_repeated(self):
        # Measurement before start_noise and after stop_noise should not be
        # repeated (only mid-circuit AND inside noise window).
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR rec[-1]
            TICK[start_noise]
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            TICK[stop_noise]
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        # Layer 1 (before start_noise): not repeated (outside window).
        # Layer 2 (inside window): repeated -> 2 copies.
        # Layer 3 (after stop_noise): final readout -> not repeated.
        # Total: 1 + 2 + 1 = 4.
        assert result.num_measurements == 4


# ---------------------------------------------------------------------------
# 10. SHIFT_COORDS — consumed by flattening (stim behaviour)
# ---------------------------------------------------------------------------

class TestShiftCoords:
    def test_shift_coords_consumed_by_flatten(self):
        # SHIFT_COORDS is consumed during flattening — coordinates are baked
        # into detector args. The transform still produces correct detectors.
        src = stim.Circuit("""
            R 0 1
            TICK
            MPP Z0*Z1
            DETECTOR(0, 0, 0) rec[-1]
            SHIFT_COORDS(0, 0, 1)
            TICK
            MPP Z0*Z1
            DETECTOR(0, 0, 0) rec[-1] rec[-2]
        """)
        result = repeat_measurements(src, 2)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0
        # The second detector should have absorbed the shift (coords = 0,0,1).
        dets = _instructions_with_name(result, "DETECTOR")
        shifted = [d for d in dets if d.gate_args_copy() == [0, 0, 1]]
        assert len(shifted) > 0


# ---------------------------------------------------------------------------
# 11. Real circuit regression
# ---------------------------------------------------------------------------

class TestRealCircuit:
    CIRCUIT_PATH = (
        pathlib.Path(__file__).parent.parent
        / "in/circuits/honeycomb_gidney"
        / "model=SDEM3,r=10,w=6,h=9,p=0.001,obs=H,sheared=False.stim"
    )

    def test_r2_deterministic(self):
        src = stim.Circuit.from_file(str(self.CIRCUIT_PATH))
        result = repeat_measurements(src, 2)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_r2_measurement_count(self):
        src = stim.Circuit.from_file(str(self.CIRCUIT_PATH))
        result = repeat_measurements(src, 2)
        # Original: 984 meas. 930 mid-circuit meas are doubled, final 54 stay.
        # 930 * 2 + 54 = 1914.
        assert result.num_measurements == 1914

    def test_r3_deterministic(self):
        src = stim.Circuit.from_file(str(self.CIRCUIT_PATH))
        result = repeat_measurements(src, 3)
        dem = result.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_pipeline_repeat_then_translate(self):
        """Full pipeline: repeat_measurements -> translate_circuit."""
        from bell_translation import translate_circuit

        src = stim.Circuit.from_file(str(self.CIRCUIT_PATH))
        repeated = repeat_measurements(src, 2)
        translated = translate_circuit(repeated)
        dem = translated.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0
