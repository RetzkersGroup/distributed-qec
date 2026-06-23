"""Tests for src/bell_translation.py."""

import pathlib

import pytest
import stim

from bell_translation import NoiseParams, translate_circuit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_instruction(circuit: stim.Circuit, name: str) -> bool:
    return any(op.name == name for op in circuit.flattened())


def _instructions_with_name(circuit: stim.Circuit, name: str):
    return [op for op in circuit.flattened() if op.name == name]


def _all_cx_targets(circuit: stim.Circuit) -> list[tuple[int, int]]:
    pairs = []
    for op in circuit.flattened():
        if op.name == "CX":
            targets = op.targets_copy()
            for i in range(0, len(targets), 2):
                pairs.append((targets[i].qubit_value, targets[i + 1].qubit_value))
    return pairs


# ---------------------------------------------------------------------------
# 1. Noiseless translation — deterministic detectors
# ---------------------------------------------------------------------------

class TestNoiselessDeterministic:
    def test_xx_mpp(self):
        src = stim.Circuit("""
            RX 0 1
            MPP X0*X1
            MPP X0*X1
            DETECTOR rec[-1] rec[-2]
        """)
        tr = translate_circuit(src)
        dem = tr.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_zz_mpp(self):
        src = stim.Circuit("""
            R 0 1
            MPP Z0*Z1
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        tr = translate_circuit(src)
        dem = tr.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_yy_mpp(self):
        src = stim.Circuit("""
            R 0 1
            H 0 1
            S 0 1
            MPP Y0*Y1
            MPP Y0*Y1
            DETECTOR rec[-1] rec[-2]
        """)
        tr = translate_circuit(src)
        dem = tr.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_mxx_instruction(self):
        src = stim.Circuit("""
            RX 0 1
            MXX 0 1
            MXX 0 1
            DETECTOR rec[-1] rec[-2]
        """)
        tr = translate_circuit(src)
        dem = tr.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_mzz_instruction(self):
        src = stim.Circuit("""
            R 0 1
            MZZ 0 1
            MZZ 0 1
            DETECTOR rec[-1] rec[-2]
        """)
        tr = translate_circuit(src)
        dem = tr.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_myy_instruction(self):
        src = stim.Circuit("""
            R 0 1
            H 0 1
            S 0 1
            MYY 0 1
            MYY 0 1
            DETECTOR rec[-1] rec[-2]
        """)
        tr = translate_circuit(src)
        dem = tr.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0


# ---------------------------------------------------------------------------
# 2. Ancilla allocation — fixed mapping
# ---------------------------------------------------------------------------

class TestAncillaAllocation:
    def test_num_qubits_doubled(self):
        src = stim.Circuit("""
            R 0 1 2 3
            MPP Z0*Z1 Z2*Z3
            DETECTOR rec[-1]
            DETECTOR rec[-2]
        """)
        tr = translate_circuit(src)
        assert tr.num_qubits == 8  # 4 data + 4 ancilla

    def test_ancilla_indices_correct(self):
        # With 4 data qubits, ancilla for q0 = 4, q1 = 5, q2 = 6, q3 = 7.
        src = stim.Circuit("""
            R 0 1 2 3
            MPP Z0*Z2
            DETECTOR rec[-1]
        """)
        tr = translate_circuit(src)
        cx_pairs = _all_cx_targets(tr)
        # Bell-pair CX: anc(0)=4 -> anc(2)=6
        assert (4, 6) in cx_pairs
        # Data-ancilla CX: data 0 -> anc(0)=4 (ZZ: data controls)
        assert (0, 4) in cx_pairs
        # Data-ancilla CX: data 2 -> anc(2)=6
        assert (2, 6) in cx_pairs


# ---------------------------------------------------------------------------
# 3. Ancilla reuse across rounds
# ---------------------------------------------------------------------------

class TestAncillaReuse:
    def test_same_ancillas_in_both_rounds(self):
        src = stim.Circuit("""
            RX 0 1
            MPP X0*X1
            MPP X0*X1
            DETECTOR rec[-1] rec[-2]
        """)
        tr = translate_circuit(src)
        # Only 4 qubits total (2 data + 2 ancilla), not 6.
        assert tr.num_qubits == 4

    def test_different_pairs_share_ancillas(self):
        # MPP X0*X1 uses ancillas (2, 3); MPP X0*X2 uses (2, 4).
        # Qubit 0's ancilla is always 3 (num_data=3, anc(0)=3).
        src = stim.Circuit("""
            RX 0 1 2
            MPP X0*X1
            MPP X0*X2
            DETECTOR rec[-1]
            DETECTOR rec[-2]
        """)
        tr = translate_circuit(src)
        assert tr.num_qubits == 6  # 3 data + 3 ancilla
        cx_pairs = _all_cx_targets(tr)
        # Ancilla for q0=3 appears in both rounds.
        q0_cx = [(c, t) for c, t in cx_pairs if c == 3 or t == 3]
        assert len(q0_cx) >= 2  # At least one CX per round involving anc(0)


# ---------------------------------------------------------------------------
# 4. Rec-ref rewriting
# ---------------------------------------------------------------------------

class TestRecRefRewriting:
    def test_mix_of_mpp_and_m(self):
        # MPP produces 2 recs, M produces 1 rec.
        src = stim.Circuit("""
            R 0 1
            MPP Z0*Z1
            M 0
            DETECTOR rec[-1] rec[-2]
        """)
        tr = translate_circuit(src)
        # Should still have deterministic detectors:
        # rec[-2] in source -> the MPP result (now 2 recs)
        # rec[-1] in source -> the M result (still 1 rec)
        # Stim would complain if the addressing is wrong.
        dem = tr.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0

    def test_observable_include_rewritten(self):
        src = stim.Circuit("""
            RX 0 1
            MPP X0*X1
            OBSERVABLE_INCLUDE(0) rec[-1]
        """)
        tr = translate_circuit(src)
        obs_ops = _instructions_with_name(tr, "OBSERVABLE_INCLUDE")
        assert len(obs_ops) == 1
        targets = obs_ops[0].targets_copy()
        # Original rec[-1] should become 2 targets.
        assert len(targets) == 2
        assert all(t.is_measurement_record_target for t in targets)


# ---------------------------------------------------------------------------
# 5. Noise injection — each channel appears
# ---------------------------------------------------------------------------

class TestNoiseInjection:
    def test_p_bell_xyz(self):
        src = stim.Circuit("R 0 1\nMPP Z0*Z1\nDETECTOR rec[-1]")
        noise = NoiseParams(p_bell_x=0.01, p_bell_y=0.02, p_bell_z=0.03)
        tr = translate_circuit(src, noise)
        assert _has_instruction(tr, "PAULI_CHANNEL_1")

    def test_p_internal(self):
        src = stim.Circuit("R 0 1\nMPP Z0*Z1\nDETECTOR rec[-1]")
        noise = NoiseParams(p_internal=0.01)
        tr = translate_circuit(src, noise)
        assert _has_instruction(tr, "DEPOLARIZE2")

    def test_p_data_z_on_data_qubits(self):
        # v2 channel: Z_ERROR(p_data_z) applied to the two participating
        # data qubits during Bell-pair prep.
        src = stim.Circuit("R 0 1\nMPP Z0*Z1\nDETECTOR rec[-1]")
        noise = NoiseParams(p_data_z=0.04)
        tr = translate_circuit(src, noise)
        z_err_ops = _instructions_with_name(tr, "Z_ERROR")
        # Exactly one Z_ERROR(0.04) on data qubits {0, 1}.
        matching = [op for op in z_err_ops if op.gate_args_copy() == [0.04]]
        assert len(matching) == 1
        qubits = {t.qubit_value for t in matching[0].targets_copy()}
        assert qubits == {0, 1}

    def test_p_m(self):
        src = stim.Circuit("R 0 1\nMPP Z0*Z1\nM 0\nDETECTOR rec[-1]")
        noise = NoiseParams(p_m=0.01)
        tr = translate_circuit(src, noise)
        # Both Bell-block measurements and single-qubit M should have (0.01).
        m_ops = _instructions_with_name(tr, "M")
        for op in m_ops:
            assert op.gate_args_copy() == [0.01]

    def test_p_1_after_h(self):
        src = stim.Circuit("R 0 1\nH 0\nMPP Z0*Z1\nDETECTOR rec[-1]")
        noise = NoiseParams(p_1=0.005)
        tr = translate_circuit(src, noise)
        dep1_ops = _instructions_with_name(tr, "DEPOLARIZE1")
        assert len(dep1_ops) > 0

    def test_p_init_data_qubit(self):
        src = stim.Circuit("R 0 1\nRX 2 3\nMPP Z0*Z1\nDETECTOR rec[-1]")
        noise = NoiseParams(p_init=0.01)
        tr = translate_circuit(src, noise)
        assert _has_instruction(tr, "X_ERROR")
        assert _has_instruction(tr, "Z_ERROR")


# ---------------------------------------------------------------------------
# 6. p_init does NOT affect ancilla resets
# ---------------------------------------------------------------------------

class TestPInitAncillaIsolation:
    def test_ancilla_resets_no_x_error(self):
        src = stim.Circuit("R 0 1\nMPP Z0*Z1\nDETECTOR rec[-1]")
        noise = NoiseParams(p_init=0.05)
        tr = translate_circuit(src, noise)
        # X_ERROR should only appear after the data-qubit R (on qubits 0, 1),
        # NOT after the ancilla RX/R inside the Bell block.
        x_err_ops = _instructions_with_name(tr, "X_ERROR")
        for op in x_err_ops:
            for t in op.targets_copy():
                assert t.qubit_value < 2  # Only data qubits 0, 1

    def test_ancilla_resets_no_z_error(self):
        src = stim.Circuit("RX 0 1\nMPP X0*X1\nDETECTOR rec[-1]")
        noise = NoiseParams(p_init=0.05)
        tr = translate_circuit(src, noise)
        z_err_ops = _instructions_with_name(tr, "Z_ERROR")
        for op in z_err_ops:
            for t in op.targets_copy():
                assert t.qubit_value < 2  # Only data qubits 0, 1


# ---------------------------------------------------------------------------
# 7. Noise stripping
# ---------------------------------------------------------------------------

class TestNoiseStripping:
    def test_input_noise_removed(self):
        src = stim.Circuit("""
            R 0 1
            DEPOLARIZE1(0.1) 0 1
            H 0
            DEPOLARIZE1(0.1) 0
            MPP Z0*Z1
            DETECTOR rec[-1]
        """)
        tr = translate_circuit(src)  # noiseless translation
        # Original DEPOLARIZE1 should be gone (noiseless output).
        assert not _has_instruction(tr, "DEPOLARIZE1")
        assert not _has_instruction(tr, "DEPOLARIZE2")
        assert not _has_instruction(tr, "PAULI_CHANNEL_1")
        assert not _has_instruction(tr, "X_ERROR")
        assert not _has_instruction(tr, "Z_ERROR")

    def test_only_biased_bell_noise_present(self):
        src = stim.Circuit("""
            R 0 1
            DEPOLARIZE1(0.99) 0 1
            MPP Z0*Z1
            DETECTOR rec[-1]
        """)
        noise = NoiseParams(p_internal=0.001)
        tr = translate_circuit(src, noise)
        dep2_ops = _instructions_with_name(tr, "DEPOLARIZE2")
        # All DEPOLARIZE2 should have arg 0.001 (not 0.99).
        for op in dep2_ops:
            assert op.gate_args_copy() == [0.001]


# ---------------------------------------------------------------------------
# 8. Allowlist rejection
# ---------------------------------------------------------------------------

class TestAllowlistRejection:
    def test_cx_rejected(self):
        src = stim.Circuit("R 0 1\nCX 0 1\nMPP X0*X1\nDETECTOR rec[-1]")
        with pytest.raises(ValueError, match="CX"):
            translate_circuit(src)

    def test_cz_rejected(self):
        src = stim.Circuit("R 0 1\nCZ 0 1\nMPP Z0*Z1\nDETECTOR rec[-1]")
        with pytest.raises(ValueError, match="CZ"):
            translate_circuit(src)

    def test_ry_rejected(self):
        src = stim.Circuit("RY 0\nMPP X0*X0\nDETECTOR rec[-1]")
        with pytest.raises(ValueError, match="RY"):
            translate_circuit(src)

    def test_3body_mpp_rejected(self):
        src = stim.Circuit("R 0 1 2\nMPP X0*X1*X2\nDETECTOR rec[-1]")
        with pytest.raises(ValueError, match="length 3"):
            translate_circuit(src)

    def test_mixed_basis_mpp_rejected(self):
        src = stim.Circuit("R 0 1\nMPP X0*Z1\nDETECTOR rec[-1]")
        with pytest.raises(ValueError, match="XZ"):
            translate_circuit(src)


# ---------------------------------------------------------------------------
# 9. REPEAT blocks — flattened correctly
# ---------------------------------------------------------------------------

class TestRepeatFlattening:
    def test_repeat_matches_unrolled(self):
        repeated = stim.Circuit("""
            R 0 1
            MPP Z0*Z1
            DETECTOR rec[-1]
            REPEAT 3 {
                MPP Z0*Z1
                DETECTOR rec[-1] rec[-2]
            }
        """)
        unrolled = stim.Circuit("""
            R 0 1
            MPP Z0*Z1
            DETECTOR rec[-1]
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
            MPP Z0*Z1
            DETECTOR rec[-1] rec[-2]
        """)
        tr_repeated = translate_circuit(repeated)
        tr_unrolled = translate_circuit(unrolled)
        # Both should have identical measurement counts and DEM.
        assert tr_repeated.num_measurements == tr_unrolled.num_measurements
        assert tr_repeated.num_detectors == tr_unrolled.num_detectors
        dem_r = tr_repeated.detector_error_model(decompose_errors=False)
        dem_u = tr_unrolled.detector_error_model(decompose_errors=False)
        assert dem_r.num_errors == dem_u.num_errors == 0


# ---------------------------------------------------------------------------
# 10. Real circuit regression
# ---------------------------------------------------------------------------

class TestRealCircuitRegression:
    def test_repeat_pairwise_xyz(self):
        path = pathlib.Path(__file__).parent.parent / "in/circuits/test/repeat_pairwise_xyz.stim"
        src = stim.Circuit.from_file(str(path))
        tr = translate_circuit(src)
        assert tr.num_qubits == 12
        assert tr.num_measurements == 12
        assert tr.num_detectors == 6
        dem = tr.detector_error_model(decompose_errors=False)
        assert dem.num_errors == 0
