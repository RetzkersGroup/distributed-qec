#!/usr/bin/env python3
"""Generate gauge honeycomb (phenomenological-noise) circuits on a (p, p_m) grid.

These are the repeated-measurement circuits of Derks et al., "Improving the
spacetime overhead of dynamic stabilizer codes by repeating measurements": the
gauge factor g repeats every check g times *with comparison detectors built in*,
so the circuit stays graphlike (plain pymatching threshold fits are meaningful),
unlike the src/repeat_measurements full-decoder route. We build them by calling
the paper's compiler (the `main` package, from OzYossarian/Kandel) once per grid
point, applying absolute phenomenological noise:

  - data qubits:   PAULI_CHANNEL_1(p/3, p/3, p/3)   (x-axis "p", the p_res analog)
  - measurements:  flip probability p_m             (y-axis, the p_z analog)

No `per`-renormalization (unlike floquet_colour_codes_numerics/.../generate_circuits.py):
p and p_m are the physical rates directly, so the two grid axes are independent —
exactly the (p_res, p_z) structure the downstream threshold/boundary stack expects.

Sidecar .json metadata is keyed so tools/decode.py + tools/_grid_scans/plot_grid.py
work UNCHANGED: p_m is stored under `p_bell_z` (the column the FSS fit sweeps),
`noise_model = peterjan_pheno_g{g}`, `expected_distance = d`, `k = 1`.

Examples:

    python tools/peterjan_tools/make_grid_circuits_peterjan.py \\
        --distances 4 8 12 --gauge-factors 1 2 3 4 5 \\
        --p-values 0.0 0.001 --pm-values 0.01 0.03 0.05 \\
        --output-dir out/circuits/peterjan_circuits/grid

    # cheap distance/DEM verification on a subset:
    python tools/peterjan_tools/make_grid_circuits_peterjan.py \\
        --distances 4 --gauge-factors 1 2 --p-values 0.001 --pm-values 0.02 \\
        --verify-distance --output-dir /tmp/peterjan_verify
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import multiprocessing as mp
import os
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

# The paper's compiler package (OzYossarian/Kandel, installed as `main`).
from main.building_blocks.detectors.Stabilizer import Stabilizer
from main.codes.tic_tac_toe.gauge.GaugeHoneycombCode import GaugeHoneycombCode
from main.compiling.compilers.AncillaPerCheckCompiler import AncillaPerCheckCompiler
from main.compiling.noise.models import PhenomenologicalNoise
from main.compiling.noise.noises import OneQubitNoise
from main.compiling.syndrome_extraction.extractors.ancilla_per_check.mixed.CxCyCzExtractor import (
    CxCyCzExtractor,
)

# PAULI_CHANNEL_1 needs px+py+pz < 1 and a measurement flip is a probability;
# values near the bound are unphysical, so skip axes above these.
_P_MAX = 0.75
_PM_MAX = 0.4

# Extended timelike-distance table (gauge factors 4, 5) produced by
# extend_timelike_distance_data.py. Kandel's shipped table only covers g<=3, so
# for g>=4 we monkeypatch GaugeHoneycombCode.get_timelike_distance to read this
# superset instead — the installed Kandel package is never mutated. Runs at
# import time so every ProcessPoolExecutor worker inherits the patch. g<=3 is
# unaffected whether or not the table exists (entries are identical to shipped).
_EXT_TABLE_PATH = pathlib.Path(__file__).resolve().parent / "data" \
    / "hcc_graphlike_td_data_phenomenological_extended.json"


def _install_extended_td_table() -> bool:
    if not _EXT_TABLE_PATH.exists():
        return False
    with _EXT_TABLE_PATH.open() as f:
        table = json.load(f)

    def _patched(self, n_rounds, pauli_letter, noise_model="phenomenological_noise"):
        # We only ever build phenomenological circuits; guard so a stray other
        # model doesn't silently read the wrong table.
        assert noise_model == "phenomenological_noise", noise_model
        return self.distance_from_timelike_distance_dict(
            n_rounds, pauli_letter, table)

    GaugeHoneycombCode.get_timelike_distance = _patched
    return True


_EXT_TABLE_LOADED = _install_extended_td_table()


def _fmt(x: float) -> str:
    return f"{x:.5g}"


def _build_settings(code: GaugeHoneycombCode, observable: str) -> dict:
    """Compiler settings for a memory experiment (mirrors the paper's generator).

    Rounds are scaled so the timelike distance matches d under phenomenological
    noise (the GaugeHoneycombCode default noise_model), the same as the paper's
    memory circuits. Only memory_x / memory_z are supported here.
    """
    rounds, d_x, d_z = code.get_number_of_rounds_for_timelike_distance(code.distance)
    assert d_x == code.distance or d_z == code.distance, (d_x, d_z, code.distance)

    if observable == "memory_x":
        observables = [code.logical_qubits[1].x]
        initial = [Stabilizer([(0, chk)], 0) for chk in code.check_schedule[0]]
    elif observable == "memory_z":
        observables = [code.logical_qubits[0].z]
        gf = code.gauge_factors
        initial = [Stabilizer([(0, chk)], 0)
                   for chk in code.check_schedule[gf[0] + gf[1]]]
    else:
        raise ValueError(f"observable must be memory_x/memory_z, got {observable!r}")

    return {"total_rounds": rounds, "initial_stabilizers": initial,
            "observables": observables}, rounds


def _compile(d: int, g: int, p: float, p_m: float, observable: str):
    """Compile one gauge-honeycomb circuit at absolute noise (p/3,p/3,p/3), p_m.

    Kandel's compiler prints an alive_progress bar to stdout; silence it so
    parallel workers don't flood the console.
    """
    code = GaugeHoneycombCode(d, [g, g, g])
    settings, rounds = _build_settings(code, observable)
    noise = PhenomenologicalNoise(OneQubitNoise(p / 3, p / 3, p / 3), p_m)
    compiler = AncillaPerCheckCompiler(noise, CxCyCzExtractor())
    with contextlib.redirect_stdout(io.StringIO()), \
            contextlib.redirect_stderr(io.StringIO()):
        circuit = compiler.compile_to_stim(code, **settings)
    return circuit, rounds


def _stem(d: int, g: int, n: int, observable: str, p: float, p_m: float) -> str:
    obs = observable.replace("memory_", "")  # x / z
    # `_d=`, `_k=`, `_n=` are the tokens plot_grid.load_grid parses as fallbacks;
    # pz=<p_m>,p=<p> mirror the bell generator so the fit stack keys line up.
    return (f"peterjan_gh_d={d}_k=1_n={n}_g{g}_obs{obs}"
            f"__peterjan_pheno_g{g},pz={_fmt(p_m)},p={_fmt(p)}")


def _process(d, g, n, p, p_m, out_dir, observable, skip_existing, verify_distance):
    if p > _P_MAX:
        return f"  d={d} g={g} p={p:.4g} pm={p_m:.4g} -> skipped (p > {_P_MAX})"
    if p_m > _PM_MAX:
        return f"  d={d} g={g} p={p:.4g} pm={p_m:.4g} -> skipped (p_m > {_PM_MAX})"

    stem = _stem(d, g, n, observable, p, p_m)
    circuit_path = out_dir / f"{stem}.stim"
    meta_path = out_dir / f"{stem}.json"
    if skip_existing and circuit_path.exists() and meta_path.exists():
        return f"  d={d} g={g} p={p:.4g} pm={p_m:.4g} -> {circuit_path.name} (skipped)"

    circuit, rounds = _compile(d, g, p, p_m, observable)

    note = ""
    if verify_distance:
        dem = circuit.detector_error_model(decompose_errors=True)
        gl = len(dem.shortest_graphlike_error())
        note = f" [graphlike={gl}{'' if gl >= d else ' < d !!'}]"

    circuit.to_file(str(circuit_path))
    metadata = {
        "code_family": "gauge_honeycomb_periodic",
        "code_name": "GaugeHoneycombCode",
        "source": "https://github.com/peter-janderks/floquet_colour_codes_numerics",
        "noise_model": f"peterjan_pheno_g{g}",
        "p": float(p),
        "p_bell_z": float(p_m),   # p_m stored here so the FSS fit sweeps it
        "p_m": float(p_m),        # honest alias
        "direction_kind": "grid",
        "direction_value": None,
        "gauge_factor": g,
        "repetitions": g,
        "distance": d,
        "expected_distance": d,
        "k": 1,
        "rounds": rounds,
        "observable": observable,
        "px": p / 3, "py": p / 3, "pz": p / 3,
        "num_qubits": circuit.num_qubits,
        "num_measurements": circuit.num_measurements,
        "num_detectors": circuit.num_detectors,
        "num_observables": circuit.num_observables,
        # plot_grid.load_grid parses n / d / k / rounds fallbacks from the
        # input_circuit filename (`_n=`, `_d=`, `_k=`), so point it at our stem.
        "input_circuit": str(circuit_path),
        "output_circuit": str(circuit_path),
    }
    with meta_path.open("w") as f:
        json.dump(metadata, f, indent=2)
    return f"  d={d} g={g} p={p:.4g} pm={p_m:.4g} -> {circuit_path.name}{note}"


def _probe_n(d: int, g: int, observable: str) -> int:
    """Physical qubit count for (d, g) — constant across p/p_m, so probe once."""
    circuit, _ = _compile(d, g, 0.0, 0.0, observable)
    return circuit.num_qubits


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--distances", nargs="+", type=int, default=[4, 8, 12])
    parser.add_argument("--gauge-factors", nargs="+", type=int,
                        default=[1, 2, 3, 4, 5], dest="gauge_factors",
                        help="Uniform gauge factor g -> [g,g,g] (g measurement reps). "
                             "g>=4 requires the extended timelike-distance table: run "
                             "tools/peterjan_tools/extend_timelike_distance_data.py once "
                             "(Kandel's shipped table only covers g<=3).")
    parser.add_argument("--p-values", nargs="+", type=float, required=True,
                        dest="p_values",
                        help="Data-noise grid: PAULI_CHANNEL_1(p/3,p/3,p/3). "
                             "0.0 allowed (pure-p_m axis).")
    parser.add_argument("--pm-values", nargs="+", type=float, required=True,
                        dest="pm_values",
                        help="Measurement-flip grid. 0.0 allowed (pure-p axis).")
    parser.add_argument("--observable", default="memory_x",
                        choices=["memory_x", "memory_z"])
    parser.add_argument("--output-dir",
                        default="out/circuits/peterjan_circuits/grid")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--verify-distance", action="store_true",
                        help="Assert shortest_graphlike_error >= d per circuit "
                             "(expensive; for the verification pass, not production).")
    parser.add_argument("--workers", type=int, default=mp.cpu_count())
    args = parser.parse_args()

    if max(args.gauge_factors) >= 4 and not _EXT_TABLE_LOADED:
        parser.error(
            f"gauge factor >= 4 requested but the extended timelike-distance "
            f"table is missing ({_EXT_TABLE_PATH}). Run "
            f"tools/peterjan_tools/extend_timelike_distance_data.py once first.")

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # n is deterministic per (d, g); probe once so filenames (and skip-existing)
    # are known without recompiling every grid point.
    print("Probing qubit counts per (d, g)...")
    n_by = {(d, g): _probe_n(d, g, args.observable)
            for d in args.distances for g in args.gauge_factors}

    jobs = [(d, g, n_by[(d, g)], p, pm)
            for d in args.distances
            for g in args.gauge_factors
            for p in args.p_values
            for pm in args.pm_values]

    n_workers = max(1, min(args.workers, len(jobs)))
    print(f"{len(args.distances)} d x {len(args.gauge_factors)} g x "
          f"{len(args.p_values)} p x {len(args.pm_values)} p_m = {len(jobs)} "
          f"circuits -> {out_dir} (workers={n_workers})")

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(_process, d, g, n, p, pm, out_dir, args.observable,
                          args.skip_existing, args.verify_distance)
                for d, g, n, p, pm in jobs]
        for done, fut in enumerate(as_completed(futs), 1):
            print(f"[{done}/{len(jobs)}]{fut.result()}")


if __name__ == "__main__":
    main()
