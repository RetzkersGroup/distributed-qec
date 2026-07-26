#!/usr/bin/env python3
"""Generate biased_bell_v3 circuits on a direct (p, p_bell_z) grid.

Unlike tools/make_threshold_circuits_bell.py -- which sweeps a single "overall
strength" p along bias rays (p_bell_z = beta * p) or fixed-p_bell_z lines -- this
scans p and p_bell_z as two *independent* axes: every (p, p_bell_z) pair on the
outer product of the two grids is emitted. (Downstream analysis is TBD.)

For each input circuit, strips existing noise and re-applies the biased_bell_v3
Bell-mediated noise model (via bell_translation) at each grid point. The noise
derivation is imported from make_threshold_circuits_bell so the two generators
share one source of truth. Metadata records p and p_bell_z directly
(direction_kind="grid"), plus bias = p_bell_z / p for convenience.

Examples:

    python tools/make_grid_circuits_bell.py \\
        --input "in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim" \\
        --output-dir out/circuits/grid_bell

    python tools/make_grid_circuits_bell.py \\
        --input "in/circuits/honeycomb_gidney/*obs=H*.stim" \\
        --p-values 0.001 0.003 0.005 \\
        --p-bell-z-values 0.0 0.001 0.01 0.1
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import stim

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
_SRC = _ROOT / "src"
_TOOLS = _ROOT / "tools"
for _p in (str(_SRC), str(_TOOLS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from bell_translation import translate_circuit  # noqa: E402
from repeat_measurements import repeat_measurements  # noqa: E402

# Reuse the noise derivation / input resolution from the bias-scan generator so
# the two share one source of truth (a circuit at a given (p, p_bell_z) is
# byte-identical regardless of which script wrote it).
from make_threshold_circuits_bell import (  # noqa: E402
    _load_metadata,
    _model_name,
    _noise_from_p,
    _resolve_inputs,
)

# Default grids: p reuses the bias scan's log grid; p_bell_z spans a slightly
# wider log range since it is now independent of p.
DEFAULT_P: tuple[float, ...] = tuple(np.logspace(-4, -2, 16))
DEFAULT_PBZ: tuple[float, ...] = tuple(np.logspace(-4, -1, 16))

# PAULI_CHANNEL_1 needs p_bell_x + p_bell_y + p_bell_z < 1; values near that
# bound are unphysical / uninteresting, so skip p_bell_z above this.
_PBZ_MAX = 0.4


def _fmt(x: float) -> str:
    return f"{x:.5g}"


def _process(
    in_path: pathlib.Path,
    p: float,
    p_bell_z: float,
    out_dir: pathlib.Path,
    skip_existing: bool,
    mzz_gadget: str,
    scheme: str,
    noiseless_start_end: bool,
    repeat: int = 1,
) -> str:
    model_name = _model_name(mzz_gadget, scheme)
    if repeat > 1:
        # Full-decoder repeated measurements: the structural repetition is
        # applied BEFORE translation, and the name gets a _fullrep{r} suffix so
        # downstream fits/plots (which key on noise_model) see a new scheme.
        model_name = f"{model_name}_fullrep{repeat}"
    out_stem = f"{in_path.stem}__{model_name},pz={_fmt(p_bell_z)},p={_fmt(p)}"
    out_circuit_path = out_dir / f"{out_stem}.stim"
    out_meta_path = out_dir / f"{out_stem}.json"

    if p_bell_z > _PBZ_MAX:
        return f"  pz={p_bell_z:.4g} p={p:.4g} -> skipped (p_bell_z > {_PBZ_MAX})"

    if skip_existing and out_circuit_path.exists() and out_meta_path.exists():
        return f"  pz={p_bell_z:.4g} p={p:.4g} -> {out_circuit_path.name} (skipped)"

    raw = stim.Circuit.from_file(str(in_path))
    if repeat > 1:
        raw = repeat_measurements(raw, repeat)
    noise = _noise_from_p(p, p_bell_z, mzz_gadget, scheme)
    translated = translate_circuit(
        raw, noise, noiseless_start_end=noiseless_start_end,
    )
    translated.to_file(str(out_circuit_path))

    src_meta = _load_metadata(in_path)
    metadata = {
        **src_meta,  # carry through code_family/k/expected_distance/rounds/...
        "p": float(p),
        "p_bell_z": float(p_bell_z),
        "direction_kind": "grid",
        "direction_value": None,
        # bias is a convenience readout, not a swept axis here.
        "bias": (float(p_bell_z / p) if p > 0 else None),
        "p_bell_z_fixed": None,
        "noise_model": model_name,
        "mzz_gadget": mzz_gadget,
        "scheme": scheme,
        "repetitions": repeat,
        "decoder_kind": "full" if repeat > 1 else "single",
        "noise_params": noise.to_dict(),
        "noiseless_start_end": noiseless_start_end,
        "translation": "bell",
        "input_circuit": str(in_path),
        "output_circuit": str(out_circuit_path),
        "num_qubits": translated.num_qubits,
        "num_measurements": translated.num_measurements,
        "num_detectors": translated.num_detectors,
        "num_observables": translated.num_observables,
    }
    with out_meta_path.open("w") as f:
        json.dump(metadata, f, indent=2)

    return f"  pz={p_bell_z:.4g} p={p:.4g} -> {out_circuit_path.name}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to .stim file, directory, or glob pattern.",
    )
    parser.add_argument(
        "--p-values", nargs="+", type=float, default=None,
        help="p grid (default: np.logspace(-4, -2, 16)).",
    )
    parser.add_argument(
        "--p-bell-z-values", nargs="+", type=float, default=None,
        dest="p_bell_z_values",
        help="p_bell_z grid (default: np.logspace(-4, -1, 16)). "
             "0.0 is allowed (Bell-Z prep noise off).",
    )
    parser.add_argument(
        "--output-dir", default="out/circuits/grid_bell",
        help="Output directory (default: out/circuits/grid_bell).",
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip (p, p_bell_z) tasks whose .stim and .json already exist.",
    )
    parser.add_argument(
        "--workers", type=int, default=mp.cpu_count(),
        help="Parallel worker processes (default: cpu_count()).",
    )
    parser.add_argument(
        "--noise-model", default="biased_bell_v3",
        choices=["biased_bell_v2", "biased_bell_v3"], dest="noise_model",
        help="Noise model name (default: biased_bell_v3). Sets the MZZ gadget "
             "shape; the noise channels are identical between v2 and v3.",
    )
    parser.add_argument(
        "--mzz-gadget", default=None, choices=["v2", "v3"], dest="mzz_gadget",
        help="Explicit MZZ gadget override. Defaults to match --noise-model.",
    )
    parser.add_argument(
        "--scheme", default="bell",
        choices=["bell", "equivalent_v3", "equivalent_v3_majority", "equivalent_v3_majority5"],
        dest="scheme",
        help="Translation scheme (default: bell = two-ancilla Bell pair).",
    )
    parser.add_argument(
        "--noiseless-start-end", action="store_true",
        dest="noiseless_start_end",
        help="If set, only inject Bell-prep noise inside the "
             "TICK[start_noise]..TICK[stop_noise] window in the input circuit.",
    )
    parser.add_argument(
        "--repeat", type=int, default=1,
        help="Full-decoder measurement repetitions: apply repeat_measurements "
             "(r copies of each mid-circuit measurement layer + comparison "
             "detectors) before translation. Default 1 = no repetition. "
             "Do NOT combine with the equivalent_v3_majority* schemes (those "
             "already fold repetition into the noise block).",
    )
    args = parser.parse_args()

    if args.repeat > 1 and "majority" in args.scheme:
        parser.error(
            "--repeat > 1 double-counts repetition when combined with a "
            "majority scheme; use --scheme bell or equivalent_v3."
        )

    if args.mzz_gadget is None:
        args.mzz_gadget = "v2" if args.noise_model == "biased_bell_v2" else "v3"

    inputs = _resolve_inputs(args.input)
    if not inputs:
        parser.error(f"No .stim files matched: {args.input}")

    p_values = list(DEFAULT_P if args.p_values is None else args.p_values)
    pbz_values = list(
        DEFAULT_PBZ if args.p_bell_z_values is None else args.p_bell_z_values
    )

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs: list[tuple[pathlib.Path, float, float]] = []
    for in_path in inputs:
        for p in p_values:
            for pbz in pbz_values:
                jobs.append((in_path, p, pbz))

    n_workers = max(1, min(args.workers, len(jobs)))
    print(
        f"{len(inputs)} circuits x {len(p_values)} p-values x "
        f"{len(pbz_values)} p_bell_z-values = {len(jobs)} circuits to write "
        f"-> {out_dir} (workers={n_workers})"
    )

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [
            ex.submit(_process, in_path, p, pbz, out_dir,
                      args.skip_existing, args.mzz_gadget, args.scheme,
                      args.noiseless_start_end, args.repeat)
            for in_path, p, pbz in jobs
        ]
        for done, fut in enumerate(as_completed(futures), 1):
            print(f"[{done}/{len(jobs)}]{fut.result()}")


if __name__ == "__main__":
    main()
