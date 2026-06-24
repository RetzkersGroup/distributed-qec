#!/usr/bin/env python3
"""Generate biased_bell_v3 circuits for a (p, p_bell_z) threshold scan.

For each input circuit, strips the existing noise and re-applies the
`biased_bell_v3` Bell-mediated noise model (via `bell_translation`) for
every (direction, p) pair. v2 is selectable via `--mzz-gadget v2`; only
the MZZ gadget shape differs between the two. A direction is either:

  * a **bias** `beta` -> p_bell_z = beta * p, or
  * a **fixed p_bell_z** value -> p_bell_z held constant while p sweeps.

Default directions cover both:
  - biases:        beta in {1, 3, 10, 30, 100}
  - fixed p_bell_z: {1e-3, 1e-2}

Examples:

    python tools/make_threshold_circuits_bell.py \\
        --input "in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim" \\
        --output-dir out/circuits/thresholds_bell

    python tools/make_threshold_circuits_bell.py \\
        --input "in/circuits/honeycomb_gidney/*.stim" \\
        --biases 0 1 10 100 \\
        --p-bell-z-values 0.001 \\
        --p-values 0.001 0.003 0.005 0.01
"""

from __future__ import annotations

import argparse
import glob
import json
import multiprocessing as mp
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import stim

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from bell_translation import NoiseParams, translate_circuit  # noqa: E402


DEFAULT_BIASES: tuple[float, ...] = (1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0)
DEFAULT_PBZ_FIXED: tuple[float, ...] = (1e-3, 1e-2)


def _resolve_inputs(input_arg: str) -> list[pathlib.Path]:
    p = pathlib.Path(input_arg)
    if p.is_dir():
        files = sorted(p.glob("*.stim"))
    elif p.exists():
        files = [p]
    else:
        files = sorted(pathlib.Path(s) for s in glob.glob(input_arg))
    return [f for f in files if f.suffix == ".stim"]


def _load_metadata(stim_path: pathlib.Path) -> dict:
    json_path = stim_path.with_suffix(".json")
    if json_path.exists():
        with open(json_path) as f:
            return json.load(f)
    return {}


def _noise_from_p(p: float, p_bell_z: float, mzz_gadget: str) -> NoiseParams:
    return NoiseParams(
        p_bell_x=p / 2,
        p_bell_y=p / 2,
        p_bell_z=p_bell_z,
        p_data_z=p,
        p_internal=0.0,
        p_m=p,
        p_1=p / 10,
        p_init=p,
        mzz_gadget=mzz_gadget,
    )


def _format_bias(b: float) -> str:
    if b == int(b):
        return str(int(b))
    return f"{b:g}".replace(".", "p")


def _process(
    in_path: pathlib.Path,
    direction_kind: str,         # "bias" or "pbz"
    direction_value: float,      # beta (for "bias") or fixed p_bell_z (for "pbz")
    p: float,
    out_dir: pathlib.Path,
    skip_existing: bool,
    mzz_gadget: str,
) -> str:
    if direction_kind == "bias":
        p_bell_z = direction_value * p
        dir_token = f"bias={_format_bias(direction_value)}"
    elif direction_kind == "pbz":
        p_bell_z = direction_value
        dir_token = f"pbz={direction_value:.5g}"
    else:
        raise ValueError(direction_kind)

    model_name = "biased_bell_v2" if mzz_gadget == "v2" else "biased_bell_v3"
    out_stem = (
        f"{in_path.stem}__{model_name},{dir_token},p={p:.5g}"
    )
    out_circuit_path = out_dir / f"{out_stem}.stim"
    out_meta_path = out_dir / f"{out_stem}.json"

    # Skip out-of-range cases: PAULI_CHANNEL_1 needs p_bell_x + p_bell_y +
    # p_bell_z < 1, and even values approaching that bound are unphysical /
    # uninteresting for threshold extraction.
    if p_bell_z > 0.4:
        return (
            f"  {dir_token:>20s}  p={p:.4g} -> skipped "
            f"(p_bell_z={p_bell_z:.4g} > 0.4)"
        )

    if skip_existing and out_circuit_path.exists() and out_meta_path.exists():
        return f"  {dir_token:>20s}  p={p:.4g} -> {out_circuit_path.name} (skipped)"

    raw = stim.Circuit.from_file(str(in_path))
    noise = _noise_from_p(p, p_bell_z, mzz_gadget)
    translated = translate_circuit(raw, noise)
    translated.to_file(str(out_circuit_path))

    src_meta = _load_metadata(in_path)
    metadata = {
        **src_meta,  # carry through code_family/k/d/rounds/w/h/observable/...
        "p": float(p),
        "p_bell_z": float(p_bell_z),
        "direction_kind": direction_kind,
        "direction_value": float(direction_value),
        "bias": float(direction_value) if direction_kind == "bias" else None,
        "p_bell_z_fixed": (
            float(direction_value) if direction_kind == "pbz" else None
        ),
        "noise_model": model_name,
        "mzz_gadget": mzz_gadget,
        "noise_params": noise.to_dict(),
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

    return f"  {dir_token:>20s}  p={p:.4g} -> {out_circuit_path.name}"


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
        help="Physical error rates (default: np.logspace(-4, -2, 16)).",
    )
    parser.add_argument(
        "--biases", nargs="*", type=float, default=None,
        help=f"Bias values beta = p_bell_z / p "
             f"(default: {list(DEFAULT_BIASES)}). Use empty list to disable.",
    )
    parser.add_argument(
        "--p-bell-z-values", nargs="*", type=float, default=None,
        dest="p_bell_z_values",
        help=f"Fixed p_bell_z values (default: {list(DEFAULT_PBZ_FIXED)}). "
             f"Use empty list to disable.",
    )
    parser.add_argument(
        "--output-dir", default="out/circuits/thresholds_bell",
        help="Output directory (default: out/circuits/thresholds_bell).",
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip (direction, p) tasks whose .stim and .json already exist.",
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
        help="Explicit MZZ gadget override. Defaults to match --noise-model "
             "(biased_bell_v2 -> v2, biased_bell_v3 -> v3).",
    )
    args = parser.parse_args()

    if args.mzz_gadget is None:
        args.mzz_gadget = "v2" if args.noise_model == "biased_bell_v2" else "v3"

    inputs = _resolve_inputs(args.input)
    if not inputs:
        parser.error(f"No .stim files matched: {args.input}")

    if args.p_values is None:
        p_values = np.logspace(-4, -2, 16).tolist()
    else:
        p_values = list(args.p_values)

    biases = list(DEFAULT_BIASES if args.biases is None else args.biases)
    pbz_fixed = list(
        DEFAULT_PBZ_FIXED if args.p_bell_z_values is None else args.p_bell_z_values
    )
    if not biases and not pbz_fixed:
        parser.error(
            "At least one of --biases or --p-bell-z-values must be non-empty."
        )

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs: list[tuple[pathlib.Path, str, float, float]] = []
    for in_path in inputs:
        for beta in biases:
            for p in p_values:
                jobs.append((in_path, "bias", beta, p))
        for pbz in pbz_fixed:
            for p in p_values:
                jobs.append((in_path, "pbz", pbz, p))

    n_workers = max(1, min(args.workers, len(jobs)))
    print(
        f"{len(inputs)} circuits x ({len(biases)} biases + "
        f"{len(pbz_fixed)} fixed pbz) x {len(p_values)} p-values "
        f"= {len(jobs)} circuits to write -> {out_dir} (workers={n_workers})"
    )

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [
            ex.submit(_process, in_path, kind, val, p, out_dir,
                      args.skip_existing, args.mzz_gadget)
            for in_path, kind, val, p in jobs
        ]
        for done, fut in enumerate(as_completed(futures), 1):
            print(f"[{done}/{len(jobs)}]{fut.result()}")


if __name__ == "__main__":
    main()
