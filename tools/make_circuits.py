#!/usr/bin/env python3
"""Generate translated (Bell-mediated) circuits for a sweep of noise levels.

Takes noiseless input circuits, translates each through the Bell
measurement scheme at multiple physical error rates, and writes the
output circuits + metadata JSONs.

The noise parameters are derived from a single physical error rate `p`:
  - p_bell_x = p_bell_y = p_bell_z = p  (Bell-pair prep noise)
  - p_internal = p                       (CX depolarizing)
  - p_m = p                              (measurement flip)
  - p_1 = p / 10                         (single-qubit gate depol)
  - p_init = p                           (data-qubit reset noise)

These defaults can be overridden with --noise-profile or individual flags.

Examples:

    # Sweep p = 0.0005, 0.001, 0.002 over all k=10 circuits:
    python tools/make_circuits.py \
        --input "in/circuits/hyperbolic_floquet/8.8.8/k=10_d3a22e5e/*.stim" \
        --p-values 0.0005 0.001 0.002 0.005 \
        --output-dir out/circuits/hyperbolic_floquet

    # Single circuit, single p:
    python tools/make_circuits.py \
        --input in/circuits/hyperbolic_floquet/8.8.8/k=10_d3a22e5e/semi_hyperbolic_floquet_l=2_n=256_k=10_d=4_basis=Z_num_sub_rounds=96_noise_model=SDEM3_p=0.0.stim \
        --p-values 0.001 \
        --output-dir out/circuits/hyperbolic_floquet
"""

from __future__ import annotations

import argparse
import glob
import json
import pathlib
import re
import sys

import stim

BELL_BIAS = 1

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from bell_translation import NoiseParams, translate_circuit  # noqa: E402


def _resolve_inputs(input_arg: str) -> list[pathlib.Path]:
    p = pathlib.Path(input_arg)
    if p.is_dir():
        files = sorted(p.glob("*.stim"))
    elif p.exists():
        files = [p]
    else:
        files = sorted(pathlib.Path(s) for s in glob.glob(input_arg))
    return [f for f in files if f.suffix == ".stim"]


def _parse_filename(filename: str) -> dict:
    """Extract n, k, d, num_sub_rounds, rounds from hyperbolic filename."""
    patterns = {
        "n": r"n=(\d+)",
        "k": r"k=(\d+)",
        "d": r"d=(\d+)",
        "num_sub_rounds": r"num_sub_rounds=(\d+)",
        "l": r"l=(\d+)",
    }
    params: dict = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, filename)
        params[key] = int(match.group(1)) if match else None

    if params["num_sub_rounds"] is not None:
        params["rounds"] = (params["num_sub_rounds"] - 6) // 6
    else:
        params["rounds"] = None

    basis_match = re.search(r"basis=([XZ])", filename)
    params["basis"] = basis_match.group(1) if basis_match else None

    return params


def _noise_from_p(p: float, bell_bias: float) -> NoiseParams:
    """All channels = p except p_bell_z = bell_bias * p and p_1 = p/10."""
    return NoiseParams(
        p_bell_x=p,
        p_bell_y=p,
        p_bell_z=bell_bias * p,
        p_internal=p,
        p_m=p,
        p_1=p / 10,
        p_init=p,
    )


def _process(
    in_path: pathlib.Path,
    p: float,
    bell_bias: float,
    out_dir: pathlib.Path,
) -> None:
    noise = _noise_from_p(p, bell_bias)
    circuit = stim.Circuit.from_file(str(in_path))
    translated = translate_circuit(circuit, noise)

    out_stem = f"{in_path.stem}_biased_bell_p={p}_bias={bell_bias}"
    out_circuit_path = out_dir / f"{out_stem}.stim"
    out_meta_path = out_dir / f"{out_stem}.json"

    translated.to_file(str(out_circuit_path))

    parsed = _parse_filename(in_path.name)
    metadata = {
        "p": p,
        "bell_bias": bell_bias,
        "noise_model": "biased_bell",
        "noise_params": noise.to_dict(),
        "translation": "bell",
        "code_family": "semi_hyperbolic_floquet",
        "n": parsed["n"],
        "k": parsed["k"],
        "expected_distance": parsed["d"],
        "rounds": parsed["rounds"],
        "num_sub_rounds": parsed["num_sub_rounds"],
        "basis": parsed["basis"],
        "l": parsed["l"],
        "input_circuit": str(in_path),
        "output_circuit": str(out_circuit_path),
        "num_qubits": translated.num_qubits,
        "num_measurements": translated.num_measurements,
        "num_detectors": translated.num_detectors,
        "num_observables": translated.num_observables,
    }

    with out_meta_path.open("w") as f:
        json.dump(metadata, f, indent=2)

    print(f"  p={p} → {out_circuit_path.name}")


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
        "--p-values", nargs="+", type=float, required=True,
        help="Physical error rates to sweep.",
    )
    parser.add_argument(
        "--bell-bias", type=float, default=BELL_BIAS,
        help=f"Multiplicative factor for p_bell_z (default: {BELL_BIAS}).",
    )
    parser.add_argument(
        "--output-dir", default="out/circuits/hyperbolic_floquet",
        help="Output directory (default: out/circuits/hyperbolic_floquet).",
    )

    args = parser.parse_args()
    inputs = _resolve_inputs(args.input)
    if not inputs:
        parser.error(f"No .stim files matched: {args.input}")

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for in_path in inputs:
        print(f"{in_path.name}:")
        for p in args.p_values:
            _process(in_path, p, args.bell_bias, out_dir)


if __name__ == "__main__":
    main()
