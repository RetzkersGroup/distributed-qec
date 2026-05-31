#!/usr/bin/env python3
"""Translate stim circuits' pairwise measurements via Bell pairs.

Reads one or more stim circuits (with optional sidecar .json metadata),
strips noise, replaces every 2-body MPP / MXX / MZZ / MYY with a
Bell-mediated equivalent, and applies the `biased_bell` noise model.
Writes each translated circuit plus a merged metadata json.

`--input` accepts:
  - a single .stim file
  - a directory (every .stim inside is processed, non-recursively)
  - a glob pattern (quote it to keep the shell from expanding it)

Examples:

    python tools/translate.py \
        --input "in/circuits/hyperbolic_floquet/circuit.stim" \
        --output-dir out/ \
        --p-init-x 0.01 --p-init-y 0.01 --p-init-z 0.01 \
        --p-internal 0.001 --p-m 0.01 --p-1 0.001

    python tools/translate.py \
        --input "in/circuits/test/*.stim" \
        --output-dir out/circuits/ \
        --noiseless
"""

from __future__ import annotations

import argparse
import glob
import json
import pathlib
import sys

import stim

# Make src/ importable when running this script directly.
_SRC = pathlib.Path(__file__).resolve().parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from bell_translation import NoiseParams, translate_circuit  # noqa: E402


def _resolve_inputs(input_arg: str) -> list[pathlib.Path]:
    """Expand --input into a list of .stim files."""
    p = pathlib.Path(input_arg)
    if p.is_dir():
        files = sorted(p.glob("*.stim"))
    elif p.exists():
        files = [p]
    else:
        # Treat as a glob pattern.
        files = sorted(pathlib.Path(s) for s in glob.glob(input_arg))
    return [f for f in files if f.suffix == ".stim"]


def _process_one(
    in_path: pathlib.Path,
    metadata_arg: str | None,
    out_dir: pathlib.Path,
    out_suffix: str,
    noise: NoiseParams | None,
    noise_model_name: str,
) -> None:
    meta_path = (
        pathlib.Path(metadata_arg) if metadata_arg else in_path.with_suffix(".json")
    )
    in_metadata: dict = {}
    if meta_path.exists():
        with meta_path.open() as f:
            in_metadata = json.load(f)
    else:
        print(
            f"[warn] metadata file not found: {meta_path} "
            f"(continuing with empty metadata)",
            file=sys.stderr,
        )

    out_stem = f"{in_path.stem}_{out_suffix}"
    out_circuit_path = out_dir / f"{out_stem}.stim"
    out_meta_path = out_dir / f"{out_stem}.json"

    circuit = stim.Circuit.from_file(str(in_path))
    translated = translate_circuit(circuit, noise)
    translated.to_file(str(out_circuit_path))

    out_metadata = {
        **in_metadata,
        "translation": "bell",
        "noise_model": noise_model_name if noise is not None else "noiseless",
        "noise_params": noise.to_dict() if noise is not None else None,
        "input_circuit": str(in_path),
        "input_metadata": str(meta_path) if meta_path.exists() else None,
        "output_circuit": str(out_circuit_path),
        "num_qubits": translated.num_qubits,
        "num_measurements": translated.num_measurements,
        "num_detectors": translated.num_detectors,
        "num_observables": translated.num_observables,
    }
    with out_meta_path.open("w") as f:
        json.dump(out_metadata, f, indent=2)

    print(f"Wrote {out_circuit_path}")
    print(f"Wrote {out_meta_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to .stim file, directory of .stim files, or glob pattern.",
    )
    parser.add_argument(
        "--metadata",
        help=(
            "Path to input metadata json. Only valid when --input resolves to a "
            "single file. By default, looks for <input-stem>.json next to each input."
        ),
    )
    parser.add_argument(
        "--output-dir", default="out/",
        help="Directory to write translated circuits and metadata into.",
    )
    parser.add_argument(
        "--output-suffix", default="biased_bell",
        help="Suffix appended to each output filename (default: biased_bell).",
    )

    parser.add_argument("--noise-model", default="biased_bell",
                        choices=["biased_bell"],
                        help="Noise model name (only biased_bell currently).")
    parser.add_argument("--p-init-x", type=float, default=0.0,
                        help="X component of Bell-pair prep noise.")
    parser.add_argument("--p-init-y", type=float, default=0.0,
                        help="Y component of Bell-pair prep noise.")
    parser.add_argument("--p-init-z", type=float, default=0.0,
                        help="Z component of Bell-pair prep noise.")
    parser.add_argument("--p-internal", type=float, default=0.0,
                        help="DEPOLARIZE2 after each ancilla<->data CX.")
    parser.add_argument("--p-m", type=float, default=0.0,
                        help="Measurement-flip probability (M(p), MX(p), ...).")
    parser.add_argument("--p-1", type=float, default=0.0, dest="p_1",
                        help="DEPOLARIZE1 after every single-qubit gate.")
    parser.add_argument("--p-init", type=float, default=0.0,
                        help="Data-qubit state-prep noise: X_ERROR after R, "
                             "Z_ERROR after RX.")

    parser.add_argument(
        "--noiseless", action="store_true",
        help="Skip noise injection (overrides all p_* flags).",
    )

    args = parser.parse_args()

    inputs = _resolve_inputs(args.input)
    if not inputs:
        parser.error(f"No .stim files matched: {args.input}")

    if args.metadata and len(inputs) > 1:
        parser.error("--metadata is only valid when --input matches a single file.")

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    noise = None
    if not args.noiseless:
        noise = NoiseParams(
            p_init_x=args.p_init_x,
            p_init_y=args.p_init_y,
            p_init_z=args.p_init_z,
            p_internal=args.p_internal,
            p_m=args.p_m,
            p_1=args.p_1,
            p_init=args.p_init,
        )

    for in_path in inputs:
        _process_one(
            in_path=in_path,
            metadata_arg=args.metadata,
            out_dir=out_dir,
            out_suffix=args.output_suffix,
            noise=noise,
            noise_model_name=args.noise_model,
        )


if __name__ == "__main__":
    main()
