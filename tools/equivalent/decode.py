#!/usr/bin/env python3
"""Decode biased_bell_v3_equivalent circuits with sinter and write a CSV.

Thin wrapper around `tools/decode.py` for circuits produced by
`tools/equivalent/make_threshold_circuits.py`. Same collection logic and CSV
output (out/collection/); only the default circuit glob and output preamble are
specialized to the equivalent model.

Examples:

    python tools/equivalent/decode.py \
        --max-shots 100_000 --max-errors 200 --print-progress

    python tools/equivalent/decode.py \
        --circuit-files out/circuits/thresholds_equivalent/*.stim \
        --decoder pymatching --output-preamble equivalent_scan
"""

from __future__ import annotations

import argparse
import glob
import pathlib
import sys

_TOOLS = pathlib.Path(__file__).resolve().parent.parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))

import decode as base_decode  # noqa: E402

DEFAULT_GLOB = "out/circuits/thresholds_equivalent/*.stim"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--circuit-files", nargs="+", default=None,
        help=f"Translated .stim files to decode (default: {DEFAULT_GLOB}).",
    )
    parser.add_argument("--output-preamble", type=str, default="equivalent_scan",
                        help="Prefix for the output CSV (default: equivalent_scan).")
    base_decode.add_decode_args(parser)
    args = parser.parse_args()

    circuit_files = args.circuit_files
    if circuit_files is None:
        circuit_files = sorted(glob.glob(DEFAULT_GLOB))
        if not circuit_files:
            parser.error(f"No .stim files matched the default glob: {DEFAULT_GLOB}")

    base_decode.run_decode(
        circuit_files,
        max_shots=args.max_shots,
        max_errors=args.max_errors,
        decoder=args.decoder,
        output_preamble=args.output_preamble,
        cores=args.cores,
        print_progress=args.print_progress,
        comment=args.comment,
    )


if __name__ == "__main__":
    main()
