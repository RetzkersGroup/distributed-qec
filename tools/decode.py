#!/usr/bin/env python3
"""Sample translated circuits using sinter and write results to CSV.

Reads one or more translated .stim circuits (each with a sidecar .json
metadata file), runs sinter.collect() with the chosen decoder, and
appends rows to a CSV file in out/collection/.

Examples:

    python tools/decode.py \
        --circuit-files out/circuits/rep_code_biased_bell.stim \
        --max-shots 100_000 --max-errors 100 --decoder pymatching

    python tools/decode.py \
        --circuit-files out/*.stim \
        --max-shots 1_000_000 --max-errors 500 \
        --output-preamble hyperbolic_scan
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import sys
import time
from datetime import date
from multiprocessing import cpu_count

import pandas as pd
import sinter
import stim

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


def get_json_metadata(circuit_file: str) -> dict:
    json_file = circuit_file.replace(".stim", ".json")
    if os.path.exists(json_file):
        with open(json_file) as f:
            return json.load(f)
    print(f"Warning: metadata file {json_file} not found. Using empty metadata.",
          file=sys.stderr)
    return {}


def get_custom_decoders(decoder_name: str) -> dict[str, sinter.Decoder] | None:
    """Return custom decoder dict if needed, None for built-in decoders."""
    if decoder_name in ("pymatching", "fusion_blossom"):
        return None
    raise ValueError(
        f"Unknown decoder: {decoder_name!r}. "
        f"Supported: pymatching, fusion_blossom."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--circuit-files", nargs="+", required=True,
        help="Translated .stim circuit files to decode.",
    )
    parser.add_argument("--max-shots", type=int, default=10_000,
                        help="Max shots per task (default: 10_000).")
    parser.add_argument("--max-errors", type=int, default=100,
                        help="Max errors per task (default: 100).")
    parser.add_argument("--decoder", type=str, default="pymatching",
                        help="Decoder name (default: pymatching).")
    parser.add_argument("--output-preamble", type=str, default="output",
                        help="Prefix for the output CSV filename.")
    parser.add_argument("--cores", type=int, default=cpu_count(),
                        help="Number of worker processes.")
    parser.add_argument("--print-progress", action="store_true",
                        help="Show sinter progress during collection.")
    parser.add_argument("--comment", type=str, default="",
                        help="Extra comment stored in metadata.")
    args = parser.parse_args()

    log_dir = _ROOT / "out" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{date.today()}_decode.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )

    circuit_files = [f for f in args.circuit_files if f.endswith(".stim")]
    if not circuit_files:
        logging.error("No .stim files found in the provided arguments.")
        sys.exit(1)

    logging.info(
        f"Starting decode: {len(circuit_files)} circuit(s), "
        f"decoder={args.decoder}, max_shots={args.max_shots}, "
        f"max_errors={args.max_errors}, cores={args.cores}"
    )
    start_time = time.time()

    tasks_list = []
    for circuit_file in circuit_files:
        circuit = stim.Circuit.from_file(circuit_file)
        metadata = get_json_metadata(circuit_file)
        metadata["circuit_file"] = circuit_file
        if args.comment:
            metadata["comment"] = args.comment
        tasks_list.append(sinter.Task(
            circuit=circuit,
            json_metadata=metadata,
        ))

    custom_decoders = get_custom_decoders(args.decoder)
    samples = sinter.collect(
        num_workers=args.cores,
        max_shots=args.max_shots,
        max_errors=args.max_errors,
        tasks=tasks_list,
        decoders=[args.decoder],
        custom_decoders=custom_decoders,
        print_progress=args.print_progress,
    )

    out_dir = _ROOT / "out" / "collection"
    out_dir.mkdir(parents=True, exist_ok=True)
    today_str = date.today().strftime("%Y-%m-%d")
    outfile = out_dir / f"{args.output_preamble},{today_str}.csv"

    df = pd.DataFrame([vars(stat) for stat in samples])
    file_exists = outfile.is_file()
    df.to_csv(outfile, mode="a", header=not file_exists, index=False)

    elapsed = time.time() - start_time
    logging.info(f"Finished. Results saved to {outfile}")
    logging.info(f"Total time: {elapsed:.2f}s")
    print(f"\nResults: {outfile}")


if __name__ == "__main__":
    main()
