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


def _write_sinter_csv(path: pathlib.Path, stats) -> None:
    """Persist sinter TaskStats in sinter's own CSV format (JSON metadata).

    sinter only writes CSVs transiently (save_resume_filepath, deleted on clean
    completion), so we serialize the returned stats ourselves using sinter's
    CSV_HEADER + TaskStats.to_csv_line(). The resulting file is what
    `--resume` reads back via existing_data_filepaths.
    """
    with open(path, "w") as f:
        f.write(sinter.CSV_HEADER + "\n")
        for s in stats:
            line = s.to_csv_line()
            f.write(line if line.endswith("\n") else line + "\n")


# Beam search decoder configs, Table I of arXiv:2512.07057. The name (used as the
# decoder name, like tesseract_short/long) encodes beam_width and total BP iters
# (= initial_iters + max_rounds * iters_per_round); "res" is num_results (retained
# solutions).
BEAMSEARCH_CONFIGS: dict[str, dict[str, int]] = {
    "beam8_230iters":        dict(beam_width=8,  max_rounds=10, initial_iters=30, iters_per_round=20, num_results=1),
    "beam32_340iters":       dict(beam_width=32, max_rounds=10, initial_iters=40, iters_per_round=30, num_results=1),
    "beam64_640iters":       dict(beam_width=64, max_rounds=20, initial_iters=40, iters_per_round=30, num_results=1),
    "beam64_32res_640iters": dict(beam_width=64, max_rounds=20, initial_iters=40, iters_per_round=30, num_results=32),
}


def get_custom_decoders(decoder_name: str) -> dict[str, sinter.Decoder] | None:
    """Return custom decoder dict if needed, None for built-in decoders."""
    # Built-in sinter decoders (no custom_decoders entry needed).
    # 'pymatching-correlated' is native correlated MWPM (pymatching >= 2.3:
    # two-pass sparse blossom that reweights edges correlated via decomposed
    # hyperedges) — handles the non-graphlike DEM better than plain pymatching.
    if decoder_name in ("pymatching", "pymatching-correlated", "fusion_blossom",
                        "hypergraph_union_find", "mw_parity_factor", "vacuous"):
        return None
    if decoder_name in ("tesseract_short", "tesseract_long"):
        from tesseract_decoder_utils import get_tesseract_decoder
        return get_tesseract_decoder(decoder_name.removeprefix("tesseract_"))
    if decoder_name == "belief_matching":
        # Correlated MWPM: belief propagation reweights the matching graph
        # (accounting for correlated / hyperedge errors that plain pymatching
        # decomposes lossily), then pymatching matches on the reweighted graph.
        from beliefmatching import BeliefMatchingSinterDecoder
        return {"belief_matching": BeliefMatchingSinterDecoder(max_bp_iters=20)}
    if decoder_name in BEAMSEARCH_CONFIGS:
        # IonQ BP + beam-search decoder for qLDPC codes (arXiv:2512.07057). Not
        # pip-installable: gitignored clone (BeamSearchDecoder/) built in place and
        # put on the venv path via a .pth file -- see README. The decoder name is
        # the Table I config (see BEAMSEARCH_CONFIGS), so it lands in the CSV
        # `decoder` column and the plot legend.
        from sinter_beamsearch import SinterDecoder_BeamSearch
        return {decoder_name: SinterDecoder_BeamSearch(**BEAMSEARCH_CONFIGS[decoder_name])}
    if decoder_name in ("relay-bp", "mem-bp", "msl-bp"):
        from relay_bp.stim import sinter_decoders
        return sinter_decoders()
    if decoder_name in ("relay-bp-light", "relay-bp-ultralight"):
        # Lighter relay ensembles: the library default (60 sets x 60 iters,
        # stop after 5 agreeing solutions) costs ~2.5 s/shot on the d=7
        # repeated honeycomb circuits near threshold; these presets are
        # ~5x / ~13x faster with no accuracy loss observed in spot checks.
        from relay_bp.stim.sinter.decoders import SinterDecoder_RelayBP
        if decoder_name == "relay-bp-light":
            dec = SinterDecoder_RelayBP(num_sets=10, stop_nconv=1)
        else:
            dec = SinterDecoder_RelayBP(
                num_sets=4, stop_nconv=1, pre_iter=30, set_max_iter=30,
            )
        return {decoder_name: dec}
    raise ValueError(
        f"Unknown decoder: {decoder_name!r}. "
        f"Supported: pymatching, pymatching-correlated, fusion_blossom, "
        f"belief_matching, tesseract_short, tesseract_long, "
        f"{', '.join(BEAMSEARCH_CONFIGS)} (beam search), "
        f"relay-bp, relay-bp-light, relay-bp-ultralight, mem-bp, msl-bp."
    )


def run_decode(
    circuit_files: list[str],
    *,
    max_shots: int = 10_000,
    max_errors: int = 100,
    decoder: str = "pymatching",
    output_preamble: str = "output",
    cores: int | None = None,
    print_progress: bool = False,
    comment: str = "",
    resume: bool = False,
) -> pathlib.Path:
    """Collect decoding statistics for `circuit_files` and write a CSV.

    Shared by the generic CLI and the per-model wrappers (e.g.
    tools/equivalent/decode.py). Returns the output CSV path.

    Default (append) mode writes a dated pandas CSV `{preamble},{date}.csv`.
    With ``resume=True`` it instead maintains a single sinter-format CSV
    `{preamble}.csv`: prior results in it are counted toward max_shots/max_errors
    (sinter's existing_data_filepaths), so a re-run only samples the shortfall.
    """
    if cores is None:
        cores = cpu_count()

    log_dir = _ROOT / "out" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{date.today()}_decode.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )

    circuit_files = [f for f in circuit_files if f.endswith(".stim")]
    if not circuit_files:
        logging.error("No .stim files found in the provided arguments.")
        sys.exit(1)

    logging.info(
        f"Starting decode: {len(circuit_files)} circuit(s), "
        f"decoder={decoder}, max_shots={max_shots}, "
        f"max_errors={max_errors}, cores={cores}"
    )
    start_time = time.time()

    tasks_list = []
    for circuit_file in circuit_files:
        circuit = stim.Circuit.from_file(circuit_file)
        metadata = get_json_metadata(circuit_file)
        metadata["circuit_file"] = circuit_file
        if comment:
            metadata["comment"] = comment
        tasks_list.append(sinter.Task(
            circuit=circuit,
            json_metadata=metadata,
        ))

    out_dir = _ROOT / "out" / "collection"
    out_dir.mkdir(parents=True, exist_ok=True)

    custom_decoders = get_custom_decoders(decoder)
    collect_kwargs = dict(
        num_workers=cores,
        max_shots=max_shots,
        max_errors=max_errors,
        tasks=tasks_list,
        decoders=[decoder],
        custom_decoders=custom_decoders,
        print_progress=print_progress,
    )

    if resume:
        # Single canonical sinter-format CSV; sinter counts its existing rows
        # toward the budget and we write the merged result back.
        outfile = out_dir / f"{output_preamble}.csv"
        if outfile.is_file():
            collect_kwargs["existing_data_filepaths"] = [str(outfile)]
            logging.info(f"Resume: counting existing results in {outfile}")
        samples = sinter.collect(**collect_kwargs)
        _write_sinter_csv(outfile, samples)
    else:
        today_str = date.today().strftime("%Y-%m-%d")
        outfile = out_dir / f"{output_preamble},{today_str}.csv"
        samples = sinter.collect(**collect_kwargs)
        df = pd.DataFrame([vars(stat) for stat in samples])
        # Serialize concurrent appends (e.g. PBS array tasks sharing one CSV)
        # with an advisory lock on a sidecar lockfile.
        import fcntl
        lockfile = outfile.with_suffix(".lock")
        with open(lockfile, "w") as lf:
            fcntl.flock(lf, fcntl.LOCK_EX)
            file_exists = outfile.is_file() and outfile.stat().st_size > 0
            df.to_csv(outfile, mode="a", header=not file_exists, index=False)
            fcntl.flock(lf, fcntl.LOCK_UN)

    elapsed = time.time() - start_time
    logging.info(f"Finished. Results saved to {outfile}")
    logging.info(f"Total time: {elapsed:.2f}s")
    print(f"\nResults: {outfile}")
    return outfile


def add_decode_args(parser: argparse.ArgumentParser) -> None:
    """Attach the shared decode options to `parser` (everything but inputs)."""
    parser.add_argument("--max-shots", type=int, default=10_000,
                        help="Max shots per task (default: 10_000).")
    parser.add_argument("--max-errors", type=int, default=100,
                        help="Max errors per task (default: 100).")
    parser.add_argument("--decoder", type=str, default="pymatching",
                        help="Decoder name (default: pymatching).")
    parser.add_argument("--cores", type=int, default=cpu_count(),
                        help="Number of worker processes.")
    parser.add_argument("--print-progress", action="store_true",
                        help="Show sinter progress during collection.")
    parser.add_argument("--comment", type=str, default="",
                        help="Extra comment stored in metadata.")
    parser.add_argument("--resume", action="store_true",
                        help="Maintain a single sinter-format {preamble}.csv and "
                             "resume from it (count existing shots/errors toward "
                             "the budget) instead of appending a dated pandas CSV.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--circuit-files", nargs="+", required=True,
        help="Translated .stim circuit files to decode.",
    )
    parser.add_argument("--output-preamble", type=str, default="output",
                        help="Prefix for the output CSV filename.")
    add_decode_args(parser)
    args = parser.parse_args()

    run_decode(
        args.circuit_files,
        max_shots=args.max_shots,
        max_errors=args.max_errors,
        decoder=args.decoder,
        output_preamble=args.output_preamble,
        cores=args.cores,
        print_progress=args.print_progress,
        comment=args.comment,
        resume=args.resume,
    )


if __name__ == "__main__":
    main()
