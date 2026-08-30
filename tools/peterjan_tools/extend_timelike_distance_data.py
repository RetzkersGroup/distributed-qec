#!/usr/bin/env python3
"""Extend Kandel's timelike-distance table to gauge factors 4 and 5.

GaugeHoneycombCode.get_number_of_rounds_for_timelike_distance reads a precomputed
graphlike-timelike-distance table shipped inside the Kandel package
(`.../timelike_distance_data/hcc_graphlike_td_data_phenomenological.json`), which
only covers uniform gauge factors 1..3 (the paper's scanned range). g >= 4 raises
a KeyError.

This one-time script computes the missing uniform-gauge-factor entries
((4,4,4), (5,5,5), ...) using Kandel's OWN generation logic
(create_timelike_distance_data.get_td_bulk_and_boundary), merges them with the
shipped table, and writes a superset to
tools/peterjan_tools/data/hcc_graphlike_td_data_phenomenological_extended.json.
The 1..3 entries are copied verbatim, so g<=3 behaviour is unchanged.

make_grid_circuits_peterjan.py loads this extended table and monkeypatches
GaugeHoneycombCode.get_timelike_distance to use it (contained: the installed
Kandel package is never mutated).

Usage:
    python tools/peterjan_tools/extend_timelike_distance_data.py            # g up to 5
    python tools/peterjan_tools/extend_timelike_distance_data.py --max-g 4
"""

from __future__ import annotations

import argparse
import json
import pathlib
import time

# Kandel (installed as `main`): the shipped table and the upstream generator.
from main.codes.tic_tac_toe.gauge import GaugeHoneycombCode as _ghc_mod
from main.codes.tic_tac_toe.gauge.create_timelike_distance_data import (
    get_td_bulk_and_boundary,
)

_DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"
_OUT = _DATA_DIR / "hcc_graphlike_td_data_phenomenological_extended.json"
# The noise-model string create_timelike_distance_data.generate_circuit expects
# (NOT the "phenomenological_noise" string get_timelike_distance keys files on).
_NOISE = "phenomenological"


def _shipped_table_path() -> pathlib.Path:
    return (pathlib.Path(_ghc_mod.__file__).parent / "timelike_distance_data"
            / "hcc_graphlike_td_data_phenomenological.json")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-g", type=int, default=5,
                    help="Highest uniform gauge factor to add (default 5).")
    args = ap.parse_args()

    shipped = _shipped_table_path()
    with shipped.open() as f:
        table = json.load(f)  # {"X": {...}, "Z": {...}}
    print(f"Loaded shipped table: {shipped}")

    new_gs = [g for g in range(4, args.max_g + 1)]
    for g in new_gs:
        gauge_factors = (g, g, g)
        key = str(gauge_factors)              # matches upstream f"({a}, {b}, {c})"
        bulk_length = 4 * sum(gauge_factors)  # GaugeHoneycombCode bulk length
        for letter in ("X", "Z"):
            if key in table.get(letter, {}):
                print(f"  [{letter}] {key} already present — skipping")
                continue
            t0 = time.time()
            td_bulk, td_boundary = get_td_bulk_and_boundary(
                _NOISE, list(gauge_factors), bulk_length, letter,
                True, "GaugeHoneycombCode")
            # JSON object keys are strings; upstream indexes td_boundary[str(i)].
            td_boundary = {str(k): v for k, v in td_boundary.items()}
            table.setdefault(letter, {})[key] = {
                "td_bulk": td_bulk, "td_boundary": td_boundary}
            print(f"  [{letter}] {key}: td_bulk={td_bulk}, "
                  f"{len(td_boundary)} boundary layers, {time.time()-t0:.1f}s")

    _DATA_DIR.mkdir(parents=True, exist_ok=True)
    with _OUT.open("w") as f:
        json.dump(table, f, indent=4)
    print(f"\nWrote extended table -> {_OUT}")


if __name__ == "__main__":
    main()
