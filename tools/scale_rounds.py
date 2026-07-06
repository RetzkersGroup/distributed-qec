#!/usr/bin/env python3
"""Rescale the REPEAT (bulk QEC round) count of memory circuits.

For each input .stim, set its single top-level REPEAT block to a fixed count
(`--rounds`) or to `--rounds-multiple * d` rounds (d parsed from the filename) --
e.g. `--rounds-multiple 3` gives 3d rounds. Writes the rescaled circuits to
--output-dir with `num_sub_rounds` updated in the filename so the downstream
metadata's round count stays correct.

Used to run the hyperbolic Floquet codes for 3*d rounds instead of the fixed 16
rounds baked into the source circuits.

Example:
    python tools/scale_rounds.py \\
        --input "in/circuits/hyperbolic_floquet/8.8.8/**/*.stim" \\
        --output-dir out/circuits/hyperbolic_3d_inputs \\
        --rounds-multiple 3
"""

from __future__ import annotations

import argparse
import glob
import pathlib
import re

import stim

_RE_D = re.compile(r"_d=(\d+)")


def set_repeat_count(circuit: stim.Circuit, count: int):
    """Return (circuit with its single REPEAT block set to `count`, n_blocks)."""
    out = stim.Circuit()
    n_rep = 0
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            out += inst.body_copy() * count
            n_rep += 1
        else:
            out.append(inst)
    return out, n_rep


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--input", required=True, help="Glob for input .stim circuits.")
    ap.add_argument("--output-dir", required=True)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--rounds", type=int, help="Fixed REPEAT (round) count.")
    grp.add_argument("--rounds-multiple", type=float,
                     help="REPEAT count = round(multiple * d), d from filename.")
    args = ap.parse_args()

    files = sorted(glob.glob(args.input))
    if not files:
        ap.error(f"No files matched: {args.input}")

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for f in files:
        name = pathlib.Path(f).name
        if args.rounds is not None:
            count = args.rounds
        else:
            m = _RE_D.search(name)
            if not m:
                raise SystemExit(f"cannot parse _d= from {name}; use --rounds instead")
            count = int(round(args.rounds_multiple * int(m.group(1))))

        circuit = stim.Circuit.from_file(f)
        new, n_rep = set_repeat_count(circuit, count)
        if n_rep != 1:
            raise SystemExit(f"{name}: expected exactly 1 REPEAT block, found {n_rep}")

        # Keep num_sub_rounds consistent so the derived round count = `count`
        # (a full round is 6 sub-rounds plus the boundary round).
        new_name = re.sub(r"num_sub_rounds=\d+",
                          f"num_sub_rounds={6 * (count + 1)}", name)
        new.to_file(str(out_dir / new_name))
        print(f"  {name}: REPEAT -> {count} rounds")


if __name__ == "__main__":
    main()
