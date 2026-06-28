#!/usr/bin/env python3
"""Generate biased_bell_v3_equivalent circuits for a (p, p_bell_z) threshold scan.

Thin wrapper around `tools/make_threshold_circuits_bell.py` that pins the
translation scheme to `equivalent_v3` (the single-ancilla XCX gadget with the
correlated `E(...)` noise block; see equivalent_bell_bias.md) and defaults the
input to the Gidney honeycomb circuits.

All other flags (`--p-values`, `--biases`, `--p-bell-z-values`, `--workers`,
`--skip-existing`, `--noiseless-start-end`, ...) are forwarded unchanged to the
underlying tool. A scan direction is either a bias `beta` (-> p_bell_z = beta*p)
or a fixed p_bell_z value.

Examples:

    # Default input (Gidney honeycomb), default directions and p-values:
    python tools/equivalent/make_threshold_circuits.py

    # Custom sweep:
    python tools/equivalent/make_threshold_circuits.py \
        --biases 1 10 100 \
        --p-bell-z-values 0.001 \
        --p-values 0.001 0.003 0.005 0.01
"""

from __future__ import annotations

import pathlib
import sys

_TOOLS = pathlib.Path(__file__).resolve().parent.parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))

import make_threshold_circuits_bell as base  # noqa: E402

DEFAULT_INPUT = "in/circuits/honeycomb_gidney/*.stim"
DEFAULT_OUTPUT_DIR = "out/circuits/thresholds_equivalent"
SCHEME = "equivalent_v3"


def main() -> None:
    argv = sys.argv[1:]

    # Equivalent-specific defaults, injected only when the user didn't pass them.
    if not any(a == "--input" or a.startswith("--input=") for a in argv):
        argv += ["--input", DEFAULT_INPUT]
    if not any(a == "--output-dir" or a.startswith("--output-dir=") for a in argv):
        argv += ["--output-dir", DEFAULT_OUTPUT_DIR]
    # The scheme is fixed for this wrapper; drop any user-supplied --scheme
    # (both "--scheme X" and "--scheme=X" forms) and append our own.
    cleaned: list[str] = []
    skip_next = False
    for a in argv:
        if skip_next:
            skip_next = False
            continue
        if a == "--scheme":
            skip_next = True  # also drop the following value token
            continue
        if a.startswith("--scheme="):
            continue
        cleaned.append(a)
    argv = cleaned + ["--scheme", SCHEME]

    sys.argv = [sys.argv[0], *argv]
    base.main()


if __name__ == "__main__":
    main()
