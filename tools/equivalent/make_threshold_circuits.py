#!/usr/bin/env python3
"""Generate equivalent-scheme circuits for a (p, p_bell_z) threshold scan.

Thin wrapper around `tools/make_threshold_circuits_bell.py` that defaults the
translation scheme to `equivalent_v3` (the single-ancilla XCX gadget with the
correlated `E(...)` noise block; see equivalent_bell_bias.md) and the input to
the Gidney honeycomb circuits. Pass `--scheme equivalent_v3_majority` for the
majority-vote variant; the default output directory follows the scheme.

All other flags (`--p-values`, `--biases`, `--p-bell-z-values`, `--workers`,
`--skip-existing`, `--noiseless-start-end`, ...) are forwarded unchanged to the
underlying tool. A scan direction is either a bias `beta` (-> p_bell_z = beta*p)
or a fixed p_bell_z value.

Examples:

    # Default input (Gidney honeycomb), default directions and p-values:
    python tools/equivalent/make_threshold_circuits.py

    # Majority-vote variant:
    python tools/equivalent/make_threshold_circuits.py --scheme equivalent_v3_majority

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
DEFAULT_SCHEME = "equivalent_v3"
# Output directory per scheme (used only when --output-dir is not supplied).
OUTPUT_DIR_BY_SCHEME = {
    "equivalent_v3": "out/circuits/thresholds_equivalent",
    "equivalent_v3_majority": "out/circuits/thresholds_equivalent_majority",
}


def _flag_value(argv: list[str], flag: str) -> str | None:
    """Return the value of `--flag X` or `--flag=X` in argv, else None."""
    for i, a in enumerate(argv):
        if a == flag and i + 1 < len(argv):
            return argv[i + 1]
        if a.startswith(flag + "="):
            return a.split("=", 1)[1]
    return None


def main() -> None:
    argv = sys.argv[1:]

    # Resolve scheme (default equivalent_v3); the underlying tool validates it.
    scheme = _flag_value(argv, "--scheme") or DEFAULT_SCHEME

    # Equivalent-specific defaults, injected only when the user didn't pass them.
    if not any(a == "--input" or a.startswith("--input=") for a in argv):
        argv += ["--input", DEFAULT_INPUT]
    if not any(a == "--output-dir" or a.startswith("--output-dir=") for a in argv):
        out_dir = OUTPUT_DIR_BY_SCHEME.get(scheme, "out/circuits/thresholds_equivalent")
        argv += ["--output-dir", out_dir]
    if _flag_value(argv, "--scheme") is None:
        argv += ["--scheme", scheme]

    sys.argv = [sys.argv[0], *argv]
    base.main()


if __name__ == "__main__":
    main()
