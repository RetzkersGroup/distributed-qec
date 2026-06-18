#!/usr/bin/env python3
"""Generate flexible_pairwise-noised circuits for a 3D threshold sweep.

For each input circuit, strips the existing noise and re-applies
`cultiv_noise.NoiseModel.flexible_pairwise` for every (direction, p) pair.
A direction is a tuple of ratios (r_m_bulk, r_spam, r_res); the
absolute noise rates at sweep value `p` are `p * r_*`.

Default directions (10) cover the named profiles plus the interior of the
ratio simplex:

    measurement-limited:   (1,    1e-6, 1e-6)
    spam-limited:          (1e-6, 1,    1e-6)
    res-limited:           (1e-6, 1e-6, 1   )
    balanced:              (1,    1,    1   )
    six interior mixes from integer compositions a_m+a_s+a_d=6 with min ratio 1,
    each normalized so max(ratio)=1.

Examples:

    python tools/make_threshold_circuits.py \\
        --input "in/circuits/honeycomb_gidney/model=SDEM3,*obs=H*.stim" \\
        --output-dir out/circuits/thresholds
"""

from __future__ import annotations

import argparse
import glob
import json
import multiprocessing as mp
import pathlib
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import stim

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import _cultiv_noise as cultiv_noise  # noqa: E402


SMALL = 1e-6


def _default_directions() -> list[tuple[str, tuple[float, float, float]]]:
    """Named corners + balanced + six interior mixes (max-normalized)."""
    dirs: list[tuple[str, tuple[float, float, float]]] = [
        ("measurement-limited", (1.0, SMALL, SMALL)),
        ("spam-limited", (SMALL, 1.0, SMALL)),
        ("res-limited", (SMALL, SMALL, 1.0)),
        ("balanced", (1.0, 1.0, 1.0)),
    ]
    # Integer compositions of 6 with all parts >= 1 and not all equal.
    # Take the (3,2,1)-type family: 6 permutations.
    mixes = [
        (3, 2, 1), (3, 1, 2),
        (2, 3, 1), (1, 3, 2),
        (2, 1, 3), (1, 2, 3),
    ]
    for a_m, a_s, a_d in mixes:
        m = max(a_m, a_s, a_d)
        ratios = (a_m / m, a_s / m, a_d / m)
        name = f"mix_{a_m}{a_s}{a_d}"
        dirs.append((name, ratios))

    # Edge directions on the three coordinate planes -- one ratio set to the
    # SMALL floor, the other two from a (2,1) composition.
    edges = [
        (2, 1, 0), (1, 2, 0),   # p_res ~= 0 plane
        (2, 0, 1), (1, 0, 2),   # p_spam ~= 0 plane
        (0, 2, 1), (0, 1, 2),   # p_m_bulk ~= 0 plane
    ]
    for a_m, a_s, a_d in edges:
        parts = (a_m, a_s, a_d)
        m = max(parts)
        ratios = tuple(a / m if a > 0 else SMALL for a in parts)
        name = f"edge_{a_m}{a_s}{a_d}"
        dirs.append((name, ratios))
    return dirs


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


def _format_ratio(r: float) -> str:
    """Filename-friendly ratio token; collapse small floor to 'e'."""
    if r <= SMALL * 10:
        return "e"
    if r == 1.0:
        return "1"
    # 2 sig figs, strip trailing zeros
    s = f"{r:.3g}"
    return s.replace(".", "p")


def _process(
    in_path: pathlib.Path,
    direction_name: str,
    ratios: tuple[float, float, float],
    p: float,
    out_dir: pathlib.Path,
    skip_existing: bool,
) -> str:
    r_m, r_s, r_d = ratios
    p_m_bulk = p * r_m
    p_spam = p * r_s
    p_res = p * r_d

    rm_tok = _format_ratio(r_m)
    rs_tok = _format_ratio(r_s)
    rd_tok = _format_ratio(r_d)
    out_stem = (
        f"{in_path.stem}__dir={direction_name},rm={rm_tok},rs={rs_tok},rd={rd_tok},p={p:.5g}"
    )
    out_circuit_path = out_dir / f"{out_stem}.stim"
    out_meta_path = out_dir / f"{out_stem}.json"

    if skip_existing and out_circuit_path.exists() and out_meta_path.exists():
        return f"  {direction_name:>22s}  p={p:.4g} -> {out_circuit_path.name} (skipped, exists)"

    raw = stim.Circuit.from_file(str(in_path)).without_noise()
    model = cultiv_noise.NoiseModel.flexible_pairwise(
        p_m_bulk=p_m_bulk,
        p_spam=p_spam,
        p_res=p_res,
    )
    noisy = model.noisy_circuit(raw)

    noisy.to_file(str(out_circuit_path))

    src_meta = _load_metadata(in_path)
    metadata = {
        "p": float(p),
        "direction": direction_name,
        "ratio_m_bulk": float(r_m),
        "ratio_spam": float(r_s),
        "ratio_res": float(r_d),
        "p_m_bulk": float(p_m_bulk),
        "p_spam": float(p_spam),
        "p_res": float(p_res),
        "noise_model": "flexible_pairwise",
        "code_family": src_meta.get("code_family"),
        "code_name": src_meta.get("code_name"),
        "k": src_meta.get("k"),
        "expected_distance": src_meta.get("expected_distance"),
        "rounds": src_meta.get("rounds"),
        "w": src_meta.get("w"),
        "h": src_meta.get("h"),
        "input_circuit": str(in_path),
        "output_circuit": str(out_circuit_path),
        "num_qubits": noisy.num_qubits,
        "num_measurements": noisy.num_measurements,
        "num_detectors": noisy.num_detectors,
        "num_observables": noisy.num_observables,
    }

    with out_meta_path.open("w") as f:
        json.dump(metadata, f, indent=2)

    return f"  {direction_name:>22s}  p={p:.4g} -> {out_circuit_path.name}"


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
        help="Physical error rates (default: np.logspace(-2.5, -1, 25)).",
    )
    parser.add_argument(
        "--output-dir", default="out/circuits/thresholds",
        help="Output directory (default: out/circuits/thresholds).",
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip (direction, p) tasks whose .stim and .json already exist.",
    )
    parser.add_argument(
        "--workers", type=int, default=mp.cpu_count(),
        help="Parallel worker processes (default: cpu_count()).",
    )
    args = parser.parse_args()

    inputs = _resolve_inputs(args.input)
    if not inputs:
        parser.error(f"No .stim files matched: {args.input}")

    if args.p_values is None:
        p_values = np.logspace(-2.5, -1, 25).tolist()
    else:
        p_values = list(args.p_values)

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    directions = _default_directions()
    jobs = [
        (in_path, direction_name, ratios, p)
        for in_path in inputs
        for direction_name, ratios in directions
        for p in p_values
    ]
    n_workers = max(1, min(args.workers, len(jobs)))
    print(
        f"{len(inputs)} circuits x {len(directions)} directions x "
        f"{len(p_values)} p-values = {len(jobs)} circuits to write "
        f"-> {out_dir} (workers={n_workers})"
    )

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [
            ex.submit(_process, in_path, direction_name, ratios, p, out_dir, args.skip_existing)
            for in_path, direction_name, ratios, p in jobs
        ]
        for done, fut in enumerate(as_completed(futures), 1):
            msg = fut.result()
            print(f"[{done}/{len(jobs)}]{msg}")


if __name__ == "__main__":
    main()
