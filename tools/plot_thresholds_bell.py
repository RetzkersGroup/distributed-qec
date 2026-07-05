#!/usr/bin/env python3
"""Fit and plot Bell-mediated threshold data, saving every figure to out/figures/.

Non-interactive counterpart to the Bell threshold notebooks. Once the decode
CSVs have been collected (via tools/run_decode_thresholds_bell.sh and
tools/equivalent/run_decode*.sh), this reproduces the notebook figures in one
command and writes them to disk with the same
`savefig(..., bbox_inches='tight')` convention the notebooks use.

Stages (each mirrors a notebook):
  1. Load `threshold_bell,*.csv`  -> fit v2/v3 thresholds (cutoff 0.30).
  2. Load `threshold_equivalent*,*.csv` -> fit the equivalents (cutoff 0.45).
  3. Concat fits, write `threshold_bell_fits.csv`.
  4. Save four figures to out/figures/.

Examples:

    python tools/plot_thresholds_bell.py
    python tools/plot_thresholds_bell.py --figure-format png --out-dir out/figures
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")  # headless: no display needed to save figures.

import pandas as pd

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from bell_threshold_analysis import (  # noqa: E402
    load_decode_csvs,
    fit_thresholds,
    plot_direction_curves,
    plot_threshold_plane,
    plot_threshold_vs_fixed_pbz,
    plot_combined_plane,
)

# Bell v2/v3 raw decode CSVs and the noise models they carry.
_BELL_GLOB = "threshold_bell,*.csv"
_BELL_MODELS = ["biased_bell_v2", "biased_bell_v3"]
_BELL_CUTOFF = 0.30

# Single-ancilla equivalent decode CSVs (fitted with a looser saturation cut).
_EQUIV_GLOBS = ["threshold_equivalent,*.csv", "threshold_equivalent_majority,*.csv"]
_EQUIV_CUTOFF = 0.45


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--collection-dir", type=pathlib.Path,
        default=_ROOT / "out" / "collection",
        help="Directory holding the decode CSVs (default: out/collection).",
    )
    parser.add_argument(
        "--out-dir", type=pathlib.Path, default=_ROOT / "out" / "figures",
        help="Directory to write figures into (default: out/figures).",
    )
    parser.add_argument(
        "--figure-format", default="pdf",
        help="Figure file extension / matplotlib format (default: pdf).",
    )
    parser.add_argument(
        "--n-per-side", type=int, default=2,
        help="Half-width of the FSS fit window (default: 2 -> 5 p-values).",
    )
    args = parser.parse_args()

    coll = args.collection_dir
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ext = args.figure_format.lstrip(".")

    def _save(fig, name: str) -> None:
        if fig is None:
            print(f"  [skip] {name}: nothing to plot.")
            return
        path = args.out_dir / f"{name}.{ext}"
        fig.savefig(path, bbox_inches="tight")
        print(f"  saved {path}")

    # --- Stage 1: v2/v3 curves + fits ------------------------------------
    bell_df = load_decode_csvs(str(coll / _BELL_GLOB), _BELL_MODELS)
    bell_fits = fit_thresholds(
        bell_df, n_per_side=args.n_per_side, logical_cutoff=_BELL_CUTOFF,
        verbose=True,
    )
    print(f"Fitted {len(bell_fits)} v2/v3 directions.")

    # --- Stage 2: equivalents (optional) ---------------------------------
    equiv_fits = pd.DataFrame()
    try:
        equiv_df = load_decode_csvs([str(coll / g) for g in _EQUIV_GLOBS])
    except FileNotFoundError:
        print("No equivalent CSVs found; plotting v2/v3 only.")
    else:
        equiv_fits = fit_thresholds(
            equiv_df, n_per_side=args.n_per_side,
            logical_cutoff=_EQUIV_CUTOFF, verbose=True,
        )
        print(f"Fitted {len(equiv_fits)} equivalent-model directions.")

    fits_df = pd.concat([bell_fits, equiv_fits], ignore_index=True)

    # --- Stage 3: persist fits (same file the notebooks share) -----------
    fits_path = coll / "threshold_bell_fits.csv"
    fits_path.parent.mkdir(parents=True, exist_ok=True)
    fits_df.to_csv(fits_path, index=False)
    print(f"Saved fits -> {fits_path}")

    # --- Stage 4: figures ------------------------------------------------
    print(f"Writing figures to {args.out_dir}:")
    _save(plot_direction_curves(bell_df, bell_fits, noise_models=_BELL_MODELS),
          "threshold_bell_direction_curves")
    _save(plot_threshold_plane(fits_df), "threshold_bell_curves")
    _save(plot_threshold_vs_fixed_pbz(fits_df), "threshold_bell_vs_fixed_pbz")
    _save(plot_combined_plane(fits_df), "threshold_bell_combined_plane")


if __name__ == "__main__":
    main()
