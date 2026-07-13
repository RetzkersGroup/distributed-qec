#!/usr/bin/env python3
"""Overlay threshold boundaries (p_z* vs p_res) from several grid CSVs.

Each --source is ``LABEL=CSV_GLOB``; the chosen scheme's boundary is fitted from
each source and drawn on one axes, so e.g. the honeycomb and hyperbolic Floquet
codes can be compared on the same (p_res, p_z) plane. Reuses load_grid /
collect_thresholds from plot_grid.py (same fit, same saturation cutoff).

Example:
    python tools/plot_grid_boundary.py \\
        --source "Honeycomb code=out/collection/grid_bell,*.csv" \\
        --source "Hyperbolic Floquet code=out/collection/grid_hyperbolic,*.csv" \\
        --out out/figures/grid_combined_threshold.png
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
for _p in (str(_ROOT / "tools" / "_grid_scans"), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import plot_grid as g  # noqa: E402

_PALETTE = ["#4e79a7", "#e15759", "#59a14f", "#f28e2b", "#b07aa1", "#76b7b2"]
_MARKERS = ["o", "s", "^", "D", "v", "P"]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--source", action="append", required=True,
                    help="LABEL=CSV_GLOB (repeatable).")
    ap.add_argument("--scheme", default="biased_bell_v3_equivalent",
                    help="noise_model to fit from each source "
                         "(default biased_bell_v3_equivalent, the 1-rep baseline).")
    ap.add_argument("--raw", action="store_true",
                    help="Fit on the raw rate instead of per-observable.")
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures" / "grid_combined_threshold.png")
    ap.add_argument("--method", choices=["minvar", "fss"], default="fss",
                    help="Threshold estimator (default fss).")
    args = ap.parse_args()
    g.THRESHOLD_METHOD = args.method

    rate_col, sigma_col = ("pL_raw", "sigma_raw") if args.raw else ("pL", "sigma")

    fig, ax = plt.subplots(figsize=(6.5, 5))
    for i, src in enumerate(args.source):
        if "=" not in src:
            ap.error(f"--source must be LABEL=CSV_GLOB; got {src!r}")
        label, csv = src.split("=", 1)
        df = g.load_grid(csv)
        sdf = df[df["noise_model"] == args.scheme]
        if sdf.empty:
            print(f"[skip] {label}: no rows with noise_model={args.scheme}")
            continue
        thr = g.collect_thresholds(sdf, rate_col=rate_col, sigma_col=sigma_col)
        if not thr:
            print(f"[skip] {label}: no thresholds fitted")
            continue
        p, pth, err = map(np.array, zip(*thr))
        ax.errorbar(p, pth, yerr=err, fmt=_MARKERS[i % len(_MARKERS)] + "-",
                    color=_PALETTE[i % len(_PALETTE)], label=label,
                    markersize=6, linewidth=1.6, capsize=3)
        print(f"{label}: {len(thr)} thresholds, "
              f"p_res {p.min():.3g}..{p.max():.3g}, p_z* {pth.min():.3g}..{pth.max():.3g}")

    ax.set_xlabel(r"$p_{\text{res}}$", fontsize=11)
    ax.set_ylabel(r"$p_z^*$ (threshold)", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {args.out}")


if __name__ == "__main__":
    main()
