#!/usr/bin/env python3
"""Tesseract vs pymatching on the original Gidney planar-honeycomb circuits.

Consumes the CSVs written by decode_array.pbs
(out/collection/honeycomb_tesseract_debug*.csv — dated append-mode files, all
globbed and aggregated) and draws one LER-vs-p panel per observable (H/V):
color = distance, solid = pymatching, dashed = tesseract_short. A summary
table prints the per-(p, d) ratio of logical error rates.

Rows are aggregated across CSVs/resubmissions by strong_id (shots and errors
sum), so rerunning the PBS array just tightens the error bars.

Usage:
    python tools/honeycomb_tesseract_debug/plot_comparison.py
    python tools/honeycomb_tesseract_debug/plot_comparison.py --per-round
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
for _p in (str(_ROOT / "tools" / "_grid_scans"), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import plot_grid as g  # noqa: E402  (load_grid: strong_id aggregation, pL)

_CSV_GLOB = str(_ROOT / "out/collection/honeycomb_tesseract_debug*.csv")

_DECODERS = [("pymatching", "-"), ("tesseract_short", "--")]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=_CSV_GLOB,
                    help="CSV glob (default: honeycomb_tesseract_debug*.csv).")
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures" / "honeycomb_tesseract_debug.pdf")
    ap.add_argument("--per-round", action="store_true",
                    help="Plot per-round logical error rate (1-(1-pL)^(1/rounds)).")
    args = ap.parse_args()

    df = g.load_grid(args.csv)
    if "observable" not in df.columns:
        df["observable"] = df["json_metadata"].apply(
            lambda m: g._parse_meta(m).get("observable", "?"))

    rate_col = "pL"
    if args.per_round:
        df["pL_round"] = 1 - (1 - df["pL"]) ** (1 / df["rounds"])
        rate_col = "pL_round"

    observables = sorted(df["observable"].dropna().unique())
    dists = sorted(int(d) for d in df["expected_distance"].dropna().unique())

    fig, axes = plt.subplots(1, len(observables),
                             figsize=(6.2 * len(observables), 5.0),
                             constrained_layout=True, squeeze=False)
    for c, obs in enumerate(observables):
        ax = axes[0][c]
        for dec, ls in _DECODERS:
            for d in dists:
                sub = df[(df["observable"] == obs) & (df["decoder"] == dec)
                         & (df["expected_distance"] == d)
                         & (df["errors"] > 0) & (df["shots"] > 0)].sort_values("p")
                if len(sub) == 0:
                    continue
                ax.errorbar(sub["p"], sub[rate_col], yerr=sub["sigma"],
                            marker=g.MARKERS.get(d, "o"), linestyle=ls,
                            color=g.COLORS.get(d), markersize=4,
                            linewidth=1.4, capsize=2, elinewidth=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"physical error rate $p$ (EM3_v3)", fontsize=13)
        ax.set_ylabel(("Per-round " if args.per_round else "")
                      + "logical error rate per observable", fontsize=12)
        ax.set_title(f"obs = {obs}", fontsize=13)
        ax.grid(True, which="both", alpha=0.3)
        handles = (
            [Line2D([0], [0], color=g.COLORS.get(d), marker=g.MARKERS.get(d, "o"),
                    linestyle="-", label=f"d={d}") for d in dists]
            + [Line2D([0], [0], color="k", linestyle=ls, label=dec)
               for dec, ls in _DECODERS]
        )
        ax.legend(handles=handles, fontsize=9, framealpha=0.9)

    fig.suptitle("Gidney planar honeycomb (EM3_v3, circuits as published): "
                 "tesseract vs pymatching", fontsize=14)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)

    # Ratio table: pymatching / tesseract per (p, d), pooled over observables.
    print("\n  LER ratio pymatching / tesseract  (>1 -> tesseract better)")
    print("  " + "-" * 58)
    pool = df[(df["errors"] > 0)].groupby(["decoder", "p", "expected_distance"]).agg(
        errors=("errors", "sum"), shots=("shots", "sum")).reset_index()
    pool["rate"] = pool["errors"] / pool["shots"]
    pm = pool[pool["decoder"] == "pymatching"].set_index(["p", "expected_distance"])
    ts = pool[pool["decoder"] == "tesseract_short"].set_index(["p", "expected_distance"])
    common = pm.index.intersection(ts.index)
    if len(common) == 0:
        print("  (no overlapping (p, d) points with errors yet)")
    for p, d in sorted(common):
        r_pm, r_ts = pm.loc[(p, d), "rate"], ts.loc[(p, d), "rate"]
        ratio = r_pm / r_ts if r_ts > 0 else np.inf
        print(f"  p={p:<8g} d={int(d):<3d} pm={r_pm:9.3e}  tess={r_ts:9.3e}"
              f"  ratio={ratio:6.2f}")
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
