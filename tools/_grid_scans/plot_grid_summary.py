#!/usr/bin/env python3
"""Summary figure comparing the honeycomb and hyperbolic Floquet codes.

  (a) top    -- threshold boundaries p_z* vs p_res for: honeycomb code,
                honeycomb + 3-repetition majority vote, and hyperbolic Floquet.
  (b) bottom -- per-observable logical error rate vs p_z at a fixed p_res (1e-3)
                for the two base codes (no repetition), the two codes drawn in
                distinct color families.

Honeycomb is read from grid_bell, hyperbolic from grid_hyperbolic. Reuses the
load/fit helpers from plot_grid.py (same saturation cutoff etc.).
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

import plot_grid as g  # noqa: E402

_HONEY = str(_ROOT / "out/collection/grid_bell.csv")
_HYP = str(_ROOT / "out/collection/grid_hyperbolic.csv")
_BELL = "biased_bell_v3"
_MAJ = "biased_bell_v3_equivalent_majority"

# Professional, distinct hues. Honeycomb = blue family, hyperbolic = orange
# family (matched between the two panels); majority vote = green.
_C_HONEY = "#4e79a7"
_C_MAJ = "#59a14f"
_C_HYP = "#f28e2b"

# (label, csv, scheme, color, marker) for the top boundary panel.
_BOUNDARIES = [
    ("Honeycomb code", _HONEY, _BELL, _C_HONEY, "o"),
    ("Honeycomb + 3-rep. majority vote", _HONEY, _MAJ, _C_MAJ, "D"),
    ("Hyperbolic Floquet code", _HYP, _BELL, _C_HYP, "s"),
]

# (name, csv, scheme, colormap) for the bottom LER panel base codes.
_CODES = [
    ("Honeycomb", _HONEY, _BELL, plt.get_cmap("Blues")),
    ("Hyperbolic Floquet", _HYP, _BELL, plt.get_cmap("Oranges")),
]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--p-example", type=float, default=1e-3,
                    help="p_res for the bottom LER panel (default 1e-3).")
    ap.add_argument("--logy-boundary", action="store_true",
                    help="Use a log y-axis for the boundary panel.")
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures" / "grid_summary.pdf")
    args = ap.parse_args()

    cache: dict[str, "object"] = {}

    def load(csv):
        if csv not in cache:
            cache[csv] = g.load_grid(csv)
        return cache[csv]

    fig = plt.figure(figsize=(6.6, 8.6), constrained_layout=True)
    ax_top = fig.add_subplot(2, 1, 1)
    ax_bot = fig.add_subplot(2, 1, 2)

    # ---- (a) boundaries ------------------------------------------------------
    for label, csv, scheme, color, marker in _BOUNDARIES:
        sdf = load(csv)
        sdf = sdf[sdf["noise_model"] == scheme]
        thr = g.collect_thresholds(sdf)
        if not thr:
            print(f"[skip boundary] {label}: no thresholds")
            continue
        p, pth, err = map(np.array, zip(*thr))
        ax_top.errorbar(p, pth, yerr=err, fmt=marker + "-", color=color,
                        label=label, markersize=5, linewidth=1.7, capsize=3)
        print(f"{label}: {len(thr)} thresholds")
    if args.logy_boundary:
        ax_top.set_yscale("log")
    ax_top.set_xlabel(r"$p_{\text{res}}$", fontsize=11)
    ax_top.set_ylabel(r"$p_z^{*}$ (threshold)", fontsize=11)
    ax_top.grid(True, alpha=0.3)
    ax_top.legend(fontsize=8.5, framealpha=0.9)

    # ---- (b) per-observable LER vs p_z at fixed p_res ------------------------
    handles = []
    for cname, csv, scheme, cmap in _CODES:
        sdf = load(csv)
        sdf = sdf[sdf["noise_model"] == scheme]
        p_vals = sorted(sdf["p"].dropna().unique())
        p = min(p_vals, key=lambda v: abs(v - args.p_example))
        dists = sorted(int(d) for d in sdf["expected_distance"].dropna().unique())
        shades = np.linspace(0.45, 0.9, len(dists))
        for d, sh in zip(dists, shades):
            color = cmap(sh)
            xs, ys, es = g._pz_curve(sdf, p, d, "pL", "sigma")
            if len(xs) == 0:
                continue
            ax_bot.errorbar(xs, ys, yerr=es, marker="o", linestyle="-",
                            color=color, markersize=4, linewidth=1.4,
                            capsize=2, elinewidth=0.7)
            handles.append(Line2D([0], [0], color=color, marker="o",
                                  linestyle="-",
                                  label=f"{cname} {g._code_label(sdf, d)}"))
        fit = g.fit_pz_threshold(sdf, p)
        if fit is not None:
            ax_bot.axvline(fit.p_th, color=cmap(0.98), linestyle="--", linewidth=1.3)

    ax_bot.set_xscale("log")
    ax_bot.set_yscale("log")
    ax_bot.set_xlabel(r"$p_z$", fontsize=11)
    ax_bot.set_ylabel("Logical error rate per observable", fontsize=11)
    ax_bot.set_title(rf"$p_{{\text{{res}}}} = {args.p_example:.2g}$", fontsize=10)
    ax_bot.grid(True, which="both", alpha=0.3)
    ax_bot.legend(handles=handles, fontsize=7.5, loc="lower right", framealpha=0.9)

    for ax, lab in [(ax_top, "(a)"), (ax_bot, "(b)")]:
        ax.text(-0.02, 1.02, lab, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {args.out}")


if __name__ == "__main__":
    main()
