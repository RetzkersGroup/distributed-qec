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

plt.rcParams.update({"xtick.labelsize": 15, "ytick.labelsize": 15})

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
for _p in (str(_ROOT / "tools" / "_grid_scans"), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import plot_grid as g  # noqa: E402
import analyze_axis_threshold as az  # noqa: E402  (reuse the check-threshold estimators)
# analyze_axis_threshold resets tick label sizes on import; restore the summary's.
plt.rcParams.update({"xtick.labelsize": 15, "ytick.labelsize": 15})

_HONEY = str(_ROOT / "out/collection/grid_bell.csv")
_PRES_CHECK = str(_ROOT / "out/collection/pres_check.csv")   # p_z = 0 slice
_PZ_CHECK = str(_ROOT / "out/collection/pz_check.csv")       # p_res ~ 0 slice
_HYP = str(_ROOT / "out/collection/grid_hyperbolic.csv")
# Non-repeated baseline is now the single-ancilla 1-rep equivalent (== two-anc).
_BELL = "biased_bell_v3_equivalent"
_MAJ = "biased_bell_v3_equivalent_majority"
_MAJ5 = "biased_bell_v3_equivalent_majority5"

# Consistent scheme: one hue per code family (honeycomb = Blues, hyperbolic =
# Oranges, matching the bottom panel), with the *shade* encoding the repetition
# level and the same gradient direction for both codes -- darker as the number
# of repetitions increases (1-rep -> 3-rep -> 5-rep).
_blues = plt.get_cmap("Blues")
_oranges = plt.get_cmap("Oranges")
_BND_HONEY = [_blues(0.45), _blues(0.68), _blues(0.92)]     # 1-rep, 3-rep, 5-rep
_BND_HYP = [_oranges(0.45), _oranges(0.68), _oranges(0.92)]
_C_HONEY = _BND_HONEY[0]   # base-code accent, reused for the bottom-panel vline
_C_HYP = _BND_HYP[0]

# (label, csv, scheme, color, marker) for the top boundary panel.
_BOUNDARIES = [
    ("Honeycomb code", _HONEY, _BELL, _BND_HONEY[0], "o"),
    ("Honeycomb + 3-rep. majority vote", _HONEY, _MAJ, _BND_HONEY[1], "D"),
    ("Honeycomb + 5-rep. majority vote", _HONEY, _MAJ5, _BND_HONEY[2], "P"),
    ("Hyperbolic Floquet code", _HYP, _BELL, _BND_HYP[0], "s"),
    ("Hyperbolic Floquet + 3-rep. majority vote", _HYP, _MAJ, _BND_HYP[1], "v"),
    ("Hyperbolic Floquet + 5-rep. majority vote", _HYP, _MAJ5, _BND_HYP[2], "X"),
]

# (name, csv, scheme, colormap, threshold-vline color) for the bottom LER panel.
# The vline uses a darker family shade so it stays visible over the distance curves.
_CODES = [
    ("Honeycomb", _HONEY, _BELL, plt.get_cmap("Blues"), _blues(0.95)),
    ("Hyperbolic Floquet", _HYP, _BELL, plt.get_cmap("Oranges"), _oranges(0.95)),
]

# Individual boundary points with unreliable FSS fits, excluded from panel (a).
# (family, scheme, p_res) -- matched with a small tolerance. The p_z=0 intercept
# already pins the corresponding end of the boundary.
_DROP_BOUNDARY_POINTS = [
    ("honeycomb", _MAJ, 0.0010715),   # p_z* fit at p_res~0.0011 is a bad, low outlier
]


def _drop_bad_points(thr, fam, scheme, tol=1e-6):
    """Remove hand-flagged bad boundary points for this (family, scheme)."""
    return [(p, pth, err) for (p, pth, err) in thr
            if not any(f == fam and s == scheme and abs(p - pp) <= tol
                       for f, s, pp in _DROP_BOUNDARY_POINTS)]


def _monotonic_decreasing(p, pth, err):
    """Keep the sub-sequence with non-increasing p_z* as p_res grows; skip breakers.

    The physical boundary p_z*(p_res) can only decrease as residual noise rises;
    a fitted point that jumps *up* relative to the last kept point is an artifact,
    so it is dropped. Scans from low p_res (the most reliable, high-p_z* end).
    """
    order = np.argsort(p)
    keep, last = [], np.inf
    for i in order:
        if pth[i] <= last * (1 + 1e-9):
            keep.append(i)
            last = pth[i]
    keep = np.array(keep, dtype=int)
    return p[keep], pth[keep], err[keep]


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
    ap.add_argument("--method", choices=["raw", "fss", "minvar"], default="fss",
                    help="Threshold estimator: fss = the grid-pipeline fit "
                         "with the scheme-aware cutoff (default); raw = "
                         "src/threshold_fit.py FSS fit with a flat 0.35 "
                         "cutoff and no pz_max cap; minvar = plateau-relative "
                         "min-variance point.")
    args = ap.parse_args()
    g.THRESHOLD_METHOD = args.method

    cache: dict[str, "object"] = {}

    def load(csv):
        if csv not in cache:
            cache[csv] = g.load_grid(csv)
        return cache[csv]

    def load_check(csv):
        """Load a check CSV and tag each row's code family from its circuit path."""
        d = g.load_grid(csv)
        d["family"] = d["input_circuit"].fillna("").apply(
            lambda s: "honeycomb" if "honeycomb" in s
            else ("hyperbolic" if "hyperbolic" in s else "?"))
        return d

    # Axis intercepts: p_z* at p_res=0 (pz_check) and p_res* at p_z=0 (pres_check).
    try:
        chk_pz = load_check(_PZ_CHECK)      # scan p_z at p_res~0
        chk_pres = load_check(_PRES_CHECK)  # scan p_res at p_z=0
    except FileNotFoundError:
        chk_pz = chk_pres = None

    fig = plt.figure(figsize=(6.6, 8.6), constrained_layout=True)
    ax_top = fig.add_subplot(2, 1, 1)
    ax_bot = fig.add_subplot(2, 1, 2)

    # ---- (a) boundaries ------------------------------------------------------
    for label, csv, scheme, color, marker in _BOUNDARIES:
        fam = "honeycomb" if csv == _HONEY else "hyperbolic"
        sdf = load(csv)
        sdf = sdf[sdf["noise_model"] == scheme]
        thr = _drop_bad_points(g.collect_thresholds(sdf), fam, scheme)
        if not thr:
            print(f"[skip boundary] {label}: no thresholds")
            continue
        p, pth, err = map(np.array, zip(*thr))
        n_all = len(p)

        # Axis intercepts from the checks, computed with the same estimator
        # (args.method): p_z* at p_res=0 (top) and p_res* at p_z=0 (right).
        pt_pz = pt_pres = None
        if chk_pz is not None:
            r = az._fit(chk_pz[(chk_pz["family"] == fam)
                               & (chk_pz["noise_model"] == scheme)],
                        "p_bell_z", scheme, args.method)
            pt_pz = r[0] if r is not None else None
        if chk_pres is not None:
            r = az._fit(chk_pres[(chk_pres["family"] == fam)
                                 & (chk_pres["noise_model"] == scheme)],
                        "p", scheme, args.method)
            pt_pres = r[0] if r is not None else None

        # Prepend the p_res=0 intercept (the true top of the boundary) *before*
        # enforcing monotonicity, so grid points that over-fit above it are
        # dropped. Then append the p_z=0 intercept as the low-p_z endpoint.
        pl, tl, el = list(p), list(pth), list(err)
        if pt_pz is not None:
            pl, tl, el = [0.0] + pl, [pt_pz] + tl, [0.0] + el
        p, pth, err = _monotonic_decreasing(np.array(pl), np.array(tl), np.array(el))
        n_mono = len(p)
        if pt_pres is not None:
            p = np.append(p, pt_pres); pth = np.append(pth, 0.0); err = np.append(err, 0.0)

        ax_top.errorbar(p, pth, yerr=err, fmt=marker + "-", color=color,
                        label=label, markersize=5, linewidth=1.7, capsize=3)
        # Highlight the two check-derived intercept points.
        for xi, yi in ((0.0, pt_pz), (pt_pres, 0.0)):
            if yi is not None and xi is not None:
                ax_top.scatter([xi], [yi], marker="*", s=130, color=color,
                               edgecolors="k", linewidths=0.5, zorder=6)
        print(f"{label}: {n_mono}/{n_all} grid thresholds (monotonic); "
              f"intercepts p_z*(p_res=0)={pt_pz}, p_res*(p_z=0)={pt_pres}")
    if args.logy_boundary:
        ax_top.set_yscale("log")
    ax_top.set_xlabel(r"$p^*_{\text{res}}$", fontsize=17)
    ax_top.set_ylabel(r"$p_z^*$ (threshold)", fontsize=17)
    ax_top.grid(True, alpha=0.3)
    ax_top.legend(fontsize=8.5, framealpha=0.9)

    # ---- (b) per-observable LER vs p_z at fixed p_res ------------------------
    handles = []
    for cname, csv, scheme, cmap, vcolor in _CODES:
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
        # Threshold line for this code, colored to match its family, with its
        # value in the legend so the two are easy to tell apart.
        res = g.estimate_pz_threshold(sdf, p)
        if res is not None:
            pth = res[0]
            ax_bot.axvline(pth, color=vcolor, linestyle="--", linewidth=1.6)
            handles.append(Line2D([0], [0], color=vcolor, linestyle="--",
                                  linewidth=1.6,
                                  label=rf"{cname} threshold: $p_z = {pth:.2g}$"))

    ax_bot.set_xscale("log")
    ax_bot.set_yscale("log")
    ax_bot.set_xlabel(r"$p_z$", fontsize=17)
    ax_bot.set_ylabel("Logical error rate per observable", fontsize=17)
    ax_bot.set_title(rf"$p_{{\text{{res}}}} = {args.p_example:.2g}$", fontsize=16)
    ax_bot.grid(True, which="both", alpha=0.3)
    ax_bot.legend(handles=handles, fontsize=7.5, loc="lower right", framealpha=0.9)

    for ax, lab in [(ax_top, "(a)"), (ax_bot, "(b)")]:
        ax.text(-0.02, 1.02, lab, transform=ax.transAxes, fontsize=19,
                fontweight="bold", va="bottom", ha="left")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {args.out}")


if __name__ == "__main__":
    main()
