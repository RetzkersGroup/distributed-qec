#!/usr/bin/env python3
"""Per-mechanism noise-probability comparison at p_z = 0, scanning p_res.

For every effective noise mechanism (a flip pattern over {M, X0, X1, Z0, Z1})
produced by any of the three deployed models -- the two-ancilla Bell gadget
(biased_bell_v3), the single-repetition equivalent, and the 3-repetition
majority-vote equivalent -- draw one subplot of its probability vs p_res, with
p_z fixed to 0. Each subplot overlays five curves:

    two-anc, 1-rep, 3-rep            (the models themselves)
    1-rep x3, two-anc x3             (the naive additive 3-round expectation)

so it is easy to see which mechanisms the 3-rep model simply triples (data
errors: 3-rep ~ 3 x 1-rep) versus quadratically suppresses (the measurement
flip {M}, where 3-rep << 3 x 1-rep thanks to majority voting).

All curves are analytic (from equivalent_noise_calculator), so no sampling is
needed.

Usage:
    python tools/_grid_scans/plot_noise_mechanisms.py
    python tools/_grid_scans/plot_noise_mechanisms.py --p-lo 1e-4 --p-hi 3e-2 -n 60
"""

from __future__ import annotations

import argparse
import collections
import math
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

import equivalent_noise_calculator as e  # noqa: E402

_C_TWO = "#4e79a7"    # two-ancilla Bell
_C_ONE = "#59a14f"    # 1-rep equivalent
_C_THREE = "#e15759"  # 3-rep majority vote
_C_FIVE = "#674ea7"   # 5-rep majority vote

# (label, function, color) -- the base models.
_MODELS = [
    ("two-anc", e.get_two_anc_noise_operations, _C_TWO),
    ("1-rep", e.get_equivalent_noise_operations_single, _C_ONE),
    ("3-rep", e.get_equivalent_noise_operations, _C_THREE),
    ("5-rep", lambda pe, pz: e.get_equivalent_noise_operations(pe, pz, 5), _C_FIVE),
]


def _mech_dict(fn, p, p_bell_z=0.0) -> dict:
    d: dict[frozenset, float] = collections.defaultdict(float)
    for pr, flips in fn(p, p_bell_z):
        d[frozenset(flips)] += pr
    return d


def _pattern_label(pat: frozenset) -> str:
    return "{" + ", ".join(sorted(pat)) + "}"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--p-lo", type=float, default=1e-4, help="min p_res (default 1e-4).")
    ap.add_argument("--p-hi", type=float, default=3e-2, help="max p_res (default 3e-2).")
    ap.add_argument("-n", "--n-points", type=int, default=60,
                    help="p_res grid points (default 60).")
    ap.add_argument("--floor", type=float, default=1e-12,
                    help="drop mechanisms whose max prob stays below this.")
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures" / "noise_mechanisms.pdf")
    args = ap.parse_args()

    ps = np.logspace(math.log10(args.p_lo), math.log10(args.p_hi), args.n_points)

    # probs[model][pattern] = probability array over ps.
    probs: dict[str, dict[frozenset, np.ndarray]] = {name: {} for name, _, _ in _MODELS}
    patterns: set[frozenset] = set()
    for name, fn, _ in _MODELS:
        series: dict[frozenset, np.ndarray] = collections.defaultdict(
            lambda: np.zeros(len(ps)))
        for i, p in enumerate(ps):
            for pat, pr in _mech_dict(fn, float(p)).items():
                series[pat][i] = pr
        probs[name] = dict(series)
        patterns |= set(series)

    def _max_over_models(pat):
        return max(probs[n].get(pat, np.zeros(1)).max() for n, _, _ in _MODELS)

    patterns = [pat for pat in patterns if _max_over_models(pat) > args.floor]
    # Order: fewest labels first, then most probable, so the leading mechanisms
    # (M, single data errors) come first.
    patterns.sort(key=lambda pat: (len(pat), -_max_over_models(pat)))

    n = len(patterns)
    ncol = min(5, n)
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.9 * nrow),
                             squeeze=False, constrained_layout=True)

    def _mask(a):  # zeros are invisible on a log axis; make them gaps
        return np.where(a > 0, a, np.nan)

    ymin = max(args.floor, 1e-12)
    for idx, pat in enumerate(patterns):
        ax = axes[idx // ncol][idx % ncol]
        for name, _, color in _MODELS:
            y = probs[name].get(pat)
            if y is not None:
                ax.plot(ps, _mask(y), color=color, linewidth=1.9, label=name)
        # Additive-expectation overlays: 1-rep scaled by the repetition count,
        # to compare against the 3-rep / 5-rep curves (data errors ~ x R; the
        # measurement flip {M} falls far below, thanks to majority voting).
        y1 = probs["1-rep"].get(pat)
        if y1 is not None:
            ax.plot(ps, _mask(3.0 * y1), color=_C_THREE, linestyle=":",
                    linewidth=1.4, label="1-rep x3")
            ax.plot(ps, _mask(5.0 * y1), color=_C_FIVE, linestyle=":",
                    linewidth=1.4, label="1-rep x5")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(bottom=ymin)
        ax.set_title(_pattern_label(pat), fontsize=11)
        ax.grid(True, which="both", alpha=0.25)
        if idx // ncol == nrow - 1:
            ax.set_xlabel(r"$p_{\mathrm{res}}$", fontsize=11)
        if idx % ncol == 0:
            ax.set_ylabel("probability", fontsize=10)

    for j in range(n, nrow * ncol):  # hide unused axes
        axes[j // ncol][j % ncol].axis("off")

    handles = [Line2D([0], [0], color=_C_TWO, lw=2, label="two-anc"),
               Line2D([0], [0], color=_C_ONE, lw=2, label="1-rep"),
               Line2D([0], [0], color=_C_THREE, lw=2, label="3-rep"),
               Line2D([0], [0], color=_C_FIVE, lw=2, label="5-rep"),
               Line2D([0], [0], color=_C_THREE, lw=1.4, ls=":", label="1-rep x3"),
               Line2D([0], [0], color=_C_FIVE, lw=1.4, ls=":", label="1-rep x5")]
    # Park the legend in the empty tail of the grid (unused bottom-row panels)
    # if there is room, else above the figure.
    if n < nrow * ncol:
        first_empty = axes[n // ncol][n % ncol]
        first_empty.legend(handles=handles, loc="center", fontsize=12,
                           framealpha=0.9, title="model")
    else:
        fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=11,
                   framealpha=0.9, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle(r"Noise-mechanism probabilities vs $p_{\mathrm{res}}$ "
                 r"(at $p_z = 0$)", fontsize=15)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)
    print(f"{n} mechanisms plotted -> {args.out}")


if __name__ == "__main__":
    main()
