#!/usr/bin/env python3
"""Full-decoder vs local (majority-vote) repeated-measurement boundaries.

One boundary panel in the style of plot_grid_summary.py panel (a): threshold
boundaries p_z* vs p_res* for the honeycomb code, comparing the two ways of
decoding repeated pairwise measurements:

  - LOCAL (majority vote): the equivalent_v3_majority / majority5 schemes,
    where the R-repetition majority vote is folded into the noise block and
    the outer decoder sees a single effective measurement. r = 1, 3, 5.
    Read from out/collection/grid_bell.csv (the existing grid runs).

  - FULL decoder: repeat_measurements applied structurally (r copies of every
    mid-circuit measurement layer + comparison detectors), then the 1-rep
    equivalent_v3 translation; pymatching decodes everything together.
    r = 2, 3, 4, 5. Read from out/collection/grid_full.csv (run_grid_full.sh).

The r=1 baseline (biased_bell_v3_equivalent) is shared by both families and
drawn once in a neutral color. Axis-intercept stars are overlaid where the
dedicated check CSVs exist (pz_check/pres_check for the local schemes, and
full_check_pz/full_check_pres if/when full-decoder check scans are run).

Usage:
    python tools/_grid_scans/plot_full_vs_local.py
    python tools/_grid_scans/plot_full_vs_local.py --method minvar --logy
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"xtick.labelsize": 15, "ytick.labelsize": 15})

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
for _p in (str(_ROOT / "tools" / "_grid_scans"), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import plot_grid as g  # noqa: E402
import analyze_axis_threshold as az  # noqa: E402  (axis-intercept estimators)
from plot_grid_summary import _monotonic_decreasing  # noqa: E402
plt.rcParams.update({"xtick.labelsize": 15, "ytick.labelsize": 15})

# Default (plain pymatching) CSVs. --correlated swaps in the `*_correlated`
# siblings written by run_grid_{bell,full}_correlated.sh (native correlated MWPM).
# Correlated rows must live in separate CSVs: sinter's strong_id encodes the
# decoder and load_grid groups by strong_id, so mixing decoders in one CSV would
# put two rows per grid point and corrupt the per-scheme threshold fits.
_GRID_LOCAL = str(_ROOT / "out/collection/grid_bell.csv")
_GRID_FULL = str(_ROOT / "out/collection/grid_full.csv")
# Axis-intercept check CSVs (optional; stars are skipped when missing).
_CHECK_PZ_LOCAL = str(_ROOT / "out/collection/pz_check.csv")
_CHECK_PRES_LOCAL = str(_ROOT / "out/collection/pres_check.csv")
_CHECK_PZ_FULL = str(_ROOT / "out/collection/full_check_pz.csv")
_CHECK_PRES_FULL = str(_ROOT / "out/collection/full_check_pres.csv")


def _correlated_csv(path: str) -> str:
    """`out/collection/foo.csv` -> `out/collection/foo_correlated.csv`."""
    p = pathlib.Path(path)
    return str(p.with_name(f"{p.stem}_correlated{p.suffix}"))

_BASE = "biased_bell_v3_equivalent"


def _fullrep(r: int) -> str:
    return f"{_BASE}_fullrep{r}"


# FSS p_z window cap per scheme. plot_grid's default (_pz_max_for) only knows
# the majority names; the fullrep thresholds sit above the majority ones
# (~7% at 3-rep, ~11% at 5-rep), so give each r matching headroom.
_PZ_MAX = {
    _BASE: 0.09,
    f"{_BASE}_majority": 0.16,
    f"{_BASE}_majority5": 0.25,
    _fullrep(2): 0.12,
    _fullrep(3): 0.18,
    _fullrep(4): 0.22,
    _fullrep(5): 0.28,
}

# Curves: (scheme, label, color, marker, linestyle). Local = Blues (matching
# the summary figure's honeycomb hue), full = Greens; shade grows with r.
_blues = plt.get_cmap("Blues")
_greens = plt.get_cmap("Greens")
_CURVES = [
    (_BASE, "no repetition (r=1)", "#444444", "o", "-"),
    (f"{_BASE}_majority", "Pymatching decoder with majority vote, r=3", _blues(0.60), "D", "--"),
    (f"{_BASE}_majority5", "Pymatching decoder with majority vote, r=5", _blues(0.90), "P", "--"),
    (_fullrep(2), "pymatching decoder without majority vote, r=2", _greens(0.45), "s", "-"),
    (_fullrep(3), "pymatching decoder without majority vote, r=3", _greens(0.62), "v", "-"),
    (_fullrep(4), "pymatching decoder without majority vote, r=4", _greens(0.78), "^", "-"),
    (_fullrep(5), "pymatching decoder without majority vote, r=5", _greens(0.94), "X", "-"),
]


def _collect_thresholds(sdf, pz_max: float):
    """collect_thresholds with an explicit FSS p_z window cap.

    Mirrors plot_grid.collect_thresholds but threads pz_max through to
    fit_pz_threshold (needed because the fullrep thresholds exceed the 0.09
    default cap). minvar/raw methods delegate to estimate_pz_threshold, which
    doesn't consume pz_max.
    """
    out = []
    for p in sorted(sdf["p"].dropna().unique()):
        if g.THRESHOLD_METHOD == "fss":
            fit = g.fit_pz_threshold(sdf, float(p), pz_max=pz_max)
            res = (fit.p_th, fit.p_th_err) if fit is not None else None
        else:
            res = g.estimate_pz_threshold(sdf, float(p))
        if res is not None and np.isfinite(res[0]):
            out.append((float(p), res[0], res[1]))
    return out


def _intercept(chk, scheme: str, xcol: str, method: str):
    """Axis-intercept threshold from a check CSV, or None."""
    if chk is None:
        return None
    sub = chk[chk["noise_model"] == scheme]
    if len(sub) == 0:
        return None
    res = az._fit(sub, xcol, scheme, method)
    return res[0] if res is not None else None


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--out", type=pathlib.Path, default=None,
                    help="Output PDF (default: full_vs_local[_correlated].pdf).")
    ap.add_argument("--method", choices=["raw", "fss", "minvar"], default="fss",
                    help="Threshold estimator (as in plot_grid_summary.py).")
    ap.add_argument("--logy", action="store_true",
                    help="Log y-axis for the boundary panel.")
    ap.add_argument("--correlated", action="store_true",
                    help="Read the *_correlated.csv data (native correlated MWPM, "
                         "--decoder pymatching-correlated) instead of plain pymatching.")
    args = ap.parse_args()
    g.THRESHOLD_METHOD = args.method

    # Resolve the CSV set and output name for the chosen decoder.
    grid_local, grid_full = _GRID_LOCAL, _GRID_FULL
    chk = [_CHECK_PZ_LOCAL, _CHECK_PRES_LOCAL, _CHECK_PZ_FULL, _CHECK_PRES_FULL]
    default_out = _ROOT / "out" / "figures" / "full_vs_local.pdf"
    if args.correlated:
        grid_local, grid_full = _correlated_csv(grid_local), _correlated_csv(grid_full)
        chk = [_correlated_csv(c) for c in chk]
        default_out = _ROOT / "out" / "figures" / "full_vs_local_correlated.pdf"
    chk_pz_local, chk_pres_local, chk_pz_full, chk_pres_full = chk
    if args.out is None:
        args.out = default_out

    def load(csv):
        try:
            df = g.load_grid(csv)
        except FileNotFoundError:
            return None
        # Honeycomb-only figure; drop other families if the CSV mixes them.
        return df[df["input_circuit"].fillna("").str.contains("honeycomb")]

    df_local = load(grid_local)
    df_full = load(grid_full)
    chk_pz = {  # scheme csv per decoder family
        "local": load(chk_pz_local), "full": load(chk_pz_full)}
    chk_pres = {
        "local": load(chk_pres_local), "full": load(chk_pres_full)}

    fig, ax = plt.subplots(figsize=(6.6, 5.2), constrained_layout=True)

    # Under correlated matching every curve uses correlated MWPM; note it once
    # in the legend title rather than lengthening each label.
    curves = _CURVES
    if args.correlated:
        curves = [(s, l.replace("Pymatching decoder", "matching")
                          .replace("pymatching decoder", "matching"),
                   c, m, ls) for s, l, c, m, ls in _CURVES]

    for scheme, label, color, marker, ls in curves:
        is_full = "fullrep" in scheme
        df = df_full if is_full else df_local
        kind = "full" if is_full else "local"
        if df is None:
            print(f"[skip] {label}: CSV not found")
            continue
        sdf = df[df["noise_model"] == scheme]
        if len(sdf) == 0:
            print(f"[skip] {label}: no rows for {scheme}")
            continue

        thr = _collect_thresholds(sdf, _PZ_MAX[scheme])
        if not thr:
            print(f"[skip] {label}: no thresholds fit")
            continue
        p, pth, err = map(np.array, zip(*thr))
        n_all = len(p)

        # Axis intercepts (stars), where check data exists. Same estimator.
        pt_pz = _intercept(chk_pz[kind], scheme, "p_bell_z", args.method)
        pt_pres = _intercept(chk_pres[kind], scheme, "p", args.method)

        # Prepend the p_res=0 intercept before enforcing monotonicity, then
        # append the p_z=0 intercept (same treatment as plot_grid_summary).
        pl, tl, el = list(p), list(pth), list(err)
        if pt_pz is not None:
            pl, tl, el = [0.0] + pl, [pt_pz] + tl, [0.0] + el
        p, pth, err = _monotonic_decreasing(np.array(pl), np.array(tl), np.array(el))
        n_mono = len(p)
        if pt_pres is not None:
            p = np.append(p, pt_pres)
            pth = np.append(pth, 0.0)
            err = np.append(err, 0.0)

        ax.errorbar(p, pth, yerr=None, fmt=marker + ls, color=color,
                    label=label, markersize=5, linewidth=1.7, capsize=3)
        for xi, yi in ((0.0, pt_pz), (pt_pres, 0.0)):
            if xi is not None and yi is not None:
                ax.scatter([xi], [yi], marker="*", s=130, color=color,
                           edgecolors="k", linewidths=0.5, zorder=6)
        print(f"{label}: {n_mono}/{n_all} grid thresholds (monotonic); "
              f"intercepts p_z*(p_res=0)={pt_pz}, p_res*(p_z=0)={pt_pres}")

    if args.logy:
        ax.set_yscale("log")
    ax.set_xlabel(r"$p^*_{\text{res}}$", fontsize=17)
    ax.set_ylabel(r"$p_z^*$ (threshold)", fontsize=17)
    ax.grid(True, alpha=0.3)
    leg_title = "correlated MWPM" if args.correlated else None
    ax.legend(fontsize=9, framealpha=0.9, title=leg_title, title_fontsize=9)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
