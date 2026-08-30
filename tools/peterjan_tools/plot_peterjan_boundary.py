#!/usr/bin/env python3
"""Boundary figure p_m* vs p* for the gauge-honeycomb phenomenological circuits.

The analog of plot_full_vs_local.py, for the native repeated-measurement route
(gauge factor g repeats each check g times with comparison detectors built in,
so the circuits stay graphlike and plain pymatching thresholds are meaningful).

One boundary curve per gauge factor g: threshold p_m* vs data-noise p*, fit from
out/collection/peterjan_grid.csv (run_grid_peterjan.sh). Axis-intercept stars
come from the 1-D check scans (peterjan_check_pz.csv for the p=0 intercept,
peterjan_check_pres.csv for the p_m=0 intercept).

Usage:
    python tools/peterjan_tools/plot_peterjan_boundary.py
    python tools/peterjan_tools/plot_peterjan_boundary.py --method minvar --logy
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
import analyze_axis_threshold as az  # noqa: E402
from plot_grid_summary import _monotonic_decreasing  # noqa: E402

_GRID = str(_ROOT / "out/collection/peterjan_grid.csv")
_CHECK_PZ = str(_ROOT / "out/collection/peterjan_check_pz.csv")
_CHECK_PRES = str(_ROOT / "out/collection/peterjan_check_pres.csv")

# Per-g FSS p_m window cap (threshold rises with g); keep in sync with
# run_grid_peterjan.sh PM_HI and analyze_peterjan_check._PM_MAX.
_PM_MAX = {1: 0.09, 2: 0.16, 3: 0.28, 4: 0.32, 5: 0.38}

_greens = plt.get_cmap("Greens")
_CURVES = [
    (1, "g=1 (no repetition)", "#444444", "o", "-"),
    (2, "g=2", _greens(0.40), "s", "-"),
    (3, "g=3", _greens(0.55), "v", "-"),
    (4, "g=4", _greens(0.72), "^", "-"),
    (5, "g=5", _greens(0.90), "X", "-"),
]


def _scheme(g_val: int) -> str:
    return f"peterjan_pheno_g{g_val}"


def _collect_thresholds(sdf, pm_max: float):
    """Grid thresholds p_m*(p) for one scheme, with an explicit p_m window cap."""
    out = []
    for p in sorted(sdf["p"].dropna().unique()):
        if g.THRESHOLD_METHOD == "fss":
            fit = g.fit_pz_threshold(sdf, float(p), pz_max=pm_max)
            res = (fit.p_th, fit.p_th_err) if fit is not None else None
        else:
            res = g.estimate_pz_threshold(sdf, float(p))
        if res is not None and np.isfinite(res[0]):
            out.append((float(p), res[0], res[1]))
    return out


def _load(csv):
    p = pathlib.Path(csv)
    try:
        df = g.load_grid([csv, str(p.with_name(f"{p.stem},*.csv"))])
    except FileNotFoundError:
        return None
    return df[df["noise_model"].fillna("").str.startswith("peterjan_pheno_")]


def _intercept(chk, scheme: str, xcol: str, method: str, pm_max: float):
    """Axis-intercept threshold from a check CSV, or None.

    For the p_m axis (xcol='p_bell_z') we bypass az._fit's scheme-name-based cap
    (0.09 for our peterjan_ schemes) and thread the per-g pm_max through, exactly
    as _collect_thresholds does for the grid; the p axis has no such cap issue.
    """
    if chk is None:
        return None
    sub = chk[chk["noise_model"] == scheme]
    if len(sub) == 0:
        return None
    if xcol == "p_bell_z" and method == "fss":
        fit = g.fit_pz_threshold(sub, float(sub["p"].median()), pz_max=pm_max)
        return fit.p_th if fit is not None else None
    res = az._fit(sub, xcol, scheme, method)
    return res[0] if res is not None else None


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures" / "peterjan_boundary.pdf")
    ap.add_argument("--method", choices=["raw", "fss", "minvar"], default="fss")
    ap.add_argument("--logy", action="store_true")
    args = ap.parse_args()
    g.THRESHOLD_METHOD = args.method

    df = _load(_GRID)
    chk_pz = _load(_CHECK_PZ)
    chk_pres = _load(_CHECK_PRES)
    if df is None:
        sys.exit(f"No grid data at {_GRID} — run run_grid_peterjan.sh first.")

    fig, ax = plt.subplots(figsize=(6.6, 5.2), constrained_layout=True)

    for g_val, label, color, marker, ls in _CURVES:
        scheme = _scheme(g_val)
        sdf = df[df["noise_model"] == scheme]
        if len(sdf) == 0:
            print(f"[skip] {label}: no rows for {scheme}")
            continue
        thr = _collect_thresholds(sdf, _PM_MAX.get(g_val, 0.28))
        if not thr:
            print(f"[skip] {label}: no thresholds fit")
            continue
        p, pth, err = map(np.array, zip(*thr))
        n_all = len(p)

        pm_max = _PM_MAX.get(g_val, 0.28)
        pt_pm = _intercept(chk_pz, scheme, "p_bell_z", args.method, pm_max)  # p≈0 intercept
        pt_p = _intercept(chk_pres, scheme, "p", args.method, pm_max)        # p_m≈0 intercept

        pl, tl, el = list(p), list(pth), list(err)
        if pt_pm is not None:
            pl, tl, el = [0.0] + pl, [pt_pm] + tl, [0.0] + el
        p, pth, err = _monotonic_decreasing(np.array(pl), np.array(tl), np.array(el))
        n_mono = len(p)
        if pt_p is not None:
            p = np.append(p, pt_p)
            pth = np.append(pth, 0.0)
            err = np.append(err, 0.0)

        ax.errorbar(p, pth, yerr=None, fmt=marker + ls, color=color,
                    label=label, markersize=5, linewidth=1.7, capsize=3)
        for xi, yi in ((0.0, pt_pm), (pt_p, 0.0)):
            if xi is not None and yi is not None:
                ax.scatter([xi], [yi], marker="*", s=130, color=color,
                           edgecolors="k", linewidths=0.5, zorder=6)
        print(f"{label}: {n_mono}/{n_all} grid thresholds (monotonic); "
              f"intercepts p_m*(p=0)={pt_pm}, p*(p_m=0)={pt_p}")

    if args.logy:
        ax.set_yscale("log")
    ax.set_xlabel(r"$p^*$ (data noise)", fontsize=17)
    ax.set_ylabel(r"$p_m^*$ (threshold)", fontsize=17)
    ax.set_title("Gauge honeycomb (phenomenological): threshold vs gauge factor",
                 fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, framealpha=0.9)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
