#!/usr/bin/env python3
"""Fit the pure-p_m thresholds for the gauge-honeycomb phenomenological circuits.

Consumes out/collection/peterjan_check_pz.csv (run_check_peterjan.sh Scan A:
p=0, p_m scanned, one scheme per gauge factor g) and fits a p_m threshold per g,
reusing plot_grid's FSS fitter (honeycomb-periodic, single memory observable).
Prints a table with the g=1 baseline gain factor and draws one LER-vs-p_m panel
per g with the fitted threshold marked.

Unlike analyze_axis_threshold._fit (which caps unknown schemes at pz_max=0.09),
this threads a per-g p_m window cap through fit_pz_threshold, since the threshold
grows with g and would otherwise be pinned against the 0.09 default.

Usage:
    python tools/peterjan_tools/analyze_peterjan_check.py
    python tools/peterjan_tools/analyze_peterjan_check.py --method minvar
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

import plot_grid as g          # noqa: E402
import analyze_axis_threshold as az  # noqa: E402

_CHECK_PZ = str(_ROOT / "out/collection/peterjan_check_pz.csv")

# Per-g p_m FSS window cap (thresholds rise with g). Widen if a fit pins against
# the edge; keep in sync with the PM_HI ceilings in run_grid_peterjan.sh.
_PM_MAX = {1: 0.09, 2: 0.16, 3: 0.28, 4: 0.32, 5: 0.38}


def _with_dated(csv: str) -> list[str]:
    p = pathlib.Path(csv)
    return [csv, str(p.with_name(f"{p.stem},*.csv"))]


def _load(csv):
    try:
        df = g.load_grid(_with_dated(csv))
    except FileNotFoundError:
        return None
    return df[df["noise_model"].fillna("").str.startswith("peterjan_pheno_")]


def _fit_pm(sub, g_val: int, method: str):
    """(p_th, p_th_err) on the p_m axis for one g, per-g window cap. None if no fit."""
    pm_max = _PM_MAX.get(g_val, 0.28)
    if method == "minvar":
        pth = az._threshold_min_variance(sub, "p_bell_z")
        return (pth, None) if pth is not None else None
    if method == "raw":
        fit = g.fit_pz_threshold_raw(sub, float(sub["p"].median()))
    else:
        fit = g.fit_pz_threshold(sub, float(sub["p"].median()), pz_max=pm_max)
    return (fit.p_th, fit.p_th_err) if fit is not None else None


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=_CHECK_PZ)
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures" / "peterjan_check_pz.pdf")
    ap.add_argument("--method", choices=["raw", "fss", "minvar"], default="fss")
    args = ap.parse_args()
    g.THRESHOLD_METHOD = args.method

    df = _load(args.csv)
    if df is None or len(df) == 0:
        sys.exit(f"No data in {args.csv} — run run_check_peterjan.sh first.")

    schemes = sorted(df["noise_model"].dropna().unique())
    reps = sorted((int(nm.rsplit("_g", 1)[1]), nm) for nm in schemes)

    thr: dict[int, float | None] = {}
    err: dict[int, float | None] = {}
    for gv, nm in reps:
        res = _fit_pm(df[df["noise_model"] == nm], gv, args.method)
        thr[gv], err[gv] = res if res is not None else (None, None)

    # ---- figure: one LER-vs-p_m panel per g -------------------------------
    ncols = len(reps)
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 4.3),
                             constrained_layout=True, squeeze=False)
    for c, (gv, nm) in enumerate(reps):
        ax = axes[0][c]
        sub = df[df["noise_model"] == nm]
        for d in sorted(int(x) for x in sub["expected_distance"].dropna().unique()):
            xs, ys, es = az._curve(sub, d, "p_bell_z")
            if len(xs) == 0:
                continue
            ax.errorbar(xs, ys, yerr=es, marker=g.MARKERS.get(d, "o"),
                        linestyle="-", color=g.COLORS.get(d), markersize=4,
                        linewidth=1.4, capsize=2, elinewidth=0.8,
                        label=g._code_label(sub, d))
        if thr[gv] is not None:
            ax.axvline(thr[gv], color="k", linestyle="--", linewidth=1.4)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$p_m$ (measurement flip)", fontsize=13)
        if c == 0:
            ax.set_ylabel("Logical error rate per observable", fontsize=12)
        thr_txt = (rf"$p_m^* = {thr[gv]:.3g}$" if thr[gv] is not None
                   else "no threshold fit")
        ax.set_title(f"gauge factor g={gv}\n{thr_txt}", fontsize=12)
        ax.grid(True, which="both", alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8, loc="lower right", framealpha=0.85)
    fig.suptitle(r"Gauge honeycomb (phenomenological): pure $p_m$ threshold "
                 r"at $p=0$", fontsize=15)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)

    # ---- table ------------------------------------------------------------
    def _s(v, w):
        return (f"{v:.4g}".rjust(w) if v is not None else "n/a".rjust(w))

    base = thr.get(1)
    print("\n  Pure p_m threshold vs gauge factor (native repeated measurements)")
    print("  " + "-" * 56)
    print(f"  {'scheme':<20s}{'p_m*':>10s}{'err':>10s}{'vs g=1':>10s}")
    for gv, _nm in reps:
        v, e = thr[gv], err[gv]
        gain = (v / base) if (v and base) else None
        print(f"  {f'gauge g={gv}':<20s}{_s(v, 10)}{_s(e, 10)}"
              f"{(f'{gain:.2f}x' if gain else 'n/a'):>10s}")
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
