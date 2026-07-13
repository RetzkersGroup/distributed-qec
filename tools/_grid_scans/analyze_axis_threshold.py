#!/usr/bin/env python3
"""Fit a single-axis threshold for the four code/scheme combos.

Consumes one of the dedicated sanity-check CSVs written by
run_pres_threshold_check.sh, each holding all four combos on a single scan line
with the *other* axis pinned to 0:

  pres_check.csv  -- p_z   = 0, p_res swept  ->  pure p_res threshold
  pz_check.csv    -- p_res = 0, p_z   swept  ->  pure p_z   threshold

The combos are honeycomb / hyperbolic x base / majority, distinguished by
input_circuit (honeycomb vs hyperbolic) and noise_model (base vs the
*_equivalent_majority scheme). The scan axis is auto-detected (the axis that
varies), and a threshold is fit along it -- reusing src/threshold_fit.py (the
same locate -> quadratic-FSS estimator the (p, p_z) grid uses). Prints each
fitted threshold +/- error and draws a 2x2 LER panel with it marked.

Usage:
    python tools/_grid_scans/analyze_axis_threshold.py --csv out/collection/pres_check.csv
    python tools/_grid_scans/analyze_axis_threshold.py --csv out/collection/pz_check.csv
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"xtick.labelsize": 13, "ytick.labelsize": 13})

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
for _p in (str(_ROOT / "tools" / "_grid_scans"), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from math import comb  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

import plot_grid as g  # noqa: E402  (load_grid, fit_pz_threshold, _pz_curve, _cutoff_for, ...)
from threshold_fit import locate_threshold, fit_threshold  # noqa: E402


def _invert(f, target):
    """Solve f(p) = target for p in (0, 0.5). f is monotonic there. None-safe."""
    if target is None:
        return None
    try:
        return brentq(lambda p: f(p) - target, 1e-9, 0.5)
    except Exception:
        return None


def _maj_transform(reps):
    """p_z transform: probability an R-repetition majority vote fails (>=R//2+1)."""
    thr = reps // 2 + 1
    return lambda p: sum(comb(reps, k) * p ** k * (1 - p) ** (reps - k)
                         for k in range(thr, reps + 1))


def _odd_transform(reps):
    """p_res transform: probability of an odd number of flips over R rounds."""
    return lambda p: sum(comb(reps, k) * p ** k * (1 - p) ** (reps - k)
                         for k in range(1, reps + 1, 2))

# Grid is families (rows) x schemes (cols). noise_model strings identify schemes.
_BASE = "biased_bell_v3"                       # two-ancilla Bell (physical)
_SINGLE = "biased_bell_v3_equivalent"          # single-rep equivalent (calculator)
_MAJ = "biased_bell_v3_equivalent_majority"    # 3-rep majority vote (calculator)
_MAJ5 = "biased_bell_v3_equivalent_majority5"  # 5-rep majority vote (calculator)
_FAMILIES = [("honeycomb", "Honeycomb"), ("hyperbolic", "Hyperbolic Floquet")]
# (noise_model, panel label, repetitions). The non-repeated baseline is the
# single-ancilla 1-rep equivalent (reproduces the two-ancilla biased_bell_v3).
_SCHEMES = [
    (_SINGLE, "1-rep equivalent", 1),
    (_MAJ, "3-rep majority vote", 3),
    (_MAJ5, "5-rep majority vote", 5),
]

# Per-axis presentation. p_z majority thresholds sit near ~7%, so the p_z fit
# needs a higher saturation cap than the 0.09 grid default to bracket them.
_AXES = {
    # `expected` is a reference curve for the 3-repetition gadget, overlaid for
    # comparison: majority-vote failure prob (>=2 of 3 wrong) on the p_z axis,
    # and odd-number-of-flips prob (1 or 3 of 3) on the p_res axis.
    # `family` selects the R-repetition transform used to predict the repeated
    # threshold from the 1-rep one: "odd" (odd # of flips) on p_res, "maj"
    # (majority failure) on p_z.
    "p": dict(col="p", fixed="p_bell_z", xlabel=r"$p_{\mathrm{res}}$",
              sym=r"p_{\mathrm{res}}^*", title=r"Pure $p_{\mathrm{res}}$ threshold at $p_z = 0$",
              name="p_res threshold at p_z = 0", family="odd"),
    "p_bell_z": dict(col="p_bell_z", fixed="p", xlabel=r"$p_z$",
                     sym=r"p_z^*", title=r"Pure $p_z$ threshold at $p_{\mathrm{res}} = 0$",
                     name="p_z threshold at p_res = 0", family="maj"),
}


def _curve(sub, d, xcol):
    """Sorted (x, pL, sigma) for one distance on the scan line, valid points only."""
    s = sub[(sub["expected_distance"] == d) & (sub["errors"] > 0)
            & (sub["shots"] > 0)].sort_values(xcol)
    return s[xcol].to_numpy(float), s["pL"].to_numpy(float), s["sigma"].to_numpy(float)


def _has_protective_regime(sub, xcol) -> bool:
    """True if the largest distance beats the smallest at the lowest scan value.

    If the biggest code is already worse than the smallest at the low-noise end,
    we are above threshold everywhere and any crossing is spurious. Mirrors
    plot_grid._has_protective_regime, generalized to either scan axis.
    """
    ds = sorted(sub["expected_distance"].dropna().unique())
    if len(ds) < 2:
        return False
    x_lo, y_lo, _ = _curve(sub, int(ds[0]), xcol)
    x_hi, y_hi, _ = _curve(sub, int(ds[-1]), xcol)
    if len(x_lo) == 0 or len(x_hi) == 0:
        return False
    j = int(np.argmin(np.abs(x_lo - x_hi[0])))
    return y_hi[0] < y_lo[j]


def _fit_p_axis(sub, n_per_side: int = 2, logical_cutoff: float | None = None):
    """Fit the p_res threshold across distances (p_z fixed at 0). None if no fit."""
    if logical_cutoff is None:
        logical_cutoff = g._cutoff_for(sub)
    s = sub[(sub["errors"] > 0) & (sub["shots"] > 0) & (sub["pL"] < logical_cutoff)]
    if s["expected_distance"].nunique() < 2 or not _has_protective_regime(sub, "p"):
        return None
    x = s["p"].to_numpy(float)
    d = s["expected_distance"].to_numpy(float)
    pL = s["pL"].to_numpy(float)
    sigma = s["sigma"].to_numpy(float)
    p_init = locate_threshold(x, d, pL)
    if p_init is None:
        return None
    try:
        fit = fit_threshold(x, d, pL, sigma=sigma, p_init=p_init, n_per_side=n_per_side)
    except Exception:
        return None
    if not (x.min() <= fit.p_th <= x.max()):
        return None
    return fit


def _threshold_min_variance(sub, xcol, rate_col="pL"):
    """Min-variance crossing with a plateau-relative saturation guard.

    At each scan value, take the mean logical error and the variance of the
    log logical error across distances. The maximal error rate is the largest of
    those means (the saturated plateau). The threshold is the scan value of
    *minimum* log-variance among the points whose mean logical error is below
    half that maximum -- i.e. the tightest crossing of the distance curves,
    restricted to the sub-saturation regime so the flat plateau (where all
    distances coincide, also low-variance) can't win. Returns the point (no fit,
    no error bar), or None.
    """
    s = sub[(sub["errors"] > 0) & (sub["shots"] > 0)]
    xs, means, varlog = [], [], []
    for x, grp in s.groupby(xcol):
        pL = grp[rate_col].to_numpy(float)
        pL = pL[pL > 0]
        if len(pL) < 2:
            continue
        xs.append(float(x))
        means.append(float(pL.mean()))
        varlog.append(float(np.var(np.log(pL))))
    if not xs:
        return None
    xs, means, varlog = np.array(xs), np.array(means), np.array(varlog)
    valid = means < 0.5 * means.max()
    if not valid.any():
        return None
    return float(xs[np.argmin(np.where(valid, varlog, np.inf))])


def _fit(sub, xcol, nm, method="minvar"):
    """Return (p_th, p_th_err) for the chosen method, or None.

    method="minvar": plateau-relative min-variance crossing point (no error bar).
    method="fss":    quadratic finite-size-scaling fit (p_th_err from the fit).
    """
    if method == "minvar":
        pth = _threshold_min_variance(sub, xcol)
        return (pth, None) if pth is not None else None
    # fss:
    if xcol == "p":
        fit = _fit_p_axis(sub)
    else:
        # Majority thresholds sit above the 0.09 grid cap (3-rep ~7%, 5-rep ~11%);
        # raise pz_max accordingly.
        p_val = float(sub["p"].median())
        pz_max = 0.25 if "majority5" in nm else (0.16 if "equivalent" in nm else 0.09)
        fit = g.fit_pz_threshold(sub, p_val, pz_max=pz_max)
    return (fit.p_th, fit.p_th_err) if fit is not None else None


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=str(_ROOT / "out/collection/pres_check.csv"),
                    help="Single-axis sanity-check CSV (pres_check.csv or pz_check.csv).")
    ap.add_argument("--out", type=pathlib.Path, default=None,
                    help="Output figure (default: out/figures/<csv-stem>.pdf).")
    ap.add_argument("--method", choices=["minvar", "fss"], default="fss",
                    help="Threshold estimator: minvar = plateau-relative "
                         "min-variance crossing point (default); fss = quadratic "
                         "finite-size-scaling fit.")
    args = ap.parse_args()

    df = g.load_grid(args.csv)
    df["family"] = df["input_circuit"].fillna("").apply(
        lambda s: "honeycomb" if "honeycomb" in s else
                  ("hyperbolic" if "hyperbolic" in s else "?"))

    # Auto-detect the scan axis: the one that varies (the other is pinned at 0).
    xcol = "p" if df["p"].nunique() >= df["p_bell_z"].nunique() else "p_bell_z"
    ax_cfg = _AXES[xcol]
    # Keep the scan line: the fixed axis is pinned at ~0 (exactly 0, or the tiny
    # p_res floor Scan B uses to keep the matching graph well-formed).
    df = df[np.isclose(df[ax_cfg["fixed"]].astype(float), 0.0, atol=1e-4)]

    out = args.out or (_ROOT / "out" / "figures" / f"{pathlib.Path(args.csv).stem}.pdf")

    fig, axes = plt.subplots(len(_FAMILIES), len(_SCHEMES),
                             figsize=(5 * len(_SCHEMES), 4.3 * len(_FAMILIES)),
                             constrained_layout=True, squeeze=False)
    thr_by: dict[tuple[str, str], float | None] = {}  # (family, scheme) -> p_th
    transform_family = {"maj": _maj_transform, "odd": _odd_transform}[ax_cfg["family"]]
    for r, (fam, fam_label) in enumerate(_FAMILIES):
        for c, (nm, s_label, reps) in enumerate(_SCHEMES):
            ax = axes[r][c]
            sub = df[(df["family"] == fam) & (df["noise_model"] == nm)]
            for d in sorted(int(x) for x in sub["expected_distance"].dropna().unique()):
                xs, ys, es = _curve(sub, d, xcol)
                if len(xs) == 0:
                    continue
                ax.errorbar(xs, ys, yerr=es, marker=g.MARKERS.get(d, "o"),
                            linestyle="-", color=g.COLORS.get(d), markersize=4,
                            linewidth=1.4, capsize=2, elinewidth=0.8,
                            label=g._code_label(sub, d))
            res = _fit(sub, xcol, nm, args.method)
            p_th = res[0] if res is not None else None
            p_th_err = res[1] if res is not None else None
            thr_by[(fam, nm)] = p_th
            if p_th is not None:
                ax.axvline(p_th, color="k", linestyle="--", linewidth=1.4,
                           label=rf"measured ${ax_cfg['sym']} = {p_th:.3g}$")
                err_txt = rf" \pm {p_th_err:.1g}" if p_th_err is not None else ""
                thr_txt = rf"${ax_cfg['sym']} = {p_th:.3g}{err_txt}$"
            else:
                thr_txt = "no threshold fit"

            # On each repeated panel, overlay the threshold predicted from the
            # paired 1-rep equivalent through the R-repetition transform:
            # p_th^Rrep = f_R^{-1}(p_th^1rep).
            if reps > 1:
                pred = _invert(transform_family(reps), thr_by.get((fam, _SINGLE)))
                if pred is not None:
                    ax.axvline(pred, color="#d62728", linestyle=":", linewidth=2.0,
                               label=rf"predicted from 1-rep: ${ax_cfg['sym']} = {pred:.3g}$")

            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlabel(ax_cfg["xlabel"], fontsize=15)
            if c == 0:
                ax.set_ylabel("Logical error rate per observable", fontsize=12)
            ax.set_title(f"{fam_label}: {s_label}\n{thr_txt}", fontsize=12)
            ax.grid(True, which="both", alpha=0.3)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=8, loc="lower right", framealpha=0.85)

    fig.suptitle(ax_cfg["title"], fontsize=16)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)

    # Comparison table: each repeated scheme's measured threshold vs the value
    # predicted from the paired 1-rep equivalent through the R-rep transform.
    def _s(v, w):
        return (f"{v:.4g}".rjust(w) if v is not None else "n/a".rjust(w))

    rep_schemes = [(nm, lbl, reps) for nm, lbl, reps in _SCHEMES if reps > 1]
    print(f"\n  {ax_cfg['name']}  (repeated thresholds predicted from 1-rep equiv)")
    print("  " + "-" * 78)
    print(f"  {'code':<14s}{'rep':>4s}{'1-rep':>11s}{'meas':>12s}"
          f"{'pred':>12s}{'meas/pred':>12s}")
    for fam, label in _FAMILIES:
        single = thr_by.get((fam, _SINGLE))
        for nm, _lbl, reps in rep_schemes:
            meas = thr_by.get((fam, nm))
            pred = _invert(transform_family(reps), single)
            ratio = (meas / pred) if (meas and pred) else None
            print(f"  {label:<14s}{reps:>4d}{_s(single, 11)}{_s(meas, 12)}"
                  f"{_s(pred, 12)}{(f'{ratio:.2f}' if ratio else 'n/a'):>12s}")
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main()
