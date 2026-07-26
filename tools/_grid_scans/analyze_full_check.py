#!/usr/bin/env python3
"""Fit pure p_z thresholds for the full-decoder repeated circuits.

Consumes out/collection/full_check_pz.csv (run_pres_threshold_check_full.sh: p_res
pinned to a floor, p_z scanned, one scheme per repetition count r) and fits a
p_z threshold per r — reusing analyze_axis_threshold's estimators, honeycomb
only. Prints a comparison table against the local-decoder references from
pz_check.csv (r=1 baseline + 3/5-rep majority vote) when that CSV exists, and
draws one LER-vs-p_z panel per repetition with the fitted threshold marked.

Usage:
    python tools/_grid_scans/analyze_full_check.py
    python tools/_grid_scans/analyze_full_check.py --method minvar
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
import analyze_axis_threshold as az  # noqa: E402

_BASE = "biased_bell_v3_equivalent"
_FULL_CHECK = str(_ROOT / "out/collection/full_check_pz.csv")
_PZ_CHECK = str(_ROOT / "out/collection/pz_check.csv")

# Local-decoder reference schemes in pz_check.csv: (noise_model, label).
_LOCAL_REFS = [
    (_BASE, "r=1 baseline"),
    (f"{_BASE}_majority", "3-rep majority"),
    (f"{_BASE}_majority5", "5-rep majority"),
]


def _load_honeycomb(csv):
    try:
        df = g.load_grid(csv)
    except FileNotFoundError:
        return None
    return df[df["input_circuit"].fillna("").str.contains("honeycomb")]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=_FULL_CHECK,
                    help="Full-decoder p_z-scan CSV (default: full_check_pz.csv).")
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures" / "full_check_pz.pdf")
    ap.add_argument("--method", choices=["raw", "fss", "minvar"], default="fss",
                    help="Threshold estimator (as in analyze_axis_threshold.py).")
    args = ap.parse_args()
    g.THRESHOLD_METHOD = args.method

    df = _load_honeycomb(args.csv)
    if df is None or len(df) == 0:
        sys.exit(f"No data in {args.csv} — run run_pres_threshold_check_full.sh first.")

    # Repetition counts present, from the noise_model suffix. A plain base
    # scheme in the same CSV (no fullrep suffix) is the r=1 same-decoder
    # baseline (e.g. the tesseract runs include it for a fair comparison).
    schemes = sorted(
        nm for nm in df["noise_model"].dropna().unique() if "fullrep" in nm)
    reps = [(int(nm.rsplit("fullrep", 1)[1]), nm) for nm in schemes]
    if _BASE in df["noise_model"].unique():
        reps.append((1, _BASE))
    reps.sort()
    if not reps:
        sys.exit(f"No fullrep schemes found in {args.csv}.")

    # Fit each repetition's p_z threshold (same estimator as the check figures).
    thr: dict[int, float | None] = {}
    err: dict[int, float | None] = {}
    for r, nm in reps:
        res = az._fit(df[df["noise_model"] == nm], "p_bell_z", nm, args.method)
        thr[r], err[r] = res if res is not None else (None, None)

    # Local-decoder references from pz_check.csv (same estimator).
    chk = _load_honeycomb(_PZ_CHECK)
    refs: dict[str, float | None] = {}
    if chk is not None:
        for nm, _label in _LOCAL_REFS:
            res = az._fit(chk[chk["noise_model"] == nm], "p_bell_z", nm, args.method)
            refs[nm] = res[0] if res is not None else None

    # ---- figure: one LER panel per r --------------------------------------
    ncols = len(reps)
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 4.3),
                             constrained_layout=True, squeeze=False)
    for c, (r, nm) in enumerate(reps):
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
        if thr[r] is not None:
            ax.axvline(thr[r], color="k", linestyle="--", linewidth=1.4,
                       label=rf"$p_z^* = {thr[r]:.3g}$")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$p_z$", fontsize=15)
        if c == 0:
            ax.set_ylabel("Logical error rate per observable", fontsize=12)
        thr_txt = (rf"$p_z^* = {thr[r]:.3g}$" if thr[r] is not None
                   else "no threshold fit")
        ax.set_title(f"Full decoder, r={r}\n{thr_txt}", fontsize=12)
        ax.grid(True, which="both", alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8, loc="lower right", framealpha=0.85)
    fig.suptitle(r"Full-decoder pure $p_z$ threshold at $p_{\mathrm{res}} \approx 0$",
                 fontsize=15)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)

    # ---- comparison table ---------------------------------------------------
    def _s(v, w):
        return (f"{v:.4g}".rjust(w) if v is not None else "n/a".rjust(w))

    # Gain baseline: the in-CSV r=1 fit when present (same decoder), else the
    # pymatching r=1 reference from pz_check.csv.
    base = thr.get(1) or refs.get(_BASE)
    print("\n  Pure p_z threshold: full decoder vs local (majority) references")
    print("  " + "-" * 64)
    print(f"  {'scheme':<24s}{'p_z*':>10s}{'err':>10s}{'vs r=1':>10s}")
    for nm, label in _LOCAL_REFS:
        v = refs.get(nm)
        gain = (v / base) if (v and base) else None
        print(f"  {label:<24s}{_s(v, 10)}{'':>10s}"
              f"{(f'{gain:.2f}x' if gain else 'n/a'):>10s}")
    for r, _nm in reps:
        v, e = thr[r], err[r]
        gain = (v / base) if (v and base) else None
        print(f"  {f'full decoder, r={r}':<24s}{_s(v, 10)}{_s(e, 10)}"
              f"{(f'{gain:.2f}x' if gain else 'n/a'):>10s}")
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
