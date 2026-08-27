#!/usr/bin/env python3
"""LER vs p_z: tesseract vs pymatching, r=1 and r=3 (plain matplotlib).

Same data and layout as plot_tesseract_vs_pymatching_ler.py, but drawn with
matplotlib directly instead of sinter.plot_error_rate. Error bars use the
Wilson score interval on the per-shot rate (so 0-error points still show a
finite upper bound, like sinter's), then mapped to the per-observable rate
pL = 1 - (1 - p_shot) ** (1/k).

Sources:
  - Tesseract: full_check_pz_tesseract{.csv, ,*.csv}
    (r=1 = biased_bell_v3_equivalent, r=3 = biased_bell_v3_equivalent_fullrep3).
  - Pymatching r=1: pz_check.csv (biased_bell_v3_equivalent).
  - Pymatching r=3: full_check_pz.csv (biased_bell_v3_equivalent_fullrep3).

Colors: r=1 tesseract reds, r=1 pymatching blues, r=3 tesseract purples,
r=3 pymatching grays; shade encodes distance.

Usage:
    python tools/_grid_scans/plot_tesseract_vs_pymatching_ler_mpl.py
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

_TESS_CSV = str(_ROOT / "out/collection/full_check_pz_tesseract.csv")
_PM_R1_CSV = str(_ROOT / "out/collection/pz_check.csv")
_MAJ_CSV = str(_ROOT / "out/collection/pz_check.csv")

_BASE = "biased_bell_v3_equivalent"
_FULLREP3 = "biased_bell_v3_equivalent_fullrep3"
_MAJORITY = "biased_bell_v3_equivalent_majority"  # 3-rep majority vote

# (label, source_csv, noise_model, decoder, colormap)
_CURVES = [
    ("r=1 tesseract", _TESS_CSV, _BASE, "tesseract_short", plt.get_cmap("Reds")),
    ("r=1 pymatching", _PM_R1_CSV, _BASE, "pymatching", plt.get_cmap("Blues")),
    ("r=3 tesseract", _TESS_CSV, _FULLREP3, "tesseract_short", plt.get_cmap("Purples")),
    ("3-rep majority+pymatching", _MAJ_CSV, _MAJORITY, "pymatching", plt.get_cmap("Greys")),
]

_Z = 1.96  # 95% Wilson interval


def _with_dated(csv: str) -> list[str]:
    p = pathlib.Path(csv)
    return [csv, str(p.with_name(f"{p.stem},*.csv"))]


def _load(csv: str):
    try:
        df = g.load_grid(_with_dated(csv))
    except FileNotFoundError:
        return None
    return df[df["input_circuit"].fillna("").str.contains("honeycomb")]


def _wilson(errors: np.ndarray, shots: np.ndarray):
    """Wilson score interval (lo, hi) for a per-shot binomial rate.

    Finite even at 0 errors (lo=0, hi>0), unlike the normal-approx sigma.
    """
    n = shots.astype(float)
    phat = errors / n
    z2 = _Z * _Z
    denom = 1 + z2 / n
    center = (phat + z2 / (2 * n)) / denom
    half = (_Z * np.sqrt(phat * (1 - phat) / n + z2 / (4 * n * n))) / denom
    return np.clip(center - half, 0, None), center + half


def _per_obs(rate: np.ndarray, k: np.ndarray):
    """Per-shot rate -> per-observable rate: 1 - (1 - r) ** (1/k)."""
    return 1 - (1 - np.clip(rate, 0, 1 - 1e-15)) ** (1.0 / k)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures"
                    / "tesseract_vs_pymatching_ler_mpl.pdf")
    ap.add_argument("--floor", type=float, default=None,
                    help="Optional y-axis lower limit (e.g. 1e-7).")
    args = ap.parse_args()

    loaded = {csv: _load(csv) for csv in {c[1] for c in _CURVES}}

    # Distances present across all sources.
    distances = sorted({
        int(d)
        for df in loaded.values() if df is not None
        for d in df["expected_distance"].dropna().unique()
    })

    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    legend_handles = []

    for label, csv, noise_model, decoder, cmap in _CURVES:
        df = loaded.get(csv)
        if df is None:
            print(f"[skip] {label}: {csv} not found")
            continue
        sub = df[(df["noise_model"] == noise_model) & (df["decoder"] == decoder)]
        if len(sub) == 0:
            print(f"[skip] {label}: no rows")
            continue

        ls = "-" if "tesseract" in label else "--"
        shades = np.linspace(0.4, 0.92, max(len(distances), 1))
        for i, d in enumerate(distances):
            color = cmap(shades[i])
            s = sub[(sub["expected_distance"] == d) & (sub["shots"] > 0)].sort_values("p_bell_z")
            if len(s) == 0:
                continue
            pz = s["p_bell_z"].to_numpy(float)
            errors = s["errors"].to_numpy(float)
            shots = s["shots"].to_numpy(float)
            k = s["k"].astype(float).clip(lower=1).to_numpy()

            rate = errors / shots
            lo_shot, hi_shot = _wilson(errors, shots)
            y = _per_obs(rate, k)
            y_lo = _per_obs(lo_shot, k)
            y_hi = _per_obs(hi_shot, k)

            # Measured points (>=1 error) only: connected line + Wilson bars.
            # Zero-error points are dropped (their censored upper bounds cluttered
            # the plot); they carry no measured rate to place on a log axis.
            meas = y > 0
            if meas.any():
                ax.errorbar(
                    pz[meas], y[meas], yerr=[(y - y_lo)[meas], (y_hi - y)[meas]],
                    marker="o", markersize=4, linewidth=1.3, capsize=2,
                    elinewidth=0.8, color=color, linestyle=ls,
                )

            # One legend entry per (method, distance) with its actual shade.
            legend_handles.append(Line2D(
                [0], [0], color=color, marker="o", linestyle=ls,
                linewidth=1.4, label=f"{label}, d={d}"))

    ax.set_xscale("log")
    ax.set_yscale("log")
    if args.floor is not None:
        ax.set_ylim(bottom=args.floor)
    ax.set_xlabel(r"$p_z$ (Bell prep Z-noise)", fontsize=14)
    ax.set_ylabel("Logical error rate per observable", fontsize=13)
    ax.set_title("Full decoder (tesseract) vs 3-rep majority vote (pymatching)",
                 fontsize=13)
    ax.grid(True, which="both", alpha=0.3)
    # Legend outside (up to 4 methods x 3 distances = 12 entries); one column
    # per method so the shade gradient per method reads top-to-bottom.
    ax.legend(handles=legend_handles, fontsize=8, framealpha=0.9,
              loc="center left", bbox_to_anchor=(1.01, 0.5),
              ncol=1, title="method, distance")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {args.out}")


if __name__ == "__main__":
    main()
