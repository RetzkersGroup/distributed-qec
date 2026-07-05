#!/usr/bin/env python3
"""Paper-style plots for the honeycomb (p, p_z) grid scan.

Consumes the shared grid CSV produced by tools/run_grid_bell.sh
(out/collection/grid_bell,*.csv), where p (= p_res) and p_z (= p_bell_z) are two
independent axes and expected_distance is the code distance. The CSV may hold
several schemes at once (bell, majority-vote, ...), separated by the
noise_model column; each is fitted and drawn as its own boundary curve, and the
LER panel overlays every scheme at a shared p_res.

Modeled on the lattice-surgery `plot_ls.py` "paper plot": sweep the Z-noise
axis (p_z) at fixed residual noise (p_res), fit a threshold per distance-family
crossing, and draw:

  grid_bell_threshold_pz_vs_p.png  fitted p_z threshold vs p_res (2D boundary), per scheme
  grid_bell_ler_vs_pz.png          LER vs p_z at one shared p_res (per d, per scheme)
  grid_bell_paper_figure.pdf       stacked: (a) boundary on top, (b) LER below
  grid_bell_ler_debug/*.png        the LER panel, one figure per sampled p_res (debug)

Loads with pandas + ast (our json_metadata is a Python-dict repr, not JSON) and
reuses src/threshold_fit.py (the same locate -> quadratic-FSS estimator) fitting
along the p_z axis. Logical rate is per-round, per-observable so the
different-rounds distances are comparable.

Usage:
    python tools/plot_grid_bell.py
    python tools/plot_grid_bell.py --p-example 1e-3 --no-debug
"""

from __future__ import annotations

import argparse
import ast
import glob
import pathlib
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from threshold_fit import locate_threshold, fit_threshold  # noqa: E402

COLORS = {3: "#e15759", 5: "#4e79a7", 7: "#59a14f",
          9: "#f28e2b", 11: "#76b7b2", 13: "#edc948", 15: "#b07aa1"}
MARKERS = {3: "o", 5: "s", 7: "^", 9: "D", 11: "v", 13: "x", 15: "P"}

# One color / linestyle / short label per scheme (noise_model). Color drives the
# boundary overlay; linestyle distinguishes schemes in the LER overlay (where
# color already encodes distance).
SCHEME_COLORS = {
    "biased_bell_v2": "#e15759",
    "biased_bell_v3": "#4e79a7",
    "biased_bell_v3_equivalent": "#59a14f",
    "biased_bell_v3_equivalent_majority": "#b07aa1",
}
SCHEME_LINESTYLES = {
    "biased_bell_v2": ":",
    "biased_bell_v3": "-",
    "biased_bell_v3_equivalent": "-.",
    "biased_bell_v3_equivalent_majority": "--",
}
SCHEME_LABELS = {
    "biased_bell_v2": "bell v2",
    "biased_bell_v3": "bell v3",
    "biased_bell_v3_equivalent": "equivalent",
    "biased_bell_v3_equivalent_majority": "majority vote",
}


def scheme_label(name: str) -> str:
    return SCHEME_LABELS.get(name, name)


def scheme_linestyle(name: str) -> str:
    return SCHEME_LINESTYLES.get(name, "-")


_AGG = {
    "shots": "sum", "errors": "sum", "discards": "sum", "seconds": "sum",
    "decoder": "first", "json_metadata": "first",
}
_META = ["p", "p_bell_z", "expected_distance", "rounds", "k",
         "direction_kind", "noise_model"]


# ── Data ──────────────────────────────────────────────────────────────────────

def load_grid(csv_globs: str | list[str]) -> pd.DataFrame:
    """Load grid decode CSV(s): merge by strong_id, unpack metadata, add pL."""
    if isinstance(csv_globs, str):
        csv_globs = [csv_globs]
    files = sorted(f for g in csv_globs for f in glob.glob(g))
    if not files:
        raise FileNotFoundError(f"No files matched: {csv_globs}")

    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.groupby("strong_id").agg(_AGG).reset_index()

    meta = df["json_metadata"].apply(
        lambda s: ast.literal_eval(s) if isinstance(s, str) else dict(s)
    )
    for c in _META:
        df[c] = meta.apply(lambda m, c=c: m.get(c))

    d = df["expected_distance"].astype(float)
    rounds = df["rounds"].astype(float)
    k = df["k"].astype(float).clip(lower=1)
    shots = df["shots"].clip(lower=1)
    pL_total = df["errors"] / shots
    # Per-round, per-observable single-event rate (comparable across distances).
    scale = d / (rounds * k)
    df["pL"] = 1 - (1 - pL_total) ** scale
    df["sigma"] = (
        np.sqrt(np.clip(pL_total * (1 - pL_total), 0, None) / shots)
        * scale * (1 - pL_total) ** (scale - 1)
    )
    return df


def _pz_curve(df: pd.DataFrame, p: float, d: int):
    """Sorted (p_z, pL, sigma) for one p_res and distance, valid points only."""
    sub = df[(np.isclose(df["p"], p)) & (df["expected_distance"] == d)
             & (df["errors"] > 0) & (df["shots"] > 0)]
    sub = sub.sort_values("p_bell_z")
    return (sub["p_bell_z"].to_numpy(), sub["pL"].to_numpy(),
            sub["sigma"].to_numpy())


def _has_protective_regime(df: pd.DataFrame, p: float) -> bool:
    """True if the largest distance beats the smallest at the lowest p_z.

    If the biggest code is already *worse* than the smallest at the low-noise
    end, we are above threshold everywhere (no p_z crossing), so any fit is
    spurious. Mirrors the `show_threshold` guard in the reference plot_ls.py.
    """
    distances = sorted(df["expected_distance"].dropna().unique())
    d_lo, d_hi = int(distances[0]), int(distances[-1])
    xs_lo, ys_lo, _ = _pz_curve(df, p, d_lo)
    xs_hi, ys_hi, _ = _pz_curve(df, p, d_hi)
    if len(xs_lo) == 0 or len(xs_hi) == 0:
        return False
    pz0 = xs_hi[0]
    j = int(np.argmin(np.abs(xs_lo - pz0)))
    return ys_hi[0] < ys_lo[j]


def fit_pz_threshold(df: pd.DataFrame, p: float):
    """Fit p_z threshold at fixed p_res across all distances. None if no fit."""
    sub = df[np.isclose(df["p"], p) & (df["errors"] > 0) & (df["shots"] > 0)]
    if sub["expected_distance"].nunique() < 2:
        return None
    if not _has_protective_regime(df, p):
        return None
    pz = sub["p_bell_z"].to_numpy(float)
    d = sub["expected_distance"].to_numpy(float)
    pL = sub["pL"].to_numpy(float)
    sigma = sub["sigma"].to_numpy(float)

    p_init = locate_threshold(pz, d, pL)
    if p_init is None:
        return None
    try:
        fit = fit_threshold(pz, d, pL, sigma=sigma, p_init=p_init, n_per_side=2)
    except Exception:
        return None
    if not (pz.min() <= fit.p_th <= pz.max()):
        return None
    return fit


def collect_thresholds(df: pd.DataFrame):
    """Fit p_z threshold at every p_res; return sorted (p, p_th, p_th_err)."""
    out = []
    for p in sorted(df["p"].dropna().unique()):
        fit = fit_pz_threshold(df, float(p))
        if fit is not None and np.isfinite(fit.p_th):
            out.append((float(p), fit.p_th, fit.p_th_err))
    return out


def _scheme_has_p(sdf: pd.DataFrame, p_res: float) -> bool:
    return bool(np.any(np.isclose(sdf["p"].to_numpy(float), p_res)))


# ── Figures ───────────────────────────────────────────────────────────────────

def plot_ler_overlay(df: pd.DataFrame, p_res: float, schemes: list[str], ax=None):
    """LER vs p_z at fixed p_res, overlaying every scheme.

    Color encodes distance, linestyle encodes scheme; each scheme's fitted p_z
    threshold is drawn as a vertical line in the scheme color. Only schemes that
    actually sampled `p_res` are drawn.
    """
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(5.2, 4))

    distances_seen: set[int] = set()
    schemes_plotted: list[str] = []
    for s in schemes:
        sdf = df[df["noise_model"] == s]
        if not _scheme_has_p(sdf, p_res):
            continue
        ls = scheme_linestyle(s)
        for d in sorted(sdf["expected_distance"].dropna().unique()):
            d = int(d)
            xs, ys, es = _pz_curve(sdf, p_res, d)
            if len(xs) == 0:
                continue
            ax.errorbar(xs, ys, yerr=es, fmt="", marker=MARKERS.get(d, "o"),
                        linestyle=ls, color=COLORS.get(d), markersize=4,
                        linewidth=1.4, capsize=3, elinewidth=0.8)
            distances_seen.add(d)
        fit = fit_pz_threshold(sdf, p_res)
        if fit is not None:
            ax.axvline(fit.p_th, color=SCHEME_COLORS.get(s, "k"),
                       linestyle=ls, linewidth=1.2)
        schemes_plotted.append(s)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$p_z$", fontsize=11)
    ax.set_ylabel("Logical error rate per round", fontsize=11)
    ax.set_title(rf"$p_{{\text{{res}}}} = {p_res:.2g}$", fontsize=10)
    ax.grid(True, which="both", alpha=0.3)

    # Two-part legend: color = distance, linestyle = scheme.
    d_handles = [Line2D([0], [0], color=COLORS.get(d, "k"),
                        marker=MARKERS.get(d, "o"), linestyle="",
                        label=rf"$d={d}$")
                 for d in sorted(distances_seen)]
    s_handles = [Line2D([0], [0], color="k", linestyle=scheme_linestyle(s),
                        label=scheme_label(s))
                 for s in schemes_plotted]
    if d_handles or s_handles:
        ax.legend(handles=d_handles + s_handles, fontsize=7,
                  loc="lower right", framealpha=0.8)
    return fig if own else None


def plot_threshold_vs_p(thr_by_scheme: dict, ax=None):
    """Fitted p_z threshold vs p_res — the 2D boundary, one curve per scheme."""
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6, 5))
    for scheme, thresholds in thr_by_scheme.items():
        if not thresholds:
            continue
        p, pth, err = map(np.array, zip(*thresholds))
        ax.errorbar(p, pth, yerr=err, fmt="o-",
                    color=SCHEME_COLORS.get(scheme, "k"),
                    label=scheme_label(scheme),
                    markersize=6, linewidth=1.6, capsize=3)
    # Linear axes (lin-lin) for the threshold boundary.
    ax.set_xlabel(r"$p_{\text{res}}$", fontsize=11)
    ax.set_ylabel(r"$p_z^*$ (threshold)", fontsize=11)
    ax.set_title(r"$p_z$ threshold vs $p_{\text{res}}$", fontsize=10)
    ax.grid(True, alpha=0.3)
    if len(thr_by_scheme) > 1:
        ax.legend(fontsize=8)
    return fig if own else None


def plot_paper_figure(df: pd.DataFrame, p_example: float, thr_by_scheme: dict,
                      schemes: list[str]):
    """Stacked paper figure: (a) threshold boundary (top), (b) LER (bottom)."""
    fig = plt.figure(figsize=(5.2, 7.4), constrained_layout=True)
    ax_a = fig.add_subplot(2, 1, 1)
    ax_b = fig.add_subplot(2, 1, 2)
    plot_threshold_vs_p(thr_by_scheme, ax=ax_a)
    plot_ler_overlay(df, p_example, schemes, ax=ax_b)
    for ax, label in [(ax_a, "(a)"), (ax_b, "(b)")]:
        ax.text(-0.02, 1.03, label, transform=ax.transAxes, fontsize=12,
                fontweight="bold", va="bottom", ha="left")
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────

def _common_p(df: pd.DataFrame, schemes: list[str]) -> list[float]:
    """p_res values present in *every* scheme (so the LER overlay shows all)."""
    sets = [set(np.round(np.log10(df[df["noise_model"] == s]["p"].dropna()
                                   .unique()), 6))
            for s in schemes]
    common = set.intersection(*sets) if sets else set()
    return sorted(10.0 ** x for x in common)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--csv", default=str(_ROOT / "out/collection/grid_bell,*.csv"),
                        help="Glob for grid decode CSV(s).")
    parser.add_argument("--out-dir", type=pathlib.Path,
                        default=_ROOT / "out" / "figures",
                        help="Directory for output figures.")
    parser.add_argument("--p-example", type=float, default=1e-3,
                        help="Shared p_res for the overlaid LER panel.")
    parser.add_argument("--no-debug", action="store_true",
                        help="Skip the per-p_res debug LER figures.")
    args = parser.parse_args()

    df = load_grid(args.csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    schemes = sorted(df["noise_model"].dropna().unique())
    print(f"Loaded {len(df)} rows | schemes: {schemes}")
    for s in schemes:
        sdf = df[df["noise_model"] == s]
        print(f"  {s}: {len(sdf)} rows, p={sdf['p'].nunique()} vals, "
              f"p_z={sdf['p_bell_z'].nunique()} vals, "
              f"d={sorted(sdf['expected_distance'].dropna().unique())}")

    thr_by_scheme = {s: collect_thresholds(df[df["noise_model"] == s])
                     for s in schemes}
    for s in schemes:
        print(f"{s}: fitted {len(thr_by_scheme[s])} thresholds")

    # Shared p_res for the overlaid LER panel: nearest common value to --p-example.
    common = _common_p(df, schemes)
    if not common:
        raise SystemExit(f"schemes share no common p_res values: {schemes}")
    p_example = min(common, key=lambda v: abs(v - args.p_example))
    print(f"Overlaid LER panel at shared p_res = {p_example:.3g} "
          f"(common values: {[f'{v:.3g}' for v in common]})")

    fig_a = plot_threshold_vs_p(thr_by_scheme)
    out_a = args.out_dir / "grid_bell_threshold_pz_vs_p.png"
    fig_a.savefig(out_a, dpi=150, bbox_inches="tight"); plt.close(fig_a)
    print(f"Saved -> {out_a}")

    fig_b = plot_ler_overlay(df, p_example, schemes)
    out_b = args.out_dir / "grid_bell_ler_vs_pz.png"
    fig_b.savefig(out_b, dpi=150, bbox_inches="tight"); plt.close(fig_b)
    print(f"Saved -> {out_b}")

    fig_p = plot_paper_figure(df, p_example, thr_by_scheme, schemes)
    out_p = args.out_dir / "grid_bell_paper_figure.pdf"
    fig_p.savefig(out_p, bbox_inches="tight"); plt.close(fig_p)
    print(f"Saved -> {out_p}")

    # Debug: the LER overlay, one independent figure per sampled p_res.
    if not args.no_debug:
        dbg_dir = args.out_dir / "grid_bell_ler_debug"
        dbg_dir.mkdir(parents=True, exist_ok=True)
        all_p = sorted(df["p"].dropna().unique())
        for i, p_res in enumerate(all_p):
            fig = plot_ler_overlay(df, float(p_res), schemes)
            fig.savefig(dbg_dir / f"{i:02d}_ler_pres={p_res:.5g}.png",
                        dpi=130, bbox_inches="tight")
            plt.close(fig)
        print(f"Saved {len(all_p)} debug LER figures -> {dbg_dir}/")


if __name__ == "__main__":
    main()
