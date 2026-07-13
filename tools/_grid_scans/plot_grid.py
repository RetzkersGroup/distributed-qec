#!/usr/bin/env python3
"""Paper-style plots for the honeycomb (p, p_z) grid scan.

Consumes the shared grid CSV produced by tools/run_grid_bell.sh
(out/collection/grid_bell.csv), where p (= p_res) and p_z (= p_bell_z) are two
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

Loads sinter-format CSVs (JSON metadata; legacy pandas repr also handled) and
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
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({"xtick.labelsize": 15, "ytick.labelsize": 15})

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from threshold_fit import locate_threshold, fit_threshold  # noqa: E402

# Distance -> color/marker, covering both honeycomb (odd d=3,5,7) and hyperbolic
# (even d=2,4,6,8,10,12). Rank by (d//2 - 1) so each successive code distance
# gets the next palette entry; honeycomb 3/5/7 keep their original red/blue/green.
_PALETTE = ["#e15759", "#4e79a7", "#59a14f", "#f28e2b", "#76b7b2", "#edc948", "#b07aa1"]
_MARKS = ["o", "s", "^", "D", "v", "P", "X"]
COLORS = {d: _PALETTE[(d // 2 - 1) % len(_PALETTE)] for d in range(2, 16)}
MARKERS = {d: _MARKS[(d // 2 - 1) % len(_MARKS)] for d in range(2, 16)}

# One color / linestyle / short label per scheme (noise_model). Color drives the
# boundary overlay; linestyle distinguishes schemes in the LER overlay (where
# color already encodes distance).
SCHEME_COLORS = {
    "biased_bell_v2": "#e15759",
    "biased_bell_v3": "#4e79a7",
    "biased_bell_v3_equivalent": "#59a14f",
    "biased_bell_v3_equivalent_majority": "#b07aa1",
    "biased_bell_v3_equivalent_majority5": "#674ea7",
}
SCHEME_LINESTYLES = {
    "biased_bell_v2": ":",
    "biased_bell_v3": "-",
    "biased_bell_v3_equivalent": "-.",
    "biased_bell_v3_equivalent_majority": "--",
    "biased_bell_v3_equivalent_majority5": (0, (3, 1, 1, 1)),  # dash-dot-dot
}
SCHEME_LABELS = {
    "biased_bell_v2": "bell v2",
    "biased_bell_v3": "Honeycomb code",
    "biased_bell_v3_equivalent": "equivalent",
    "biased_bell_v3_equivalent_majority": "Honeycomb with 3 repetitions and majority vote",
    "biased_bell_v3_equivalent_majority5": "5-rep majority vote",
}


def scheme_label(name: str) -> str:
    return SCHEME_LABELS.get(name, name)


def scheme_linestyle(name: str) -> str:
    return SCHEME_LINESTYLES.get(name, "-")


# Short filename tag per scheme for the per-scheme debug figures.
SCHEME_TAGS = {
    "biased_bell_v2": "norepeat_v2",
    "biased_bell_v3": "norepeat",
    "biased_bell_v3_equivalent": "norepeat_equiv",
    "biased_bell_v3_equivalent_majority": "repeat",
    "biased_bell_v3_equivalent_majority5": "repeat5",
}


def scheme_tag(name: str) -> str:
    return SCHEME_TAGS.get(name, name)


_AGG = {
    "shots": "sum", "errors": "sum", "discards": "sum", "seconds": "sum",
    "decoder": "first", "json_metadata": "first",
}
_META = ["p", "p_bell_z", "expected_distance", "rounds", "k", "w", "h",
         "direction_kind", "noise_model", "input_circuit"]

# Fallbacks parsed from the input-circuit filename when the sidecar metadata
# doesn't carry d / k / n / rounds (e.g. the hyperbolic Floquet circuits, which
# have no .json). Honeycomb rows keep their explicit metadata (and get n=w*h).
_RE_D = re.compile(r"_d=(\d+)")
_RE_K = re.compile(r"_k=(\d+)")
_RE_N = re.compile(r"_n=(\d+)")
_RE_NSR = re.compile(r"num_sub_rounds=(\d+)")


def _fn_int(text, pat):
    m = pat.search(text or "")
    return int(m.group(1)) if m else None


def _rounds_from_fn(text):
    # Floquet: a QEC round is 6 sub-rounds, minus the boundary round.
    nsr = _fn_int(text, _RE_NSR)
    return (nsr - 6) // 6 if nsr else None


# ── Data ──────────────────────────────────────────────────────────────────────

def _parse_meta(s):
    """Parse a json_metadata cell: sinter-format (JSON) or legacy pandas (repr)."""
    if not isinstance(s, str):
        return dict(s)
    try:
        return json.loads(s)              # sinter's CSV (JSON)
    except (json.JSONDecodeError, ValueError):
        return ast.literal_eval(s)        # legacy pandas CSV (Python repr)


def load_grid(csv_globs: str | list[str]) -> pd.DataFrame:
    """Load grid decode CSV(s): merge by strong_id, unpack metadata, add pL."""
    if isinstance(csv_globs, str):
        csv_globs = [csv_globs]
    files = sorted(f for g in csv_globs for f in glob.glob(g))
    if not files:
        raise FileNotFoundError(f"No files matched: {csv_globs}")

    # skipinitialspace: sinter's CSV_HEADER pads column names with spaces.
    df = pd.concat([pd.read_csv(f, skipinitialspace=True) for f in files],
                   ignore_index=True)
    df = df.groupby("strong_id").agg(_AGG).reset_index()

    meta = df["json_metadata"].apply(_parse_meta)
    for c in _META:
        df[c] = meta.apply(lambda m, c=c: m.get(c))

    # Fill d / k / rounds from the input-circuit filename where missing.
    ic = df["input_circuit"].fillna("")
    for col, getter in (
        ("expected_distance", lambda s: _fn_int(s, _RE_D)),
        ("k", lambda s: _fn_int(s, _RE_K)),
        ("rounds", _rounds_from_fn),
    ):
        miss = df[col].isna()
        if miss.any():
            df.loc[miss, col] = ic[miss].apply(getter)

    # n (physical qubits) for the [[n, k, d]] labels: filename n= (hyperbolic)
    # or w*h (honeycomb).
    df["n"] = ic.apply(lambda s: _fn_int(s, _RE_N))
    miss = df["n"].isna()
    if miss.any():
        wh = df["w"].astype(float) * df["h"].astype(float)
        df.loc[miss, "n"] = wh[miss]

    k = df["k"].astype(float).clip(lower=1)
    shots = df["shots"].clip(lower=1)
    pL_total = df["errors"] / shots
    # Raw per-shot logical error rate (used by the debug figures).
    df["pL_raw"] = pL_total
    df["sigma_raw"] = np.sqrt(np.clip(pL_total * (1 - pL_total), 0, None) / shots)
    # Per-observable single-event rate: split the any-observable failure evenly
    # among the k observables (no round normalisation). Used everywhere else.
    scale = 1.0 / k
    df["pL"] = 1 - (1 - pL_total) ** scale
    df["sigma"] = (
        np.sqrt(np.clip(pL_total * (1 - pL_total), 0, None) / shots)
        * scale * (1 - pL_total) ** (scale - 1)
    )
    return df


def _code_label(df: pd.DataFrame, d: int) -> str:
    """[[n, k, d]] label for the code at distance d (n, k read from the data)."""
    row = df[df["expected_distance"] == d]
    n = row["n"].iloc[0] if not row.empty else None
    k = row["k"].iloc[0] if not row.empty else None
    ns = str(int(n)) if pd.notna(n) else "?"
    ks = str(int(k)) if pd.notna(k) else "?"
    return f"[[{ns}, {ks}, {int(d)}]]"


def _pz_curve(df: pd.DataFrame, p: float, d: int,
              rate_col: str = "pL", sigma_col: str = "sigma"):
    """Sorted (p_z, rate, sigma) for one p_res and distance, valid points only."""
    sub = df[(np.isclose(df["p"], p)) & (df["expected_distance"] == d)
             & (df["errors"] > 0) & (df["shots"] > 0)]
    sub = sub.sort_values("p_bell_z")
    return (sub["p_bell_z"].to_numpy(), sub[rate_col].to_numpy(),
            sub[sigma_col].to_numpy())


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


def _cutoff_for(df: pd.DataFrame) -> float:
    """Per-observable saturation cutoff for the fit, by scheme.

    The equivalent/majority schemes cross at a higher per-observable rate (their
    p_z threshold is large, so the k=10 rate is already sizeable there), so a 0.3
    cutoff would clip points at the crossing; use 0.45 for them, 0.3 otherwise.
    Assumes `df` is a single scheme (fit_pz_threshold is always called per-scheme).

    """
    nm = df["noise_model"].dropna()
    return 0.45 if (len(nm) and "equivalent" in str(nm.iloc[0])) else 0.3


def _pz_cutoff_for(df: pd.DataFrame) -> float:
    """Saturation cutoff for the *p_z*-axis FSS fit specifically.

    Same as `_cutoff_for` except for k=1 (honeycomb) 5-rep: there the per-obs
    rate *is* the raw rate and saturates near 0.5, and the 5-rep p_z crossing
    (~11%) sits right below the ceiling, so 0.45 strips almost every
    above-crossing point and the quadratic overshoots. Keeping points up to 0.495
    anchors the fit to the crossing (~0.112, matching the majority-5 prediction).
    Scoped to the p_z fit so it never perturbs the p_res-axis fit (whose crossing
    is nowhere near saturation).
    """
    nm = df["noise_model"].dropna()
    s = str(nm.iloc[0]) if len(nm) else ""
    if "equivalent" not in s:
        return 0.3
    k = df["k"].dropna()
    if "majority5" in s and len(k) and int(k.iloc[0]) == 1:
        return 0.495
    return 0.45


def _pz_max_for(df: pd.DataFrame) -> float:
    """FSS p_z window cap, by scheme: the 5-rep majority threshold sits ~11-13%,
    above the 0.09 default, so it needs more headroom (others keep 0.09)."""
    nm = df["noise_model"].dropna()
    return 0.25 if (len(nm) and "majority5" in str(nm.iloc[0])) else 0.09


def fit_pz_threshold(df: pd.DataFrame, p: float, n_per_side: int = 2,
                     pz_max: float | None = None, logical_cutoff: float | None = None,
                     rate_col: str = "pL", sigma_col: str = "sigma"):
    """Fit p_z threshold at fixed p_res across all distances. None if no fit.

    ``n_per_side`` sets the FSS fit window: n_per_side p_z values on each side of
    the coarse threshold plus the center (default 2 -> 5 p_z values). Points with
    p_z above ``pz_max`` are excluded so the crossing is never pinned against the
    grid ceiling. Points whose per-observable rate exceeds ``logical_cutoff``
    (default: scheme-aware, see `_cutoff_for`) are dropped so the min-variance
    locate isn't seeded by the saturated high-p_z plateau.
    """
    if logical_cutoff is None:
        logical_cutoff = _pz_cutoff_for(df)
    if pz_max is None:
        pz_max = _pz_max_for(df)
    sub = df[np.isclose(df["p"], p) & (df["errors"] > 0) & (df["shots"] > 0)
             & (df["p_bell_z"] <= pz_max) & (df["pL"] < logical_cutoff)]
    if sub["expected_distance"].nunique() < 2:
        return None
    if not _has_protective_regime(df, p):
        return None
    pz = sub["p_bell_z"].to_numpy(float)
    d = sub["expected_distance"].to_numpy(float)
    pL = sub[rate_col].to_numpy(float)
    sigma = sub[sigma_col].to_numpy(float)

    p_init = locate_threshold(pz, d, pL)
    if p_init is None:
        return None
    try:
        fit = fit_threshold(pz, d, pL, sigma=sigma, p_init=p_init,
                            n_per_side=n_per_side)
    except Exception:
        return None
    if not (pz.min() <= fit.p_th <= pz.max()):
        return None
    return fit


# Threshold estimator selected across the whole module: "fss" (quadratic
# finite-size scaling, default) or "minvar" (plateau-relative min-variance point).
THRESHOLD_METHOD = "fss"


def pz_threshold_min_variance(df: pd.DataFrame, p: float, rate_col: str = "pL"):
    """Min-variance p_z crossing at fixed p_res, with a plateau-relative guard.

    At each p_z, take the mean logical error and the variance of log(rate) across
    distances. The maximal error rate is the largest mean (the saturated plateau);
    the threshold is the p_z of minimum log-variance among points whose mean is
    below half that maximum. Returns the point (no error bar), or None. The
    protective-regime guard makes the boundary end where there is no crossing.
    """
    if not _has_protective_regime(df, p):
        return None
    sub = df[np.isclose(df["p"], p) & (df["errors"] > 0) & (df["shots"] > 0)]
    xs, means, varlog = [], [], []
    for x, grp in sub.groupby("p_bell_z"):
        r = grp[rate_col].to_numpy(float)
        r = r[r > 0]
        if len(r) < 2:
            continue
        xs.append(float(x))
        means.append(float(r.mean()))
        varlog.append(float(np.var(np.log(r))))
    if not xs:
        return None
    xs, means, varlog = np.array(xs), np.array(means), np.array(varlog)
    valid = means < 0.5 * means.max()
    if not valid.any():
        return None
    return float(xs[np.argmin(np.where(valid, varlog, np.inf))])


def estimate_pz_threshold(df: pd.DataFrame, p: float, rate_col: str = "pL",
                          sigma_col: str = "sigma"):
    """Dispatch the active threshold estimator; return (p_th, p_th_err) or None."""
    if THRESHOLD_METHOD == "minvar":
        pth = pz_threshold_min_variance(df, p, rate_col=rate_col)
        return (pth, 0.0) if pth is not None else None
    fit = fit_pz_threshold(df, p, rate_col=rate_col, sigma_col=sigma_col)
    return (fit.p_th, fit.p_th_err) if fit is not None else None


def collect_thresholds(df: pd.DataFrame, rate_col: str = "pL",
                       sigma_col: str = "sigma"):
    """Estimate p_z threshold at every p_res; return sorted (p, p_th, p_th_err)."""
    out = []
    for p in sorted(df["p"].dropna().unique()):
        res = estimate_pz_threshold(df, float(p), rate_col=rate_col, sigma_col=sigma_col)
        if res is not None and np.isfinite(res[0]):
            out.append((float(p), res[0], res[1]))
    return out


def _scheme_has_p(sdf: pd.DataFrame, p_res: float) -> bool:
    return bool(np.any(np.isclose(sdf["p"].to_numpy(float), p_res)))


# ── Figures ───────────────────────────────────────────────────────────────────

def plot_ler_overlay(df: pd.DataFrame, p_res: float, schemes: list[str], ax=None,
                     xlim=None, ylim=None, show_title=False, per_observable=True):
    """LER vs p_z at fixed p_res, overlaying every scheme.

    Color encodes distance, linestyle encodes scheme; each scheme's fitted p_z
    threshold is drawn as a vertical line in the scheme color. Only schemes that
    actually sampled `p_res` are drawn. With ``per_observable=False`` the raw
    per-shot logical error rate is plotted instead (used by the debug figures);
    the threshold vline is unchanged (it's a p_z position, y-axis independent).
    """
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(5.2, 4))
    rate_col, sigma_col = ("pL", "sigma") if per_observable else ("pL_raw", "sigma_raw")

    distances_seen: set[int] = set()
    schemes_plotted: list[str] = []
    for s in schemes:
        sdf = df[df["noise_model"] == s]
        if not _scheme_has_p(sdf, p_res):
            continue
        ls = scheme_linestyle(s)
        for d in sorted(sdf["expected_distance"].dropna().unique()):
            d = int(d)
            xs, ys, es = _pz_curve(sdf, p_res, d, rate_col, sigma_col)
            if len(xs) == 0:
                continue
            ax.errorbar(xs, ys, yerr=es, fmt="", marker=MARKERS.get(d, "o"),
                        linestyle=ls, color=COLORS.get(d), markersize=4,
                        linewidth=1.4, capsize=3, elinewidth=0.8)
            distances_seen.add(d)
        res = estimate_pz_threshold(sdf, p_res, rate_col=rate_col, sigma_col=sigma_col)
        if res is not None:
            ax.axvline(res[0], color=SCHEME_COLORS.get(s, "k"),
                       linestyle=ls, linewidth=1.2)
        schemes_plotted.append(s)

    ax.set_xscale("log")
    ax.set_yscale("log")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xlabel(r"$p_z$", fontsize=17)
    ax.set_ylabel("Logical error rate per observable" if per_observable
                  else "Logical error rate", fontsize=17)
    if show_title:
        title = rf"$p_{{\text{{res}}}} = {p_res:.2g}$"
        if len(schemes_plotted) == 1:
            title += f"  [{scheme_tag(schemes_plotted[0])}]"
        ax.set_title(title, fontsize=16)
    ax.grid(True, which="both", alpha=0.3)

    # Two-part legend: color = distance, linestyle = scheme.
    d_handles = [Line2D([0], [0], color=COLORS.get(d, "k"),
                        marker=MARKERS.get(d, "o"), linestyle="",
                        label=_code_label(df, d))
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
        ls = "--" if scheme == "biased_bell_v3" else "-"
        ax.errorbar(p, pth, yerr=err, marker="o", linestyle=ls,
                    color=SCHEME_COLORS.get(scheme, "k"),
                    label=scheme_label(scheme),
                    markersize=6, linewidth=1.6, capsize=3)
    # Linear axes (lin-lin) for the threshold boundary.
    ax.set_xlabel(r"$p_{\text{res}}$", fontsize=17)
    ax.set_ylabel(r"$p_z^*$ (threshold)", fontsize=17)
    ax.grid(True, alpha=0.3)
    if len(thr_by_scheme) > 1:
        ax.legend(fontsize=8)
    return fig if own else None


def plot_paper_figure(df: pd.DataFrame, p_example: float, thr_by_scheme: dict,
                      ler_schemes: list[str], per_observable: bool = True):
    """Stacked paper figure: (a) threshold boundary (top), (b) LER (bottom).

    The boundary panel shows every fitted scheme; the LER panel is restricted to
    ``ler_schemes``.
    """
    fig = plt.figure(figsize=(5.2, 7.4), constrained_layout=True)
    ax_a = fig.add_subplot(2, 1, 1)
    ax_b = fig.add_subplot(2, 1, 2)
    plot_threshold_vs_p(thr_by_scheme, ax=ax_a)
    plot_ler_overlay(df, p_example, ler_schemes, ax=ax_b, per_observable=per_observable)
    for ax, label in [(ax_a, "(a)"), (ax_b, "(b)")]:
        ax.text(-0.02, 1.03, label, transform=ax.transAxes, fontsize=18,
                fontweight="bold", va="bottom", ha="left")
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────

def _select_ler_schemes(schemes: list[str], requested: list[str] | None) -> list[str]:
    """Which schemes to draw in the LER panel.

    With ``requested`` (substrings), keep matching schemes; otherwise default to
    every scheme except the majority-vote one.
    """
    if requested:
        sel = [s for s in schemes if any(r in s for r in requested)]
        return sel or schemes
    non_maj = [s for s in schemes if "majority" not in s]
    return non_maj or schemes


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
    parser.add_argument("--csv", default=str(_ROOT / "out/collection/grid_bell.csv"),
                        help="Glob for grid decode CSV(s).")
    parser.add_argument("--prefix", default="grid_bell",
                        help="Basename prefix for the output figures "
                             "(e.g. grid_hyperbolic). Default: grid_bell.")
    parser.add_argument("--code-label", default="Honeycomb code",
                        help="Legend label for the base biased_bell_v3 code "
                             "(e.g. 'Hyperbolic Floquet code').")
    parser.add_argument("--majority-label",
                        default="Honeycomb with 3 repetitions and majority vote",
                        help="Legend label for the majority-vote scheme.")
    parser.add_argument("--out-dir", type=pathlib.Path,
                        default=_ROOT / "out" / "figures",
                        help="Directory for output figures.")
    parser.add_argument("--p-example", type=float, default=1e-3,
                        help="Shared p_res for the overlaid LER panel.")
    parser.add_argument("--ler-schemes", nargs="*", default=None,
                        help="Scheme-name substrings to show in the LER panel "
                             "(default: all schemes except majority vote). The "
                             "boundary panel and debug figures always show all.")
    parser.add_argument("--no-debug", action="store_true",
                        help="Skip the per-p_res debug LER figures.")
    parser.add_argument("--raw", action="store_true",
                        help="Plot and fit on the raw (P any-observable-fails) "
                             "rate instead of the per-observable rate.")
    parser.add_argument("--method", choices=["minvar", "fss"], default="fss",
                        help="Threshold estimator: minvar = plateau-relative "
                             "min-variance crossing (default); fss = quadratic FSS.")
    args = parser.parse_args()
    global THRESHOLD_METHOD
    THRESHOLD_METHOD = args.method
    # The non-repeated baseline is now the single-ancilla 1-rep equivalent.
    SCHEME_LABELS["biased_bell_v3_equivalent"] = args.code_label
    SCHEME_LABELS["biased_bell_v3_equivalent_majority"] = args.majority_label

    # Raw vs per-observable: which rate column feeds both plotting and fitting.
    per_obs = not args.raw
    rate_col, sigma_col = ("pL", "sigma") if per_obs else ("pL_raw", "sigma_raw")

    df = load_grid(args.csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    schemes = sorted(df["noise_model"].dropna().unique())
    print(f"Loaded {len(df)} rows | schemes: {schemes}")
    for s in schemes:
        sdf = df[df["noise_model"] == s]
        print(f"  {s}: {len(sdf)} rows, p={sdf['p'].nunique()} vals, "
              f"p_z={sdf['p_bell_z'].nunique()} vals, "
              f"d={sorted(sdf['expected_distance'].dropna().unique())}")

    thr_by_scheme = {s: collect_thresholds(df[df["noise_model"] == s],
                                           rate_col=rate_col, sigma_col=sigma_col)
                     for s in schemes}
    for s in schemes:
        print(f"{s}: fitted {len(thr_by_scheme[s])} thresholds")

    # Schemes shown in the LER panel (default: all but majority vote).
    ler_schemes = _select_ler_schemes(schemes, args.ler_schemes)

    # Shared p_res for the LER panel: nearest common value to --p-example.
    common = _common_p(df, ler_schemes)
    if not common:
        raise SystemExit(f"LER schemes share no common p_res values: {ler_schemes}")
    p_example = min(common, key=lambda v: abs(v - args.p_example))
    print(f"LER panel schemes: {ler_schemes} | shared p_res = {p_example:.3g}")

    fig_a = plot_threshold_vs_p(thr_by_scheme)
    out_a = args.out_dir / f"{args.prefix}_threshold_pz_vs_p.png"
    fig_a.savefig(out_a, dpi=150, bbox_inches="tight"); plt.close(fig_a)
    print(f"Saved -> {out_a}")

    fig_b = plot_ler_overlay(df, p_example, ler_schemes, per_observable=per_obs)
    out_b = args.out_dir / f"{args.prefix}_ler_vs_pz.png"
    fig_b.savefig(out_b, dpi=150, bbox_inches="tight"); plt.close(fig_b)
    print(f"Saved -> {out_b}")

    fig_p = plot_paper_figure(df, p_example, thr_by_scheme, ler_schemes,
                              per_observable=per_obs)
    out_p = args.out_dir / f"{args.prefix}_paper_figure.pdf"
    fig_p.savefig(out_p, bbox_inches="tight"); plt.close(fig_p)
    print(f"Saved -> {out_p}")

    # Debug: per-scheme LER figures, one per sampled p_res, with fixed axes so
    # they flip through / compare cleanly. Repeat and non-repeat schemes get a
    # distinct filename prefix and their own index.
    if not args.no_debug:
        dbg_dir = args.out_dir / f"{args.prefix}_ler_debug"
        dbg_dir.mkdir(parents=True, exist_ok=True)
        for old in dbg_dir.glob("*.png"):
            old.unlink()
        xlim = (float(df["p_bell_z"].min()), float(df["p_bell_z"].max()))
        ylim = (1e-5, 1.0)
        total = 0
        for s in schemes:
            for i, p_res in enumerate(sorted(df[df["noise_model"] == s]["p"]
                                             .dropna().unique())):
                # Two panels: per-observable (left, where the crossing is clear)
                # and overall raw P(any-fail) (right, which saturates at 1).
                fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.4, 4.2),
                                               constrained_layout=True)
                plot_ler_overlay(df, float(p_res), [s], ax=axL,
                                 xlim=xlim, ylim=ylim, per_observable=True)
                plot_ler_overlay(df, float(p_res), [s], ax=axR,
                                 xlim=xlim, ylim=ylim, per_observable=False)
                fig.suptitle(rf"$p_{{\text{{res}}}} = {p_res:.2g}$  [{scheme_tag(s)}]",
                             fontsize=17)
                fig.savefig(dbg_dir / f"{scheme_tag(s)}_{i:02d}_ler_pres={p_res:.5g}.png",
                            dpi=130, bbox_inches="tight")
                plt.close(fig)
                total += 1
        print(f"Saved {total} debug LER figures (per-observable | overall) "
              f"-> {dbg_dir}/")


if __name__ == "__main__":
    main()
