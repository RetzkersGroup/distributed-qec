"""Load, fit, and plot Bell-mediated threshold-scan data.

Shared analysis backend for the Bell threshold notebooks
(`notebooks/threshold_comparison_bell.ipynb`,
`notebooks/threshold_vs_bias_bell.ipynb`) and the standalone
`tools/plot_thresholds_bell.py` driver. Keeping the CSV-loading,
threshold-fitting, and plotting logic in one place means the notebooks and
the CLI produce byte-for-byte identical figures.

Pipeline (per experiment family):

    decode CSVs (out/collection/threshold_*,*.csv)
        -> load_decode_csvs()   # merge by strong_id, parse metadata
        -> fit_thresholds()     # locate -> quadratic-FSS fit per direction
        -> plot_* ()            # matplotlib Figures, saved by the caller

The figures mirror the notebook cells exactly:

    plot_direction_curves    -> per-(model, direction) p_L-vs-p grid
    plot_threshold_plane     -> headline (p, p_bell_z) bias-scan curve
    plot_threshold_vs_fixed_pbz -> p* vs fixed p_bell_z
    plot_combined_plane      -> both scans on the (p_bell_z, p) plane
"""

from __future__ import annotations

import ast
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from threshold_fit import locate_threshold, fit_threshold

# The two Bell noise models compared side by side in the comparison notebook.
DEFAULT_NOISE_MODELS = ["biased_bell_v2", "biased_bell_v3"]

# Stable color per noise model (shared by every figure). Union of the palettes
# used across the two notebooks so a single dict drives all plots.
NOISE_MODEL_COLORS = {
    "biased_bell_v2": "#e15759",  # red-ish
    "biased_bell_v3": "#4e79a7",  # blue-ish
    "biased_bell_v3_equivalent": "#59a14f",  # green-ish
    "biased_bell_v3_equivalent_majority": "#b07aa1",  # purple-ish
}

# Marker per code distance for the per-direction curves.
DISTANCE_MARKERS = {2: "o", 3: "o", 4: "s", 5: "s", 6: "^", 7: "^"}

# Metadata columns pulled out of each row's json_metadata blob.
_META_COLUMNS = [
    "p", "p_bell_z", "direction_kind", "direction_value",
    "bias", "p_bell_z_fixed", "expected_distance", "rounds", "k",
    "noise_model",
]

# sinter stat columns aggregated when merging repeated runs of one circuit.
_AGG = {
    "shots": "sum",
    "errors": "sum",
    "discards": "sum",
    "seconds": "sum",
    "decoder": "first",
    "json_metadata": "first",
}


def _parse_metadata(meta):
    return ast.literal_eval(meta) if isinstance(meta, str) else dict(meta)


def _direction_label(row) -> str:
    """One string per scan direction, used for grouping, titles, and filenames."""
    if row["direction_kind"] == "bias":
        b = row["direction_value"]
        return f"bias={int(b) if b == int(b) else b:g}"
    return f'pbz={row["direction_value"]:.4g}'


def load_decode_csvs(
    csv_globs: str | list[str],
    noise_models: list[str] | None = None,
) -> pd.DataFrame:
    """Load one or more decode CSVs into a tidy per-circuit DataFrame.

    Merges repeated runs of the same circuit by ``strong_id`` (summing shots
    and errors), unpacks ``json_metadata`` into columns, and attaches a
    ``direction`` label. When ``noise_models`` is given, rows are filtered to
    those models; otherwise every model present is kept (each CSV's metadata
    already carries the right ``noise_model``, so mixed families separate
    themselves downstream).

    Raises FileNotFoundError if no file matches any glob.
    """
    if isinstance(csv_globs, str):
        csv_globs = [csv_globs]
    files = sorted(f for g in csv_globs for f in glob.glob(g))
    if not files:
        raise FileNotFoundError(f"No files matched: {csv_globs}")

    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.groupby("strong_id").agg(_AGG).reset_index()

    meta = df["json_metadata"].apply(_parse_metadata)
    for col in _META_COLUMNS:
        df[col] = meta.apply(lambda m, c=col: m.get(c))

    if noise_models is not None:
        df = df[df["noise_model"].isin(noise_models)].reset_index(drop=True)

    df["direction"] = df.apply(_direction_label, axis=1)
    return df


def _per_round_logical_rate(grp: pd.DataFrame, logical_cutoff: float):
    """Convert per-circuit (shots, errors) to per-round, per-observable p_L.

    Returns (p_arr, d_arr, pL, sigma) with saturated / empty points masked to
    NaN. Matches the conversion used in both Bell notebooks.
    """
    p_arr = grp["p"].to_numpy(dtype=float)
    d_arr = grp["expected_distance"].to_numpy(dtype=float)
    rounds_arr = grp["rounds"].to_numpy(dtype=float)
    k_arr = grp["k"].to_numpy(dtype=float)
    shots_arr = grp["shots"].to_numpy(dtype=float)
    errors_arr = grp["errors"].to_numpy(dtype=float)

    valid = (
        (shots_arr > 0)
        & (errors_arr > 0)
        & (errors_arr < shots_arr * logical_cutoff)
    )
    pL_total = np.where(valid, errors_arr / np.maximum(shots_arr, 1), np.nan)
    sigma_total = np.where(
        valid,
        np.sqrt(
            np.clip(pL_total * (1 - pL_total), 0, None)
            / np.maximum(shots_arr, 1)
        ),
        np.nan,
    )

    # Per-round, per-observable single-event rate.
    scale = d_arr / (rounds_arr * np.maximum(k_arr, 1))
    pL = 1 - (1 - pL_total) ** scale
    sigma = sigma_total * scale * (1 - pL_total) ** (scale - 1)
    return p_arr, d_arr, pL, sigma


def fit_thresholds(
    df: pd.DataFrame,
    *,
    n_per_side: int = 2,
    logical_cutoff: float = 0.3,
    sort: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """Fit a threshold per ``(noise_model, direction)`` group.

    Same pipeline as ``threshold_3d_surface.ipynb``: convert to a per-round,
    per-observable logical rate, coarse-locate the crossing, then run the
    quadratic finite-size-scaling refinement (`threshold_fit.fit_threshold`).

    Parameters mirror the notebook knobs: ``n_per_side`` sets the fit window
    (``2`` -> 5 p-values), ``logical_cutoff`` drops saturated points (the
    comparison notebook uses ``0.3`` for v2/v3; the equivalents use ``0.45``).

    Returns a DataFrame with one row per fitted direction.
    """
    fit_rows = []
    for (noise_model, direction_name), grp in df.groupby(["noise_model", "direction"]):
        p_arr, d_arr, pL, sigma = _per_round_logical_rate(grp, logical_cutoff)

        p_init = locate_threshold(p_arr, d_arr, pL)
        if p_init is None:
            continue
        try:
            fit = fit_threshold(
                p_arr, d_arr, pL,
                sigma=sigma,
                p_init=p_init,
                n_per_side=n_per_side,
            )
        except Exception as exc:
            if verbose:
                print(f"[skip] {noise_model} / {direction_name}: {exc}")
            continue
        if fit is None:
            continue

        fit_rows.append({
            "noise_model": noise_model,
            "direction": direction_name,
            "direction_kind": grp["direction_kind"].iloc[0],
            "bias": grp["bias"].iloc[0],
            "p_bell_z_fixed": grp["p_bell_z_fixed"].iloc[0],
            "p_th": fit.p_th,
            "p_th_err": fit.p_th_err,
            "mu": fit.mu,
            "mu_err": fit.mu_err,
            "p_init": fit.p_init,
            "n_points": fit.n_points,
        })

    fits_df = pd.DataFrame(fit_rows)
    if sort and not fits_df.empty:
        # noise_model, then biases first (numeric), then fixed pbz (numeric).
        fits_df["_sort"] = (
            fits_df["direction_kind"].map({"bias": 0, "pbz": 1}) * 1e9
            + fits_df["bias"].fillna(0).astype(float)
            + fits_df["p_bell_z_fixed"].fillna(0).astype(float) * 1e3
        )
        fits_df = (
            fits_df.sort_values(["noise_model", "_sort"])
            .drop(columns="_sort")
            .reset_index(drop=True)
        )
    return fits_df


def _direction_sort_key(df: pd.DataFrame, direction: str):
    sub = df[df["direction"] == direction].iloc[0]
    if sub["direction_kind"] == "bias":
        return (0, float(sub["bias"]))
    return (1, float(sub["p_bell_z_fixed"]))


def plot_direction_curves(
    df: pd.DataFrame,
    fits_df: pd.DataFrame,
    *,
    noise_models: list[str] | None = None,
    ncols: int = 3,
) -> plt.Figure:
    """Per-(model, direction) logical-error-rate-vs-p grid.

    Color encodes the noise model, marker encodes the distance, and each
    fitted threshold is drawn as a vertical dashed line. Mirrors the plotting
    cell of ``threshold_comparison_bell.ipynb``.
    """
    if noise_models is None:
        noise_models = DEFAULT_NOISE_MODELS

    direction_names = sorted(
        df["direction"].unique(),
        key=lambda d: _direction_sort_key(df, d),
    )
    present_models = [m for m in noise_models if m in df["noise_model"].unique()]

    fits_lookup = {}
    if not fits_df.empty:
        fits_lookup = {
            (row["noise_model"], row["direction"]): row.to_dict()
            for _, row in fits_df.iterrows()
        }

    n = len(direction_names)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False
    )

    for idx, direction_name in enumerate(direction_names):
        ax = axes[idx // ncols][idx % ncols]

        for model in present_models:
            sub = df[(df["noise_model"] == model) & (df["direction"] == direction_name)]
            if sub.empty:
                continue
            color = NOISE_MODEL_COLORS.get(model, "k")

            scale = sub["expected_distance"] / (
                sub["rounds"] * sub["k"].clip(lower=1)
            )
            pL_total = sub["errors"] / sub["shots"].clip(lower=1)
            pL = 1 - (1 - pL_total) ** scale
            sigma_total = np.sqrt(
                np.clip(pL_total * (1 - pL_total), 0, None)
                / sub["shots"].clip(lower=1)
            )
            sigma = sigma_total * scale * (1 - pL_total) ** (scale - 1)

            for d, sub_d in sub.groupby("expected_distance"):
                mask = sub_d.index
                ax.errorbar(
                    sub_d["p"], pL.loc[mask], yerr=sigma.loc[mask],
                    fmt=DISTANCE_MARKERS.get(int(d), "o") + "-",
                    color=color, alpha=0.85, markersize=4,
                    label=f"{model} d={int(d)}",
                )

            fit = fits_lookup.get((model, direction_name))
            if fit is not None and np.isfinite(fit["p_th"]):
                ax.axvline(fit["p_th"], color=color, linestyle="--",
                           linewidth=1, alpha=0.7)

        parts = [direction_name]
        for model in present_models:
            fit = fits_lookup.get((model, direction_name))
            if fit is not None and np.isfinite(fit["p_th"]):
                parts.append(
                    f"{model.replace('biased_bell_', '')}: "
                    f"$p_{{th}}={fit['p_th']:.3g}$"
                )
        ax.set_title("\n".join(parts))

        ax.loglog()
        ax.set_xlabel("Physical Error Rate (p)")
        ax.set_ylabel("Logical Error Rate per round")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    for idx in range(len(direction_names), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.tight_layout()
    return fig


def _bias_pbz_by_model(fits_df: pd.DataFrame):
    models = sorted(
        fits_df["noise_model"].unique(),
        key=lambda m: list(NOISE_MODEL_COLORS).index(m)
        if m in NOISE_MODEL_COLORS else 99,
    )
    bias_by_model = {
        m: fits_df[(fits_df["noise_model"] == m)
                   & (fits_df["direction_kind"] == "bias")]
        .sort_values("bias").reset_index(drop=True)
        for m in models
    }
    pbz_by_model = {
        m: fits_df[(fits_df["noise_model"] == m)
                   & (fits_df["direction_kind"] == "pbz")]
        .sort_values("p_bell_z_fixed").reset_index(drop=True)
        for m in models
    }
    return models, bias_by_model, pbz_by_model


def plot_threshold_plane(fits_df: pd.DataFrame) -> plt.Figure:
    """Headline bias-scan figure: each fit -> point (p*, beta * p*).

    One curve per noise model on the (p, p_bell_z) plane. v3, its
    single-ancilla equivalent, and the majority variant should overlap.
    Mirrors the saved cell of ``threshold_vs_bias_bell.ipynb``.
    """
    models, bias_by_model, _ = _bias_pbz_by_model(fits_df)

    fig, ax = plt.subplots(figsize=(7, 5))
    for model in models:
        bdf = bias_by_model[model]
        if bdf.empty:
            continue
        p = bdf["p_th"].to_numpy(dtype=float)
        pbz = bdf["bias"].to_numpy(dtype=float) * p
        ax.plot(p, pbz, "o-", color=NOISE_MODEL_COLORS.get(model, "k"),
                label=model)

    ax.set_xlabel(r"$p_\text{res}$")
    ax.set_ylabel(r"$p_{\mathrm{bell},z}$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_threshold_vs_fixed_pbz(fits_df: pd.DataFrame) -> plt.Figure | None:
    """p* vs fixed p_bell_z (independent hardware scan), one curve per model.

    Returns None if no fixed-pbz fits are present. Mirrors the corresponding
    cell of ``threshold_vs_bias_bell.ipynb``.
    """
    models, _, pbz_by_model = _bias_pbz_by_model(fits_df)

    fig, ax = plt.subplots(figsize=(7, 5))
    plotted = False
    for model in models:
        pdf = pbz_by_model[model]
        if pdf.empty:
            continue
        plotted = True
        ax.errorbar(
            pdf["p_bell_z_fixed"], pdf["p_th"], yerr=pdf["p_th_err"],
            fmt="s-", color=NOISE_MODEL_COLORS.get(model, "k"),
            label=f"{model}",
        )
    if not plotted:
        plt.close(fig)
        return None

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$p_{\mathrm{bell},z}$ (fixed)")
    ax.set_ylabel(r"Threshold $p^*$")
    ax.set_title(r"Threshold vs fixed $p_{\mathrm{bell},z}$")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_combined_plane(fits_df: pd.DataFrame) -> plt.Figure:
    """Both scans on one (p_bell_z, p) plane, one set of curves per model.

    bias scan -> (p_th * beta, p_th); fixed-pbz scan -> (p_bell_z_fixed, p_th).
    Mirrors the final cell of ``threshold_vs_bias_bell.ipynb``.
    """
    models, bias_by_model, pbz_by_model = _bias_pbz_by_model(fits_df)

    fig, ax = plt.subplots(figsize=(7, 5))
    for model in models:
        color = NOISE_MODEL_COLORS.get(model, "k")
        bdf = bias_by_model[model]
        if not bdf.empty:
            bx = bdf["bias"] * bdf["p_th"]
            ax.errorbar(
                bx, bdf["p_th"], yerr=bdf["p_th_err"],
                fmt="o-", color=color, label=f"{model} bias scan",
            )
        pdf = pbz_by_model[model]
        if not pdf.empty:
            ax.errorbar(
                pdf["p_bell_z_fixed"], pdf["p_th"], yerr=pdf["p_th_err"],
                fmt="s--", color=color, label=f"{model} fixed-pbz",
            )

    ax.plot(
        [1e-6, 1.0], [1e-6, 1.0],
        color="gray", linestyle=":", alpha=0.5,
        label=r"$p^* = p_{\mathrm{bell},z}$",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$p_{\mathrm{bell},z}$ at threshold")
    ax.set_ylabel(r"Threshold $p^*$")
    ax.set_title(r"Threshold trajectories on the $(p_{\mathrm{bell},z}, p)$ plane")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig
