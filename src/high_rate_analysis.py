"""Load and plot high-rate code decode data (hyperbolic Floquet + stairway).

Shared backend for `notebooks/high_rate.ipynb`,
`notebooks/high_rate_stairway.ipynb`, and `tools/plot_high_rate.py`. Reads the
raw decode CSVs written by `tools/high_rate/run_decode_*.sh`, merges repeated
runs, resolves each `[[n, k, d]]` code, and draws logical-error-rate-vs-p panels
(via `sinter.plot_error_rate`) with an optional secondary `p_bell_z` axis.

Unlike `bell_threshold_analysis` (which fits thresholds), this family just plots
the raw error-rate curves at a fixed Z-bias `p_bell_z = bias * p`.
"""

from __future__ import annotations

import ast
import json
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sinter

# n / k / d parsed from a hyperbolic input-circuit filename, used as a fallback
# when the sidecar metadata doesn't carry them directly.
_NKD_FILE = re.compile(r"n=(\d+)_k=(\d+)_d=(\d+)")

# sinter stat columns aggregated when merging repeated runs of one circuit.
_AGG = {
    "shots": "sum",
    "errors": "sum",
    "discards": "sum",
    "seconds": "sum",
    "low_confidence_shots": "sum",
    "decoder": "first",
    "json_metadata": "first",
}


def _parse_custom_counts(s) -> dict:
    """Parse a `custom_counts` CSV cell into a plain dict.

    Handles both formats decode.py emits: sinter's native JSON
    (``{"low_confidence_shots":3}``, from ``--resume``) and pandas' stringified
    ``collections.Counter(...)`` repr (from the default append path).
    """
    if not isinstance(s, str) or not s.strip() or s == "Counter()":
        return {}
    # sinter native (resume) format: JSON object.
    try:
        return {k: int(v) for k, v in json.loads(s).items()}
    except (json.JSONDecodeError, AttributeError, ValueError):
        pass
    # pandas default (append) format: "Counter({'k': v, ...})".
    m = re.search(r"Counter\(\{(.+)\}\)", s)
    if not m:
        return {}
    return {k: int(v) for k, v in re.findall(r"'(\w+)':\s*(\d+)", m.group(1))}


def _nkd_from_meta(m) -> tuple:
    """Resolve (n, k, d): prefer sidecar metadata, fall back to the filename.

    Stairway sidecar metadata carries n / k / expected_distance directly
    (propagated by make_threshold_circuits_bell.py); the hyperbolic circuits
    only encode them in the input-circuit filename.
    """
    n = m.get("n")
    k = m.get("k")
    d = m.get("expected_distance", m.get("d"))
    if n is not None and k is not None and d is not None:
        return (int(n), int(k), int(d))
    mm = _NKD_FILE.search(m.get("input_circuit", ""))
    if not mm:
        return (None, None, None)
    return tuple(int(x) for x in mm.groups())


def load_collection(csv_glob: str, bias: float) -> tuple[pd.DataFrame, list[str]]:
    """Load a high-rate decode collection, filtered to one Z-bias.

    Merges repeated runs by ``strong_id``, unpacks ``custom_counts`` into a
    ``low_confidence_shots`` column, resolves the ``[[n, k, d]]`` code, and
    keeps only rows whose metadata ``bias`` matches ``bias``.

    Returns ``(df, code_order)`` where ``code_order`` lists the codes sorted by
    ``n`` (used to color curves consistently across panels).
    """
    import glob

    files = sorted(glob.glob(csv_glob))
    if not files:
        raise FileNotFoundError(f"No files matched: {csv_glob}")

    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

    if "custom_counts" in df.columns:
        counts = df["custom_counts"].apply(_parse_custom_counts)
        df["low_confidence_shots"] = counts.apply(
            lambda c: c.get("low_confidence_shots", 0)
        )
    else:
        df["low_confidence_shots"] = 0

    df = df.groupby("strong_id").agg(_AGG).reset_index()

    meta = df["json_metadata"].apply(
        lambda s: ast.literal_eval(s) if isinstance(s, str) else dict(s)
    )
    df["p"] = meta.apply(lambda m: m.get("p"))
    df["p_bell_z"] = meta.apply(lambda m: m.get("p_bell_z"))
    df["bias"] = meta.apply(lambda m: m.get("bias"))

    df[["n", "k", "d"]] = pd.DataFrame(
        meta.apply(_nkd_from_meta).tolist(), index=df.index
    )
    df["code"] = df.apply(lambda r: f"[[{r['n']},{r['k']},{r['d']}]]", axis=1)

    df["pL"] = df["errors"] / df["shots"].clip(lower=1)
    df["pL_err"] = np.sqrt(
        np.clip(df["pL"] * (1 - df["pL"]), 0, None) / df["shots"].clip(lower=1)
    )

    available_biases = sorted(df["bias"].dropna().unique())
    df = df[np.isclose(df["bias"].astype(float), float(bias))].reset_index(drop=True)
    if df.empty:
        raise ValueError(
            f"no rows for bias={bias} in {csv_glob}; "
            f"available biases: {available_biases}"
        )

    code_order = df.drop_duplicates("code").sort_values("n")["code"].tolist()
    return df, code_order


def df_to_tasks(df: pd.DataFrame) -> list:
    """Rebuild `sinter.TaskStats` from merged rows, injecting n / k / d.

    `sinter.plot_error_rate` consumes `TaskStats`, so we round-trip the merged
    DataFrame back into stats with n / k / expected_distance available in the
    metadata for the group / color / failure-unit callbacks.
    """
    tasks = []
    for _, row in df.iterrows():
        m = (
            ast.literal_eval(row["json_metadata"])
            if isinstance(row["json_metadata"], str)
            else dict(row["json_metadata"])
        )
        m["n"] = int(row["n"])
        m["k"] = int(row["k"])
        m["expected_distance"] = int(row["d"])
        tasks.append(
            sinter.TaskStats(
                strong_id=row["strong_id"],
                decoder=row["decoder"],
                json_metadata=m,
                shots=int(row["shots"]),
                errors=int(row["errors"]),
                discards=int(row["discards"]),
                seconds=float(row["seconds"]),
            )
        )
    return tasks


def _default_failure_units(stat) -> float:
    """Failure units per shot: 1 for pymatching, d/16 for the (tesseract) rest.

    Matches `high_rate.ipynb`'s two-panel convention, where pymatching-decoded
    hyperbolic curves are per-observable and tesseract-decoded stairway curves
    are rescaled by distance.
    """
    if stat.decoder == "pymatching":
        return 1
    return stat.json_metadata["expected_distance"] / 16


PAPER_RCPARAMS = {
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
}


def set_paper_style() -> None:
    """Apply the shared paper rcParams used by the high-rate figures."""
    plt.rcParams.update(PAPER_RCPARAMS)


def plot_panel(
    ax,
    df: pd.DataFrame,
    code_order: list[str],
    *,
    label_suffix: str = "",
    failure_units_per_shot_func=None,
) -> None:
    """Draw one `[[n, k, d]]` family on `ax`.

    Viridis color per code (ordered by ``code_order``), sorted legend, ylim,
    and major+minor grid. The secondary ``p_bell_z`` axis is added separately by
    the caller via :func:`add_pbz_axis` (it depends on layout). Pass
    ``failure_units_per_shot_func`` to override the default per-decoder scaling.
    """
    if failure_units_per_shot_func is None:
        failure_units_per_shot_func = _default_failure_units

    cmap = plt.get_cmap("viridis")
    colors = {
        code: cmap(i / max(1, len(code_order) - 1))
        for i, code in enumerate(code_order)
    }

    def group_func(stat):
        n = stat.json_metadata["n"]
        k = stat.json_metadata["k"]
        d = stat.json_metadata["expected_distance"]
        base = f"[[{n}, {k}, {d}]]"
        return f"{base} {label_suffix}" if label_suffix else base

    def plot_args_func(index, group_key, group_stats):
        # `code_order` keys have no spaces / no suffix — strip both.
        key = group_key.split("]]")[0].replace(" ", "") + "]]"
        return {
            "color": colors.get(key, "k"),
            "marker": "o", "markersize": 4, "linewidth": 1.2,
        }

    sinter.plot_error_rate(
        ax=ax,
        stats=df_to_tasks(df),
        x_func=lambda s: s.json_metadata["p"],
        group_func=group_func,
        plot_args_func=plot_args_func,
        failure_units_per_shot_func=failure_units_per_shot_func,
        failure_values_func=lambda s: s.json_metadata["k"],
    )
    ax.loglog()

    handles, labels = ax.get_legend_handles_labels()

    def _extract_n(lbl):
        m = re.search(r"\[\[(\d+),", lbl)
        return int(m.group(1)) if m else 0

    pairs = sorted(zip(handles, labels), key=lambda hl: _extract_n(hl[1]))
    if pairs:
        h, l = zip(*pairs)
        ax.legend(h, l, loc="lower right",
                  frameon=True, edgecolor="black", fancybox=False)

    ax.set_ylim(1e-5, 4)
    ax.set_xlabel(r"Residual error rate $p_{\text{res}}$")
    ax.set_ylabel(r"Logical error rate per observable")
    ax.grid(which="major", linestyle="-", linewidth=0.5, alpha=0.7)
    ax.grid(which="minor", linestyle=":", linewidth=0.3, alpha=0.5)


def add_pbz_axis(ax, bias: float):
    """Attach a secondary x-axis showing p_bell_z = bias * p below `ax`."""
    ax_pbz = ax.secondary_xaxis(
        -0.20,
        functions=(lambda p: p * bias, lambda pbz: pbz / bias),
    )
    ax_pbz.set_xlabel(
        rf"Bell-pair $Z$ error rate $p_{{bell,z}} = {bias:g}\,p_{{\text{{res}}}}$"
    )
    return ax_pbz
