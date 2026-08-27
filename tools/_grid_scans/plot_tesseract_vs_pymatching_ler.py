#!/usr/bin/env python3
"""LER vs p_z: tesseract vs pymatching, r=1 and r=3, sinter error bars.

Combines:
  - Tesseract data from full_check_pz_tesseract{,.csv, ,soca.csv, ,2026-07-28.csv}
    (r=1 = biased_bell_v3_equivalent, r=3 = biased_bell_v3_equivalent_fullrep3).
  - Pymatching r=1 from pz_check.csv (biased_bell_v3_equivalent).
  - Pymatching r=3 from full_check_pz.csv (biased_bell_v3_equivalent_fullrep3).

Uses sinter.plot_error_rate for Wilson-interval error bars (handles 0-error
points with upper confidence bounds).

Colors:
  - r=1 tesseract: reds (shade by d)
  - r=1 pymatching: blues (shade by d)
  - r=3 tesseract: purples (shade by d)
  - r=3 pymatching: grays (shade by d)

Usage:
    python tools/_grid_scans/plot_tesseract_vs_pymatching_ler.py
"""

from __future__ import annotations

import argparse
import glob
import pathlib
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sinter
import stim

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
for _p in (str(_ROOT / "tools" / "_grid_scans"), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import plot_grid as g  # noqa: E402


# ─── Data loading ─────────────────────────────────────────────────────────────

def _load_as_task_stats(csv_globs: list[str]) -> list[sinter.TaskStats]:
    """Load CSV(s) into TaskStats, handling both sinter-format and pandas-append.

    Uses load_grid (which handles both formats and aggregates by strong_id),
    then reconstructs TaskStats objects for sinter's plotting API.
    """
    try:
        df = g.load_grid(csv_globs)
    except FileNotFoundError:
        return []
    # Filter honeycomb only.
    df = df[df["input_circuit"].fillna("").str.contains("honeycomb")]
    stats = []
    for _, row in df.iterrows():
        meta = g._parse_meta(row["json_metadata"])
        stats.append(sinter.TaskStats(
            decoder=row["decoder"],
            json_metadata=meta,
            strong_id=row["strong_id"],
            shots=int(row["shots"]),
            errors=int(row["errors"]),
            discards=int(row.get("discards", 0)),
            seconds=float(row.get("seconds", 0)),
        ))
    return stats


def _with_dated(csv: str) -> list[str]:
    """Expand a resume-mode CSV path to include its dated append siblings."""
    p = pathlib.Path(csv)
    return [csv, str(p.with_name(f"{p.stem},*.csv"))]


# ─── Curve definitions ────────────────────────────────────────────────────────

_TESS_CSV = str(_ROOT / "out/collection/full_check_pz_tesseract.csv")
_PM_R1_CSV = str(_ROOT / "out/collection/pz_check.csv")
_PM_R3_CSV = str(_ROOT / "out/collection/full_check_pz.csv")

_BASE = "biased_bell_v3_equivalent"
_FULLREP3 = "biased_bell_v3_equivalent_fullrep3"

# (label, noise_model, decoder_filter, color_map)
_CURVES = [
    ("r=1 tesseract", _BASE, "tesseract_short", plt.get_cmap("Reds")),
    ("r=1 pymatching", _BASE, "pymatching", plt.get_cmap("Blues")),
    ("r=3 tesseract", _FULLREP3, "tesseract_short", plt.get_cmap("Purples")),
    ("r=3 pymatching", _FULLREP3, "pymatching", plt.get_cmap("Greys")),
]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=pathlib.Path,
                    default=_ROOT / "out" / "figures" / "tesseract_vs_pymatching_ler.pdf")
    args = ap.parse_args()

    # Load all data sources.
    tess_stats = _load_as_task_stats(_with_dated(_TESS_CSV))
    pm_r1_stats = _load_as_task_stats(_with_dated(_PM_R1_CSV))
    pm_r3_stats = _load_as_task_stats(_with_dated(_PM_R3_CSV))

    # Pool: tag each stat with decoder for filtering.
    all_stats = tess_stats + pm_r1_stats + pm_r3_stats

    if not all_stats:
        sys.exit("No data found.")

    distances = sorted({
        s.json_metadata.get("expected_distance")
        for s in all_stats if s.json_metadata.get("expected_distance") is not None
    })

    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)

    legend_handles = []
    for label, noise_model, dec_filter, cmap in _CURVES:
        subset = [
            s for s in all_stats
            if s.json_metadata.get("noise_model") == noise_model
            and s.decoder == dec_filter
        ]
        if not subset:
            continue

        shades = np.linspace(0.4, 0.9, max(len(distances), 1))
        for i, d in enumerate(distances):
            color = cmap(shades[i])
            d_subset = [s for s in subset if s.json_metadata.get("expected_distance") == d]
            if not d_subset:
                continue

            sinter.plot_error_rate(
                ax=ax,
                stats=d_subset,
                x_func=lambda s: s.json_metadata.get("p_bell_z"),
                failure_units_per_shot_func=lambda s: s.json_metadata.get("k", 1),
                group_func=lambda s: (s.json_metadata.get("noise_model"), s.decoder,
                                      s.json_metadata.get("expected_distance")),
                plot_args_func=lambda _idx, _key, _pts, c=color: dict(
                    color=c, marker="o", markersize=4, linewidth=1.4),
            )

        # Legend entry per curve family.
        legend_handles.append(
            Line2D([0], [0], color=cmap(0.7), marker="o", linestyle="-",
                   linewidth=1.4, label=label))

    # Distance legend entries.
    for i, d in enumerate(distances):
        legend_handles.append(
            Line2D([0], [0], color="k", marker="o", linestyle="none",
                   markersize=4 + i * 1.5, label=f"d={d}"))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$p_z$ (Bell prep Z-noise)", fontsize=14)
    ax.set_ylabel("Logical error rate per observable", fontsize=13)
    ax.set_title("Tesseract vs pymatching: full-decoder honeycomb", fontsize=14)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(handles=legend_handles, fontsize=9, framealpha=0.9, loc="lower right")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {args.out}")


if __name__ == "__main__":
    main()
