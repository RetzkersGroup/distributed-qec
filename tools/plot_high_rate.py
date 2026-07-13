#!/usr/bin/env python3
"""Plot high-rate code decode data, saving figures to out/figures/.

Non-interactive counterpart to notebooks/high_rate.ipynb and
notebooks/high_rate_stairway.ipynb. Once the decode CSVs are collected (via
tools/high_rate/run_decode_*.sh), this reproduces both notebook figures in one
command, at a chosen Z-bias `p_bell_z = BIAS * p`:

  - {date}_high_rate_bias{BIAS}.pdf          two panels (hyperbolic + stairway)
  - {date}_high_rate_stairway_bias{BIAS}.pdf single stairway panel

Figures are written with the notebooks' `savefig(..., bbox_inches='tight')`
convention. A missing collection is skipped with a message.

Examples:

    python tools/plot_high_rate.py
    python tools/plot_high_rate.py --bias 10 --figure-format png
"""

from __future__ import annotations

import argparse
import datetime
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")  # headless: no display needed to save figures.

import matplotlib.pyplot as plt

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from high_rate_analysis import (  # noqa: E402
    load_collection,
    plot_panel,
    add_pbz_axis,
    set_paper_style,
)

_HYP_GLOB = "threshold_hyperbolic_floquet,*.csv"
_STAIRWAY_GLOB = "threshold_stairway,*.csv"


def _two_panel(df_hyp, codes_hyp, df_st, codes_st, bias):
    """Reproduce high_rate.ipynb: hyperbolic (a) over stairway (b)."""
    fig, (ax_hyp, ax_st) = plt.subplots(2, 1, figsize=(5, 8))

    plot_panel(ax_hyp, df_hyp, codes_hyp)
    ax_hyp.set_title("Hyperbolic Floquet codes")
    ax_hyp.set_xlim(float(df_hyp["p"].min()) * 0.7, 3e-3)

    plot_panel(ax_st, df_st, codes_st)
    ax_st.set_title("Stairway codes")
    ax_st.set_xlim(float(df_st["p"].min()) * 0.7, float(df_st["p"].max()) * 1.4)

    for ax, tag in [(ax_hyp, "a)"), (ax_st, "b)")]:
        ax.text(-0.16, 1.02, tag, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="bottom", ha="left")

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.55)
    for ax in (ax_hyp, ax_st):
        add_pbz_axis(ax, bias)
    fig.set_dpi(150)
    return fig


def _stairway_panel(df_st, codes_st, bias):
    """Reproduce high_rate_stairway.ipynb: single stairway panel."""
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    # The standalone stairway figure counts a per-observable failure (unit 1),
    # unlike the two-panel figure which rescales the stairway curves by d/16.
    plot_panel(ax, df_st, codes_st, label_suffix="stairway code",
               failure_units_per_shot_func=lambda s: 1)

    ax.set_xlim(float(df_st["p"].min()) * 0.7, float(df_st["p"].max()) * 1.4)
    ax.set_xlabel(r"Physical error rate $p$")

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.24)
    add_pbz_axis(ax, bias)
    fig.set_dpi(150)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bias", type=float, default=10,
                        help="Z-bias p_bell_z / p to plot (default: 10).")
    parser.add_argument(
        "--collection-dir", type=pathlib.Path,
        default=_ROOT / "out" / "collection",
        help="Directory holding the decode CSVs (default: out/collection).",
    )
    parser.add_argument(
        "--out-dir", type=pathlib.Path, default=_ROOT / "out" / "figures",
        help="Directory to write figures into (default: out/figures).",
    )
    parser.add_argument("--figure-format", default="pdf",
                        help="Figure format / extension (default: pdf).")
    args = parser.parse_args()

    set_paper_style()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ext = args.figure_format.lstrip(".")
    date = datetime.date.today().strftime("%Y%m%d")
    coll = args.collection_dir
    bias = args.bias

    def _load(glob_name, label):
        try:
            return load_collection(str(coll / glob_name), bias)
        except (FileNotFoundError, ValueError) as exc:
            print(f"  [skip] {label}: {exc}")
            return None, None

    df_hyp, codes_hyp = _load(_HYP_GLOB, "hyperbolic")
    df_st, codes_st = _load(_STAIRWAY_GLOB, "stairway")

    def _save(fig, name):
        path = args.out_dir / f"{name}.{ext}"
        fig.savefig(path, bbox_inches="tight")
        print(f"  saved {path}")

    print(f"Writing figures to {args.out_dir} (bias={bias:g}):")
    if df_hyp is not None and df_st is not None:
        _save(_two_panel(df_hyp, codes_hyp, df_st, codes_st, bias),
              f"{date}_high_rate_bias{int(bias)}")
    else:
        print("  [skip] two-panel figure: needs both hyperbolic and stairway data.")

    if df_st is not None:
        _save(_stairway_panel(df_st, codes_st, bias),
              f"{date}_high_rate_stairway_bias{int(bias)}")


if __name__ == "__main__":
    main()
