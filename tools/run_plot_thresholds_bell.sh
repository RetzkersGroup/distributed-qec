#!/bin/bash
set -e

# Fit and plot the Bell threshold data collected by
# tools/run_decode_thresholds_bell.sh (and tools/equivalent/run_decode*.sh).
# Writes fits to out/collection/threshold_bell_fits.csv and figures to
# out/figures/.
python tools/plot_thresholds_bell.py "$@"
