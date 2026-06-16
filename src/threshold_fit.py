"""Estimate the threshold of a code family from sinter-style sweep data.

Two-stage estimator (matches the approach in
~/Projects/Titus/2025/Distributed/Lattice_Surgery/lattice_surgery/merge_plot.py):

1. `locate_threshold` — coarse estimate. For each sampled p, take the variance
   of log(p_L) across distances; the p with the smallest variance is where the
   curves are closest, i.e. the rough crossing point.

2. `fit_threshold` — quadratic finite-size-scaling refinement, restricted to a
   small window of p values around the coarse estimate:

       p_L(p, d) = A + B * x + C * x**2,    x = (p - p_th) * d**mu

   with binomial-sigma-weighted curve_fit and bounded parameters.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit


@dataclass
class ThresholdFit:
    p_th: float
    p_th_err: float
    mu: float
    mu_err: float
    A: float
    B: float
    C: float
    p_init: float
    p_window: tuple[float, ...]
    n_points: int

    def predict(self, p: np.ndarray, d: np.ndarray) -> np.ndarray:
        x = (np.asarray(p) - self.p_th) * np.asarray(d) ** self.mu
        return self.A + self.B * x + self.C * x * x


def locate_threshold(
    p: np.ndarray,
    d: np.ndarray,
    p_L: np.ndarray,
) -> float | None:
    """Return the p minimizing variance of log(p_L) across distances.

    This is a coarse, fit-free estimate of where the distance curves cross —
    use it to pick a window for `fit_threshold`. Returns None if no p has
    samples at >=2 distances.
    """
    p = np.asarray(p, dtype=float)
    d = np.asarray(d, dtype=float)
    p_L = np.asarray(p_L, dtype=float)

    log_pL_by_p: dict[float, list[float]] = {}
    for pi, di, pli in zip(p, d, p_L):
        if not (np.isfinite(pi) and np.isfinite(di) and np.isfinite(pli) and pli > 0):
            continue
        log_pL_by_p.setdefault(float(pi), []).append(np.log(pli))

    best_p, min_var = None, np.inf
    for pi, log_rates in sorted(log_pL_by_p.items()):
        if len(log_rates) < 2:
            continue
        var = float(np.var(log_rates))
        if var < min_var:
            min_var, best_p = var, pi
    return best_p


def fit_threshold(
    p: np.ndarray,
    d: np.ndarray,
    p_L: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
    p_init: float | None = None,
    n_per_side: int = 2,
) -> ThresholdFit:
    """Fit the quadratic FSS ansatz around a coarse threshold estimate.

    Parameters
    ----------
    p, d, p_L : arrays of physical error rate, code distance, and logical
        error rate (same length).
    sigma : optional 1-sigma uncertainty on each p_L (e.g. binomial). If
        provided, used for weighted least-squares with absolute_sigma=True.
    p_init : optional coarse threshold estimate. If None, computed via
        `locate_threshold`.
    n_per_side : window half-width in units of unique-p indices around
        p_init. Default 2 → 5 p-values total.
    """
    p = np.asarray(p, dtype=float)
    d = np.asarray(d, dtype=float)
    p_L = np.asarray(p_L, dtype=float)

    mask = np.isfinite(p) & np.isfinite(d) & np.isfinite(p_L) & (p_L > 0) & (p > 0)
    if sigma is not None:
        sigma = np.asarray(sigma, dtype=float)
        mask &= np.isfinite(sigma) & (sigma > 0)
    p, d, p_L = p[mask], d[mask], p_L[mask]
    if sigma is not None:
        sigma = sigma[mask]

    if p_init is None:
        p_init = locate_threshold(p, d, p_L)
        if p_init is None:
            raise ValueError("could not locate threshold (need >=2 distances per p)")

    # Window: n_per_side neighbors on each side of p_init in the unique-p list.
    all_p = np.array(sorted(set(p.tolist())))
    idx = int(np.argmin(np.abs(all_p - p_init)))
    lo, hi = max(0, idx - n_per_side), min(len(all_p), idx + n_per_side + 1)
    window = set(all_p[lo:hi].tolist())
    in_window = np.array([float(pi) in window for pi in p])
    p, d, p_L = p[in_window], d[in_window], p_L[in_window]
    if sigma is not None:
        sigma = sigma[in_window]

    if len(p) < 5:
        raise ValueError(
            f"need at least 5 points in the fit window; got {len(p)} "
            f"(window={sorted(window)})"
        )
    if len(np.unique(d)) < 2:
        raise ValueError("need at least 2 distinct distances in the window")

    A0 = float(np.mean(p_L[np.isclose(p, all_p[idx])])) if (
        np.any(np.isclose(p, all_p[idx]))
    ) else float(np.median(p_L))

    def model(pd, A, B, C, p_th, mu):
        pp, dd = pd
        x = (pp - p_th) * dd ** mu
        return A + B * x + C * x * x

    p0 = [A0, 1.0, 0.0, p_init, 1.0]
    bounds = (
        [0.0, -np.inf, -np.inf, max(0.0, p_init * 0.5), 0.1],
        [1.0, np.inf, np.inf, 1.0, 5.0],
    )

    popt, pcov = curve_fit(
        model,
        (p, d),
        p_L,
        p0=p0,
        sigma=sigma,
        absolute_sigma=sigma is not None,
        bounds=bounds,
        maxfev=50_000,
    )
    perr = np.sqrt(np.diag(pcov))
    A, B, C, p_th, mu = popt
    return ThresholdFit(
        p_th=float(p_th),
        p_th_err=float(perr[3]),
        mu=float(mu),
        mu_err=float(perr[4]),
        A=float(A),
        B=float(B),
        C=float(C),
        p_init=float(p_init),
        p_window=tuple(sorted(window)),
        n_points=int(len(p)),
    )
