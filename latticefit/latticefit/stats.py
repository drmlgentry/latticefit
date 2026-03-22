"""
latticefit.stats
================
Null hypothesis tests for lattice alignment.
"""

from __future__ import annotations
import numpy as np
from .core import fit, LatticeFitResult
from typing import Sequence
from dataclasses import dataclass


@dataclass
class NullTestResult:
    observed_rms:  float
    null_mean:     float
    null_std:      float
    p_value:       float
    n_trials:      int
    test_name:     str

    def summary(self) -> str:
        z = (self.null_mean - self.observed_rms) / max(self.null_std, 1e-12)
        return (
            f"Null test: {self.test_name}\n"
            f"  Observed RMS = {self.observed_rms:.5f}\n"
            f"  Null mean    = {self.null_mean:.5f} ± {self.null_std:.5f}\n"
            f"  Z-score      = {z:.2f}\n"
            f"  p-value      = {self.p_value:.4f}  (n={self.n_trials})"
        )


def log_uniform_null(
    result:   LatticeFitResult,
    n_trials: int = 10_000,
    seed:     int = 42,
) -> NullTestResult:
    """
    Null: n masses drawn uniformly in [log x_min, log x_max].
    Tests whether observed RMS could arise by chance from any
    log-uniform distribution over the same range.
    """
    rng = np.random.default_rng(seed)
    x = result.data
    lo, hi = np.log(x.min()), np.log(x.max())
    n = len(x)

    null_rms = np.empty(n_trials)
    for i in range(n_trials):
        sample = np.exp(rng.uniform(lo, hi, n))
        null_rms[i] = fit(
            sample, anchor=result.anchor,
            base=result.base, denom=result.denom
        ).rms

    p = float(np.mean(null_rms <= result.rms))
    return NullTestResult(
        observed_rms = result.rms,
        null_mean    = float(null_rms.mean()),
        null_std     = float(null_rms.std()),
        p_value      = p,
        n_trials     = n_trials,
        test_name    = "Log-uniform",
    )


def sector_anchor_null(
    result:      LatticeFitResult,
    sector_ids:  Sequence[Sequence[int]],
    n_trials:    int = 10_000,
    seed:        int = 42,
) -> NullTestResult:
    """
    Structure-preserving null: keep within-sector ratios fixed,
    randomise sector anchor positions.

    sector_ids : list of index groups, e.g. [[0,1,2],[3,4,5],[6,7,8]]
                 for three sectors of three particles each.
    """
    rng = np.random.default_rng(seed)
    x = result.data
    log_x = np.log(x)
    lo, hi = log_x.min(), log_x.max()

    # Within-sector log-differences (relative to sector minimum)
    sectors = []
    for idx in sector_ids:
        idx = list(idx)
        lx = log_x[idx]
        sectors.append(lx - lx.min())  # offsets from anchor

    null_rms = np.empty(n_trials)
    for i in range(n_trials):
        sample = np.empty(len(x))
        for idx, offsets in zip(sector_ids, sectors):
            span = offsets.max()
            anchor = rng.uniform(lo, hi - span)
            sample[list(idx)] = np.exp(anchor + offsets)
        null_rms[i] = fit(
            sample, anchor=result.anchor,
            base=result.base, denom=result.denom
        ).rms

    p = float(np.mean(null_rms <= result.rms))
    return NullTestResult(
        observed_rms = result.rms,
        null_mean    = float(null_rms.mean()),
        null_std     = float(null_rms.std()),
        p_value      = p,
        n_trials     = n_trials,
        test_name    = "Sector-anchor",
    )
