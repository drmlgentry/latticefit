"""
latticefit.core
===============
Deterministic engine for fitting positive data to a geometric lattice

    x_i ≈ A · r^(k_i / d),   k_i ∈ Z

All operations are reproducible and parameter-free given (A, r, d).
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class LatticeFitResult:
    """Container for a single lattice fit."""
    data:       np.ndarray          # input data (positive reals)
    labels:     np.ndarray          # best-fit integer labels k_i
    predicted:  np.ndarray          # A · r^(k/d) for each point
    residuals:  np.ndarray          # |log_r(x/A) - k/d| for each point
    anchor:     float               # A
    base:       float               # r
    denom:      int                 # d
    rms:        float               # RMS residual
    names:      list[str] | None = None  # optional particle/item names

    def summary(self) -> str:
        lines = [
            f"LatticeFit  r={self.base:.6g}  d={self.denom}  A={self.anchor:.6g}",
            f"RMS residual = {self.rms:.5f}  (max possible = {0.5/self.denom:.5f})",
            "",
            f"{'Name':<12} {'x_obs':>14} {'k':>6} {'x_pred':>14} {'δ':>8} {'|Δ|%':>7}",
            "-" * 62,
        ]
        names = self.names or [f"[{i}]" for i in range(len(self.data))]
        for n, xo, k, xp, res in zip(names, self.data,
                                      self.labels, self.predicted,
                                      self.residuals):
            pct = abs(xo - xp) / xo * 100
            lines.append(f"{n:<12} {xo:>14.6g} {k:>6d} {xp:>14.6g} "
                         f"{res:>8.4f} {pct:>7.1f}%")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"<LatticeFitResult rms={self.rms:.4f} n={len(self.data)}>"


def fit(
    data:   Sequence[float],
    anchor: float,
    base:   float,
    denom:  int = 1,
    names:  Sequence[str] | None = None,
    q_min:  int = -500,
    q_max:  int = 500,
) -> LatticeFitResult:
    """
    Fit positive data to the lattice {A · r^(k/d) : k ∈ Z}.

    Parameters
    ----------
    data   : positive real measurements
    anchor : reference value A (fixes origin of lattice)
    base   : multiplicative base r (e.g. phi, 2, e)
    denom  : denominator d (lattice spacing = log(r)/d)
    names  : optional labels for each data point
    q_min  : lower bound on integer scan (safety)
    q_max  : upper bound on integer scan (safety)

    Returns
    -------
    LatticeFitResult
    """
    x = np.asarray(data, dtype=float)
    if np.any(x <= 0):
        raise ValueError("All data values must be strictly positive.")
    if base <= 0 or base == 1:
        raise ValueError("base must be a positive number ≠ 1.")
    if anchor <= 0:
        raise ValueError("anchor must be strictly positive.")
    if denom < 1:
        raise ValueError("denom must be a positive integer.")

    log_r = np.log(base)
    y = np.log(x / anchor) / log_r  # log_r(x/A)

    k = np.round(denom * y).astype(int)
    k = np.clip(k, q_min, q_max)

    y_pred = k / denom
    residuals = np.abs(y - y_pred)
    predicted = anchor * base ** y_pred
    rms = float(np.sqrt(np.mean(residuals ** 2)))

    return LatticeFitResult(
        data=x,
        labels=k,
        predicted=predicted,
        residuals=residuals,
        anchor=anchor,
        base=base,
        denom=denom,
        rms=rms,
        names=list(names) if names is not None else None,
    )
