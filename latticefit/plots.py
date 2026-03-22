"""
latticefit.plots
================
Publication-quality visualisation of lattice fits.
"""

from __future__ import annotations
import numpy as np
from .core import LatticeFitResult
from .stats import NullTestResult


def plot_fit(
    result:   LatticeFitResult,
    outfile:  str | None = None,
    title:    str | None = None,
    show:     bool = False,
):
    """
    Log-scale scatter plot of observed vs predicted values,
    with lattice grid lines.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required for plotting.")

    x    = result.data
    xp   = result.predicted
    k    = result.labels
    A    = result.anchor
    r    = result.base
    d    = result.denom
    names = result.names or [str(i) for i in range(len(x))]

    fig, ax = plt.subplots(figsize=(8, 5))

    # Lattice grid
    k_lo, k_hi = k.min() - 2, k.max() + 2
    for ki in range(k_lo, k_hi + 1):
        yl = A * r ** (ki / d)
        ax.axhline(yl, color="#dddddd", linewidth=0.4, zorder=0)

    idx = np.arange(len(x))
    ax.scatter(idx, x,  color="#1f77b4", s=70, zorder=5,
               label="Observed", edgecolors="black", linewidths=0.5)
    ax.scatter(idx, xp, color="#1f77b4", s=70, zorder=4,
               facecolors="none", edgecolors="#1f77b4",
               linewidths=1.5, label="Lattice prediction")
    for i, (xo, xpr) in enumerate(zip(x, xp)):
        ax.plot([i, i], [min(xo, xpr), max(xo, xpr)],
                color="#1f77b4", linewidth=0.8, alpha=0.5)

    ax.set_yscale("log")
    ax.set_xticks(idx)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel(f"Value (log scale)", fontsize=11)
    ax.set_title(title or
                 fr"Lattice fit  $r={r:.4g}$, $d={d}$, $A={A:.4g}$"
                 f"\nRMS residual = {result.rms:.4f}", fontsize=11)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    if outfile:
        plt.savefig(outfile, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    plt.close()
    return fig


def plot_null(
    null:    NullTestResult,
    outfile: str | None = None,
    show:    bool = False,
):
    """Histogram of null RMS distribution vs observed."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required for plotting.")

    # Regenerate null distribution for plotting
    # (stored as summary stats only — recompute if needed)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.axvline(null.observed_rms, color="red", linewidth=2,
               label=f"Observed RMS = {null.observed_rms:.4f}")
    ax.axvline(null.null_mean, color="gray", linewidth=1.5,
               linestyle="--", label=f"Null mean = {null.null_mean:.4f}")
    ax.set_xlabel("RMS residual", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title(f"Null test: {null.test_name}  (p = {null.p_value:.3f})",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    plt.close()
    return fig
