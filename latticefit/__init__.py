"""
LatticeFit
==========
Deterministic engine for discovering discrete multiplicative
structure in positive real data.

    x_i ≈ A · r^(k_i / d),   k_i ∈ Z

Quick start
-----------
>>> import latticefit
>>> result = latticefit.fit(data, anchor=data[0], base=1.618, denom=4)
>>> print(result.summary())
>>> latticefit.plots.plot_fit(result, outfile="fit.png")
"""

from .core     import fit, LatticeFitResult
from .optimize import discover, PHI, KNOWN
from .stats    import log_uniform_null, sector_anchor_null, NullTestResult
from . import plots

__version__ = "0.2.0"
__author__  = "Marvin L. Gentry"
__all__     = [
    "fit", "LatticeFitResult",
    "discover", "PHI", "KNOWN",
    "log_uniform_null", "sector_anchor_null", "NullTestResult",
    "plots",
]

from .lucas import fit_lucas, lucas, is_prime_lucas, LucasResult
