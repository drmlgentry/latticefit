"""
latticefit.optimize
===================
Automatic lattice discovery: search over anchor, base, and denom.
"""

from __future__ import annotations
import numpy as np
from itertools import product as iproduct
from .core import fit, LatticeFitResult
from typing import Sequence


PHI   = (1 + 5 ** 0.5) / 2
KNOWN = {"phi": PHI, "e": np.e, "2": 2.0, "sqrt2": 2**0.5, "10": 10.0}


def discover(
    data:       Sequence[float],
    bases:      Sequence[float] | None = None,
    denoms:     Sequence[int]   | None = None,
    anchors:    Sequence[float] | None = None,
    names:      Sequence[str]   | None = None,
    top_n:      int = 5,
) -> list[LatticeFitResult]:
    """
    Search for the best-fitting lattice over a grid of parameters.

    Parameters
    ----------
    bases   : candidate bases r (default: phi, e, 2, sqrt2, 10)
    denoms  : candidate denominators d (default: 1,2,3,4,6)
    anchors : candidate anchors A (default: each data point + geometric mean)
    top_n   : return the top_n results sorted by RMS residual

    Returns
    -------
    list of LatticeFitResult sorted by ascending RMS
    """
    x = np.asarray(data, dtype=float)

    if bases   is None: bases   = list(KNOWN.values())
    if denoms  is None: denoms  = [1, 2, 3, 4, 6]
    if anchors is None:
        # Use each data point and the geometric mean as candidate anchors
        gm = float(np.exp(np.mean(np.log(x))))
        anchors = list(x) + [gm]

    results = []
    for r, d, A in iproduct(bases, denoms, anchors):
        try:
            res = fit(x, anchor=A, base=r, denom=d, names=names)
            results.append(res)
        except Exception:
            pass

    results.sort(key=lambda r: r.rms)
    return results[:top_n]
