"""
latticefit/lucas.py
===================
Lucas mode: fixed-base phi fit for arithmetic hyperbolic 3-manifold data.

Lucas-geodesic bridge theorem (Gentry 2026):
    ell = k * log(phi)  iff  |tr(gamma)| = L_k = phi^k + phi^(-k)

This is exact. Usage:
    from latticefit.lucas import fit_lucas, lucas
    result = fit_lucas([5.11e-4, 0.1057, 1.777, ...])
    print(result.summary())
"""
import math
import numpy as np
from dataclasses import dataclass
from typing import Optional, List

PHI = (1 + math.sqrt(5)) / 2
LOG_PHI = math.log(PHI)

_LUCAS = [2,1,3,4,7,11,18,29,47,76,123,199,322,521,843,
          1364,2207,3571,5778,9349,15127,24476,39603,64079]
PRIME_LUCAS = {2, 3, 7, 11, 29, 47, 199, 521, 2207, 3571}


def lucas(k: float) -> float:
    """L_k = phi^k + phi^(-k). Integer values: 2,1,3,4,7,11,18,29,47,..."""
    return PHI**k + PHI**(-k)


def is_prime_lucas(n: int) -> bool:
    return n in PRIME_LUCAS


@dataclass
class LucasAssignment:
    value: float
    q: float
    residual: float
    predicted: float
    delta_pct: float
    k: float
    lk: float
    is_integer_lucas: bool
    lucas_integer: Optional[int]
    is_prime_lucas_: bool


@dataclass
class LucasResult:
    anchor: float
    q_grid: int
    rms: float
    p_value: float
    null_mean: float
    assignments: List[LucasAssignment]
    lucas_integer_fraction: float
    prime_lucas_hits: List[LucasAssignment]

    def summary(self) -> str:
        lines = [
            "LucasFit  base=phi  (Lucas-geodesic bridge theorem)",
            f"  anchor = {self.anchor:.6g}",
            f"  RMS    = {self.rms:.5f}",
            f"  p      = {self.p_value:.4f}  (null mean={self.null_mean:.5f})",
            f"  Lucas-integer fraction = {self.lucas_integer_fraction:.3f}",
            "",
            f"  {'value':>12}  {'q':>6}  {'L_k':>8}  {'delta%':>8}  notes",
            f"  {'-'*58}",
        ]
        for a in self.assignments:
            lk_s = str(a.lucas_integer) if a.is_integer_lucas else f"{a.lk:.3f}"
            note = ("*** prime Lucas" if a.is_prime_lucas_
                    else "** Lucas int"  if a.is_integer_lucas
                    else "***"           if a.residual < 0.02
                    else "**"            if a.residual < 0.05
                    else "*"             if a.residual < 0.09
                    else "")
            lines.append(
                f"  {a.value:>12.5g}  {a.q:>+6.2f}  {lk_s:>8}  "
                f"{a.delta_pct:>7.3f}%  {note}")
        return "\n".join(lines)


def fit_lucas(data, anchor=None, q_grid=4, n_null=5000, rng_seed=42):
    """
    Fit data to phi-lattice {anchor * phi^(q/q_grid) : q in Z}.
    Base fixed to phi by the Lucas-geodesic bridge theorem.
    """
    vals = [float(x) for x in data if x > 0]
    if len(vals) < 2:
        raise ValueError("Need at least 2 positive values")
    if anchor is None:
        anchor = min(vals)

    def _q_res(x):
        ratio = math.log(x / anchor) / LOG_PHI
        q = round(ratio * q_grid) / q_grid
        return q, abs(ratio - q)

    assignments = []
    residuals = []
    for x in vals:
        q, res = _q_res(x)
        k = q / 2
        lk = lucas(k)
        k_int = abs(k - round(k)) < 1e-9
        k_r = round(k)
        lk_int = _LUCAS[k_r] if (k_int and 0 <= k_r < len(_LUCAS)) else None
        predicted = anchor * PHI**q
        a = LucasAssignment(
            value=x, q=q, residual=res,
            predicted=predicted,
            delta_pct=100*abs(x-predicted)/predicted,
            k=k, lk=lk,
            is_integer_lucas=k_int,
            lucas_integer=lk_int,
            is_prime_lucas_=(lk_int in PRIME_LUCAS) if lk_int else False,
        )
        assignments.append(a)
        residuals.append(res**2)

    rms = math.sqrt(sum(residuals)/len(residuals))

    rng = np.random.default_rng(rng_seed)
    log_lo, log_hi = math.log(min(vals)), math.log(max(vals))
    null_rms = [
        math.sqrt(sum(_q_res(x)[1]**2
                       for x in np.exp(rng.uniform(log_lo,log_hi,len(vals))))
                  /len(vals))
        for _ in range(n_null)
    ]
    p_value = float(np.mean(np.array(null_rms) <= rms))

    lucas_int_frac = sum(1 for a in assignments if a.is_integer_lucas)/len(assignments)

    return LucasResult(
        anchor=anchor, q_grid=q_grid, rms=rms,
        p_value=p_value, null_mean=float(np.mean(null_rms)),
        assignments=sorted(assignments, key=lambda a: a.q),
        lucas_integer_fraction=lucas_int_frac,
        prime_lucas_hits=[a for a in assignments if a.is_prime_lucas_],
    )
