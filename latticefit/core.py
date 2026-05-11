"""
latticefit.core
===============
Deterministic engine for fitting positive data to a geometric lattice

    x_i ≈ A · r^(k_i / d),   k_i ∈ Z

All operations are reproducible and parameter-free given (A, r, d).

v0.3.0 additions
----------------
- LatticeFitResult.chi2 and .p_value (weighted by experimental uncertainties)
- weighted_fit() for uncertainty-aware fitting
- pdg2024_masses() built-in Standard Model mass table (PDG 2024)
- neutrino_prediction() lattice prediction for neutrino masses
- gauge_bosons() W/Z/H lattice analysis with EW symmetry caveat
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class LatticeFitResult:
    """Container for a single lattice fit."""
    data:        np.ndarray           # input data (positive reals)
    labels:      np.ndarray           # best-fit integer labels k_i
    predicted:   np.ndarray           # A · r^(k/d) for each point
    residuals:   np.ndarray           # |log_r(x/A) - k/d| for each point
    anchor:      float                # A
    base:        float                # r
    denom:       int                  # d
    rms:         float                # RMS residual (unweighted)
    names:       list[str] | None = None   # optional particle/item names
    uncertainties: np.ndarray | None = None  # experimental uncertainties
    chi2:        float | None = None   # chi-squared (if uncertainties given)
    chi2_dof:    float | None = None   # chi2 per degree of freedom
    p_value:     float | None = None   # p-value (if uncertainties given)

    def summary(self) -> str:
        lines = [
            f"LatticeFit  r={self.base:.6g}  d={self.denom}  "
            f"A={self.anchor:.6g}",
            f"RMS residual = {self.rms:.5f}  "
            f"(max possible = {0.5/self.denom:.5f})",
        ]
        if self.chi2 is not None:
            lines.append(
                f"chi2/dof = {self.chi2_dof:.3f}  "
                f"p-value = {self.p_value:.4f}"
            )
        lines += [
            "",
            f"{'Name':<12} {'x_obs':>14} {'k':>6} {'x_pred':>14} "
            f"{'δ':>8} {'|Δ|%':>7}",
            "-" * 62,
        ]
        names = self.names or [f"[{i}]" for i in range(len(self.data))]
        for n, xo, k, xp, res in zip(
            names, self.data, self.labels, self.predicted, self.residuals
        ):
            pct = abs(xo - xp) / xo * 100
            lines.append(
                f"{n:<12} {xo:>14.6g} {k:>6d} {xp:>14.6g} "
                f"{res:>8.4f} {pct:>7.1f}%"
            )
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
        raise ValueError("base must be a positive number != 1.")
    if anchor <= 0:
        raise ValueError("anchor must be strictly positive.")
    if denom < 1:
        raise ValueError("denom must be a positive integer.")

    log_r = np.log(base)
    y = np.log(x / anchor) / log_r   # log_r(x/A)

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


def weighted_fit(
    data:          Sequence[float],
    uncertainties: Sequence[float],
    anchor:        float,
    base:          float,
    denom:         int = 1,
    names:         Sequence[str] | None = None,
    q_min:         int = -500,
    q_max:         int = 500,
) -> LatticeFitResult:
    """
    Fit with experimental uncertainties, computing chi-squared and p-value.

    Parameters
    ----------
    data          : positive real measurements
    uncertainties : 1-sigma absolute uncertainties (same units as data).
                    Use (upper + lower) / 2 for asymmetric errors.
                    Points with uncertainty=0 are treated as exact.
    anchor        : reference value A
    base          : multiplicative base r
    denom         : denominator d
    names         : optional labels
    q_min, q_max  : integer search bounds

    Returns
    -------
    LatticeFitResult with .chi2, .chi2_dof, .p_value populated
    """
    from scipy.stats import chi2 as chi2_dist

    x = np.asarray(data, dtype=float)
    sigma = np.asarray(uncertainties, dtype=float)

    if x.shape != sigma.shape:
        raise ValueError("data and uncertainties must have the same length.")

    # Basic fit first
    result = fit(data, anchor, base, denom, names, q_min, q_max)

    # Propagate experimental uncertainty to q-space residual space
    # sigma_q = (sigma_m / m) / log(base)   [by error propagation of q = log(m/A)/log(r)]
    log_r = np.log(base)
    sigma_q = np.where(sigma > 0, sigma / x / log_r, np.inf)

    # chi2: sum of (residual_q / sigma_q)^2 for constrained points only
    constrained = sigma > 0
    if constrained.sum() == 0:
        chi2_val = None
        chi2_dof_val = None
        p_val = None
    else:
        pulls = result.residuals[constrained] / sigma_q[constrained]
        chi2_val = float(np.sum(pulls ** 2))
        dof = int(constrained.sum()) - 1   # -1 for the anchor degree of freedom
        chi2_dof_val = chi2_val / max(dof, 1)
        p_val = float(1 - chi2_dist.cdf(chi2_val, df=max(dof, 1)))

    result.uncertainties = sigma
    result.chi2 = chi2_val
    result.chi2_dof = chi2_dof_val
    result.p_value = p_val
    return result


# ---------------------------------------------------------------------------
# Built-in Standard Model data (PDG 2024)
# ---------------------------------------------------------------------------

def pdg2024_masses() -> dict:
    """
    Return PDG 2024 fermion masses in GeV with 1-sigma uncertainties.

    Masses use MS-bar scheme at standard reference scales:
      - light quarks (u, d, s): MS-bar at 2 GeV
      - c, b: MS-bar at m_c, m_b respectively
      - t: from direct reconstruction

    Returns
    -------
    dict with keys:
      'names'         : list of particle names
      'masses'        : np.ndarray of central values (GeV)
      'uncertainties' : np.ndarray of 1-sigma uncertainties (GeV)
                        (arithmetic mean of upper/lower for asymmetric errors)
    """
    data = [
        # (name,   mass_GeV,     unc_GeV)
        ('e',      0.51099895e-3, 0.0),            # exact to this precision
        ('mu',     0.1056583755,  2.3e-9),
        ('tau',    1.77686,       0.00012),
        ('u',      2.16e-3,       (0.49e-3 + 0.26e-3) / 2),
        ('d',      4.67e-3,       (0.48e-3 + 0.17e-3) / 2),
        ('s',      93.4e-3,       (8.6e-3  + 3.4e-3)  / 2),
        ('c',      1.273,         0.046),
        ('b',      4.183,         0.085),
        ('t',      172.57,        0.29),
    ]
    names         = [d[0] for d in data]
    masses        = np.array([d[1] for d in data])
    uncertainties = np.array([d[2] for d in data])
    return {
        'names': names,
        'masses': masses,
        'uncertainties': uncertainties,
    }


def neutrino_prediction(
    base:   float | None = None,
    anchor: float | None = None,
    denom:  int = 4,
) -> dict:
    """
    Compute the golden-ratio lattice prediction for neutrino masses.

    Uses the unique integer triple (q1, q2, q3) satisfying both
    PDG 2024 mass-squared difference constraints simultaneously within 1σ:
      Delta_m21^2 = (7.53 +/- 0.18) x 10^{-5} eV^2
      Delta_m31^2 = (2.51 +/- 0.03) x 10^{-3} eV^2

    Parameters
    ----------
    base   : lattice base (default: golden ratio phi)
    anchor : reference mass in GeV (default: electron mass)
    denom  : lattice denominator (default: 4)

    Returns
    -------
    dict with predicted masses, q-labels, and oscillation parameter check
    """
    phi = (1 + 5**0.5) / 2
    m_e = 0.51099895e-3   # GeV

    base   = base   if base   is not None else phi
    anchor = anchor if anchor is not None else m_e

    # Unique solution (normal ordering)
    q_labels = np.array([-149, -146, -134])
    log_r = np.log(base)
    masses_GeV = anchor * base ** (q_labels / denom)
    masses_eV  = masses_GeV * 1e9          # GeV -> eV
    masses_meV = masses_eV  * 1e3          # eV -> meV

    # Compute oscillation parameters
    dm21_sq = masses_eV[1]**2 - masses_eV[0]**2   # eV^2
    dm31_sq = masses_eV[2]**2 - masses_eV[0]**2   # eV^2
    sum_nu  = float(np.sum(masses_meV))

    # PDG 2024 constraints
    dm21_pdg = 7.53e-5;  dm21_unc = 0.18e-5   # eV^2
    dm31_pdg = 2.51e-3;  dm31_unc = 0.03e-3   # eV^2

    dm21_pull = (dm21_sq - dm21_pdg) / dm21_unc
    dm31_pull = (dm31_sq - dm31_pdg) / dm31_unc

    return {
        'q_labels':   q_labels.tolist(),
        'masses_meV': masses_meV.tolist(),
        'sum_meV':    sum_nu,
        'dm21_sq':    dm21_sq,
        'dm31_sq':    dm31_sq,
        'dm21_pull':  dm21_pull,   # sigma deviation from PDG
        'dm31_pull':  dm31_pull,
        'cosmology_bound_meV': 120.0,  # Planck 2018 upper bound
        'cmbs4_sensitivity_meV': 30.0, # CMB-S4 projected sensitivity
        'ordering': 'normal',
        'base':   base,
        'anchor': anchor,
        'denom':  denom,
    }


def gauge_bosons(
    base:   float | None = None,
    anchor: float | None = None,
    denom:  int = 4,
) -> dict:
    """
    Compute golden-ratio lattice positions for electroweak gauge bosons.

    Note: W and Z masses are derived quantities (not independent Yukawa
    parameters), so individual fits are expected to be poor. The geometric
    mean sqrt(M_W * M_Z) represents the pre-symmetry-breaking EW scale
    and fits much better.

    Returns
    -------
    dict with q values, residuals, and interpretation for W, Z, H,
    and their geometric mean.
    """
    phi = (1 + 5**0.5) / 2
    m_e = 0.51099895e-3   # GeV

    base   = base   if base   is not None else phi
    anchor = anchor if anchor is not None else m_e
    log_r  = np.log(base)

    # PDG 2024
    M_W = 80.3692   # GeV (PDG 2024 combination)
    M_Z = 91.1876   # GeV
    M_H = 125.20    # GeV

    results = {}
    for name, mass in [('W', M_W), ('Z', M_Z), ('H', M_H)]:
        q_exact = denom * np.log(mass / anchor) / log_r
        q_int   = round(q_exact)
        res     = abs(q_exact - q_int) / denom
        results[name] = {
            'mass_GeV': mass,
            'q_exact':  float(q_exact),
            'q_int':    int(q_int),
            'residual': float(res),
            'predicted_GeV': float(anchor * base ** (q_int / denom)),
        }

    # Geometric mean (pre-EW-symmetry-breaking scale)
    M_geom = float(np.sqrt(M_W * M_Z))
    q_geom_exact = denom * np.log(M_geom / anchor) / log_r
    q_geom_int   = round(q_geom_exact)
    res_geom     = abs(q_geom_exact - q_geom_int) / denom
    results['geom_mean_WZ'] = {
        'mass_GeV':     M_geom,
        'q_exact':      float(q_geom_exact),
        'q_int':        int(q_geom_int),
        'residual':     float(res_geom),
        'predicted_GeV': float(anchor * base ** (q_geom_int / denom)),
        'note': ('sqrt(M_W*M_Z) is the pre-EW-symmetry-breaking scale. '
                 'W and Z individual fits are poor because their masses '
                 'are derived quantities (M_W = M_Z * cos(theta_W)), '
                 'not independent Yukawa parameters.'),
    }
    results['Weinberg_angle_deg'] = float(
        np.degrees(np.arccos(M_W / M_Z))
    )
    return results
