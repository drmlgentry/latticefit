import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHI = (1+5**0.5)/2
np.random.seed(42)

def latticefit(values, base=PHI, spacing=0.25, n_null=5000):
    values = np.asarray(values, dtype=float)
    values = values[values > 0]
    if len(values) < 10: return None
    log_v = np.log(values)/np.log(base)
    res = np.abs(log_v - np.round(log_v/spacing)*spacing)
    rms = np.sqrt(np.mean(res**2))
    lo, hi = log_v.min(), log_v.max()
    null = np.array([
        np.sqrt(np.mean(np.abs(
            u-np.round(u/spacing)*spacing)**2))
        for u in (np.random.uniform(lo,hi,len(log_v))
                  for _ in range(n_null))
    ])
    p = (null<=rms).mean()
    z = (null.mean()-rms)/null.std()
    return dict(n=len(values), rms=rms,
                null_mean=null.mean(), p=p, z=z)

# ── Dataset 1: Enzyme Km values ───────────────────────────────────
# Based on Bar-Even et al 2011 PNAS
# log10(Km/M) ~ Normal(-3.5, 1.8) across 2000+ enzymes
# Range ~0.001 uM to 100 mM = 8 orders of magnitude
print("=== Simulated: Enzyme Km values (Bar-Even 2011 parameters) ===")
log10_km = np.random.normal(-3.5, 1.8, 2000)
km_M = 10**log10_km  # in Molar
km_uM = km_M * 1e6   # convert to uM
km_uM = km_uM[(km_uM > 0.0001) & (km_uM < 1e8)]
print(f"n={len(km_uM)}, range {km_uM.min():.4f}-{km_uM.max():.0f} uM")
print(f"Orders: {np.log10(km_uM.max()/km_uM.min()):.1f}")

# ── Dataset 2: mRNA half-lives ────────────────────────────────────
# Tani et al 2012: median ~3h, range 0.1-30h, log-normal
# Yang et al 2003 mouse: mean log(t1/2) = 1.2h, sd = 0.8
print("\n=== Simulated: mRNA half-lives (Tani 2012 parameters) ===")
log_hl = np.random.normal(np.log(3.0), 0.9, 5000)
hl = np.exp(log_hl)
hl = hl[(hl > 0.05) & (hl < 200)]
print(f"n={len(hl)}, range {hl.min():.2f}-{hl.max():.1f} hours")
print(f"Orders: {np.log10(hl.max()/hl.min()):.1f}")

# ── Dataset 3: Metabolite plasma concentrations ───────────────────
# HMDB normal serum: range ~1nM to 10mM = 7 orders
# log-normal distribution
print("\n=== Simulated: Plasma metabolite concentrations (HMDB) ===")
log10_conc = np.random.normal(-4.0, 2.0, 1000)
conc_M = 10**log10_conc
conc_uM = conc_M * 1e6
conc_uM = conc_uM[(conc_uM > 0.0001) & (conc_uM < 1e7)]
print(f"n={len(conc_uM)}, range {conc_uM.min():.4f}-{conc_uM.max():.0f} uM")
print(f"Orders: {np.log10(conc_uM.max()/conc_uM.min()):.1f}")

# ── LatticeFit all three ──────────────────────────────────────────
print(f"\n{'Dataset':<35} {'n':>6} {'RMS':>8} {'Null':>8} "
      f"{'p':>8} {'z':>8} {'Sig':>5}")
print("-"*75)

for label, vals in [
    ("Enzyme Km (Bar-Even params)", km_uM),
    ("mRNA half-lives (Tani params)", hl),
    ("Plasma metabolites (HMDB params)", conc_uM),
]:
    r = latticefit(vals)
    sig = "✓" if r['p'] < 0.05 else ""
    print(f"{label:<35} {r['n']:>6} {r['rms']:>8.4f} "
          f"{r['null_mean']:>8.4f} {r['p']:>8.4f} "
          f"{r['z']:>+8.2f} {sig:>5}")

print("\nNOTE: These are log-normal simulations using published")
print("distributional parameters, not measured values.")
print("Results show expected behavior for log-normal data.")
print("Significant results would require real measured values.")
