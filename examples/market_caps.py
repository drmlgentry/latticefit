"""
Example 3: S&P 500 Sector Market Capitalisations
=================================================
Tests whether US equity sector market caps (approximate 2024 values)
exhibit discrete multiplicative structure.

This is a real-world financial dataset spanning about 1.5 orders of
magnitude. The auto-discovery mode is used since the correct lattice
is not known a priori.

Note: market caps change daily — this uses approximate 2024 figures
for demonstration. Replace with current data for live analysis.
"""
import sys
sys.path.insert(0, "../")
import numpy as np
import latticefit

# Approximate S&P 500 sector market caps, USD billions, early 2024
sector_names = [
    "Energy", "Materials", "Utilities", "Real_Estate",
    "Consumer_Staples", "Industrials", "Healthcare",
    "Consumer_Disc", "Comm_Services", "Financials", "Technology"
]
market_caps = np.array([
    800, 600, 900, 800, 1300, 2200, 2800, 3500, 4000, 4500, 14000
])

print("=" * 60)
print("S&P 500 Sector Market Caps — LatticeFit Analysis")
print("=" * 60)
print("(Approximate 2024 figures, USD billions)\n")

# Fit with phi-lattice
result_phi = latticefit.fit(market_caps, anchor=market_caps[0],
                             base=latticefit.PHI, denom=4,
                             names=sector_names)
print("φ-lattice fit:")
print(result_phi.summary())

# Auto-discover best lattice
print("\nAuto-discovery (top 3 lattices):")
top = latticefit.discover(market_caps, names=sector_names, top_n=3)
best = top[0]
for i, r in enumerate(top):
    print(f"  #{i+1}: base={r.base:.5g}  d={r.denom}  "
          f"anchor={r.anchor:.4g}  RMS={r.rms:.4f}")

print(f"\nBest fit uses base={best.base:.5g}, d={best.denom}:")
print(best.summary())

# Null tests on best fit
print()
null_lu = latticefit.log_uniform_null(best, n_trials=10_000)
print(null_lu.summary())

# Plot both
latticefit.plots.plot_fit(result_phi,
    outfile="market_caps_phi_fit.png",
    title="S&P 500 sector market caps on the φ-lattice")
latticefit.plots.plot_fit(best,
    outfile="market_caps_best_fit.png",
    title=f"S&P 500 sector market caps — best fit\n"
          f"base={best.base:.4g}, d={best.denom}, "
          f"RMS={best.rms:.4f}")
print("\nPlots saved.")
