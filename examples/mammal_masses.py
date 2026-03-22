"""
Example 1: Mammal Body Masses
==============================
Tests whether mammal body masses cluster near a geometric lattice.
Data spans 8 orders of magnitude (shrew to blue whale).

Source: Calder (1984), Schmidt-Nielsen (1984).
"""
import sys
sys.path.insert(0, "../")
import numpy as np
import latticefit

# Body masses in kg (approximate species averages)
names = ["shrew","mouse","rat","rabbit","cat","dog",
         "human","lion","horse","polar_bear","elephant","blue_whale"]
masses = np.array([0.003, 0.020, 0.200, 2.0, 4.0, 30.0,
                   70.0, 180.0, 500.0, 450.0, 5000.0, 150000.0])

print("=" * 60)
print("Mammal Body Masses — LatticeFit Analysis")
print("=" * 60)

# Fit with phi-lattice anchored at shrew
result = latticefit.fit(masses, anchor=masses[0],
                        base=latticefit.PHI, denom=4, names=names)
print(result.summary())

# Auto-discover best lattice
print("\nAuto-discovery (top 3 lattices):")
top = latticefit.discover(masses, names=names, top_n=3)
for i, r in enumerate(top):
    base_name = {latticefit.PHI: "φ", 2.0: "2",
                 2**0.5: "√2", 10.0: "10"}.get(
                 round(r.base, 4), f"{r.base:.4g}")
    print(f"  #{i+1}: base={base_name} d={r.denom} "
          f"anchor={r.anchor:.3g}  RMS={r.rms:.4f}")

# Null test
print()
null = latticefit.log_uniform_null(result, n_trials=10_000)
print(null.summary())

# Save outputs
latticefit.plots.plot_fit(result, outfile="mammal_masses_fit.png",
    title="Mammal body masses on the φ-lattice")
print("\nPlot saved: mammal_masses_fit.png")
