"""
Standard Model fermion masses — reproduces Table 1 of Gentry (2026).
"""
import sys
sys.path.insert(0, "..")
import latticefit

PHI = latticefit.PHI

masses = [5.10999e-4, 0.105658, 1.77686,   # leptons
          0.00216,    1.275,    172.76,     # up-type quarks
          0.00467,    0.0934,   4.18,       # down-type quarks
          80.379,     91.1876,  125.25]     # EW bosons

names = ["e", "mu", "tau", "u", "c", "t", "d", "s", "b",
         "W", "Z", "H"]

result = latticefit.fit(masses, anchor=5.10999e-4,
                        base=PHI, denom=4, names=names)
print(result.summary())

# Null tests
print()
null_lu = latticefit.log_uniform_null(result, n_trials=50_000)
print(null_lu.summary())

# Sector-anchor null (leptons=0-2, up=3-5, down=6-8)
null_sa = latticefit.sector_anchor_null(
    result,
    sector_ids=[[0,1,2],[3,4,5],[6,7,8]],
    n_trials=50_000,
)
print()
print(null_sa.summary())

# Plot
latticefit.plots.plot_fit(result, outfile="sm_masses_fit.png",
                          title=r"SM masses on the $\varphi$-lattice")
print("\nPlot saved to sm_masses_fit.png")
