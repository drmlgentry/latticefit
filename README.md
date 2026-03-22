# LatticeFit

**Deterministic engine for discovering discrete multiplicative structure in positive real data.**

Given measurements xâ‚, xâ‚‚, â€¦, xâ‚™, LatticeFit tests whether they cluster near a geometric lattice:

    xáµ¢ â‰ˆ A Â· r^(káµ¢/d),   káµ¢ âˆˆ Z

and quantifies whether the alignment is statistically non-accidental.

## Installation

```bash
pip install latticefit            # core only
pip install "latticefit[plots]"   # with matplotlib
```

## Quick start

```python
import latticefit
import numpy as np

# Standard Model fermion masses (GeV)
masses = [5.11e-4, 0.1057, 1.777, 0.00216, 1.275, 172.76,
          0.00467, 0.0934, 4.18]
names  = ["e", "mu", "tau", "u", "c", "t", "d", "s", "b"]

result = latticefit.fit(masses, anchor=5.11e-4,
                        base=latticefit.PHI, denom=4, names=names)
print(result.summary())

# Statistical validation
null = latticefit.log_uniform_null(result, n_trials=10_000)
print(null.summary())

# Plot
latticefit.plots.plot_fit(result, outfile="fit.png")
```

## Command-line

```bash
latticefit masses.csv --anchor 5.11e-4 --base phi --denom 4 --null 10000 --plot fit.png
latticefit data.csv --auto
```

## Applications

- Particle physics mass spectra
- Financial return distributions
- Biological scaling laws
- Engineering failure-rate hierarchies
- Signal amplitude spectra
- Any dataset spanning multiple orders of magnitude

## Patentable method

The bounded integer scan + structure-preserving null test combination
is a novel software method. See `PATENT_NOTES.md` for provisional
patent guidance.

## Citation

If you use LatticeFit in research, please cite:

> M. L. Gentry, "Geometric Unification of Flavor: Masses, Mixing,
> and CP from the Golden Ratio Lattice," submitted (2026).

## License

MIT

---
> **Patent pending.** US Provisional Application No. 64/013,306 (filed March 22, 2026).

