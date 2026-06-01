# LatticeFit

**Deterministic engine for discovering discrete multiplicative structure in positive real data.**

Given measurements x_1, x_2, ..., x_n, LatticeFit tests whether they cluster near a geometric lattice:

```
x_i ~ A * r^(k_i/d),   k_i in Z
```

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
```

## Features

- **Deterministic scan** -- no random initialization, reproducible results
- **Anchor-shift invariance** -- fits relative ratios, not absolute scale
- **Structure-preserving null hypothesis** -- permutation test that preserves
  the multiplicative spread of the data
- **Auto-discovery** -- scans over base r and denominator d automatically
- **CLI** -- `latticefit data.csv --base phi --denom 4`
- **Publication-quality plots** -- ladder diagrams, residual histograms

## Background

LatticeFit was developed as the core fitting engine for the
[Hyperbolic Flavour Geometry](https://drmlgentry.github.io/hyperbolic-flavor-geometry/)
programme, which derives Standard Model flavor parameters from the arithmetic
geometry of compact hyperbolic 3-manifolds.

The golden ratio base `r = phi = (1+sqrt(5))/2` is motivated by the
Lucas trace identity: closed geodesics of length `4m*log(phi)` in
hyperbolic 3-manifolds have integer holonomy traces equal to Lucas numbers.

## Lucas mode

```python
# Scan with golden ratio base (motivated by hyperbolic geometry)
result = latticefit.fit(masses, base=latticefit.PHI, lucas=True)
```

## Citation

If you use LatticeFit in your research, please cite:

```
Gentry, M.L. (2026). LatticeFit v0.3.0.
GitHub: https://github.com/drmlgentry/latticefit
PyPI: https://pypi.org/project/latticefit/
```

## License

MIT
