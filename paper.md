---
title: 'LatticeFit: A Python package for detecting discrete multiplicative structure in empirical data'
tags:
  - Python
  - statistics
  - data analysis
  - multiplicative structure
  - geometric lattice
  - null hypothesis testing
authors:
  - name: Marvin L. Gentry
    orcid: 0009-0006-4550-2663
    affiliation: 1
affiliations:
  - name: Independent Researcher, Seattle, WA, USA
    index: 1
date: 24 March 2026
bibliography: paper.bib
---

# Summary

Many empirical datasets across scientific domains exhibit structure
where observations cluster near a discrete geometric lattice of the
form $x \approx A \cdot r^{k/d}$, where $A$ is an anchor value,
$r$ is a base, $d$ is a denominator, and $k$ ranges over integers.
Examples include particle masses clustering near a golden-ratio lattice
[@Gentry2026masses], earthquake energy release following a
Gutenberg-Richter power law [@GutenbergRichter1944], equal-tempered
musical frequencies forming an exact geometric series, and mammalian
body masses spanning orders of magnitude.

`LatticeFit` is a Python package that detects, quantifies, and
statistically validates such discrete multiplicative structure in
arbitrary positive-valued datasets. It provides a rigorous
methodological framework including fixed anchoring, bounded integer
scans, and multiple null tests — allowing researchers to distinguish
genuine lattice clustering from random coincidence.

# Statement of Need

Despite the prevalence of multiplicative scaling in nature, no
general-purpose tool exists for detecting and validating discrete
geometric lattice structure across scientific domains. Researchers
typically either notice such patterns informally (without statistical
validation) or construct bespoke tests for specific applications.

`LatticeFit` fills this gap by providing:

1. A unified fitting framework applicable to any positive-valued dataset
2. Statistical null tests (log-uniform and sector-anchor) to assess significance
3. Bootstrap confidence intervals for parameter uncertainty
4. AIC/BIC model selection across candidate bases and denominators
5. Publication-ready output including LaTeX tables, figures, and methods paragraphs

The package is domain-agnostic and has been validated on datasets
from particle physics, geophysics, acoustics, and finance.

# Implementation

`LatticeFit` is implemented in pure Python with dependencies on
`numpy`, `pandas`, `scipy`, and optionally `matplotlib`. The core
algorithm is:

1. **Label assignment**: for each observation $x_i$, find the integer
   $k_i = \text{round}(\log_r(x_i / A) \cdot d)$ minimising the
   logarithmic residual $\delta_i = |\log_r(x_i/A) - k_i/d|$.

2. **RMS computation**: $\text{RMS} = \sqrt{\frac{1}{n}\sum_i \delta_i^2}$.

3. **Null testing**: compare the observed RMS against the distribution
   of RMS values from datasets drawn from a log-uniform distribution
   over the same range (log-uniform null), or from datasets preserving
   within-group ratios while randomising group anchors (sector-anchor null).

4. **Model selection**: grid search over candidate bases $r \in
   \{\varphi, e, 2, \sqrt{2}, 10, \pi, \ldots\}$ and denominators
   $d \in \{1, 2, 3, 4, 6, 8\}$, ranked by AIC/BIC.

The auto-discovery mode (`discover()`) exhaustively searches the
parameter space and returns the top-$n$ models by RMS, enabling
hypothesis-free exploration of lattice structure.

# Usage

```python
import numpy as np
from latticefit import fit, discover
from latticefit.stats import log_uniform_null

# Standard Model fermion masses (GeV)
masses = np.array([5.11e-4, 0.106, 1.777, 0.00216, 1.275,
                   172.76, 0.00467, 0.0934, 4.18])

# Fit to golden-ratio lattice
result = fit(masses, anchor=5.11e-4, base=1.6180339887, denom=4)
print(f"RMS = {result.rms:.4f}")

# Statistical validation
null = log_uniform_null(result, n_trials=50000)
print(f"p = {null.p_value:.3f}")

# Auto-discover best lattice
top = discover(masses, top_n=5)
```

The command-line interface supports common workflows:

```bash
# Fit and report
latticefit masses.csv --base 1.618 --denom 4

# Auto-discover
latticefit masses.csv --auto

# Publication bundle
latticefit masses.csv --bundle output/ --bootstrap --latex-table

# Built-in demos
latticefit --demo sm_masses
latticefit --demo earthquakes
latticefit --demo musical_notes
```

# Validation

`LatticeFit` has been applied to four cross-domain benchmark datasets:

| Dataset | Base $r$ | $d$ | RMS | $p$-value |
|---------|-----------|-----|-----|-----------|
| SM fermion masses | $\varphi$ | 4 | 0.069 | 0.07 |
| Mammalian body masses | $e$ | 6 | 0.043 | 0.08 |
| Equal-tempered notes | $2^{1/12}$ | 1 | $\approx 0$ | $< 0.001$ |
| Earthquake energies (M4.5+) | $\sqrt{2}$ | 2 | 0.011 | $< 0.001$ |

The musical notes example achieves exact recovery (RMS $\approx 0$)
by construction, validating the algorithm. The earthquake energy
example recovers the Gutenberg-Richter scaling law automatically.

# Related Software

To the authors' knowledge, no existing Python package provides
equivalent functionality. Related tools include:

- `scipy.stats` [@scipy]: provides statistical distributions and
  tests but no lattice-fitting capability.
- `lmfit` [@lmfit]: general curve fitting but not discrete lattice
  structure detection.
- `powerlaw` [@powerlaw]: fits continuous power-law distributions
  but not discrete geometric lattices.

# Patent Status

A US provisional patent application (No. 64/013,306, filed
22 March 2026) covers the sector-anchor null test method, which
is the primary novel statistical contribution of this work.
The core fitting and log-uniform null test are released under
the MIT License.

# Acknowledgements

The author thanks the SnapPy development team [@SnapPy] for
hyperbolic geometry tools used in companion research that motivated
this package.

# References
