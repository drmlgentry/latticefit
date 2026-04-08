---
title: 'LatticeFit: A Python package for detecting and validating discrete multiplicative structure in empirical data'
tags:
  - Python
  - statistics
  - multiplicative scaling
  - geometric lattice
  - null hypothesis testing
  - model selection
  - data analysis
authors:
  - name: Marvin L. Gentry
    orcid: 0009-0006-4550-2663
    affiliation: 1
affiliations:
  - name: Independent Researcher, Seattle, WA, USA
    index: 1
date: 7 April 2026
bibliography: paper.bib
---

# Summary

Discrete multiplicative structure — where positive-valued observations
cluster near a geometric lattice $x \approx A \cdot r^{k/d}$ for
integer $k$ — arises in diverse empirical domains. Known examples
include the equal-tempered musical scale ($r = 2^{1/12}$, exact by
construction), earthquake energy release following the Gutenberg-Richter
law [@GutenbergRichter1944], and Standard Model fermion masses clustering
near a golden-ratio lattice [@Gentry2026masses]. Despite the breadth of
phenomena exhibiting this structure, no general-purpose statistical tool
has existed for detecting, quantifying, and validating such patterns
across arbitrary datasets.

`LatticeFit` is a Python package that fills this gap. Given a
positive-valued dataset, it fits the best-matching geometric lattice,
computes a root-mean-square (RMS) residual, and assesses statistical
significance via null hypothesis testing. An auto-discovery mode
searches over candidate bases and denominators, returning models ranked
by information criterion. A publication bundle generator produces
LaTeX tables, figures, methods paragraphs, and self-contained
reproduction scripts in a single command. An interactive Streamlit web
application (v0.3.0) provides a validity-checked analysis interface for
non-programmer users.

# Statement of Need

Multiplicative scaling is ubiquitous in natural and social systems.
Power laws, log-normal distributions, and fractal structures all
exhibit multiplicative self-similarity at the continuous level.
*Discrete* multiplicative structure — clustering at specific
rational powers of a fixed base — is less studied but empirically
common. Examples encountered in this work include:

- Social media engagement counts (YouTube likes, $\sqrt{2}$-spaced, $p < 0.003$)
- Retail prices (Amazon India actual prices, $\varphi$-spaced, $p < 0.001$)
- Currency exchange rates ($\sqrt{2}$ and base-2, $p < 0.001$)
- Biological morphology (fungal root lengths, $\varphi$-spaced, $p = 0.040$)
- Demographic data (country populations, $\varphi$-spaced, $p = 0.025$)
- Allele frequency spectra (rice 44K MAF panel, $\varphi$-spaced, $p < 0.001$)
- Cetacean body masses ($\varphi$-spaced, $p = 0.008$)

Existing tools do not address this problem:

- `scipy.stats` [@scipy] provides continuous distribution fitting and
  hypothesis tests but no lattice-fitting capability.
- `lmfit` [@lmfit] performs general nonlinear curve fitting but
  operates in linear (not log-integer) parameter space.
- `powerlaw` [@powerlaw] fits continuous power-law distributions to
  complementary cumulative distribution functions, not discrete
  lattice clustering.
- `statsmodels` [@statsmodels] provides regression and time-series
  tools without multiplicative-structure detection.

`LatticeFit` provides the first unified framework for this class of
problem, with statistical rigour comparable to established tools.

# Algorithm

## Core fitting

Given $n$ positive observations $\{x_i\}$ and parameters
$(A, r, d)$ — anchor, base, and denominator — the algorithm:

1. Assigns each observation its nearest lattice label:
$$k_i = \operatorname{round}\!\left(d \cdot \log_r(x_i / A)\right)$$

2. Computes the logarithmic residual:
$$\delta_i = \left|\log_r(x_i / A) - k_i / d\right|$$

3. Returns the root-mean-square residual:
$$\text{RMS} = \sqrt{\frac{1}{n} \sum_{i=1}^n \delta_i^2}$$

The maximum possible RMS for denominator $d$ is $0.5/d$, occurring
when all observations fall exactly midway between lattice points.
Reporting $\text{RMS} / (0.5/d)$ as a percentage of maximum provides
a scale-free quality measure. The fractional RMS reduction
$({\rm RMS}_{\rm null} - {\rm RMS}_{\rm obs}) / {\rm RMS}_{\rm null}$
serves as an effect size independent of dataset size.

## Statistical validation

Two null tests assess whether observed clustering is consistent with
random chance:

**Log-uniform null.** Generate $N$ datasets of size $n$ by drawing
uniformly in $[\log x_{\min}, \log x_{\max}]$. The $p$-value is the
fraction of random datasets achieving RMS $\leq$ observed RMS.
This tests against the weakest possible null — a featureless
log-uniform distribution over the same range.

**Sector-anchor null.** When data have a natural group structure
(e.g.\ lepton/quark/boson sectors in particle physics), this test
holds within-group ratios fixed while randomising the group anchor
values. It isolates inter-sector alignment with the lattice as the
signal of interest. This test is covered by US provisional patent
application No.\ 64/013,306.

## Model selection

Auto-discovery searches over candidate bases
$r \in \{\varphi, \sqrt{2}, 2, e, \pi, 10, \varphi^2, \ldots\}$
and denominators $d \in \{1, 2, 3, 4, 6, 8, 12\}$, ranking models
by Akaike Information Criterion (AIC) with a complexity penalty for
larger $d$. This enables hypothesis-free exploration while
penalising over-fitting.

## Bootstrap confidence intervals

The `bootstrap_ci()` function resamples data with replacement
($n_{\text{boot}}$ times) and re-fits the lattice at each resample,
returning confidence intervals on RMS, anchor, and base, together
with per-observation label stability — the fraction of bootstrap
runs in which each observation retains its original lattice label.

# Validity Criteria

Systematic application across 32 datasets has identified four
criteria that must be satisfied before a LatticeFit result should
be interpreted as evidence of genuine lattice structure:

**1. Minimum dynamic range ($\geq 3$ orders of magnitude).**
When data span fewer than $\sim 2$ lattice bins in log-space, any
lattice will show apparent clustering regardless of base, because
the log-uniform null distribution becomes non-uniform at narrow
ranges. Concretely, nuclear binding energies per nucleon (AME2020)
span only 0.97 orders of magnitude; 97\% of the 3,553 nuclei fall
in two adjacent $\varphi$-lattice bins near 8~MeV/nucleon, producing
a spuriously large $z = +7.2$. Shell-correction residuals
(binding energy minus Bethe-Weizsäcker prediction) span less than
one order and give $z = -0.5$ (null), correctly identifying the
absence of $\varphi$-structure in nuclear physics.
**Rule:** Datasets with fewer than 3 orders of magnitude should be
excluded from lattice analysis, or analysed only after
physics-motivated normalisation (e.g.\ dividing by a theoretical
prediction).

**2. Data must not be binned or discretised.**
Fixed bin boundaries can align with lattice points by coincidence,
generating spurious signal. A patent-derived IC50 dataset
(CHEMBL3706050) reports only seven distinct values
(1, 5.5, 20, 65, 200, 650, 1000~nM) and gives $z = +9.8$,
entirely attributable to the bin-spacing structure of the assay
protocol. Screen for binned data with the diagnostic
$n_{\rm unique} / n_{\rm total} < 0.3$ before analysis.

**3. Use a physics-motivated null where available.**
The log-uniform null assumes no prior knowledge of the data
distribution. Where a theoretical model predicts the distribution
(Gutenberg-Richter for earthquakes, Bethe-Weizsäcker for nuclear
masses, random walk for financial returns), residuals from that
model provide a more discriminating test.

**4. Report effect size, not only $p$-value.**
At large $n$, even negligible deviations from the null achieve
statistical significance. The fractional RMS reduction provides
a $n$-independent effect size; values below 2\% should be treated
as practically null regardless of $p$.

These criteria are implemented as automatic diagnostics in the
v0.3.0 Streamlit application and are documented in the
`VALIDITY_CRITERIA.md` file in the repository.

# Implementation

`LatticeFit` is implemented in pure Python ($\geq$ 3.10) with core
dependencies `numpy` [@numpy], `pandas` [@pandas], and `scipy`
[@scipy]. Visualisation uses `matplotlib` [@matplotlib] and is
optional. The package is installable via `pip`:

```bash
pip install latticefit
```

Key modules:

- `latticefit.core`: `fit()`, `discover()`, `FitResult`
- `latticefit.stats`: `log_uniform_null()`, `sector_anchor_null()`,
  `NullTestResult`
- `latticefit.bootstrap`: `bootstrap_ci()`, `propagate_uncertainty()`
- `latticefit.models`: `select_model()` with AIC/BIC ranking
- `latticefit.bundle`: `generate_bundle()`, publication output

A Streamlit web application (`latticefit_app.py`) provides an
interactive interface with automatic validity checking, multi-base
comparison, z-score and effect-size display, a cross-domain survey
tab, and an optional AI assistant for result interpretation.

# Usage

## Python API

```python
import numpy as np
from latticefit import fit, discover
from latticefit.stats import log_uniform_null
from latticefit.bootstrap import bootstrap_ci

# Standard Model fermion masses (GeV, PDG 2024)
masses = np.array([5.11e-4, 0.106, 1.777, 0.00216, 1.275,
                   172.76, 0.00467, 0.0934, 4.18])
names  = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']

# Fit to golden-ratio lattice
result = fit(masses, anchor=5.11e-4, base=1.6180339887,
             denom=4, names=names)
print(f"RMS = {result.rms:.4f}")  # 0.0688

# Statistical validation
null = log_uniform_null(result, n_trials=50000)
print(f"p = {null.p_value:.3f}")  # 0.074

# Bootstrap confidence intervals
boot = bootstrap_ci(result, n_bootstrap=2000)
print(boot.summary())

# Auto-discover best lattice
top5 = discover(masses, top_n=5)

# Full publication bundle
from latticefit.bundle import generate_bundle
generate_bundle(result, null, boot, outdir='sm_masses_bundle/')
```

## Command-line interface

```bash
# Fit with specified parameters
latticefit masses.csv --base 1.618 --denom 4 --anchor 5.11e-4

# Auto-discover best lattice
latticefit masses.csv --auto --null 50000

# Publication bundle (table, figure, methods text, reproduce.py)
latticefit masses.csv --bundle output/ --bootstrap 2000

# Built-in demonstrations
latticefit --demo sm_masses        # Standard Model masses
latticefit --demo earthquakes      # USGS M4.5+ weekly catalog
latticefit --demo musical_notes    # Equal temperament (exact recovery)
latticefit --demo mammal_masses    # AnAge body mass database

# Output formats
latticefit data.csv --latex-table  # LaTeX-ready table
latticefit data.csv --json         # Machine-readable output
latticefit data.csv --cite         # BibTeX citation
```

# Validation

`LatticeFit` has been systematically applied to 32 real-world datasets
spanning physics, biology, geophysics, pharmacology, genomics, finance,
social media, and e-commerce. Datasets satisfying all four validity
criteria are listed in Table 1; datasets excluded by validity criteria
are discussed in the Artifact Identification section below.

| Dataset | $n$ | Dec. | Best $r$ | $d$ | RMS | $p$ | Effect |
|---------|-----|------|----------|-----|-----|-----|--------|
| Equal-tempered notes | 13 | 0.3 | $2^{1/12}$ | 1 | $\approx 0$ | $< 0.001$ | 100% |
| YouTube Brazil comments | 15,722 | 5.2 | $\sqrt{2}$ | 12 | 0.0226 | $< 0.001$ | 31% |
| YouTube Brazil likes | 16,121 | 6.3 | $\sqrt{2}$ | 12 | 0.0238 | 0.002 | 28% |
| Rice 44K allele frequencies | 36,901 | 2.5 | $\varphi$ | 4 | 0.0712 | $< 0.001$ | 1.2% |
| COVID-19 variant frequencies | 103,348 | 3.0 | $\varphi$ | 4 | 0.0712 | $< 0.001$ | 1.2%* |
| Earthquake energies (M4.5+) | 115 | 3.1 | $\sqrt{2}$ | 2 | 0.0107 | $< 0.001$ | 67% |
| Global FX rates (all years) | 5,538 | 11.2 | $\sqrt{2}$ | 12 | 0.0230 | $< 0.001$ | 30% |
| Amazon India actual prices | 1,465 | 3.6 | $\varphi$ | 12 | 0.0230 | $< 0.001$ | 30% |
| Amazon India discounted | 1,465 | 3.3 | $2$ | 12 | 0.0218 | $< 0.001$ | 33% |
| Intel daily volume | 6,559 | 2.7 | $\varphi$ | 12 | 0.0235 | $< 0.001$ | 29% |
| Country populations 2024 | 199 | 6.5 | $\varphi$ | 12 | 0.0225 | 0.025 | 31% |
| Export values 2024 | 213 | 5.4 | $\varphi$ | 12 | 0.0225 | 0.019 | 31% |
| Cetacean body masses | 76 | 3.7 | $\varphi$ | 4 | 0.0631 | 0.008 | 12% |
| Fungal root lengths | 66 | 1.1 | $\varphi$ | 4 | 0.0656 | 0.040 | 9% |
| EGFR SAR (CHEMBL1064829) | 32 | 3.4 | $\varphi$ | 4 | 0.0600 | 0.021 | 17% |
| SM fermion masses | 9 | 5.5 | $\varphi$ | 4 | 0.0688 | 0.074 | 4% |
| AnAge mammal body masses | 627 | 6.7 | $e$ | 8 | 0.0353 | 0.083 | 14% |
| S\&P 500 daily returns | 11,345 | 4.6 | $\varphi$ | 4 | 0.0726 | 0.924 | $-$1% |
| NIST ionisation energies | 1,631 | 4.4 | $\varphi$ | 4 | 0.0721 | 0.378 | 0% |
| COD crystal unit cell volumes | 525,224 | 3.3 | $\varphi$ | 4 | 0.0721 | 0.850 | $-$1% |
| HIV RT IC50 (mixed assays) | 10,041 | 9.1 | $\varphi$ | 4 | 0.0733 | 0.9996 | $-$2% |
| HYG stellar luminosities | 119,626 | 14.5 | $\sqrt{2}$ | 12 | 0.0240 | 0.270 | 27% |
| GDP growth rates | 4,307 | 4.3 | $10$ | 12 | 0.0239 | 0.184 | 27% |

Table: Cross-domain validation results. Dec.\ = $\log_{10}(x_{\max}/x_{\min})$;
Effect = $({\rm RMS}_{\rm null} - {\rm RMS}_{\rm obs})/{\rm RMS}_{\rm null}$;
$p$ = log-uniform null test p-value ($N = 5{,}000$ trials).
Rows above the horizontal rule have $p < 0.05$.
*COVID-19 variant frequencies require a compositional null (see text).

Several results merit specific comment:

**Exact recovery.** Equal-tempered musical notes achieve RMS $\approx 0$
by mathematical construction ($r = 2^{1/12}$, $d = 1$). This
validates the algorithm's correctness.

**Known law recovery.** The USGS earthquake energy dataset
automatically recovers the Gutenberg-Richter scaling law
($E \propto 10^{1.5M}$, equivalent to base-$\sqrt{2}$, $d=2$) via
auto-discovery, without prior knowledge of the law. Testing the
$\varphi$-lattice on earthquake energies gives $z = -3.6$ (null),
demonstrating correct identification of the wrong base.

**Genuine $\varphi$-signals in biological data.** The rice 44K
allele frequency panel ($n = 36{,}901$, $z = +4.25$, $p < 0.001$)
and cetacean body masses ($n = 76$, $z = +2.43$, $p = 0.008$)
show significant $\varphi$-lattice structure. The allele frequency
signal is consistent with evolutionary constraints on multiplicative
population-genetic processes; the cetacean signal is consistent
with hydrodynamic constraints on body size scaling.

**Pharmacological SAR.** The EGFR inhibitor series CHEMBL1064829
($n = 32$ continuous IC50 values, $z = +2.12$, $p = 0.021$) shows
genuine $\varphi$-lattice structure after removal of assay ceiling
values. This is the first report of $\varphi$-spacing in a
single-target structure-activity relationship series.

**Correct negatives.** S\&P 500 daily returns ($z = -1.4$, $p = 0.92$),
crystal unit cell volumes ($z = -1.0$, $p = 0.85$), and NIST
ionisation energies ($z = +0.3$, $p = 0.38$) correctly return
null results. Financial returns are governed by a random walk;
crystal volumes by geometric packing constraints incommensurate
with $\varphi$; ionisation energies by quantum shell structure
scaling as $Z^2$. The method correctly identifies the absence
of $\varphi$-structure in all three cases.

**Kleiber conjugacy.** For the AnAge mammal dataset, the golden
ratio $\varphi$ fits body mass while $\sqrt{2}$ fits metabolic rate.
Since $\varphi^{0.71} \approx \sqrt{2}$ and Kleiber's law gives
metabolic rate $\propto$ mass$^{0.71}$, these are not independent
signals — the two best-fit bases are conjugate under the biological
scaling law. `LatticeFit` correctly identifies both and the
relationship between them.

**COVID-19 variant frequencies** ($n = 103{,}348$, $z = +9.38$,
$p < 0.001$) show a very large apparent signal, but variant
frequencies are compositional data (summing to 1 per time point)
and require a Dirichlet null rather than a log-uniform null.
This result is flagged as requiring a physics-motivated null
before interpretation.

## Artifact Identification

The cross-domain survey identified two systematic artifact classes
that produce false positive signals:

**Narrow-range artifacts.** Nuclear binding energies per nucleon
(AME2020, $n = 3{,}553$) span only 0.97 orders of magnitude
(940--8795 keV/nucleon). The $\varphi$-lattice gives $z = +7.2$,
but 97\% of all nuclei fall in just two adjacent lattice bins
(q=18.50: $n = 962$; q=18.75: $n = 2{,}473$) because binding
energies saturate near 8~MeV/nucleon. Shell-correction residuals
(binding energy minus Bethe-Weizsäcker prediction) span less than
one order and give $z = -0.5$ (null), confirming the absence of
$\varphi$-structure in nuclear physics.

**Binned-data artifacts.** Patent IC50 datasets from ChEMBL
(e.g.\ CHEMBL3706050, CHEMBL3707998) report compounds as belonging
to one of seven concentration classes
(1, 5.5, 20, 65, 200, 650, 1000~nM). These bins happen to be
close to $\varphi$-lattice points ($z = +9.8$), but the signal
is entirely attributable to the assay protocol's fixed
concentration series, not to biological structure-activity
relationships. Continuous IC50 datasets from the same target
(CHEMBL1064829, $n = 32$, $z = +2.1$, $p = 0.021$) show genuine
but modest signal.

# Reproducibility

All validation datasets and analysis scripts are available at
\url{https://github.com/drmlgentry/latticefit/examples/}.
Each example includes:

- Raw data (or download instructions for licensed data)
- A `reproduce.py` script generated by `generate_bundle()`
- Expected outputs for verification

Results are deterministic given a fixed random seed (passed via
`--seed` or the `random_seed` parameter).

# Patent Status

US provisional patent application No.\ 64/013,306 (filed
22 March 2026) covers the sector-anchor null test method described
above. The core fitting algorithm, log-uniform null test, and all
other functionality are released under the MIT License without
restriction.

# Related Software

No existing Python package provides equivalent discrete
multiplicative lattice detection. The closest related tools and
their distinctions are:

- `scipy.stats` [@scipy]: continuous distribution fitting; no lattice
  detection.
- `lmfit` [@lmfit]: nonlinear least squares; no log-integer
  parameterisation.
- `powerlaw` [@powerlaw]: continuous power-law distribution fitting;
  tests against continuous null, not discrete lattice null.
- `statsmodels` [@statsmodels]: regression and time series; no
  multiplicative structure detection.

# Acknowledgements

The author thanks the SnapPy development team [@SnapPy] for
hyperbolic geometry tools used in companion research that motivated
this package, and acknowledges the USGS Earthquake Hazards Program,
the AnAge database [@AnAge], and Kaggle contributors for the
publicly available datasets used in validation.

# References
