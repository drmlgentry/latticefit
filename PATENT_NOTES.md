# LatticeFit — Provisional Patent Notes
## Systems and Methods for Detecting Discrete Multiplicative Structure in Data

**Inventor:** Marvin L. Gentry  
**Contact:** drmlgentry@protonmail.com  
**ORCID:** 0009-0006-4550-2663  
**Date of conception:** March 2026

---

## Title

**Systems and Methods for Fitting Positive Real-Valued Data to a
Parameterized Geometric Lattice with Statistical Validation**

---

## Field of the Invention

The invention relates to computational methods for data analysis,
specifically to the detection and quantification of discrete
multiplicative structure in sets of positive real-valued measurements.

---

## Background

Many physical, biological, financial, and engineering systems produce
measurements that span multiple orders of magnitude. Existing methods
for analyzing such data include log-linear regression, power-law
fitting, and spectral analysis. These methods detect continuous
parametric relationships but do not test whether data clusters near
a *discrete* set of the form {A · r^(k/d) : k ∈ Z}.

No prior art exists for the specific combination of:
1. A bounded integer scan over lattice labels
2. An anchor-selection protocol with invariance guarantees
3. Structure-preserving null hypothesis tests

---

## Summary of the Invention

The invention provides a computer-implemented method for determining
whether a set of positive real measurements {x₁, …, xₙ} clusters
near a geometric lattice parameterized by (A, r, d):

    x_i ≈ A · r^(k_i / d),   k_i ∈ Z

The method produces:
- Integer lattice labels k_i for each measurement
- Residuals δ_i = |log_r(x_i/A) − k_i/d|
- RMS residual as a global fit quality metric
- p-values from two novel null hypothesis tests

---

## Detailed Description of Preferred Embodiments

### Step 1 — Log Encoding

Given input data {x_i} and parameters (A, r, d), compute:

    y_i = log(x_i / A) / log(r)

This maps the multiplicative structure to an additive one.

### Step 2 — Bounded Integer Scan

For each observation, determine the nearest lattice label:

    k_i = round(d · y_i)

constrained to a pre-declared range [k_min, k_max] chosen to
enclose all expected data values. The range is declared before
data inspection.

### Step 3 — Residual Computation

    δ_i = |y_i − k_i / d|

The RMS residual E = sqrt(mean(δ_i²)) summarizes fit quality.
The theoretical maximum for a uniform distribution is 1/(2d).

### Step 4 — Anchor-Shift Invariance Check

Replacing anchor A with any other value A' produces labels k'_i
satisfying:

    k'_i − k'_j = k_i − k_j   for all i, j

This verifies that relative differences are anchor-independent,
a necessary condition for the structure to be physically meaningful.

### Step 5 — Log-Uniform Null Test

Generate N synthetic datasets by drawing n values uniformly in
[log(x_min), log(x_max)]. For each synthetic dataset, compute
the RMS residual using the same (A, r, d). The p-value is the
fraction of synthetic datasets achieving RMS ≤ observed.

### Step 6 — Structure-Preserving Null Test (Novel)

When the data has known grouping structure (sectors S₁, …, Sₘ),
preserve the within-sector log-ratios and randomise only the
sector anchor positions. This is strictly more powerful than
the log-uniform null when sector structure is independently known.

Specifically:
- For each sector Sⱼ, compute within-sector offsets
  Δᵢⱼ = log(xᵢ) − min_k log(xₖ)  for i ∈ Sⱼ
- In each null trial, draw a random anchor aⱼ ∈ [log x_min,
  log x_max − max(Δᵢⱼ)] for each sector
- Reconstruct null data as exp(aⱼ + Δᵢⱼ)
- Compute null RMS and compare to observed

### Step 7 — Automatic Discovery (Optional)

When (A, r, d) are unknown, the method performs a grid search
over candidate values:
- Anchors: each data point, geometric mean
- Bases: {φ, e, 2, √2, 10, …}
- Denominators: {1, 2, 3, 4, 6, …}

The search returns the top-N parameter sets by ascending RMS.

---

## Claims (Draft for Provisional)

**Claim 1.** A computer-implemented method comprising:
(a) receiving a plurality of positive real-valued measurements;
(b) computing log-encoded coordinates with respect to a declared
    anchor and base;
(c) assigning integer lattice labels by rounding within a
    pre-declared bounded range;
(d) computing residuals between encoded coordinates and assigned labels;
(e) outputting a fit quality metric.

**Claim 2.** The method of claim 1, further comprising verifying
anchor-shift invariance of relative label differences.

**Claim 3.** The method of claim 1, further comprising a log-uniform
null hypothesis test producing a p-value.

**Claim 4.** The method of claim 1, further comprising a
structure-preserving null test that randomises sector anchor positions
while preserving within-sector ratios.

**Claim 5.** The method of claim 1, further comprising an automatic
discovery step searching over a grid of (anchor, base, denominator)
parameters and returning the top-N results by fit quality.

**Claim 6.** A system comprising a processor configured to execute
the method of any preceding claim.

**Claim 7.** A non-transitory computer-readable medium storing
instructions that, when executed, perform the method of any of
claims 1-5.

---

## Novelty Arguments

**vs. log-linear regression:** Linear regression detects continuous
trends; this invention detects discrete structure with no continuous
parameters beyond (A, r, d).

**vs. Fourier analysis:** Fourier methods detect periodic structure
in the linear domain; this invention detects periodic structure in
the log domain with integer-constrained labels.

**vs. power-law fitting:** Power-law fitting optimizes continuous
exponents; this invention constrains exponents to be rational numbers
k/d with k ∈ Z.

**Key novel element:** The structure-preserving null test (Claim 4)
has no prior art. It is the first null test that conditions on
empirically known grouping structure to isolate the specific
sub-hypothesis being tested.

---

## Filing Instructions

1. File a Provisional Patent Application (PPA) with the USPTO.
   Form: SB/16 (available at USPTO.gov)
   Fee: $75-320 depending on entity size (micro entity = $75)
   
2. Include this document + the source code listing from `core.py`
   and `stats.py` as the technical disclosure.

3. You have 12 months from filing date to convert to a non-provisional
   application. Use this time to:
   - Publish the companion academic paper
   - Identify potential licensees
   - Engage a patent attorney if warranted

4. The provisional establishes priority date immediately upon filing.
   You may say "Patent Pending" after filing.

---

## Prior Art Search Notes (Preliminary)

- "Geometric series fitting" — finds continuous parameters, not discrete
- "Lattice-based cryptography" — different domain (security)  
- "Log-periodic oscillations" — continuous oscillatory fits, not discrete labels
- "Froggatt-Nielsen mechanism" — physics model, not a data analysis method
- "Discrete scale invariance" — related concept but no bounded-scan implementation

None of the above disclose the specific combination of bounded integer
scan + structure-preserving null test + anchor-shift invariance verification.
