"""
run_all_examples.py
===================
Run all LatticeFit examples and reproduce the JOSS paper Table 1.

Usage:
    cd C:/dev/latticefit
    python examples/scripts/run_all_examples.py

Skips datasets whose files are not present — run download_data.py first.
"""

import sys, os, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

try:
    from latticefit import fit
    from latticefit.stats import log_uniform_null
except ImportError:
    print("ERROR: latticefit not installed. Run: pip install -e .")
    sys.exit(1)

PHI          = (1 + 5**0.5) / 2
EXAMPLES_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR     = os.path.join(EXAMPLES_DIR, 'data')
results      = []

def run(label, vals, base, denom, anchor=None, n_null=10000):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if len(vals) < 5:
        print(f"  SKIP {label}: too few values")
        return
    if anchor is None:
        anchor = vals.min()
    t0     = time.time()
    result = fit(vals, anchor=anchor, base=base, denom=denom)
    null   = log_uniform_null(result, n_trials=n_null, random_seed=42)
    p      = null.p_value
    sig    = "✓✓" if p<0.01 else ("✓" if p<0.05 else ("~" if p<0.10 else "✗"))
    decades = np.log10(vals.max()/vals.min())
    print(f"  {sig}  {label:<42} n={len(vals):>7,}  "
          f"p={p:.3f}  RMS={result.rms:.4f}  ({time.time()-t0:.0f}s)")
    results.append(dict(dataset=label, n=len(vals),
                        decades=round(decades,1), base=base,
                        denom=denom, rms=round(result.rms,5),
                        p_value=round(p,3)))

def load(fname, col=None, cleaner=None):
    """Load a column from examples/data/fname. Returns None if missing."""
    path = os.path.join(DATA_DIR, fname)
    if not os.path.exists(path):
        # Also check examples/ root (existing files)
        path2 = os.path.join(EXAMPLES_DIR, fname)
        if not os.path.exists(path2):
            return None
        path = path2
    df = pd.read_csv(path, on_bad_lines='skip')
    if col is None:
        col = df.select_dtypes(include=[np.number]).columns[0]
    s = df[col]
    if cleaner:
        s = cleaner(s)
    return pd.to_numeric(s, errors='coerce').dropna().values


print("=" * 72)
print("LatticeFit — JOSS paper Table 1 reproduction")
print("=" * 72)
print()

# ── Built-in (always available) ───────────────────────────────────
print("── Built-in datasets")

# Musical notes
notes = load('musical_notes.csv', 'frequency_hz')
if notes is not None:
    run("Equal-tempered notes A4–A5", notes,
        base=2**(1/12), denom=1, anchor=440.0)

# SM masses
sm = load('sm_masses.csv', 'mass_gev')
if sm is not None:
    run("SM fermion + boson masses", sm,
        base=PHI, denom=4, anchor=0.000511)

# Mammal masses (from existing examples/mammal_masses.csv)
def clean_mass(s):
    return pd.to_numeric(s.astype(str).str.replace(',',''), errors='coerce')

for fname in ['mammal_masses.csv', 'data/anage_bodymass.csv']:
    mm = load(fname)
    if mm is not None:
        run("Mammal body masses", mm[mm>0],
            base=PHI**2, denom=8)
        break

# Market caps (from existing examples/market_caps.csv)
mc = load('market_caps.csv')
if mc is not None:
    run("S&P 500 sector market caps", mc[mc>0],
        base=PHI, denom=4)

# ── Downloaded datasets ───────────────────────────────────────────
print("\n── Downloaded datasets (run download_data.py if missing)")

eq = load('earthquake_energies.csv', 'energy_joules_proxy')
if eq is not None:
    run("USGS earthquake energies M4.5+", eq,
        base=2**0.5, denom=2)
else:
    print("  SKIP: earthquake_energies.csv (run download_data.py)")

pop = load('populations_2024.csv', 'population_2024')
if pop is not None:
    run("Country populations 2024", pop,
        base=PHI, denom=12)
else:
    print("  SKIP: populations_2024.csv (run download_data.py)")

# ── Licensed datasets ─────────────────────────────────────────────
print("\n── Licensed datasets (manual download required)")

def clean_rupee(s):
    return s.astype(str).str.replace('₹','').str.replace(',','').str.strip()

amz = load('amazon_prices.csv', 'actual_price', cleaner=clean_rupee)
if amz is None:
    amz = load('amazon.csv', 'actual_price', cleaner=clean_rupee)
if amz is not None:
    run("Amazon India actual prices (₹)", amz[amz>0],
        base=PHI, denom=12)
else:
    print("  SKIP: amazon_prices.csv (see DOWNLOAD_INSTRUCTIONS.txt)")

yt = load('youtube_brazil.csv', 'likes')
if yt is None:
    yt = load('BR_Trending.csv', 'likes')
if yt is not None:
    run("YouTube Brazil likes", yt[yt>0], base=2**0.5, denom=12)
else:
    print("  SKIP: youtube_brazil.csv (see DOWNLOAD_INSTRUCTIONS.txt)")

intc_vol = load('intel_stock.csv', 'Volume')
if intc_vol is None:
    intc_vol = load('INTC.csv', 'Volume')
if intc_vol is not None:
    run("Intel daily trading volume", intc_vol[intc_vol>0],
        base=PHI, denom=12)
else:
    print("  SKIP: intel_stock.csv (see DOWNLOAD_INSTRUCTIONS.txt)")

fx = load('fx_rates.csv', 'Exchange_Rate')
if fx is not None:
    run("Global FX rates vs USD", fx[fx>0], base=2**0.5, denom=12)
else:
    print("  SKIP: fx_rates.csv (see DOWNLOAD_INSTRUCTIONS.txt)")

# ── Summary ───────────────────────────────────────────────────────
if results:
    print()
    print("=" * 72)
    print("RESULTS TABLE")
    print("=" * 72)
    print(f"  {'Dataset':<42} {'n':>7} {'base':>8} "
          f"{'d':>3} {'RMS':>8} {'p':>7} {'sig':>4}")
    print("  " + "-" * 70)
    for r in sorted(results, key=lambda x: x['p_value']):
        sig = ("✓✓" if r['p_value']<0.01 else
               "✓"  if r['p_value']<0.05 else
               "~"  if r['p_value']<0.10 else "✗")
        print(f"  {r['dataset']:<42} {r['n']:>7,} "
              f"{r['base']:>8.5g} {r['denom']:>3} "
              f"{r['rms']:>8.5f} {r['p_value']:>7.3f} {sig:>4}")
    n_sig = sum(1 for r in results if r['p_value'] < 0.05)
    n_neg = sum(1 for r in results if r['p_value'] > 0.10)
    print(f"\n  Significant (p<0.05): {n_sig}/{len(results)}")
    print(f"  Clean negatives (p>0.10): {n_neg}/{len(results)}")
