"""
prepare_benchmark_data.py
=========================
Prepares all data files needed by validate_latticefit.py.

Sources:
  - Freely generated: musical notes, SM masses
  - Live download: USGS earthquakes, HYG stars
  - From uploaded files: all others (place in SOURCE_DIR first)

Run:
    python prepare_benchmark_data.py

Output: C:\dev\latticefit\examples\data\*.csv (name,value format)
"""

import numpy as np
import pandas as pd
from pathlib import Path
import urllib.request, io, os

DATA_DIR   = Path(r"C:\dev\latticefit\examples\data")
SOURCE_DIR = Path(r"C:\Users") / os.environ.get('USERNAME','') / "Downloads"

DATA_DIR.mkdir(parents=True, exist_ok=True)
PHI = (1 + 5**0.5) / 2

generated = []
skipped   = []

def save(fname, names, values):
    path = DATA_DIR / fname
    pd.DataFrame({'name': names, 'value': values}).to_csv(path, index=False)
    print(f"  ✓ {fname}  ({len(values)} rows)")
    generated.append(fname)

# ── 1. Musical notes (generated) ──────────────────────────────────
print("Generating musical notes...")
note_names = ['A4','A#4','B4','C5','C#5','D5','D#5','E5',
              'F5','F#5','G5','G#5','A5']
freqs = [440.0 * 2**(n/12) for n in range(13)]
save('musical_notes.csv', note_names, freqs)

# ── 2. SM masses (generated from PDG 2024) ─────────────────────────
print("Generating SM masses...")
sm_names  = ['e','mu','tau','u','c','t','d','s','b','W','Z','H']
sm_masses = [5.11e-4,0.10566,1.77686,0.00216,1.275,172.76,
             0.00467,0.0934,4.18,80.379,91.1876,125.25]
save('sm_masses.csv', sm_names, sm_masses)

# ── 3. USGS earthquakes (live download) ───────────────────────────
print("Downloading USGS earthquake data...")
try:
    url = ("https://earthquake.usgs.gov/fdsnws/event/1/query"
           "?format=csv&minmagnitude=4.5&orderby=time&limit=500"
           "&starttime=2025-01-01")
    with urllib.request.urlopen(url, timeout=15) as r:
        eq_df = pd.read_csv(io.StringIO(r.read().decode()))
    mags = pd.to_numeric(eq_df['mag'], errors='coerce').dropna().values
    mags = mags[mags > 0]
    energies = 10**(1.5 * mags)
    places = eq_df['place'].fillna('unknown').values[:len(energies)]
    save('earthquakes.csv', places, energies)
except Exception as e:
    print(f"  ✗ USGS download failed: {e}")
    # Fallback: hardcoded sample
    mags_s = np.array([4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,
                       5.5,5.7,5.9,6.1,6.3,6.5])
    save('earthquakes.csv',
         [f'eq_{i}' for i in range(len(mags_s))],
         10**(1.5*mags_s))
    skipped.append('earthquakes (fallback sample)')

# ── 4. From uploaded files (check Downloads folder) ───────────────
print("\nLooking for uploaded source files...")

def from_csv(src_path, fname, col_value, col_name=None,
             transform=None, filter_fn=None, year_col=None, year_val=None):
    """Load a source CSV, extract a column, save to benchmark format."""
    if not Path(src_path).exists():
        print(f"  ✗ {fname}: source not found at {src_path}")
        skipped.append(fname)
        return False
    try:
        df = pd.read_csv(src_path, on_bad_lines='skip', low_memory=False)
        if year_col and year_val:
            df = df[df[year_col].isin([year_val, year_val-1])]
            if len(df) == 0:
                df = pd.read_csv(src_path, on_bad_lines='skip')
        vals = pd.to_numeric(df[col_value], errors='coerce').dropna().values
        if filter_fn:
            vals = vals[filter_fn(vals)]
        if transform:
            vals = transform(vals)
        vals = vals[vals > 0]
        names = (df[col_name].astype(str).values[:len(vals)]
                 if col_name else [f'obs_{i}' for i in range(len(vals))])
        save(fname, names[:len(vals)], vals)
        return True
    except Exception as e:
        print(f"  ✗ {fname}: {e}")
        skipped.append(fname)
        return False

def clean_inr(series):
    return pd.to_numeric(
        series.astype(str).str.replace('₹','').str.replace(',','').str.strip(),
        errors='coerce').values

# AnAge body masses
anage_paths = [
    SOURCE_DIR / 'anage_data.txt',
    DATA_DIR.parent / 'anage_data.txt',
    Path(r'C:\dev\latticefit\examples\data\anage_data.txt'),
]
anage_found = False
for ap in anage_paths:
    if ap.exists():
        try:
            df = pd.read_csv(ap, sep='\t', low_memory=False)
            col = 'Adult weight (g)' if 'Adult weight (g)' in df.columns else 'Body mass (g)'
            vals = pd.to_numeric(df[col], errors='coerce').dropna().values
            vals = vals[vals > 0]
            names = (df['Common name'].fillna('unknown').astype(str).values
                     if 'Common name' in df.columns
                     else [f'sp_{i}' for i in range(len(vals))])
            # Mammals only
            if 'Class' in df.columns:
                mask = df['Class'].isin(['Mammalia'])
                vals_m = pd.to_numeric(df.loc[mask, col], errors='coerce').dropna().values
                vals_m = vals_m[vals_m > 0]
                if len(vals_m) > 50:
                    vals = vals_m
                    names_m = df.loc[mask,'Common name'].fillna('unknown').astype(str).values
                    names = names_m[:len(vals)]
            save('mammal_masses.csv', names[:len(vals)], vals)
            anage_found = True
            break
        except Exception as e:
            print(f"  ✗ AnAge from {ap}: {e}")

if not anage_found:
    print("  ✗ mammal_masses.csv: place anage_data.txt in examples/data/")
    skipped.append('mammal_masses')

# Fuel/crude oil
from_csv(SOURCE_DIR / 'fuel_prices_1970_2026.csv',
         'crude_oil.csv', 'Crude_Oil_Price', 'Date')

# Intel stock
from_csv(SOURCE_DIR / 'INTC.csv',
         'intel_close.csv', 'Close', 'Date')
from_csv(SOURCE_DIR / 'INTC.csv',
         'intel_volume.csv', 'Volume', 'Date',
         filter_fn=lambda v: v > 0)

# Amazon India
amz_path = SOURCE_DIR / 'amazon.csv'
if not amz_path.exists():
    amz_path = Path(r'C:\Users') / os.environ.get('USERNAME','') / 'Downloads' / 'archive' / 'amazon.csv'
if amz_path.exists():
    try:
        df = pd.read_csv(amz_path, on_bad_lines='skip')
        dp = clean_inr(df['discounted_price'])
        ap = clean_inr(df['actual_price'])
        names = df['product_name'].astype(str).values
        dp_v = dp[np.isfinite(dp) & (dp > 0)]
        ap_v = ap[np.isfinite(ap) & (ap > 0)]
        save('amazon_discounted.csv', names[:len(dp_v)], dp_v)
        save('amazon_actual.csv',     names[:len(ap_v)], ap_v)
    except Exception as e:
        print(f"  ✗ amazon: {e}")
        skipped.extend(['amazon_discounted', 'amazon_actual'])
else:
    print(f"  ✗ amazon_discounted/actual: place amazon.csv in Downloads")
    skipped.extend(['amazon_discounted', 'amazon_actual'])

# YouTube Brazil
br_path = SOURCE_DIR / 'BR_Trending.csv'
if br_path.exists():
    try:
        df = pd.read_csv(br_path, on_bad_lines='skip')
        df_d = df.groupby('video_id').agg(
            likes=('likes','max'), comments=('comments','max')).reset_index()
        lk = df_d['likes'][df_d['likes']>0].values
        cm = df_d['comments'][df_d['comments']>0].values
        vid = df_d['video_id'].values
        save('youtube_likes.csv',    vid[:len(lk)], lk)
        save('youtube_comments.csv', vid[:len(cm)], cm)
    except Exception as e:
        print(f"  ✗ youtube: {e}")
        skipped.extend(['youtube_likes','youtube_comments'])
else:
    print(f"  ✗ youtube: place BR_Trending.csv in Downloads")
    skipped.extend(['youtube_likes','youtube_comments'])

# FX rates
trade_path = SOURCE_DIR / 'Global_Trade_Exchange_Rates_2000_2025_BALANCED.csv'
if trade_path.exists():
    try:
        df = pd.read_csv(trade_path)
        fx = df['Exchange_Rate'][df['Exchange_Rate']>0].values
        exp = df['Exports_Value'][df['Exports_Value']>0].values
        ctry = df['Country'].values
        save('fx_rates.csv', ctry[:len(fx)], fx)
        # 2024 exports
        df24 = df[df['Year'].isin([2024,2023])]
        exp24 = df24['Exports_Value'][df24['Exports_Value']>0].values
        ctry24 = df24['Country'].values[:len(exp24)]
        save('export_values.csv', ctry24, exp24)
    except Exception as e:
        print(f"  ✗ fx/trade: {e}")
        skipped.extend(['fx_rates','export_values'])
else:
    print(f"  ✗ fx_rates: place Global_Trade_Exchange_Rates... in Downloads")
    skipped.extend(['fx_rates','export_values'])

# Populations
pop_path = SOURCE_DIR / 'popolazione-globale-per-paese-1950-2024.csv'
if pop_path.exists():
    try:
        df = pd.read_csv(pop_path)
        yr = 2024 if 2024 in df['year'].values else df['year'].max()
        df24 = df[df['year']==yr]
        pop = df24['population'][df24['population']>0].values
        ctry = df24['country'].values[:len(pop)]
        save('country_populations.csv', ctry, pop)
    except Exception as e:
        print(f"  ✗ populations: {e}")
        skipped.append('country_populations')
else:
    print(f"  ✗ populations: place popolazione-globale... in Downloads")
    skipped.append('country_populations')

# GDP
gdp_path = SOURCE_DIR / 'global_gdp_inflation_2000_2024.csv'
from_csv(gdp_path, 'gdp_growth.csv', 'GDP_Growth_Percent', 'Country',
         filter_fn=lambda v: v > 0)

# Fungal roots
fung_path = SOURCE_DIR / 'fungal_mycorrhizal_tropical_new.csv'
if fung_path.exists():
    try:
        df = pd.read_csv(fung_path)
        col = [c for c in df.columns if 'RootLength' in c and 'Specific' not in c][0]
        vals = pd.to_numeric(
            df[col].astype(str).str.replace(',','.'), errors='coerce').dropna().values
        vals = vals[vals > 0]
        names = df['Gen_sp（属_种）(40属，66种）'].astype(str).values[:len(vals)]
        save('fungal_roots.csv', names, vals)
    except Exception as e:
        print(f"  ✗ fungal: {e}")
        skipped.append('fungal_roots')
else:
    print(f"  ✗ fungal: place fungal_mycorrhizal_tropical_new.csv in Downloads")
    skipped.append('fungal_roots')

# HYG stellar luminosities
hyg_path = SOURCE_DIR / 'hygdata_v42_csv.gz'
if hyg_path.exists():
    try:
        df = pd.read_csv(hyg_path, compression='gzip', low_memory=False)
        lum = pd.to_numeric(df['lum'], errors='coerce').dropna().values
        lum = lum[lum > 0]
        names = df['proper'].fillna('').astype(str).values[:len(lum)]
        # subsample to 10k for speed
        idx = np.random.default_rng(42).choice(len(lum), size=min(10000,len(lum)),
                                                replace=False)
        save('hyg_luminosities.csv', names[idx], lum[idx])
    except Exception as e:
        print(f"  ✗ HYG: {e}")
        skipped.append('hyg_luminosities')
else:
    print(f"  ✗ hyg: place hygdata_v42_csv.gz in Downloads")
    skipped.append('hyg_luminosities')

# E-commerce prices
ecom_path = SOURCE_DIR / 'global_ecommerce_sales.csv'
from_csv(ecom_path, 'ecommerce_prices.csv', 'Unit_Price', 'Product_Name',
         filter_fn=lambda v: v > 0)

# ── Summary ────────────────────────────────────────────────────────
print(f"\n{'='*50}")
print(f"Generated: {len(generated)} files")
if skipped:
    print(f"Skipped:   {len(skipped)} files")
    for s in skipped:
        print(f"  - {s}")
print(f"\nData directory: {DATA_DIR}")
print(f"Now run: python validate_latticefit.py")
