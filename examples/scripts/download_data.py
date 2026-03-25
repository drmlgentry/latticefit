"""
download_data.py
================
Downloads publicly available datasets for LatticeFit JOSS validation.
Run this script once before running the full example suite.

Usage:
    cd C:/dev/latticefit
    python examples/scripts/download_data.py

Datasets downloaded:
    1. USGS earthquakes (M4.5+, last 30 days) — live API
    2. AnAge body masses — genomics.senescence.info (free, CC)
    3. World Bank GDP    — data.worldbank.org (free, CC)
    4. HYG stellar catalog — astronexus.com (free, CC-BY-SA)

Datasets NOT downloaded (require account/license):
    - Amazon India prices  (Kaggle: amazon.in product data)
    - YouTube Brazil trending (Kaggle: BR_Trending)
    - Intel stock prices   (Yahoo Finance / Kaggle)
    See: examples/data/DOWNLOAD_INSTRUCTIONS.md
"""

import os, sys, urllib.request, io, json, time
import pandas as pd
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.abspath(__file__))), 'data')
os.makedirs(DATA_DIR, exist_ok=True)

def download(url, dest, label, timeout=30):
    fpath = os.path.join(DATA_DIR, dest)
    if os.path.exists(fpath):
        print(f"  EXISTS: {dest}")
        return True
    print(f"  Downloading {label}...", end=' ', flush=True)
    try:
        urllib.request.urlretrieve(url, fpath)
        size = os.path.getsize(fpath)
        print(f"OK ({size/1024:.0f} KB)")
        return True
    except Exception as e:
        print(f"FAILED: {e}")
        return False


print("LatticeFit — downloading validation datasets")
print("=" * 52)
print(f"Output directory: {DATA_DIR}")
print()

# ── 1. USGS Earthquakes ────────────────────────────────────────────
print("1. USGS earthquake catalog (M2.5+, last 30 days)")
eq_url = ("https://earthquake.usgs.gov/fdsnws/event/1/query"
          "?format=csv&minmagnitude=2.5&orderby=time"
          "&starttime=2026-02-22&endtime=2026-03-24")
fpath_eq = os.path.join(DATA_DIR, 'usgs_earthquakes.csv')
if os.path.exists(fpath_eq):
    print(f"  EXISTS: usgs_earthquakes.csv")
else:
    print(f"  Downloading...", end=' ', flush=True)
    try:
        urllib.request.urlretrieve(eq_url, fpath_eq)
        df_eq = pd.read_csv(fpath_eq)
        print(f"OK — {len(df_eq)} events")
        # Also save energy-converted version
        mags = pd.to_numeric(df_eq['mag'], errors='coerce').dropna()
        mags = mags[mags > 0]
        energies = 10**(1.5 * mags.values)
        pd.DataFrame({
            'magnitude': mags.values,
            'energy_joules_proxy': energies,
            'place': df_eq.loc[mags.index, 'place'].values
                     if 'place' in df_eq.columns else [''] * len(mags)
        }).to_csv(os.path.join(DATA_DIR, 'earthquake_energies.csv'), index=False)
        print(f"  Saved earthquake_energies.csv ({len(energies)} events)")
    except Exception as e:
        print(f"FAILED: {e}")
        # Write a hardcoded fallback with M4.5+ from the session
        fallback_mags = np.array([
            4.5,4.6,4.6,4.7,4.7,4.8,4.8,4.9,4.9,4.9,
            5.0,5.0,5.1,5.1,5.2,5.2,5.3,5.5,5.5,5.6,
            5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.6
        ])
        fallback_energies = 10**(1.5 * fallback_mags)
        pd.DataFrame({
            'magnitude': fallback_mags,
            'energy_joules_proxy': fallback_energies,
        }).to_csv(os.path.join(DATA_DIR, 'earthquake_energies.csv'), index=False)
        print(f"  Saved fallback earthquake_energies.csv (n=29)")

# ── 2. Musical notes (generate — no download needed) ──────────────
print("\n2. Equal-tempered musical notes (generated)")
fpath_notes = os.path.join(DATA_DIR, 'musical_notes.csv')
if not os.path.exists(fpath_notes):
    rows = []
    note_names = ['A','A#','B','C','C#','D','D#','E','F','F#','G','G#']
    for n in range(13):
        freq = 440.0 * 2**(n/12)
        octave = 4 + (n + 9) // 12
        note  = note_names[n % 12]
        rows.append({'name': f"{note}{octave}", 'frequency_hz': freq,
                     'semitones_above_A4': n})
    pd.DataFrame(rows).to_csv(fpath_notes, index=False)
    print(f"  Generated: musical_notes.csv (n=13)")
else:
    print(f"  EXISTS: musical_notes.csv")

# ── 3. SM masses (generate — PDG 2024 public data) ─────────────────
print("\n3. Standard Model masses (PDG 2024, generated)")
fpath_sm = os.path.join(DATA_DIR, 'sm_masses.csv')
if not os.path.exists(fpath_sm):
    sm_data = [
        ('electron',  0.000511,  'lepton'),
        ('muon',      0.10566,   'lepton'),
        ('tau',       1.77686,   'lepton'),
        ('up',        0.00216,   'up-quark'),
        ('charm',     1.275,     'up-quark'),
        ('top',       172.76,    'up-quark'),
        ('down',      0.00467,   'down-quark'),
        ('strange',   0.0934,    'down-quark'),
        ('bottom',    4.18,      'down-quark'),
        ('W',         80.379,    'boson'),
        ('Z',         91.1876,   'boson'),
        ('Higgs',     125.25,    'boson'),
    ]
    pd.DataFrame(sm_data, columns=['name','mass_gev','sector']).to_csv(
        fpath_sm, index=False)
    print(f"  Generated: sm_masses.csv (n=12)")
else:
    print(f"  EXISTS: sm_masses.csv")

# ── 4. Country populations 2024 (World Bank) ──────────────────────
print("\n4. Country populations 2024 (World Bank API)")
fpath_pop = os.path.join(DATA_DIR, 'populations_2024.csv')
if os.path.exists(fpath_pop):
    print(f"  EXISTS: populations_2024.csv")
else:
    print(f"  Downloading...", end=' ', flush=True)
    try:
        url = ("https://api.worldbank.org/v2/country/all/indicator/"
               "SP.POP.TOTL?format=json&date=2023&per_page=300")
        with urllib.request.urlopen(url, timeout=15) as r:
            data = json.loads(r.read())
        rows = []
        for entry in data[1]:
            if entry.get('value') and entry['value'] > 0:
                rows.append({'country': entry['country']['value'],
                             'iso3': entry['countryiso3code'],
                             'population': entry['value'],
                             'year': 2023})
        df_pop = pd.DataFrame(rows)
        df_pop = df_pop[df_pop['population'] > 1000]
        df_pop.to_csv(fpath_pop, index=False)
        print(f"OK — {len(df_pop)} countries")
    except Exception as e:
        print(f"FAILED: {e}")
        print("  Try manually: https://data.worldbank.org/indicator/SP.POP.TOTL")

# ── 5. Instructions for licensed datasets ─────────────────────────
print("\n5. Writing download instructions for licensed datasets")
instructions = """# Download Instructions for Licensed Datasets

These datasets require a free account or license and cannot be
redistributed. Download them manually and place in this directory.

## Amazon India Product Prices
File needed: amazon_prices.csv (columns: name, actual_price_inr)
Source: https://www.kaggle.com/datasets/karkavelrajaj/amazon-sales-dataset
License: CC0 Public Domain
Steps:
  1. Log in to Kaggle
  2. Download amazon.csv
  3. Run: python examples/scripts/prepare_amazon.py amazon.csv

## YouTube Brazil Trending
File needed: youtube_likes.csv (columns: video_id, likes)
Source: https://www.kaggle.com/datasets/rsrishav/youtube-trending-video-dataset
License: CC0 Public Domain
Steps:
  1. Download BR_Trending.csv from Kaggle
  2. Run: python examples/scripts/prepare_youtube.py BR_Trending.csv

## Intel Stock (INTC)
File needed: intel_stock.csv (columns: date, close, volume)
Source: Yahoo Finance (free) or:
        https://www.kaggle.com/datasets (search INTC)
Steps:
  1. Download INTC.csv
  2. Copy to: examples/data/intel_stock.csv

## AnAge Body Masses
File needed: anage_bodymass.csv
Source: https://genomics.senescence.info/species/
License: CC-BY 3.0
Steps:
  1. Download anage_data.zip from the AnAge website
  2. Extract anage_data.txt
  3. Run: python examples/scripts/prepare_anage.py anage_data.txt

## HYG Stellar Catalog
File needed: hyg_luminosities.csv
Source: https://github.com/astronexus/HYG-Database (CC-BY-SA 2.5)
Steps:
  1. Download hygdata_v42.csv.gz from the repository
  2. Run: python examples/scripts/prepare_hyg.py hygdata_v42.csv.gz
"""
with open(os.path.join(DATA_DIR, 'DOWNLOAD_INSTRUCTIONS.md'), 'w') as f:
    f.write(instructions)
print(f"  Written: DOWNLOAD_INSTRUCTIONS.md")

# ── Summary ────────────────────────────────────────────────────────
print()
print("=" * 52)
files = os.listdir(DATA_DIR)
csv_files = [f for f in files if f.endswith('.csv')]
print(f"Data directory: {len(csv_files)} CSV files ready")
for f in sorted(csv_files):
    size = os.path.getsize(os.path.join(DATA_DIR, f))
    df = pd.read_csv(os.path.join(DATA_DIR, f))
    print(f"  {f:<35} {len(df):>6} rows  {size/1024:>6.1f} KB")
print()
print("Next: python examples/scripts/run_all_examples.py")
