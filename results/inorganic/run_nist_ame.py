import re, numpy as np, pandas as pd

PHI = (1+5**0.5)/2

def latticefit(values, base=PHI, spacing=0.25, n_null=5000):
    values = np.asarray(values, dtype=float)
    values = values[values > 0]
    if len(values) < 10: return None
    log_v = np.log(values)/np.log(base)
    res = np.abs(log_v - np.round(log_v/spacing)*spacing)
    rms = np.sqrt(np.mean(res**2))
    lo, hi = log_v.min(), log_v.max()
    null = np.array([
        np.sqrt(np.mean(np.abs(
            u-np.round(u/spacing)*spacing)**2))
        for u in (np.random.uniform(lo,hi,len(log_v))
                  for _ in range(n_null))
    ])
    p = (null<=rms).mean()
    z = (null.mean()-rms)/null.std()
    return dict(n=len(values), rms=rms,
                null_mean=null.mean(), p=p, z=z)

# ── 1. NIST Ionization Energies ───────────────────────────────────
print("=== NIST Ionization Energies ===")

# File is CSV with ="value" Excel quoting
with open("nist_ionization.html", encoding="utf-8", errors="ignore") as f:
    content = f.read()

# Extract all numbers from ="number" pattern
ie_raw = re.findall(r'=""([\d]+\.[\d]+)""', content)
ie = []
for v in ie_raw:
    try:
        f = float(v)
        if 1.0 <= f <= 100000:
            ie.append(f)
    except: pass

ie = np.array(ie)
print(f"Extracted: {len(ie)} ionization energies")
print(f"Range: {ie.min():.3f} - {ie.max():.1f} eV")
print(f"Orders: {np.log10(ie.max()/ie.min()):.1f}")
print(f"Sample: {ie[:8]}")

# Run latticefit all bases
print(f"\n{'Base':>6} {'RMS':>8} {'Null':>8} {'p':>8} {'z':>8} {'Sig':>5}")
print("-"*45)
for name, base, sp in [('phi',PHI,0.25),('2',2.0,1.0),
                        ('e',np.e,1.0),('10',10.0,1.0)]:
    r = latticefit(ie, base, sp)
    sig = "✓" if r['p'] < 0.05 else ""
    print(f"{name:>6} {r['rms']:>8.4f} {r['null_mean']:>8.4f} "
          f"{r['p']:>8.4f} {r['z']:>+8.2f} {sig:>5}")

# ── 2. AME2020 Binding Energies ───────────────────────────────────
print("\n=== AME2020 Nuclear Binding Energies per Nucleon ===")
with open("ame2020_mass.txt", encoding="utf-8", errors="ignore") as f:
    lines = f.readlines()

# Show a data line to understand format
data_lines = [l for l in lines if not l.startswith('0') 
              and not l.startswith('#')
              and len(l) > 60]
print(f"Data lines: {len(data_lines)}")
print(f"Sample line:\n  '{data_lines[5].rstrip()}'")
print(f"  Length: {len(data_lines[5])}")

# Parse AME2020 fixed-width format
# The format is documented in the file header
# Try different column positions
be_values = []
for line in data_lines:
    if len(line) < 60:
        continue
    # Try columns 55-65 (binding energy/A in keV)
    for col_start in [54, 55, 56]:
        try:
            be_str = line[col_start:col_start+11].strip()
            be_str = re.sub(r'[*#]', '', be_str)
            if be_str:
                be = float(be_str)
                if 100 < be < 9000:
                    be_values.append(be)
                    break
        except:
            pass

be = np.array(be_values)
print(f"\nExtracted: {len(be)} binding energies")
if len(be) > 10:
    print(f"Range: {be.min():.2f} - {be.max():.2f} keV/nucleon")
    print(f"Orders: {np.log10(be.max()/be.min()):.2f}")
    print(f"\n{'Base':>6} {'RMS':>8} {'Null':>8} {'p':>8} {'z':>8} {'Sig':>5}")
    print("-"*45)
    for name, base, sp in [('phi',PHI,0.25),('2',2.0,1.0),
                            ('e',np.e,1.0),('10',10.0,1.0)]:
        r = latticefit(be, base, sp)
        sig = "✓" if r['p'] < 0.05 else ""
        print(f"{name:>6} {r['rms']:>8.4f} {r['null_mean']:>8.4f} "
              f"{r['p']:>8.4f} {r['z']:>+8.2f} {sig:>5}")
