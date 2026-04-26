import re, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
                null_mean=null.mean(), p=p, z=z,
                log_v=log_v, null_arr=null)

with open("ame2020_mass.txt", encoding="utf-8", errors="ignore") as f:
    lines = f.readlines()

# Parse fixed-width format
# Format: a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,f10.5
# Columns (0-indexed):
#  1-3:   NZ
#  4-7:   N
#  8-12:  Z  
#  13-17: A
#  19-21: El
#  23-26: O
#  29-42: mass excess (keV)
#  43-54: mass excess unc
#  55-67: binding energy/A (keV)
#  68-77: binding energy/A unc

records = []
for line in lines:
    if len(line) < 70:
        continue
    if line[0] not in ' 0 1':
        continue
    try:
        # Extract element info
        z_str = line[9:14].strip()
        a_str = line[14:19].strip()
        el_str = line[20:23].strip()
        # Binding energy/A
        be_str = line[55:68].strip().replace('#','').replace('*','')
        be_unc_str = line[68:78].strip().replace('#','').replace('*','')
        
        if not be_str or not z_str or not a_str:
            continue
            
        z = int(z_str)
        a = int(a_str)
        be = float(be_str)
        
        if be > 100:  # skip unbound nuclei
            records.append({'Z':z, 'A':a, 'el':el_str, 'BE_A':be})
    except:
        continue

import pandas as pd
df = pd.DataFrame(records)
print(f"Parsed {len(df)} nuclei")
print(f"Z range: {df['Z'].min()} - {df['Z'].max()}")
print(f"A range: {df['A'].min()} - {df['A'].max()}")
print(f"BE/A range: {df['BE_A'].min():.2f} - {df['BE_A'].max():.2f} keV")
print(f"Orders of magnitude: {np.log10(df['BE_A'].max()/df['BE_A'].min()):.2f}")

# Check for duplicates / verify parsing
print(f"\nSample (stable nuclei):")
stable = df[df['Z'].isin([1,2,6,8,26,82])]
for _,r in stable.head(10).iterrows():
    print(f"  Z={r['Z']:3d} A={r['A']:3d} {r['el']:3s}  BE/A={r['BE_A']:.4f} keV")

# ── Key question: is the range too narrow for a meaningful test? ──
be = df['BE_A'].values
print(f"\nLog-phi span: {np.log(be.max()/be.min())/np.log(PHI):.2f} phi-units")
print(f"Lattice points in range: {(np.log(be.max()/be.min())/np.log(PHI))/0.25:.1f}")

# ── LatticeFit with proper null ───────────────────────────────────
# Standard null: log-uniform in same range
r = latticefit(be)
print(f"\n=== Standard null ===")
print(f"phi: n={r['n']}, RMS={r['rms']:.4f}, "
      f"null={r['null_mean']:.4f}, p={r['p']:.4f}, z={r['z']:+.2f}σ")

# ── Physics-motivated subsets ─────────────────────────────────────
print(f"\n=== By mass region ===")
print(f"{'Region':<20} {'n':>5} {'BE range':>15} {'z':>8} {'p':>8} {'Sig':>5}")
print("-"*60)
for label, mask in [
    ("Light (A<=40)",    df['A']<=40),
    ("Medium (40<A<=100)", (df['A']>40)&(df['A']<=100)),
    ("Heavy (A>100)",    df['A']>100),
    ("Even-even",        (df['Z']%2==0)&((df['A']-df['Z'])%2==0)),
    ("Odd-A",            df['A']%2==1),
]:
    sub = df[mask]['BE_A'].values
    if len(sub) < 20: continue
    r2 = latticefit(sub)
    sig = "✓" if r2['p'] < 0.05 else ""
    be_range = f"{sub.min():.0f}-{sub.max():.0f}"
    print(f"{label:<20} {r2['n']:>5} {be_range:>15} "
          f"{r2['z']:>+8.2f} {r2['p']:>8.4f} {sig:>5}")

# ── Plot ──────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

ax = axes[0]
ax.scatter(df['A'], df['BE_A'], s=1, alpha=0.3, c=df['Z'], cmap='viridis')
ax.set_xlabel("Mass number A")
ax.set_ylabel("Binding energy/A (keV)")
ax.set_title(f"Nuclear Binding Energies\nAME2020, n={len(df)}")

ax = axes[1]
log_v = r['log_v']
ax.hist(log_v, bins=60, color='steelblue', alpha=0.7, density=True)
lo, hi = log_v.min(), log_v.max()
for lp in np.arange(np.floor(lo/0.25)*0.25,
                     np.ceil(hi/0.25)*0.25+0.25, 0.25):
    ax.axvline(lp, color='red', alpha=0.3, linewidth=1.0)
ax.set_xlabel("log_φ(BE/A)")
ax.set_title("Log-φ transform\nred = lattice points")

ax = axes[2]
ax.hist(r['null_arr'], bins=40, color='gray', alpha=0.7,
        density=True, label='Null')
ax.axvline(r['rms'], color='red', linewidth=2,
           label=f"Obs {r['rms']:.4f}")
ax.axvline(r['null_mean'], color='blue', linewidth=1.5,
           linestyle='--', label=f"Null {r['null_mean']:.4f}")
ax.set_title(f"φ-Lattice Null Test\np={r['p']:.4f}, z={r['z']:+.2f}σ")
ax.legend(fontsize=8)

plt.suptitle("LatticeFit: AME2020 Nuclear Binding Energies",
             fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig("ame2020_latticefit.pdf", bbox_inches='tight', dpi=150)
print("\nSaved ame2020_latticefit.pdf")
