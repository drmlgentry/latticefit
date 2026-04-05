import json, numpy as np
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
    return dict(n=len(values), rms=rms, null_mean=null.mean(),
                p=p, z=z, log_v=log_v, null_arr=null)

print("Parsing COD crystal data (streaming)...")
# Stream parse - file is 795MB, load in chunks
vols = []
space_groups = []
n_elements = []

import re
# Use regex to extract vol and sg fields efficiently
# rather than loading full JSON into memory
with open("cod_cells.json", encoding="utf-8",
          errors="ignore") as f:
    content = f.read()

print("Extracting unit cell volumes...")
# Extract vol values
vol_matches = re.findall(r'"vol":"([\d.]+)"', content)
sg_matches  = re.findall(r'"sgNumber":"(\d+)"', content)
nel_matches = re.findall(r'"nel":"(\d+)"', content)

vols = np.array([float(v) for v in vol_matches if v])
sgs  = np.array([int(s) for s in sg_matches if s])
nels = np.array([int(n) for n in nel_matches if n])

print(f"Extracted: {len(vols)} unit cell volumes")
print(f"Vol range: {vols.min():.1f} - {vols.max():.1f} Å³")
print(f"Orders:    {np.log10(vols.max()/vols.min()):.1f}")
print(f"Space groups: {len(np.unique(sgs))} unique")

# Remove outliers (>99.9th percentile - some entries are erroneous)
p999 = np.percentile(vols, 99.9)
p001 = np.percentile(vols, 0.1)
vols_clean = vols[(vols >= p001) & (vols <= p999)]
print(f"\nAfter outlier removal (0.1-99.9%): n={len(vols_clean)}")
print(f"Clean range: {vols_clean.min():.1f} - {vols_clean.max():.1f} Å³")
print(f"Orders: {np.log10(vols_clean.max()/vols_clean.min()):.1f}")

# ── LatticeFit all bases ──────────────────────────────────────────
print(f"\n=== LatticeFit: COD Unit Cell Volumes ===")
print(f"{'Base':>6} {'RMS':>8} {'Null':>8} {'p':>8} {'z':>8} {'Sig':>5}")
print("-"*45)
results = {}
for name, base, sp in [('phi',PHI,0.25),('2',2.0,1.0),
                        ('e',np.e,1.0),('10',10.0,1.0)]:
    r = latticefit(vols_clean, base, sp)
    results[name] = r
    sig = "✓" if r['p'] < 0.05 else ""
    print(f"{name:>6} {r['rms']:>8.4f} {r['null_mean']:>8.4f} "
          f"{r['p']:>8.4f} {r['z']:>+8.2f} {sig:>5}")

# ── By crystal system (space group ranges) ───────────────────────
# Space group ranges: triclinic 1-2, monoclinic 3-15,
# orthorhombic 16-74, tetragonal 75-142,
# trigonal 143-167, hexagonal 168-194, cubic 195-230
if len(sgs) == len(vols):
    sg_aligned = sgs[(vols >= p001) & (vols <= p999)]
    crystal_systems = {
        'Triclinic (1-2)':    (sg_aligned >= 1)   & (sg_aligned <= 2),
        'Monoclinic (3-15)':  (sg_aligned >= 3)   & (sg_aligned <= 15),
        'Orthorhombic(16-74)':(sg_aligned >= 16)  & (sg_aligned <= 74),
        'Tetragonal(75-142)': (sg_aligned >= 75)  & (sg_aligned <= 142),
        'Trigonal(143-167)':  (sg_aligned >= 143) & (sg_aligned <= 167),
        'Hexagonal(168-194)': (sg_aligned >= 168) & (sg_aligned <= 194),
        'Cubic(195-230)':     (sg_aligned >= 195) & (sg_aligned <= 230),
    }
    print(f"\n=== By Crystal System ===")
    print(f"{'System':<22} {'n':>7} {'z':>8} {'p':>8} {'Sig':>5}")
    print("-"*50)
    for label, mask in crystal_systems.items():
        sub = vols_clean[mask]
        if len(sub) < 100: continue
        r2 = latticefit(sub)
        sig = "✓" if r2['p'] < 0.05 else ""
        print(f"{label:<22} {r2['n']:>7} {r2['z']:>+8.2f} "
              f"{r2['p']:>8.4f} {sig:>5}")

# ── Plot ──────────────────────────────────────────────────────────
r = results['phi']
fig, axes = plt.subplots(1, 3, figsize=(14,4))

ax = axes[0]
ax.hist(np.log10(vols_clean), bins=80,
        color='steelblue', alpha=0.7, density=True)
ax.set_xlabel("log10(volume / Å³)")
ax.set_title(f"COD Unit Cell Volumes\nn={len(vols_clean):,}")

ax = axes[1]
log_v = r['log_v']
ax.hist(log_v, bins=80, color='steelblue', alpha=0.7, density=True)
lo, hi = log_v.min(), log_v.max()
for lp in np.arange(np.floor(lo/0.25)*0.25,
                     np.ceil(hi/0.25)*0.25+0.25, 0.25):
    ax.axvline(lp, color='red', alpha=0.1, linewidth=0.5)
ax.set_xlabel("log_φ(volume)")
ax.set_title("Log-φ transform\nred = lattice points")

ax = axes[2]
ax.hist(r['null_arr'], bins=40, color='gray',
        alpha=0.7, density=True, label='Null')
ax.axvline(r['rms'], color='red', linewidth=2,
           label=f"Obs {r['rms']:.4f}")
ax.axvline(r['null_mean'], color='blue', linewidth=1.5,
           linestyle='--', label=f"Null {r['null_mean']:.4f}")
ax.set_title(f"φ-Lattice Null Test\np={r['p']:.4f}, z={r['z']:+.2f}σ")
ax.legend(fontsize=8)

plt.suptitle("LatticeFit: Mineral Crystal Unit Cell Volumes (COD)",
             fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("cod_latticefit.pdf", bbox_inches="tight", dpi=150)
print("\nSaved cod_latticefit.pdf")
