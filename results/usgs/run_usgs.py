import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHI = (1+5**0.5)/2

print("Loading USGS earthquake data...")
df = pd.read_csv("usgs_earthquakes.csv")
print(f"Shape: {df.shape}")
print(f"Magnitude range: {df['mag'].min():.1f} - {df['mag'].max():.1f}")
print(f"Date range: {df['time'].min()[:10]} to {df['time'].max()[:10]}")

# Convert magnitude to energy (Gutenberg-Richter)
# E ~ 10^(1.5*M + 4.8) joules
mag = df['mag'].dropna().values
energy = 10**(1.5 * mag + 4.8)
print(f"\nEnergy range: {energy.min():.2e} to {energy.max():.2e} J")
print(f"Orders of magnitude: {np.log10(energy.max()/energy.min()):.1f}")

def latticefit(values, base=PHI, spacing=0.25, n_null=5000):
    log_v = np.log(values) / np.log(base)
    res = np.abs(log_v - np.round(log_v/spacing)*spacing)
    rms = np.sqrt(np.mean(res**2))
    lo, hi = log_v.min(), log_v.max()
    null = np.array([
        np.sqrt(np.mean(np.abs(
            u - np.round(u/spacing)*spacing)**2))
        for u in (np.random.uniform(lo,hi,len(log_v))
                  for _ in range(n_null))
    ])
    p = (null <= rms).mean()
    z = (null.mean()-rms)/null.std()
    return rms, null.mean(), null.std(), p, z, log_v, null

print("\n=== LatticeFit: Earthquake Energy ===")
print(f"{'Base':>6} {'RMS':>8} {'Null':>8} {'p':>8} {'z':>8} {'Sig':>5}")
print("-"*45)

BASES = {'phi':(PHI,0.25), '2':(2.0,1.0), 'e':(np.e,1.0), '10':(10.0,1.0)}
results = {}
for name,(base,sp) in BASES.items():
    rms,nm,ns,p,z,lv,null = latticefit(energy, base, sp)
    results[name] = (rms,nm,ns,p,z,lv,null)
    sig = "✓" if p < 0.05 else ""
    print(f"{name:>6} {rms:>8.4f} {nm:>8.4f} {p:>8.4f} {z:>+8.2f} {sig:>5}")

# Also test magnitude directly
print("\n=== LatticeFit: Magnitude (direct) ===")
print(f"{'Base':>6} {'RMS':>8} {'Null':>8} {'p':>8} {'z':>8} {'Sig':>5}")
print("-"*45)
for name,(base,sp) in BASES.items():
    rms,nm,ns,p,z,lv,null = latticefit(mag+10, base, sp)  # shift to positive
    sig = "✓" if p < 0.05 else ""
    print(f"{name:>6} {rms:>8.4f} {nm:>8.4f} {p:>8.4f} {z:>+8.2f} {sig:>5}")

# Plot for phi
rms,nm,ns,p,z,log_v,null = results['phi']

fig, axes = plt.subplots(1, 3, figsize=(14,4))

ax = axes[0]
ax.hist(mag, bins=60, color='steelblue', alpha=0.7, density=True)
ax.set_xlabel("Magnitude")
ax.set_title(f"USGS Earthquakes 2024\nn={len(mag)}, M≥3.0")

ax = axes[1]
ax.hist(log_v, bins=60, color='steelblue', alpha=0.7, density=True)
lo, hi = log_v.min(), log_v.max()
for lp in np.arange(np.floor(lo/0.25)*0.25,
                     np.ceil(hi/0.25)*0.25+0.25, 0.25):
    ax.axvline(lp, color='red', alpha=0.15, linewidth=0.5)
ax.set_xlabel("log_φ(energy)")
ax.set_title("Log-φ transform\nred = lattice points")

ax = axes[2]
ax.hist(null, bins=40, color='gray', alpha=0.7, density=True, label="Null")
ax.axvline(rms, color='red', linewidth=2, label=f"Obs {rms:.4f}")
ax.axvline(nm, color='blue', linewidth=1.5,
           linestyle='--', label=f"Null {nm:.4f}")
ax.set_title(f"φ-Lattice Null Test\np={p:.4f}, z={z:+.2f}σ")
ax.legend(fontsize=8)

plt.suptitle("LatticeFit: USGS Earthquake Energies 2024",
             fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("usgs_latticefit.pdf", bbox_inches="tight", dpi=150)
print("\nSaved usgs_latticefit.pdf")

print("\n=== CROSS-DOMAIN SUMMARY ===")
rms_e,nm_e,ns_e,p_e,z_e = results['phi'][:5]
print(f"{'Domain':<32} {'n':>7} {'RMS':>8} {'z':>8} {'p':>8}")
print("-"*62)
for row in [
    ("SM fermion masses",       12,    0.055,  +4.50, "<0.001"),
    ("Rice 44K MAF",         36901,   0.0712,  +4.25, "<0.001"),
    ("COVID variant freq",  103348,   0.0712,  +9.38, "<0.001"),
    ("Cetacea body masses",    76,    0.0631,  +2.43,  "0.008"),
    (f"Earthquakes 2024 (E)", len(energy), rms_e, z_e, f"{p_e:.4f}"),
    ("All mammals",          3542,    0.0714,  +1.45,  "0.074"),
    ("HIV RT IC50",         10041,    0.0733,  -3.73, "0.9996"),
]:
    print(f"{row[0]:<32} {row[1]:>7} {row[2]:>8.4f} "
          f"{row[3]:>+8.2f} {row[4]:>8}")
