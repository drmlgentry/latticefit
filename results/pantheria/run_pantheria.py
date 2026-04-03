import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHI = (1 + 5**0.5) / 2

print("Loading PanTHERIA...")
df = pd.read_csv("pantheria.txt", sep="\t", low_memory=False)

# Extract adult body mass, drop missing (-999 is the null sentinel)
masses = df["5-1_AdultBodyMass_g"].replace(-999, np.nan).dropna()
masses = masses[masses > 0]
print(f"Species with valid adult body mass: {len(masses)}")
print(f"Range: {masses.min():.2f} g ({masses.min()/1000:.2f} kg) "
      f"to {masses.max():.2f} g ({masses.max()/1000:.0f} kg)")
print(f"Orders of magnitude: {np.log10(masses.max()/masses.min()):.1f}")

# Show some familiar species for sanity check
idx_map = df.set_index("MSW05_Binomial")["5-1_AdultBodyMass_g"].replace(-999, np.nan)
for sp in ["Mus musculus","Homo sapiens","Loxodonta africana",
           "Balaenoptera musculus","Canis lupus"]:
    if sp in idx_map.index:
        v = idx_map[sp]
        print(f"  {sp}: {v:.0f} g")

# ── LatticeFit phi analysis ───────────────────────────────────────
log_v = np.log(masses.values) / np.log(PHI)
nearest = np.round(log_v / 0.25) * 0.25
residuals = log_v - nearest
rms = np.sqrt(np.mean(residuals**2))

print(f"\nRunning null test (5000 permutations)...")
lo, hi = log_v.min(), log_v.max()
null_rms = np.array([
    np.sqrt(np.mean(np.abs(
        u - np.round(u/0.25)*0.25
    )**2))
    for u in (np.random.uniform(lo, hi, len(log_v))
              for _ in range(5000))
])
p = (null_rms <= rms).mean()
z = (null_rms.mean() - rms) / null_rms.std()

print(f"\n=== RESULTS: Mammal Body Masses (PanTHERIA) ===")
print(f"n species:     {len(masses)}")
print(f"Mass range:    {masses.min():.1f} g to {masses.max():.0f} g")
print(f"Log10 span:    {np.log10(masses.max()/masses.min()):.1f} orders")
print(f"Lattice:       phi = {PHI:.6f}, spacing = 1/4")
print(f"Observed RMS:  {rms:.4f} log-phi units")
print(f"Null mean:     {null_rms.mean():.4f}")
print(f"Null std:      {null_rms.std():.4f}")
print(f"z-score:       +{z:.2f} sigma")
print(f"p-value:       {p:.4f}")
print(f"Significant:   {'YES p<0.05' if p < 0.05 else 'NO'}")

# ── Plots ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

# 1. Mass histogram (log scale)
ax = axes[0]
ax.hist(np.log10(masses), bins=60, color="steelblue",
        alpha=0.7, density=True)
ax.set_xlabel("log10(body mass / g)")
ax.set_ylabel("Density")
ax.set_title(f"Mammal Body Masses\n{len(masses)} species, PanTHERIA")

# 2. Log-phi with lattice lines
ax = axes[1]
ax.hist(log_v, bins=80, color="steelblue", alpha=0.7, density=True)
lattice_pts = np.arange(np.floor(lo/0.25)*0.25,
                         np.ceil(hi/0.25)*0.25+0.25, 0.25)
for lp in lattice_pts:
    ax.axvline(lp, color="red", alpha=0.15, linewidth=0.5)
ax.set_xlabel("log_φ(mass)")
ax.set_ylabel("Density")
ax.set_title("Log-φ Transform\nred = φ-lattice points")

# 3. Null test
ax = axes[2]
ax.hist(null_rms, bins=40, color="gray", alpha=0.7,
        density=True, label="Null")
ax.axvline(rms, color="red", linewidth=2,
           label=f"Observed {rms:.4f}")
ax.axvline(null_rms.mean(), color="blue", linewidth=1.5,
           linestyle="--", label=f"Null mean {null_rms.mean():.4f}")
ax.set_xlabel("RMS residual (log-φ units)")
ax.set_ylabel("Density")
ax.set_title(f"φ-Lattice Null Test\np={p:.4f}, z=+{z:.2f}σ")
ax.legend(fontsize=8)

plt.suptitle("LatticeFit: Mammal Body Masses — φ-Lattice Analysis",
             fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("pantheria_latticefit.pdf", bbox_inches="tight", dpi=150)
plt.savefig("pantheria_latticefit.png", bbox_inches="tight", dpi=150)
print("\nSaved pantheria_latticefit.pdf/.png")
