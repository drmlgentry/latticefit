import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHI = (1+5**0.5)/2

df = pd.read_csv("pantheria.txt", sep="\t", low_memory=False)
mass_col  = "5-1_AdultBodyMass_g"
order_col = "MSW05_Order"
name_col  = "MSW05_Binomial"

cet = df[df[order_col]=="Cetacea"][[name_col, mass_col]].copy()
cet[mass_col] = cet[mass_col].replace(-999, np.nan)
cet = cet.dropna()
cet = cet[cet[mass_col] > 0].sort_values(mass_col)

print(f"Cetacea: {len(cet)} species")
print(f"Range: {cet[mass_col].min():.0f} g to {cet[mass_col].max():.0f} g")
print(f"Orders of magnitude: {np.log10(cet[mass_col].max()/cet[mass_col].min()):.1f}")

# Compute log-phi values and nearest lattice points
log_v   = np.log(cet[mass_col].values) / np.log(PHI)
nearest = np.round(log_v / 0.25) * 0.25
residuals = log_v - nearest
cet = cet.copy()
cet["log_phi"] = log_v
cet["nearest"] = nearest
cet["residual"] = residuals
cet["predicted_g"] = PHI**(nearest) 

# Show species near lattice points (residual < 0.05)
print(f"\nSpecies within 0.05 log-phi units of a lattice point:")
close = cet[np.abs(cet["residual"]) < 0.05].sort_values("log_phi")
for _, row in close.iterrows():
    print(f"  {row[name_col]:<35} {row[mass_col]:>12.0f} g  "
          f"q={row['nearest']:.2f}  resid={row['residual']:+.3f}")

# Null test
lo, hi = log_v.min(), log_v.max()
null_rms = np.array([
    np.sqrt(np.mean(np.abs(
        u - np.round(u/0.25)*0.25)**2))
    for u in (np.random.uniform(lo,hi,len(log_v))
              for _ in range(5000))
])
rms = np.sqrt(np.mean(residuals**2))
p   = (null_rms <= rms).mean()
z   = (null_rms.mean()-rms)/null_rms.std()

print(f"\n=== CETACEA LATTICEFIT RESULT ===")
print(f"n:            {len(cet)}")
print(f"RMS:          {rms:.4f} log-phi units")
print(f"Null mean:    {null_rms.mean():.4f}")
print(f"z-score:      +{z:.2f} sigma")
print(f"p-value:      {p:.4f}")

# Compare to SM masses result
print(f"\n=== CROSS-DOMAIN COMPARISON ===")
print(f"{'Domain':<30} {'n':>6} {'RMS':>8} {'z':>8} {'p':>8}")
print("-"*55)
for domain, n, rms_v, z_v, p_v in [
    ("SM fermion masses",        12,  0.055, 4.5,  "<0.001"),
    ("Rice 44K MAF",          36901,  0.0712, 4.25, "<0.001"),
    ("Cetacea body masses",      76,  rms,   z,    f"{p:.4f}"),
    ("All mammals",            3542,  0.0714, 1.45, "0.074"),
]:
    print(f"{domain:<30} {n:>6} {rms_v:>8.4f} {z_v:>+8.2f} {p_v:>8}")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

# 1. Mass distribution
ax = axes[0]
ax.barh(range(len(cet)), np.log10(cet[mass_col].values),
        color="steelblue", alpha=0.7, height=0.8)
ax.set_xlabel("log10(mass / g)")
ax.set_title(f"Cetacea body masses\n{len(cet)} species")
ax.set_yticks([])

# 2. Log-phi with lattice
ax = axes[1]
ax.hist(log_v, bins=25, color="steelblue", alpha=0.7, density=True)
lattice_pts = np.arange(np.floor(lo/0.25)*0.25,
                         np.ceil(hi/0.25)*0.25+0.25, 0.25)
for lp in lattice_pts:
    ax.axvline(lp, color="red", alpha=0.25, linewidth=1.0)
ax.set_xlabel("log_φ(mass)")
ax.set_title("Log-φ with lattice lines\n(red = φ-lattice points)")

# 3. Null test
ax = axes[2]
ax.hist(null_rms, bins=40, color="gray", alpha=0.7,
        density=True, label="Null")
ax.axvline(rms, color="red", linewidth=2,
           label=f"Observed {rms:.4f}")
ax.axvline(null_rms.mean(), color="blue", linewidth=1.5,
           linestyle="--", label=f"Null {null_rms.mean():.4f}")
ax.set_xlabel("RMS residual (log-φ units)")
ax.set_title(f"Null test: Cetacea\np={p:.4f}, z=+{z:.2f}σ")
ax.legend(fontsize=8)

plt.suptitle("LatticeFit: Cetacea Body Masses — φ-Lattice Analysis",
             fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("cetacea_latticefit.pdf", bbox_inches="tight", dpi=150)
plt.savefig("cetacea_latticefit.png", bbox_inches="tight", dpi=150)
print("\nSaved cetacea_latticefit.pdf/.png")
