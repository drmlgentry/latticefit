import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHI = (1+5**0.5)/2

print("Loading CDC variant data...")
df = pd.read_csv("cdc_variants.csv")
print(f"Shape: {df.shape}")
print(f"Variants: {df['variant'].nunique()}")
print(f"Date range: {df['week_ending'].min()} to {df['week_ending'].max()}")
print(f"Variants present:\n{df['variant'].value_counts().head(20)}")

# Focus on USA national, weekly, nowcast model
usa = df[(df['usa_or_hhsregion']=='USA') & 
         (df['share'] > 0) &
         (df['share'] < 1)].copy()

print(f"\nUSA national rows: {len(usa)}")

# ── Analysis 1: all variant frequency observations ────────────────
shares = usa['share'].values
shares = shares[(shares > 0.001) & (shares < 0.999)]
print(f"Valid share values: {len(shares)}")
print(f"Share range: {shares.min():.4f} to {shares.max():.4f}")

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
    return rms, null.mean(), null.std(), p, z

print("\n=== LatticeFit: All variant frequency observations ===")
print(f"{'Base':>6} {'RMS':>8} {'Null':>8} {'p':>8} {'z':>8} {'Sig':>5}")
print("-"*45)
for name, base, sp in [('phi',PHI,0.25),('2',2.0,1.0),
                        ('e',np.e,1.0),('10',10.0,1.0)]:
    rms,nm,ns,p,z = latticefit(shares, base, sp)
    sig = "✓" if p < 0.05 else ""
    print(f"{name:>6} {rms:>8.4f} {nm:>8.4f} {p:>8.4f} {z:>+8.2f} {sig:>5}")

# ── Analysis 2: peak frequency per variant ────────────────────────
peak = usa.groupby('variant')['share'].max()
peak = peak[peak > 0.01]  # variants that reached >1%
print(f"\n=== Peak frequency per variant (n={len(peak)}) ===")
print(peak.sort_values(ascending=False).to_string())

rms,nm,ns,p,z = latticefit(peak.values)
print(f"\nLatticeFit on peak frequencies:")
print(f"  RMS={rms:.4f}, null={nm:.4f}, p={p:.4f}, z={z:+.2f}σ")

# ── Analysis 3: growth rates (week-over-week ratios) ──────────────
usa2 = usa.sort_values(['variant','week_ending'])
ratios = []
for variant, grp in usa2.groupby('variant'):
    grp = grp.sort_values('week_ending')
    shares_v = grp['share'].values
    for i in range(1, len(shares_v)):
        if shares_v[i-1] > 0.005 and shares_v[i] > 0.005:
            r = shares_v[i] / shares_v[i-1]
            if 0.1 < r < 10:  # reasonable growth ratio
                ratios.append(r)

ratios = np.array(ratios)
print(f"\n=== Week-over-week growth ratios (n={len(ratios)}) ===")
print(f"Range: {ratios.min():.3f} to {ratios.max():.3f}")

if len(ratios) > 20:
    rms,nm,ns,p,z = latticefit(ratios)
    print(f"LatticeFit on growth ratios:")
    print(f"  RMS={rms:.4f}, null={nm:.4f}, p={p:.4f}, z={z:+.2f}σ")
    sig = "YES" if p < 0.05 else "NO"
    print(f"  Significant: {sig}")

# ── Summary plot ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

# 1. Variant share over time for top variants
ax = axes[0]
top_vars = usa.groupby('variant')['share'].max().nlargest(8).index
for v in top_vars:
    vdata = usa[usa['variant']==v].sort_values('week_ending')
    ax.plot(pd.to_datetime(vdata['week_ending']),
            vdata['share'], label=v, linewidth=1.5)
ax.set_xlabel("Date")
ax.set_ylabel("Frequency")
ax.set_title("SARS-CoV-2 Variant Frequencies\n(USA, CDC)")
ax.legend(fontsize=6, ncol=2)
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

# 2. Log-phi of all share values with lattice
log_v = np.log(shares)/np.log(PHI)
ax = axes[1]
ax.hist(log_v, bins=60, color="steelblue", alpha=0.7, density=True)
lo, hi = log_v.min(), log_v.max()
for lp in np.arange(np.floor(lo/0.25)*0.25,
                     np.ceil(hi/0.25)*0.25+0.25, 0.25):
    ax.axvline(lp, color="red", alpha=0.15, linewidth=0.5)
ax.set_xlabel("log_φ(variant frequency)")
ax.set_title("Log-φ Transform\nred = φ-lattice points")

# 3. Null test for shares
rms_s,nm_s,ns_s,p_s,z_s = latticefit(shares)
null_plot = np.array([
    np.sqrt(np.mean(np.abs(
        u-np.round(u/0.25)*0.25)**2))
    for u in (np.random.uniform(lo,hi,len(log_v))
              for _ in range(1000))
])
ax = axes[2]
ax.hist(null_plot, bins=40, color="gray", alpha=0.7,
        density=True, label="Null")
ax.axvline(rms_s, color="red", linewidth=2,
           label=f"Observed {rms_s:.4f}")
ax.axvline(null_plot.mean(), color="blue", linewidth=1.5,
           linestyle="--", label=f"Null {null_plot.mean():.4f}")
ax.set_xlabel("RMS residual (log-φ units)")
ax.set_title(f"φ-Lattice Null Test\np={p_s:.4f}, z={z_s:+.2f}σ")
ax.legend(fontsize=8)

plt.suptitle("LatticeFit: SARS-CoV-2 Variant Frequencies (CDC)",
             fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("covid_variants_latticefit.pdf", bbox_inches="tight", dpi=150)
plt.savefig("covid_variants_latticefit.png", bbox_inches="tight", dpi=150)
print("\nSaved covid_variants_latticefit.pdf/.png")

# ── Cross-domain table ────────────────────────────────────────────
print("\n=== CROSS-DOMAIN SUMMARY ===")
print(f"{'Domain':<32} {'n':>7} {'RMS':>8} {'z':>8} {'p':>8}")
print("-"*60)
for row in [
    ("SM fermion masses",        12,    0.055, +4.50, "<0.001"),
    ("Rice 44K MAF",          36901,   0.0712, +4.25, "<0.001"),
    ("Cetacea body masses",      76,   0.0631, +2.43, "0.008"),
    (f"COVID variant freq ({len(shares)})", len(shares), rms_s, z_s, f"{p_s:.4f}"),
    ("All mammals",            3542,   0.0714, +1.45, "0.074"),
]:
    print(f"{row[0]:<32} {row[1]:>7} {row[2]:>8.4f} {row[3]:>+8.2f} {row[4]:>8}")
