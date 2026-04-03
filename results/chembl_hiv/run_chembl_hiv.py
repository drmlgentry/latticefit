import json, urllib.request, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHI = (1+5**0.5)/2

# Load first page already downloaded
print("Loading page 1...")
with open("chembl_rt_ic50.json") as f:
    data = json.load(f)

total = data["page_meta"]["total_count"]
print(f"Total records: {total}")

records = data["activities"]

# Fetch remaining pages
offset = 1000
while offset < total:
    url = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
           f"?target_chembl_id=CHEMBL247&standard_type=IC50"
           f"&standard_units=nM&limit=1000&offset={offset}&format=json")
    print(f"  Fetching offset {offset}/{total}...")
    try:
        req = urllib.request.Request(url, headers={"Accept":"application/json"})
        with urllib.request.urlopen(req, timeout=30) as r:
            page = json.load(r)
        records.extend(page["activities"])
        offset += 1000
        time.sleep(0.5)
    except Exception as e:
        print(f"  Error at offset {offset}: {e}")
        break

print(f"Total fetched: {len(records)}")

# Extract IC50 values
ic50_vals = []
for r in records:
    v = r.get("standard_value")
    if v is not None:
        try:
            val = float(v)
            if val > 0:
                ic50_vals.append(val)
        except:
            pass

ic50 = np.array(ic50_vals)
print(f"\nValid IC50 values: {len(ic50)}")
print(f"Range: {ic50.min():.3f} - {ic50.max():.0f} nM")
print(f"Orders of magnitude: {np.log10(ic50.max()/ic50.min()):.1f}")

# Save
pd.Series(ic50, name="IC50_nM").to_csv("chembl_hiv_ic50.csv", index=False)
print("Saved chembl_hiv_ic50.csv")

# ── LatticeFit ────────────────────────────────────────────────────
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

print("\nRunning LatticeFit...")
rms,nm,ns,p,z,log_v,null = latticefit(ic50)

print(f"\n=== ChEMBL HIV RT IC50 LatticeFit ===")
print(f"n:           {len(ic50)}")
print(f"Range:       {ic50.min():.3f} - {ic50.max():.0f} nM")
print(f"Log10 span:  {np.log10(ic50.max()/ic50.min()):.1f} orders")
print(f"RMS:         {rms:.4f} log-phi units")
print(f"Null mean:   {nm:.4f}")
print(f"z-score:     {z:+.2f} sigma")
print(f"p-value:     {p:.4f}")
print(f"Significant: {'YES' if p < 0.05 else 'NO'}")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(14,4))

ax = axes[0]
ax.hist(np.log10(ic50), bins=60, color="steelblue", alpha=0.7, density=True)
ax.set_xlabel("log10(IC50 / nM)")
ax.set_title(f"HIV RT IC50 distribution\nn={len(ic50)}, ChEMBL")

ax = axes[1]
ax.hist(log_v, bins=60, color="steelblue", alpha=0.7, density=True)
lo, hi = log_v.min(), log_v.max()
for lp in np.arange(np.floor(lo/0.25)*0.25,
                     np.ceil(hi/0.25)*0.25+0.25, 0.25):
    ax.axvline(lp, color="red", alpha=0.15, linewidth=0.5)
ax.set_xlabel("log_φ(IC50)")
ax.set_title("Log-φ transform\nred = lattice points")

ax = axes[2]
ax.hist(null, bins=40, color="gray", alpha=0.7, density=True, label="Null")
ax.axvline(rms, color="red", linewidth=2, label=f"Obs {rms:.4f}")
ax.axvline(nm, color="blue", linewidth=1.5,
           linestyle="--", label=f"Null {nm:.4f}")
ax.set_title(f"φ-Lattice Null Test\np={p:.4f}, z={z:+.2f}σ")
ax.legend(fontsize=8)

plt.suptitle("LatticeFit: HIV RT IC50 Values (ChEMBL)",
             fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("chembl_hiv_latticefit.pdf", bbox_inches="tight", dpi=150)
print("\nSaved chembl_hiv_latticefit.pdf")

# Cross-domain table
print("\n=== CROSS-DOMAIN SUMMARY ===")
print(f"{'Domain':<32} {'n':>7} {'RMS':>8} {'z':>8} {'p':>8}")
print("-"*62)
for row in [
    ("SM fermion masses",       12,    0.055,  +4.50, "<0.001"),
    ("Rice 44K MAF",         36901,   0.0712,  +4.25, "<0.001"),
    ("COVID variant freq",  103348,   0.0712,  +9.38, "<0.001"),
    ("Cetacea body masses",    76,   0.0631,  +2.43,  "0.008"),
    (f"HIV RT IC50 (ChEMBL)", len(ic50), rms, z, f"{p:.4f}"),
    ("All mammals",          3542,   0.0714,  +1.45,  "0.074"),
]:
    print(f"{row[0]:<32} {row[1]:>7} {row[2]:>8.4f} "
          f"{row[3]:>+8.2f} {row[4]:>8}")
