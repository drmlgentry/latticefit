import numpy as np
import pandas as pd

PHI = (1+5**0.5)/2

# Parse AME2020
with open("ame2020_mass.txt", encoding="utf-8", errors="ignore") as f:
    lines = f.readlines()

records = []
for line in lines:
    if len(line) < 70: continue
    try:
        Z  = int(line[9:14].strip())
        N  = int(line[4:9].strip())
        A  = int(line[14:19].strip())
        el = line[20:23].strip()
        be_str = line[55:68].strip().replace("#","").replace("*","")
        if not be_str: continue
        be = float(be_str)
        if be > 100:
            records.append({"Z":Z,"N":N,"A":A,"el":el,"BE_A":be})
    except: continue

df = pd.DataFrame(records)
MAGIC = [2,8,20,28,50,82,126]

def dist_from_magic(x):
    return min(abs(x-m) for m in MAGIC)

df["dZ"] = df["Z"].apply(dist_from_magic)
df["dN"] = df["N"].apply(dist_from_magic)
df["far_from_magic"] = (df["dZ"]>5) & (df["dN"]>5)

# ── Physics-motivated null for "far from magic" ───────────────────
# These nuclei have BE/A predicted by the Bethe-Weizsacker formula:
# BE/A = aV - aS*A^(-1/3) - aC*Z^2*A^(-4/3) - aA*(A-2Z)^2/A^2
# Parameters: aV=15.85, aS=18.34, aC=0.71, aA=23.21 MeV
aV = 15.85; aS = 18.34; aC = 0.71; aA = 23.21  # MeV

def bethe_weizsacker(Z, A):
    N = A - Z
    if A <= 0: return 0
    bw = (aV - aS*A**(-1/3)
          - aC*Z**2*A**(-4/3)
          - aA*(A-2*Z)**2/A**2) * 1000  # keV
    return max(bw, 0)

far = df[df["far_from_magic"]].copy()
far["BW"] = far.apply(lambda r: bethe_weizsacker(r["Z"],r["A"]), axis=1)
far["residual_BW"] = far["BE_A"] - far["BW"]  # deviation from BW formula
far["ratio_BW"] = far["BE_A"] / far["BW"]     # ratio to BW prediction

print(f"Far-from-magic nuclei: {len(far)}")
print(f"BE/A range: {far['BE_A'].min():.1f} - {far['BE_A'].max():.1f} keV")
print(f"BW prediction range: {far['BW'].min():.1f} - {far['BW'].max():.1f} keV")
print(f"Ratio BE/BW: {far['ratio_BW'].min():.4f} - {far['ratio_BW'].max():.4f}")
print(f"Orders of ratio: {np.log10(far['ratio_BW'].max()/far['ratio_BW'].min()):.3f}")

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

print(f"\n=== LatticeFit: far-from-magic nuclei ===")
print(f"{'Quantity':<30} {'n':>5} {'RMS':>8} {'Null':>8} "
      f"{'p':>8} {'z':>8} {'Sig':>5}")
print("-"*70)

for label, vals in [
    ("BE/A (raw)",          far["BE_A"].values),
    ("BW prediction",       far["BW"].values),
    ("BE/BW ratio",         far["ratio_BW"].values),
    ("BE/A residual+offset",far["BE_A"].values - far["BW"].values + 1000),
]:
    vals = vals[vals > 0]
    r = latticefit(vals)
    if r:
        sig = "✓" if r["p"] < 0.05 else ""
        print(f"{label:<30} {r['n']:>5} {r['rms']:>8.4f} "
              f"{r['null_mean']:>8.4f} {r['p']:>8.4f} "
              f"{r['z']:>+8.2f} {sig:>5}")

# ── Key check: is the far-from-magic signal a range artifact? ────
print(f"\n=== Range analysis ===")
print(f"Far-from-magic log_phi span: "
      f"{np.log(far['BE_A'].max()/far['BE_A'].min())/np.log(PHI):.2f} phi-units")
print(f"Lattice bins occupied: "
      f"{(np.log(far['BE_A'].max()/far['BE_A'].min())/np.log(PHI))/0.25:.1f}")

log_be = np.log(far["BE_A"].values)/np.log(PHI)
nearest = np.round(log_be/0.25)*0.25
from collections import Counter
print(f"\nLattice bin occupancy for far-from-magic:")
for q,n in sorted(Counter(nearest).items()):
    print(f"  q={q:.2f}  BE={PHI**q:.1f} keV  n={n:4d} {'█'*(n//20)}")

# ── The BW ratio: is the deviation from liquid drop structured? ───
print(f"\n=== Shell correction structure ===")
# Shell corrections (BE - BW) should reflect nuclear shell structure
# If phi-lattice appears in shell corrections, that's the interesting result
shell_corr = far["BE_A"].values - far["BW"].values
print(f"Shell corrections range: {shell_corr.min():.1f} to {shell_corr.max():.1f} keV")
print(f"Mean shell correction: {shell_corr.mean():.2f} keV")
print(f"Std: {shell_corr.std():.2f} keV")

# Shift to positive for latticefit
sc_pos = shell_corr - shell_corr.min() + 1
r_sc = latticefit(sc_pos)
if r_sc:
    print(f"\nLatticeFit on |shell corrections| (shifted positive):")
    print(f"  n={r_sc['n']}, RMS={r_sc['rms']:.4f}, "
          f"null={r_sc['null_mean']:.4f}, "
          f"p={r_sc['p']:.4f}, z={r_sc['z']:+.2f}σ")
