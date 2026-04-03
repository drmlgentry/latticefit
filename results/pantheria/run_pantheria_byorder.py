import pandas as pd
import numpy as np

PHI = (1+5**0.5)/2

df = pd.read_csv("pantheria.txt", sep="\t", low_memory=False)
masses_col = "5-1_AdultBodyMass_g"
order_col  = [c for c in df.columns if "Order" in c][0]
print(f"Order column: {order_col}")

df = df[[order_col, masses_col]].copy()
df[masses_col] = df[masses_col].replace(-999, np.nan)
df = df.dropna()
df = df[df[masses_col] > 0]

def latticefit(values, n_null=2000):
    log_v = np.log(values) / np.log(PHI)
    res = np.abs(log_v - np.round(log_v/0.25)*0.25)
    rms = np.sqrt(np.mean(res**2))
    lo, hi = log_v.min(), log_v.max()
    null = np.array([
        np.sqrt(np.mean(np.abs(
            u - np.round(u/0.25)*0.25)**2))
        for u in (np.random.uniform(lo,hi,len(log_v))
                  for _ in range(n_null))
    ])
    p = (null <= rms).mean()
    z = (null.mean()-rms)/null.std()
    return rms, null.mean(), p, z

print(f"\n{'Order':<22} {'n':>5} {'RMS':>7} {'Null':>7} {'p':>7} {'z':>7} {'Sig':>5}")
print("-"*60)

results = []
for order, grp in df.groupby(order_col):
    vals = grp[masses_col].values
    if len(vals) < 30:
        continue
    rms, null_m, p, z = latticefit(vals)
    sig = "✓" if p < 0.05 else ""
    print(f"{order:<22} {len(vals):>5} {rms:>7.4f} {null_m:>7.4f} "
          f"{p:>7.4f} {z:>+7.2f} {sig:>5}")
    results.append((order, len(vals), rms, null_m, p, z))

results.sort(key=lambda x: x[4])
print(f"\nTop 5 by p-value:")
for r in results[:5]:
    print(f"  {r[0]}: n={r[1]}, p={r[4]:.4f}, z={r[5]:+.2f}σ")
