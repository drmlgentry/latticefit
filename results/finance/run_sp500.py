import numpy as np
import pandas as pd
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

# ── Load data ─────────────────────────────────────────────────────
prices = pd.read_csv("sp500_prices.csv", index_col=0,
                     parse_dates=True).squeeze()
abs_ret = pd.read_csv("sp500_returns.csv", index_col=0,
                      parse_dates=True).squeeze()

print(f"S&P 500: {len(prices)} days, {prices.index[0].date()} - "
      f"{prices.index[-1].date()}")
print(f"Price range: {prices.min():.1f} - {prices.max():.1f} "
      f"({np.log10(prices.max()/prices.min()):.1f} orders)")
print(f"Abs returns: {abs_ret.min():.6f} - {abs_ret.max():.4f} "
      f"({np.log10(abs_ret.max()/abs_ret.min()):.1f} orders)")

# ── Test 1: absolute daily returns ────────────────────────────────
print("\n=== Absolute daily returns ===")
print(f"{'Base':>6} {'RMS':>8} {'Null':>8} {'p':>8} {'z':>8} {'Sig':>5}")
print("-"*45)
ret_results = {}
for name, base, sp in [('phi',PHI,0.25),('2',2.,1.),
                        ('e',np.e,1.),('10',10.,1.)]:
    r = latticefit(abs_ret.values, base, sp)
    ret_results[name] = r
    sig = "✓" if r['p'] < 0.05 else ""
    print(f"{name:>6} {r['rms']:>8.4f} {r['null_mean']:>8.4f} "
          f"{r['p']:>8.4f} {r['z']:>+8.2f} {sig:>5}")

# ── Test 2: price levels ──────────────────────────────────────────
print("\n=== Price levels (1.8 orders - marginal range) ===")
for name, base, sp in [('phi',PHI,0.25),('2',2.,1.),
                        ('e',np.e,1.),('10',10.,1.)]:
    r = latticefit(prices.values, base, sp)
    sig = "✓" if r['p'] < 0.05 else ""
    print(f"{name:>6} {r['rms']:>8.4f} {r['null_mean']:>8.4f} "
          f"{r['p']:>8.4f} {r['z']:>+8.2f} {sig:>5}")

# ── Test 3: by decade ─────────────────────────────────────────────
print("\n=== Absolute returns by decade ===")
print(f"{'Decade':<12} {'n':>5} {'z':>8} {'p':>8} {'Sig':>5}")
print("-"*38)
for decade in [1980, 1990, 2000, 2010, 2020]:
    mask = (abs_ret.index.year >= decade) & \
           (abs_ret.index.year < decade+10)
    sub = abs_ret[mask].values
    if len(sub) < 100: continue
    r = latticefit(sub)
    sig = "✓" if r['p'] < 0.05 else ""
    print(f"{decade}s       {r['n']:>5} {r['z']:>+8.2f} "
          f"{r['p']:>8.4f} {sig:>5}")

# ── Test 4: crash periods vs calm periods ────────────────────────
print("\n=== Crash vs calm periods ===")
threshold = abs_ret.quantile(0.95)  # top 5% = volatile days
volatile = abs_ret[abs_ret >= threshold].values
calm     = abs_ret[abs_ret <  threshold].values
for label, vals in [("Volatile (top 5%)", volatile),
                     ("Calm (bottom 95%)", calm)]:
    r = latticefit(vals)
    sig = "✓" if r['p'] < 0.05 else ""
    print(f"{label:<22} n={r['n']:5d} z={r['z']:+.2f} "
          f"p={r['p']:.4f} {sig}")

# ── Plot ──────────────────────────────────────────────────────────
r = ret_results['phi']
fig, axes = plt.subplots(1, 3, figsize=(14,4))

ax = axes[0]
ax.semilogy(prices.index, prices.values, linewidth=0.5, color='steelblue')
ax.set_title(f"S&P 500 Price (1980-2024)\n{len(prices)} trading days")
ax.set_ylabel("Index level")

ax = axes[1]
log_v = r['log_v']
ax.hist(log_v, bins=80, color='steelblue', alpha=0.7, density=True)
lo, hi = log_v.min(), log_v.max()
for lp in np.arange(np.floor(lo/0.25)*0.25,
                     np.ceil(hi/0.25)*0.25+0.25, 0.25):
    ax.axvline(lp, color='red', alpha=0.15, linewidth=0.7)
ax.set_xlabel("log_φ(|return|)")
ax.set_title("Absolute returns: log-φ\nred = φ-lattice points")

ax = axes[2]
ax.hist(r['null_arr'], bins=40, color='gray',
        alpha=0.7, density=True, label='Null')
ax.axvline(r['rms'], color='red', linewidth=2,
           label=f"Obs {r['rms']:.4f}")
ax.axvline(r['null_mean'], color='blue', linewidth=1.5,
           linestyle='--', label=f"Null {r['null_mean']:.4f}")
ax.set_title(f"φ-Lattice Null Test\np={r['p']:.4f}, z={r['z']:+.2f}σ")
ax.legend(fontsize=8)

plt.suptitle("LatticeFit: S&P 500 Daily Returns (1980-2024)",
             fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig("sp500_latticefit.pdf", bbox_inches='tight', dpi=150)
print("\nSaved sp500_latticefit.pdf")
