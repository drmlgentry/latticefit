with open("latticefit_app.py", encoding="utf-8") as f:
    src = f.read()

old = '''def run_latticefit_full(vals, n_null=5000, fixed_anchor=None):
    """Run latticefit across all standard bases, return dict of results."""
    rng = np.random.default_rng(42)
    lo, hi = np.log(vals.min()), np.log(vals.max())
    results = {}

    for bname, base in BASES.items():
        for d in [2, 3, 4, 6, 8]:
            anchor = fixed_anchor if fixed_anchor is not None else vals.min()
            ks   = np.round(d * np.log(vals/anchor) / np.log(base))
            pred = anchor * base**(ks/d)
            res  = np.abs(np.log(vals/pred) / np.log(base))
            rms  = float(np.sqrt(np.mean(res**2)))
            max_rms = 0.5 / d

            null_rms = np.array([
                np.sqrt(np.mean(np.abs(
                    np.log(np.exp(rng.uniform(lo, hi, len(vals)))/
                           (anchor * base**(np.round(d*np.log(
                               np.exp(rng.uniform(lo, hi, len(vals)))/anchor)/np.log(base))/d))
                ) / np.log(base))**2))
                for _ in range(n_null)
            ])
            # Vectorised null
            rand_mat = np.exp(rng.uniform(lo, hi, (n_null, len(vals))))
            ks_mat   = np.round(d * np.log(rand_mat/anchor) / np.log(base))
            pred_mat = anchor * base**(ks_mat/d)
            res_mat  = np.abs(np.log(rand_mat/pred_mat) / np.log(base))
            null_rms = np.sqrt(np.mean(res_mat**2, axis=1))'''

new = '''def run_latticefit_full(vals, n_null=5000, fixed_anchor=None):
    """Run latticefit across all standard bases, return dict of results."""
    rng = np.random.default_rng(42)
    results = {}

    for bname, base in BASES.items():
        for d in [2, 3, 4, 6, 8]:
            anchor = fixed_anchor if fixed_anchor is not None else vals.min()
            ks   = np.round(d * np.log(vals/anchor) / np.log(base))
            pred = anchor * base**(ks/d)
            res  = np.abs(np.log(vals/pred) / np.log(base))
            rms  = float(np.sqrt(np.mean(res**2)))
            max_rms = 0.5 / d

            # Null: draw log-uniform in same range AS SEEN BY THE LATTICE
            # i.e. uniform in log-space between min and max label
            k_min = np.floor(d * np.log(vals.min()/anchor) / np.log(base))
            k_max = np.ceil( d * np.log(vals.max()/anchor) / np.log(base))
            lo_k = k_min / d * np.log(base) + np.log(anchor)
            hi_k = k_max / d * np.log(base) + np.log(anchor)
            # Vectorised null using same anchor
            rand_mat = np.exp(rng.uniform(lo_k, hi_k, (n_null, len(vals))))
            ks_mat   = np.round(d * np.log(rand_mat/anchor) / np.log(base))
            pred_mat = anchor * base**(ks_mat/d)
            res_mat  = np.abs(np.log(rand_mat/pred_mat) / np.log(base))
            null_rms = np.sqrt(np.mean(res_mat**2, axis=1))'''

if old in src:
    src = src.replace(old, new)
    print("Replaced successfully.")
else:
    print("Pattern not found - checking...")
    idx = src.find("def run_latticefit_full")
    print(src[idx:idx+1200])

with open("latticefit_app.py", "w", encoding="utf-8") as f:
    f.write(src)
