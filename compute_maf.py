import pandas as pd
import numpy as np

print("Loading...")
df = pd.read_csv("RiceDiversity.44K.MSU6.Genotypes.csv.gz", index_col=0)
print(f"Shape: {df.shape}")

# Drop the chr and position columns (first two)
# They appear as integer values mixed into the genotype data
geno = df.iloc[:, 2:]  # skip chr, position columns
print(f"Genotype matrix: {geno.shape}")

ALLELES = set('ACGT')
MISSING  = {'0', 'N'}

def maf(row):
    counts = {}
    for v in row:
        if v in ALLELES:
            counts[v] = counts.get(v, 0) + 1
    if len(counts) < 2:
        return np.nan
    total = sum(counts.values())
    minor = min(counts.values())
    return minor / total

print("Computing MAF...")
mafs = geno.apply(maf, axis=1).dropna()
mafs = mafs[(mafs > 0) & (mafs <= 0.5)]
print(f"\nSNPs with valid MAF: {len(mafs)}")
print(f"MAF range: {mafs.min():.4f} - {mafs.max():.4f}")
print(f"\nMAF distribution:\n{mafs.describe()}")
print(f"\nMAF histogram (10 bins):")
hist, edges = np.histogram(mafs, bins=10, range=(0, 0.5))
for i, (lo, hi, n) in enumerate(zip(edges, edges[1:], hist)):
    bar = '#' * (n // 200)
    print(f"  {lo:.2f}-{hi:.2f}: {n:6d} {bar}")

mafs.to_csv("rice_44k_maf.csv", header=True)
print("\nSaved rice_44k_maf.csv")
