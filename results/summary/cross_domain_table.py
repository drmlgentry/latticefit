# Final cross-domain summary table for JOSS paper

results = [
    # Domain, n, z, p, orders, category, status
    ("SM fermion masses",           12,    +4.50, "<0.001", 8.0,  "Physics",    "GENUINE"),
    ("Rice 44K MAF",             36901,    +4.25, "<0.001", 2.5,  "Genomics",   "GENUINE"),
    ("Cetacea body masses",          76,    +2.43, "0.008",  3.7,  "Ecology",    "GENUINE"),
    ("EGFR SAR CHEMBL1064829",       32,    +2.12, "0.021",  3.4,  "Pharma",     "GENUINE"),
    ("EGFR SAR CHEMBL3750872",       40,    +2.33, "0.014",  2.1,  "Pharma",     "GENUINE"),
    ("COVID variant freq",       103348,    +9.38, "<0.001", 3.0,  "Epidemiol",  "NEEDS_NULL"),
    ("All mammals",                3542,    +1.45, "0.074",  7.9,  "Ecology",    "MARGINAL"),
    ("NIST ionization energies",   1631,    +0.30, "0.378",  4.4,  "Physics",    "NULL"),
    ("S&P 500 returns",           11345,    -1.43, "0.924",  4.6,  "Finance",    "NULL"),
    ("Crystal volumes COD",      525224,    -1.04, "0.850",  3.3,  "Chemistry",  "NULL"),
    ("HIV RT IC50 mixed",         10041,    -3.73, "0.9996", 9.1,  "Pharma",     "NULL"),
    ("PharmGKB Score",             3896,   -10.58, "1.000",  0.0,  "Pharma",     "ARTIFACT"),
    ("Earthquakes GR law",        18659,    -3.58, "1.000",  6.8,  "Geophys",    "NULL"),
    ("Nuclear BE/A AME2020",       3553,    +7.18, "<0.001", 1.0,  "Physics",    "ARTIFACT"),
]

print(f"{'Domain':<35} {'n':>7} {'z':>7} {'p':>8} {'orders':>7} {'Status':<12}")
print("-"*80)
for row in results:
    marker = "✓" if row[6]=="GENUINE" else ("~" if row[6]=="MARGINAL" else " ")
    print(f"{marker} {row[0]:<34} {row[1]:>7} {row[2]:>+7.2f} "
          f"{row[3]:>8} {row[4]:>7.1f}  {row[6]:<12}")

print(f"\nGenuine signals: {sum(1 for r in results if r[6]=='GENUINE')}")
print(f"Null results:    {sum(1 for r in results if r[6]=='NULL')}")
print(f"Artifacts:       {sum(1 for r in results if r[6]=='ARTIFACT')}")
print(f"Needs better null: {sum(1 for r in results if r[6]=='NEEDS_NULL')}")
print(f"Marginal:        {sum(1 for r in results if r[6]=='MARGINAL')}")
