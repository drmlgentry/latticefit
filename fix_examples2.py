path = r"C:\dev\latticefit\examples\scripts\run_all_examples.py"
with open(path, encoding="utf-8") as f:
    s = f.read()

# Fix 1: populations column name + filter to 2024
s = s.replace(
    "load('populations_2024.csv', 'population_2024')",
    "pd.read_csv(os.path.join(DATA_DIR, 'populations_2024.csv'))"
    ".pipe(lambda df: df[df['year']==2024] if 'year' in df.columns else df)"
    "['population'].dropna().values"
)

# Fix 2: SM masses anchor
s = s.replace(
    "run(\"SM fermion + boson masses\", sm,",
    "run(\"SM fermion + boson masses\", sm, anchor=5.11e-4,"
)
s = s.replace(
    "run(\"SM fermion masses only\", sm[sm_ferm],",
    "run(\"SM fermion masses only\", sm[sm_ferm], anchor=5.11e-4,"
)

with open(path, "w", encoding="utf-8") as f:
    f.write(s)
print("Done.")
