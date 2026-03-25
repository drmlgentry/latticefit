path = r"C:\dev\latticefit\examples\scripts\run_all_examples.py"
with open(path, encoding="utf-8") as f:
    lines = f.readlines()

# Fix line 87-88: duplicate anchor keyword
for i, line in enumerate(lines):
    if 'run("SM fermion + boson masses"' in line:
        lines[i] = '    run("SM fermion + boson masses", sm,\n'
        lines[i+1] = '        base=PHI, denom=4, anchor=5.11e-4)\n'
        print(f"Fixed SM masses at line {i+1}")

# Fix populations column: line 117 truncated — check full content
for i, line in enumerate(lines):
    if "populations_2024.csv" in line and "pipe" in line:
        lines[i] = ("    pop = pd.read_csv(os.path.join(DATA_DIR, 'populations_2024.csv'))\n")
        lines.insert(i+1, "    pop = pop[pop['year']==2024]['population'].dropna().values\n")
        print(f"Fixed populations at line {i+1}")
        break

with open(path, "w", encoding="utf-8") as f:
    f.writelines(lines)
print("Done.")
