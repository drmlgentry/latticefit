path = r"C:\dev\latticefit\examples\scripts\run_all_examples.py"
with open(path, encoding="utf-8") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    # Fix SM masses anchor - ensure electron mass is used
    if 'run("SM fermion + boson masses"' in line:
        lines[i] = '    run("SM fermion + boson masses", sm,\n'
        lines[i+1] = '        base=PHI, denom=4, anchor=5.11e-4)\n'
        print(f"Fixed SM anchor at line {i+1}")

    # Fix populations year filter - try 2024 then 2023
    if "pop[pop['year']==2024]" in line:
        lines[i] = ("pop = pd.read_csv(os.path.join(DATA_DIR, 'populations_2024.csv'))\n"
                    "pop_df = pop[pop['year'].isin([2024,2023])]\n"
                    "if len(pop_df) == 0: pop_df = pop\n"
                    "pop = pop_df.groupby('country')['population']"
                    ".last().values\n"
                    "pop = pop[pop > 0]\n")
        print(f"Fixed populations filter at line {i+1}")
        break

with open(path, "w", encoding="utf-8") as f:
    f.writelines(lines)
print("Done.")
