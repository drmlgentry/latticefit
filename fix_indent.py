path = r"C:\dev\latticefit\examples\scripts\run_all_examples.py"
with open(path, encoding="utf-8") as f:
    lines = f.readlines()

# Fix: lines 117-118 are incorrectly indented inside else block
# Remove the 4-space indent from lines 117 and 118
for i, line in enumerate(lines):
    if "pop = pd.read_csv" in line and line.startswith("    "):
        lines[i] = line.lstrip("    ") # remove one level of indent
        print(f"Fixed indent line {i+1}: {lines[i].rstrip()}")
    if "pop = pop[pop" in line and line.startswith("    "):
        lines[i] = line.lstrip("    ")
        print(f"Fixed indent line {i+1}: {lines[i].rstrip()}")

with open(path, "w", encoding="utf-8") as f:
    f.writelines(lines)
print("Done.")
