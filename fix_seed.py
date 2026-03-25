path = r"C:\dev\latticefit\examples\scripts\run_all_examples.py"
with open(path, encoding="utf-8") as f:
    s = f.read()
s = s.replace("n_trials=n_null, random_seed=42", "n_trials=n_null, seed=42")
with open(path, "w", encoding="utf-8") as f:
    f.write(s)
print("Fixed:", "random_seed" not in s)
