with open("latticefit_app.py", encoding="utf-8") as f:
    src = f.read()

# Fix the note that incorrectly claimed z=+4.50
src = src.replace(
    "Note: The published SM mass result (z=+4.50) uses `sector_anchor_null`",
    "Note: The documented SM mass result (p=0.074, z≈+1.5) uses `log_uniform_null`"
)
src = src.replace(
    "which conditions on the electron mass anchor being fixed by theory (stricter test). "
    "The log-uniform null shown here is more conservative and gives z~+1.5 for n=9 "
    "— both are valid depending on prior knowledge.",
    "as shown in the JOSS paper (paper.md). "
    "The sector-anchor null (structure-preserving) gives p≈0.38 for n=9 fermions. "
    "Neither test is significant at p<0.05 for the 9-fermion dataset alone."
)

# Fix the sector null success message - was incorrectly claiming z=+4.50
src = src.replace(
    'f"— matches published result (z≈+4.50, p<0.001)."',
    'f"— matches documented result (p≈0.38 sector-anchor, p=0.074 log-uniform)."'
)

with open("latticefit_app.py", "w", encoding="utf-8") as f:
    f.write(src)
print("Fixed.")
