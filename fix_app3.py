with open("latticefit_app.py", encoding="utf-8") as f:
    src = f.read()

# Fix 1: SM masses demo should recommend fixed anchor + phi d=4
old = '''        "SM fermion masses (PDG 2024)": {
            "names": ["e","mu","tau","u","c","t","d","s","b"],
            "values": [5.11e-4,0.10566,1.77686,0.00216,1.275,172.76,
                       0.00467,0.0934,4.18],
            "unit": "GeV",
            "note": "9 fermion masses spanning 8 orders. phi-lattice, d=4.",
        },'''

new = '''        "SM fermion masses (PDG 2024)": {
            "names": ["e","mu","tau","u","c","t","d","s","b"],
            "values": [5.11e-4,0.10566,1.77686,0.00216,1.275,172.76,
                       0.00467,0.0934,4.18],
            "unit": "GeV",
            "note": ("9 fermion masses spanning 8 orders. "
                     "To reproduce the published result, set: "
                     "Base=phi, d=4, Anchor=Minimum (electron mass). "
                     "Auto-discovery with free anchor may find spurious "
                     "near-equivalent lattices at n=9."),
        },'''

# Fix 2: Add domain context to AI interpretation
old2 = '''def ai_interpret(result, dataset_info, orders, validity_issues,
                 chat_history=None):
    system = """You are a scientific assistant specialising in multiplicative 
scaling and geometric lattice structure. Interpret LatticeFit results clearly:
- What the best base means physically
- Whether p-value indicates genuine structure
- Possible mechanisms (evolutionary, quantum, geometric)
- Cross-domain context
- Validity caveats if flagged
Be specific, accurate, and cautious. 3-5 sentences unless asked for more."""'''

new2 = '''def ai_interpret(result, dataset_info, orders, validity_issues,
                 chat_history=None):
    system = """You are a scientific assistant specialising in multiplicative 
scaling and geometric lattice structure. Interpret LatticeFit results clearly:
- What the best base means physically
- Whether p-value indicates genuine structure
- Possible mechanisms (evolutionary, quantum, geometric)
- Cross-domain context from validated datasets:
  GENUINE phi-signals: SM fermion masses (phi d=4 with FIXED electron mass anchor,
  not free anchor), rice allele frequencies (z=+4.25), cetacean body masses (z=+2.43),
  single EGFR SAR series (z=+2.12)
  NULL results: S&P 500 returns, crystal volumes, earthquakes (Gutenberg-Richter law),
  NIST ionization energies
  ARTIFACTS: nuclear BE/A (narrow range 0.97 orders), binned patent IC50 data
- CRITICAL: with n<20 and free anchor, auto-discovery routinely finds spurious
  near-equivalent lattices. The SM fermion mass result requires fixed anchor=electron
  mass and phi d=4. Free-anchor results on small n datasets are unreliable.
- Validity caveats if flagged
Be specific, accurate, and cautious. 3-5 sentences unless asked for more."""'''

src = src.replace(old, new).replace(old2, new2)

with open("latticefit_app.py", "w", encoding="utf-8") as f:
    f.write(src)
print("Fixed.")
