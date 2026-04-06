with open("latticefit_app.py", encoding="utf-8") as f:
    src = f.read()

# Fix neutrino demo note to clarify it is a prediction not for fitting
old = '''        "Neutrino masses (predicted)": {
            "names": ["nu1","nu2","nu3"],
            "values": [8.387e-3, 12.033e-3, 50.972e-3],
            "unit": "eV",
            "note": "phi-lattice prediction: q=(-149,-146,-134), sum=71.4 meV. "
                    "Testable by CMB-S4.",
        },'''

new = '''        "Neutrino masses (predicted)": {
            "names": ["nu_1 (8.4 meV)","nu_2 (12.0 meV)","nu_3 (51.0 meV)"],
            "values": [8.387e-3, 12.033e-3, 50.972e-3],
            "unit": "eV",
            "note": ("PREDICTION (not for fitting): phi-lattice with anchor=m_e "
                     "and q=(-149,-146,-134) gives sum=71.4 meV. "
                     "This dataset is too small (n=3, 0.78 orders) for LatticeFit. "
                     "It is shown to illustrate the prediction; use the Python API "
                     "with fixed anchor=m_e to reproduce the full derivation."),
        },'''

# Also fix ceiling value detector - should not flag the max of a 3-point dataset
old2 = '''    max_count = np.sum(vals == vals.max())
    if max_count > len(vals) * 0.1:
        warnings.append(f"Ceiling values: {max_count} observations at maximum "
                        f"({vals.max():.3g}) -- may be assay limit sentinels.")'''

new2 = '''    max_count = np.sum(vals == vals.max())
    if max_count > len(vals) * 0.1 and len(vals) >= 10:
        warnings.append(f"Ceiling values: {max_count} observations at maximum "
                        f"({vals.max():.3g}) -- may be assay limit sentinels.")'''

src = src.replace(old, new).replace(old2, new2)

with open("latticefit_app.py", "w", encoding="utf-8") as f:
    f.write(src)
print("Fixed.")
