with open("latticefit_app.py", encoding="utf-8") as f:
    src = f.read()

# Fix run_latticefit_full to accept anchor parameter
old = "def run_latticefit_full(vals, n_null=5000):"
new = "def run_latticefit_full(vals, n_null=5000, fixed_anchor=None):"

old2 = """    for bname, base in BASES.items():
        for d in [2, 3, 4, 6, 8]:
            anchor = vals.min()"""
new2 = """    for bname, base in BASES.items():
        for d in [2, 3, 4, 6, 8]:
            anchor = fixed_anchor if fixed_anchor is not None else vals.min()"""

src = src.replace(old, new).replace(old2, new2)

# Pass anchor_val into run_latticefit_full call
old3 = "        all_results, best_key = run_latticefit_full(x, n_null=n_null)"
new3 = ("        fixed = anchor_val if not auto_mode and "
        "anchor_choice == 'Custom' else None\n"
        "        all_results, best_key = run_latticefit_full("
        "x, n_null=n_null, fixed_anchor=fixed)")

src = src.replace(old3, new3)

# Also fix the note to say "fixed anchor" when one is set
old4 = ('    st.caption(\n'
        '        "**Note:** Log-uniform null with free anchor. "\n'
        '        "For fixed-anchor or sector-aware null, use the Python API directly."\n'
        '    )')
new4 = ('    anchor_note = (\n'
        '        f"**Note:** Fixed anchor = {anchor_val:.4g} (electron mass preset). "\n'
        '        "Null test uses same anchor."\n'
        '        if not auto_mode and anchor_choice == "Custom"\n'
        '        else "**Note:** Log-uniform null with free anchor. "\n'
        '             "For fixed-anchor null, set Anchor=Custom and choose a preset."\n'
        '    )\n'
        '    st.caption(anchor_note)')

src = src.replace(old4, new4)

with open("latticefit_app.py", "w", encoding="utf-8") as f:
    f.write(src)
print("Fixed.")
