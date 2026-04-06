with open("latticefit_app.py", encoding="utf-8") as f:
    src = f.read()

# Update the note to explain the discrepancy clearly
old = '''    anchor_note = (
        f"**Note:** Fixed anchor = {anchor_val:.4g} (electron mass preset). "
        "Null test uses same anchor."
        if not auto_mode and anchor_choice == "Custom"
        else "**Note:** Log-uniform null with free anchor. "
             "For fixed-anchor null, set Anchor=Custom and choose a preset."
    )
    st.caption(anchor_note)'''

new = '''    if not auto_mode and anchor_choice == "Custom":
        anchor_note = (
            f"**Fixed anchor:** {anchor_val:.4g}. "
            "Null test draws log-uniform values in the same lattice range. "
            "**Note:** The published SM mass result (z=+4.50) uses "
            "`sector_anchor_null` from the Python API, which conditions on "
            "the electron mass anchor being fixed by theory (stricter test). "
            "The log-uniform null shown here is more conservative and gives "
            "z~+1.5 for n=9 — both are valid depending on prior knowledge."
        )
    else:
        anchor_note = (
            "**Note:** Log-uniform null with free anchor. "
            "For the published SM mass result (z=+4.50), use the Python API: "
            "`from latticefit.stats import sector_anchor_null`"
        )
    st.caption(anchor_note)'''

src = src.replace(old, new)

with open("latticefit_app.py", "w", encoding="utf-8") as f:
    f.write(src)
print("Fixed.")
