with open("latticefit_app.py", encoding="utf-8") as f:
    src = f.read()

# Pass filename/demo_key into dataset_info so AI knows what it is
old = '''        st.session_state.update({
        'fit_results': all_results,
        'best_key': best_key,
        'vals': x,
        'dataset_info': {
            "filename": filename or data_source,
            "column": names[0] if len(names)==1 else "values",
            "n": len(x),
            "range": f"{x.min():.3g}--{x.max():.3g}",
            "decades": orders,
        },'''

new = '''        st.session_state.update({
        'fit_results': all_results,
        'best_key': best_key,
        'vals': x,
        'dataset_info': {
            "filename": (demo_key if data_source == "Built-in demo"
                         else filename or data_source),
            "column": names[0] if len(names)==1 else "values",
            "n": len(x),
            "range": f"{x.min():.3g}--{x.max():.3g}",
            "decades": orders,
        },'''

# Add electron mass preset to custom anchor
old2 = '''        if anchor_choice == "Custom":
                anchor_val = st.number_input("Custom anchor", value=float(x[0]),
                                             min_value=1e-30, format="%.6g")'''

new2 = '''        if anchor_choice == "Custom":
                preset = st.selectbox("Preset anchors",
                    ["(none)", "Electron mass (5.11e-4 GeV)",
                     "Electron mass (0.511 MeV)",
                     "Proton mass (938.3 MeV)"])
                preset_vals = {
                    "Electron mass (5.11e-4 GeV)": 5.11e-4,
                    "Electron mass (0.511 MeV)":   0.511,
                    "Proton mass (938.3 MeV)":     938.3,
                }
                default_anchor = preset_vals.get(preset, float(x[0]))
                anchor_val = st.number_input("Custom anchor",
                                             value=default_anchor,
                                             min_value=1e-30, format="%.6g")'''

src = src.replace(old, new).replace(old2, new2)

with open("latticefit_app.py", "w", encoding="utf-8") as f:
    f.write(src)
print("Fixed.")
