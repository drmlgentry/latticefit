"""
LatticeFit — Interactive Web Demo
==================================
Detect discrete multiplicative structure in your data.

Deploy:
    streamlit run app.py

Or: https://streamlit.io/cloud (free hosting)
"""

import streamlit as st
import numpy as np
import pandas as pd
import io
import json
import sys
import os

# Add latticefit to path if running from repo
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import latticefit
from latticefit import fit, discover, PHI, KNOWN
from latticefit.stats import log_uniform_null, sector_anchor_null

# ── Page config ────────────────────────────────────────────────────
def _is_float(s):
    try: float(s); return True
    except: return False

st.set_page_config(
    page_title="LatticeFit",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── CSS ────────────────────────────────────────────────────────────
st.markdown("""
<style>
  .metric-box {
    background: #f0f4ff;
    border-left: 4px solid #4a90d9;
    border-radius: 6px;
    padding: 12px 16px;
    margin: 6px 0;
  }
  .rms-good  { color: #1a9850; font-weight: bold; font-size: 1.3em }
  .rms-ok    { color: #d68910; font-weight: bold; font-size: 1.3em }
  .rms-poor  { color: #c0392b; font-weight: bold; font-size: 1.3em }
  .sig-yes   { color: #1a9850 }
  .sig-no    { color: #888 }
  footer { visibility: hidden }
</style>
""", unsafe_allow_html=True)

# ── Built-in datasets ──────────────────────────────────────────────
DEMOS = {
    "Standard Model masses (physics)": {
        "names":  ["e","mu","tau","u","c","t","d","s","b","W","Z","H"],
        "values": [5.11e-4,0.1057,1.777,0.00216,1.275,172.76,
                   0.00467,0.0934,4.18,80.38,91.19,125.25],
        "anchor": 5.11e-4, "base": PHI, "denom": 4,
        "unit": "GeV",
        "note": "Fermion and boson masses from PDG 2024."
    },
    "Mammal body masses (biology)": {
        "names":  ["shrew","mouse","rat","rabbit","cat","dog",
                   "human","lion","horse","elephant","blue_whale"],
        "values": [0.003,0.02,0.2,2.0,4.0,30.0,
                   70.0,180.0,500.0,5000.0,150000.0],
        "anchor": None, "base": np.e, "denom": 6,
        "unit": "kg",
        "note": "Body masses spanning 8 orders of magnitude (Calder 1984)."
    },
    "Musical notes — equal temperament (acoustics)": {
        "names":  ["A4","A#4","B4","C5","C#5","D5","D#5","E5",
                   "F5","F#5","G5","G#5","A5"],
        "values": [440.0 * 2**(n/12) for n in range(13)],
        "anchor": 440.0, "base": 2**(1/12), "denom": 1,
        "unit": "Hz",
        "note": "Exact geometric lattice — demonstrates perfect recovery (RMS ≈ 0)."
    },
    "S&P 500 sector market caps (finance)": {
        "names":  ["Energy","Materials","Utilities","Real_Estate",
                   "Consumer_Staples","Industrials","Healthcare",
                   "Consumer_Disc","Comm_Services","Financials","Technology"],
        "values": [800,600,900,800,1300,2200,2800,3500,4000,4500,14000],
        "anchor": None, "base": PHI, "denom": 4,
        "unit": "$B",
        "note": "Approximate 2024 S&P 500 sector market capitalisations."
    },
}

# ── Sidebar ────────────────────────────────────────────────────────
with st.sidebar:
    st.image("https://raw.githubusercontent.com/drmlgentry/latticefit/main/fig_mass_lattice.png",
             use_container_width=True) if False else None
    st.title("🔬 LatticeFit")
    st.caption("Discrete lattice fitting with statistical validation")
    st.markdown("---")

    data_source = st.radio("Data source",
                           ["Built-in demo", "Upload CSV", "Paste data"])

    names, values = [], []

    if data_source == "Built-in demo":
        demo_key = st.selectbox("Select dataset", list(DEMOS.keys()))
        demo = DEMOS[demo_key]
        names  = demo["names"]
        values = demo["values"]
        st.info(demo["note"])
        unit = demo.get("unit", "")

    elif data_source == "Upload CSV":
        uploaded = st.file_uploader(
            "CSV file — any format (multi-column, USGS, custom)",
            type=["csv","txt"])
        unit = st.text_input("Unit label (optional)", "")

        if uploaded:
            raw = uploaded.read().decode("utf-8", errors="replace")
            lines = raw.splitlines()

            # Detect if multi-column (has header)
            first = lines[0].strip() if lines else ""
            is_multicolumn = "," in first and not _is_float(first.split(",")[0])

            if is_multicolumn:
                import io as _io
                df_preview = pd.read_csv(_io.StringIO(raw), nrows=3)
                st.markdown("**Multi-column CSV detected.** Select columns:")
                all_cols = list(df_preview.columns)

                col_name = st.selectbox("Name column (optional)",
                                        ["(none)"] + all_cols)
                col_value = st.selectbox("Value column", all_cols,
                    index=min(1, len(all_cols)-1))

                # Transformation
                transform = st.selectbox(
                    "Transform value",
                    ["None", "10^(1.5x) — earthquake energy from magnitude",
                     "log10(x)", "abs(x)", "x^2"])

                # Filter
                filter_positive = st.checkbox("Exclude zero/negative values", value=True)
                max_rows = st.slider("Max rows to use", 10, 1000, 200)

                # Load full data
                df_full = pd.read_csv(_io.StringIO(raw))
                col_data = pd.to_numeric(df_full[col_value], errors="coerce")
                df_full = df_full[col_data.notna()].copy()
                col_data = col_data[col_data.notna()]

                if filter_positive:
                    mask = col_data > 0
                    df_full = df_full[mask]
                    col_data = col_data[mask]

                col_data = col_data.iloc[:max_rows]
                df_full  = df_full.iloc[:max_rows]

                # Apply transform
                if "10^(1.5x)" in transform:
                    col_data = 10 ** (1.5 * col_data)
                    unit = unit or "energy"
                elif "log10" in transform:
                    col_data = np.log10(col_data.clip(lower=1e-30))
                    unit = unit or "log10"
                elif "x^2" in transform:
                    col_data = col_data ** 2

                values = list(col_data.values)
                if col_name != "(none)":
                    names = [str(v)[:20] for v in df_full[col_name].values]
                else:
                    names = [f"[{i}]" for i in range(len(values))]

                st.success(f"Loaded {len(values)} values from column '{col_value}'")
                with st.expander("Preview"):
                    preview_df = pd.DataFrame({"name": names[:10],
                                               "value": values[:10]})
                    st.dataframe(preview_df, hide_index=True)
            else:
                # Simple format
                for line in lines:
                    line = line.strip()
                    if not line or line.startswith("#"): continue
                    parts = [p.strip() for p in line.split(",")]
                    try:
                        if len(parts) >= 2 and not _is_float(parts[0]):
                            values.append(float(parts[-1]))
                            names.append(",".join(parts[:-1]))
                        elif len(parts) >= 2:
                            values.append(float(parts[-1]))
                            names.append(parts[0])
                        else:
                            values.append(float(parts[0]))
                    except ValueError:
                        pass
                if values:
                    st.success(f"Loaded {len(values)} values")

    else:  # Paste data
        unit = st.text_input("Unit label (optional)", "")
        pasted = st.text_area("Paste data (name,value or value per line)",
                              height=160,
                              placeholder="electron,0.000511\nmuon,0.1057\ntau,1.777")
        for line in pasted.splitlines():
            line = line.strip()
            if not line or line.startswith("#"): continue
            parts = [p.strip() for p in line.split(",")]
            try:
                if len(parts) >= 2:
                    values.append(float(parts[-1]))
                    names.append(",".join(parts[:-1]))
                else:
                    values.append(float(parts[0]))
            except ValueError:
                pass
        if values:
            st.success(f"Parsed {len(values)} values")

    st.markdown("---")
    st.subheader("Lattice parameters")

    auto_mode = st.checkbox("Auto-discover best lattice", value=False)

    if not auto_mode:
        base_choice = st.selectbox("Base r",
            ["φ (golden ratio)", "e (Euler)", "2", "√2", "10", "Custom"])
        base_map = {
            "φ (golden ratio)": PHI,
            "e (Euler)": np.e,
            "2": 2.0,
            "√2": 2**0.5,
            "10": 10.0,
        }
        if base_choice == "Custom":
            base_val = st.number_input("Custom base", value=PHI,
                                       min_value=1.001, format="%.6f")
        else:
            base_val = base_map[base_choice]

        denom = st.slider("Denominator d", 1, 8, 4)

        anchor_choice = st.selectbox("Anchor",
            ["First value", "Geometric mean", "Minimum", "Custom"])

    null_n = st.slider("Null test trials", 1000, 50000, 10000, step=1000)

    run_btn = st.button("▶  Run LatticeFit", type="primary",
                        use_container_width=True)

    st.markdown("---")
    st.caption("v0.1.1 · [GitHub](https://github.com/drmlgentry/latticefit) · "
               "Patent pending US 64/013,306")

# ── Main panel ─────────────────────────────────────────────────────
st.title("LatticeFit — Discrete Scaling Analysis")
st.markdown(
    "Detect whether your data clusters near a geometric lattice "
    r"$x \approx A \cdot r^{k/d}$ and quantify statistical significance."
)

if not values:
    st.info("👈 Select a demo dataset or upload your own data to get started.")

    # Show example output
    with st.expander("What does LatticeFit do?"):
        st.markdown("""
**LatticeFit** tests whether a dataset of positive measurements clusters near
a discrete geometric lattice:

$$x_i \\approx A \\cdot r^{k_i/d}, \\quad k_i \\in \\mathbb{Z}$$

Given anchor $A$, base $r$, and denominator $d$, it:

1. Assigns the nearest integer label $k_i$ to each observation
2. Computes logarithmic residuals $\\delta_i = |\\log_r(x_i/A) - k_i/d|$
3. Tests whether the RMS residual is smaller than random chance

**Example**: Standard Model fermion masses fit the $\\varphi$-lattice 
($r=\\varphi$, $d=4$) with RMS = 0.069 (out of max 0.125).

**Validated on**: particle physics, mammal biology, acoustics, finance.
        """)
    st.stop()

x = np.array(values, dtype=float)

if np.any(x <= 0):
    st.error("All values must be strictly positive.")
    st.stop()

if len(x) < 2:
    st.error("Need at least 2 data points.")
    st.stop()

if not run_btn and data_source != "Built-in demo":
    st.info("Configure parameters in the sidebar and click **Run LatticeFit**.")
    st.stop()

# ── Resolve anchor ────────────────────────────────────────────────
if not auto_mode:
    if anchor_choice == "First value":
        anchor_val = float(x[0])
    elif anchor_choice == "Geometric mean":
        anchor_val = float(np.exp(np.mean(np.log(x))))
    elif anchor_choice == "Minimum":
        anchor_val = float(x.min())
    else:
        anchor_val = st.number_input("Custom anchor", value=float(x[0]),
                                     min_value=1e-30, format="%.6g")

# ── Run fit ───────────────────────────────────────────────────────
with st.spinner("Fitting lattice…"):
    if auto_mode:
        results = discover(x, names=names if names else None, top_n=5)
        result = results[0]
        st.success(f"Auto-discovered: base = {result.base:.5g}, "
                   f"d = {result.denom}, anchor = {result.anchor:.4g}")
    else:
        result = fit(x, anchor=anchor_val, base=base_val,
                     denom=denom, names=names if names else None)

# ── Run null tests ────────────────────────────────────────────────
with st.spinner(f"Running null tests ({null_n:,} trials)…"):
    null = log_uniform_null(result, n_trials=null_n)

# ── Layout ────────────────────────────────────────────────────────
col1, col2, col3, col4 = st.columns(4)

max_rms = 0.5 / result.denom
frac = result.rms / max_rms
rms_class = "rms-good" if frac < 0.5 else "rms-ok" if frac < 0.75 else "rms-poor"
p_class = "sig-yes" if null.p_value < 0.05 else "sig-no"

with col1:
    st.metric("RMS residual", f"{result.rms:.4f}",
              delta=f"{frac*100:.0f}% of max")
with col2:
    st.metric("Max possible RMS", f"{max_rms:.4f}")
with col3:
    st.metric("p-value (log-uniform)", f"{null.p_value:.3f}",
              delta="significant" if null.p_value < 0.05 else "marginal")
with col4:
    st.metric("N observations", len(x))

st.markdown("---")

# ── Results table ─────────────────────────────────────────────────
tab1, tab2, tab3 = st.tabs(["📊 Results table", "📈 Plot", "📋 Export"])

with tab1:
    display_names = result.names or [str(i) for i in range(len(x))]
    df = pd.DataFrame({
        "Name":      display_names,
        f"Observed ({unit})": [f"{v:.5g}" for v in result.data],
        "k":         result.labels,
        f"Predicted ({unit})": [f"{v:.5g}" for v in result.predicted],
        "Residual δ": [f"{v:.4f}" for v in result.residuals],
        "|Δ|/m (%)": [f"{abs(o-p)/o*100:.1f}%"
                       for o, p in zip(result.data, result.predicted)],
    })
    st.dataframe(df, use_container_width=True, hide_index=True)

    # Residual bar chart
    st.markdown("**Residuals by observation** (lower = better lattice alignment)")
    res_df = pd.DataFrame({
        "Name": display_names,
        "Residual": result.residuals,
        "Max possible": [max_rms] * len(x),
    }).set_index("Name")
    st.bar_chart(res_df[["Residual"]])

with tab2:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 5))

        # Lattice lines
        k_lo, k_hi = result.labels.min()-2, result.labels.max()+2
        for ki in range(k_lo, k_hi+1):
            yl = result.anchor * result.base**(ki/result.denom)
            ax.axhline(yl, color="#eeeeee", linewidth=0.5, zorder=0)

        idx = np.arange(len(x))
        ax.scatter(idx, result.data, color="#2166ac", s=80, zorder=5,
                   edgecolors="black", linewidths=0.5, label="Observed")
        ax.scatter(idx, result.predicted, color="#2166ac", s=80, zorder=4,
                   facecolors="none", edgecolors="#2166ac",
                   linewidths=1.8, label="Lattice prediction")
        for i, (o, p) in enumerate(zip(result.data, result.predicted)):
            ax.plot([i,i], [min(o,p), max(o,p)],
                    color="#2166ac", alpha=0.4, linewidth=1)

        ax.set_yscale("log")
        ax.set_xticks(idx)
        ax.set_xticklabels(display_names, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel(f"Value ({unit})" if unit else "Value", fontsize=11)
        ax.set_title(
            f"LatticeFit  r={result.base:.5g}  d={result.denom}  "
            f"A={result.anchor:.4g}\nRMS={result.rms:.4f}  "
            f"p={null.p_value:.3f}", fontsize=11)
        ax.legend(fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.tight_layout()

        buf = io.BytesIO()
        plt.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        st.image(buf, use_container_width=True)
        plt.close()

        # Download plot
        st.download_button("⬇ Download plot (PNG)", buf.getvalue(),
                           "latticefit_plot.png", "image/png")
    except Exception as e:
        st.warning(f"Plot unavailable: {e}")

with tab3:
    # CSV download
    csv_lines = ["name,observed,k,predicted,residual,pct_error"]
    for n, o, k, p, r in zip(display_names, result.data,
                               result.labels, result.predicted,
                               result.residuals):
        csv_lines.append(f"{n},{o:.6g},{k},{p:.6g},{r:.5f},"
                         f"{abs(o-p)/o*100:.2f}")
    csv_str = "\n".join(csv_lines)
    st.download_button("⬇ Download results CSV", csv_str,
                       "latticefit_results.csv", "text/csv")

    # JSON download
    json_out = {
        "latticefit_version": "0.1.1",
        "parameters": {"anchor": result.anchor, "base": result.base,
                       "denom": result.denom},
        "fit": {"rms": result.rms, "n": len(result.data),
                "max_possible_rms": max_rms},
        "null_test": {"test": null.test_name, "p_value": null.p_value,
                      "null_mean": null.null_mean, "n_trials": null.n_trials},
        "observations": [
            {"name": n, "observed": float(o), "k": int(k),
             "predicted": float(p), "residual": float(r)}
            for n, o, k, p, r in zip(display_names, result.data,
                                      result.labels, result.predicted,
                                      result.residuals)
        ],
    }
    st.download_button("⬇ Download JSON", json.dumps(json_out, indent=2),
                       "latticefit_result.json", "application/json")

    # Citation
    st.markdown("---")
    st.markdown("**Cite LatticeFit:**")
    st.code("""Gentry, M. L. (2026). LatticeFit: Discrete lattice fitting
with statistical validation (v0.1.1). GitHub.
https://github.com/drmlgentry/latticefit

@software{latticefit2026,
  author  = {Gentry, Marvin L.},
  title   = {{LatticeFit}: Discrete lattice fitting with statistical validation},
  year    = {2026},
  version = {0.1.1},
  url     = {https://github.com/drmlgentry/latticefit}
}""", language="bibtex")

# ── Auto-discover results ─────────────────────────────────────────
if auto_mode:
    st.markdown("---")
    st.subheader("Auto-discovery: top 5 lattices")
    auto_df = pd.DataFrame([
        {"Rank": i+1,
         "Base r": f"{r.base:.5g}",
         "Denom d": r.denom,
         "Anchor A": f"{r.anchor:.4g}",
         "RMS": f"{r.rms:.5f}",
         "% of max": f"{r.rms/(0.5/r.denom)*100:.1f}%"}
        for i, r in enumerate(results)
    ])
    st.dataframe(auto_df, use_container_width=True, hide_index=True)

# ── Null test detail ──────────────────────────────────────────────
st.markdown("---")
with st.expander("📊 Null test details"):
    c1, c2 = st.columns(2)
    with c1:
        st.metric("Observed RMS", f"{null.observed_rms:.5f}")
        st.metric("Null mean ± std",
                  f"{null.null_mean:.5f} ± {null.null_std:.5f}")
    with c2:
        z = (null.null_mean - null.observed_rms) / max(null.null_std, 1e-12)
        st.metric("Z-score", f"{z:.2f}")
        st.metric("p-value", f"{null.p_value:.4f}")

    if null.p_value < 0.05:
        st.success("Significant (p < 0.05): lattice alignment unlikely by chance.")
    elif null.p_value < 0.10:
        st.warning("Marginal (0.05 ≤ p < 0.10): weak evidence of lattice structure.")
    else:
        st.info("Not significant (p ≥ 0.10): lattice alignment consistent with chance.")

    st.caption(f"Log-uniform null, {null.n_trials:,} trials. "
               "p-value = fraction of random spectra achieving RMS ≤ observed.")
