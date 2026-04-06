"""
latticefit_app.py  v0.3.0
=========================
LatticeFit interactive demo — upgraded with:
  - Validity checker (range, binning, effect size)
  - Cross-domain survey results tab
  - Z-score + effect size display
  - AI assistant (column identification + interpretation + chat)
  - Multi-base comparison
  - Normal/inverted ordering neutrino prediction

Usage:
    streamlit run latticefit_app.py

Requirements:
    pip install streamlit pandas numpy matplotlib latticefit anthropic
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io, json, urllib.request

# ── Page config ────────────────────────────────────────────────────
st.set_page_config(
    page_title="LatticeFit",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

PHI = (1 + 5**0.5) / 2

# ── Imports ────────────────────────────────────────────────────────
try:
    from latticefit import fit, discover, PHI as LF_PHI, KNOWN
    from latticefit.stats import log_uniform_null, sector_anchor_null
    HAS_LATTICEFIT = True
except ImportError:
    HAS_LATTICEFIT = False

try:
    import anthropic
    HAS_ANTHROPIC = True
except ImportError:
    HAS_ANTHROPIC = False

# ── CSS ────────────────────────────────────────────────────────────
st.markdown("""
<style>
  .validity-ok   {background:#e8f5e9;border-left:4px solid #2e7d32;
                  border-radius:6px;padding:10px 14px;margin:4px 0;}
  .validity-warn {background:#fff8e1;border-left:4px solid #f57f17;
                  border-radius:6px;padding:10px 14px;margin:4px 0;}
  .validity-fail {background:#ffebee;border-left:4px solid #c62828;
                  border-radius:6px;padding:10px 14px;margin:4px 0;}
  .domain-genuine{background:#e8f5e9;border-radius:4px;
                  padding:2px 8px;color:#1b5e20;font-weight:600;}
  .domain-null   {background:#fafafa;border-radius:4px;
                  padding:2px 8px;color:#888;}
  .domain-artifact{background:#fff3e0;border-radius:4px;
                   padding:2px 8px;color:#e65100;}
  footer {visibility:hidden}
</style>
""", unsafe_allow_html=True)

# ── Constants ──────────────────────────────────────────────────────
BASES = {"φ (golden ratio)": PHI, "√2": 2**0.5, "2": 2.0,
         "e": np.e, "10": 10.0}

CROSS_DOMAIN = [
    # domain, n, rms, null_rms, z, p_str, orders, category, status
    ("SM fermion masses",          12,  0.0550, 0.0720, +4.50, "<0.001", 8.0, "Physics",    "GENUINE"),
    ("Rice 44K MAF",            36901,  0.0712, 0.0721, +4.25, "<0.001", 2.5, "Genomics",   "GENUINE"),
    ("Cetacea body masses",         76,  0.0631, 0.0720, +2.43, "0.008",  3.7, "Ecology",    "GENUINE"),
    ("EGFR SAR (CHEMBL1064829)",    32,  0.0600, 0.0721, +2.12, "0.021",  3.4, "Pharma",     "GENUINE"),
    ("EGFR SAR (CHEMBL3750872)",    40,  0.0600, 0.0721, +2.33, "0.014",  2.1, "Pharma",     "GENUINE"),
    ("COVID variant freq",      103348,  0.0712, 0.0721, +9.38, "<0.001", 3.0, "Epidemiol.", "NEEDS NULL"),
    ("All mammals",               3542,  0.0714, 0.0721, +1.45, "0.074",  7.9, "Ecology",    "MARGINAL"),
    ("NIST ionization energies",  1631,  0.0721, 0.0723, +0.30, "0.378",  4.4, "Physics",    "NULL"),
    ("S&P 500 returns",          11345,  0.0726, 0.0722, -1.43, "0.924",  4.6, "Finance",    "NULL"),
    ("Crystal volumes (COD)",   525224,  0.0721, 0.0721, -1.04, "0.850",  3.3, "Chemistry",  "NULL"),
    ("HIV RT IC50 (mixed)",      10041,  0.0733, 0.0722, -3.73, "0.9996", 9.1, "Pharma",     "NULL"),
    ("Earthquakes (GR law)",     18659,  0.0730, 0.0721, -3.58, "1.000",  6.8, "Geophys.",   "NULL"),
    ("Nuclear BE/A (AME2020)",    3553,  0.0681, 0.0720, +7.18, "<0.001", 1.0, "Physics",    "ARTIFACT"),
    ("PharmGKB Score",            3896,  0.0776, 0.0722,-10.58, "1.000",  0.0, "Pharma",     "ARTIFACT"),
]

# ── Validity checker ───────────────────────────────────────────────
def check_validity(vals):
    issues = []
    warnings = []
    passes = []

    # 1. Minimum range
    orders = np.log10(vals.max() / vals.min()) if vals.min() > 0 else 0
    if orders < 1.0:
        issues.append(f"Range too narrow: {orders:.2f} orders of magnitude (need ≥ 3). "
                      "Result will be artifactually significant for any base.")
    elif orders < 3.0:
        warnings.append(f"Range marginal: {orders:.2f} orders (recommend ≥ 3 for reliable test).")
    else:
        passes.append(f"Range: {orders:.1f} orders of magnitude ✓")

    # 2. Binned data check
    n_unique = len(np.unique(vals))
    frac_unique = n_unique / len(vals)
    if frac_unique < 0.3:
        issues.append(f"Data appears binned: only {n_unique} unique values "
                      f"out of {len(vals)} ({frac_unique:.0%}). "
                      "Bin boundaries may create spurious lattice structure.")
    elif frac_unique < 0.6:
        warnings.append(f"Possible binning: {n_unique}/{len(vals)} unique values ({frac_unique:.0%}).")
    else:
        passes.append(f"Unique values: {n_unique}/{len(vals)} ({frac_unique:.0%}) ✓")

    # 3. Minimum sample size
    if len(vals) < 20:
        warnings.append(f"Small sample: n={len(vals)}. Results may be unreliable (recommend n ≥ 30).")
    else:
        passes.append(f"Sample size: n={len(vals)} ✓")

    # 4. Ceiling values (e.g. IC50 assay limits)
    mode_val = np.bincount(
        np.searchsorted(np.unique(vals), vals)
    ).max()
    max_count = np.sum(vals == vals.max())
    if max_count > len(vals) * 0.1:
        warnings.append(f"Ceiling values: {max_count} observations at maximum "
                        f"({vals.max():.3g}) — may be assay limit sentinels.")

    return passes, warnings, issues, orders


# ── Core latticefit ────────────────────────────────────────────────
def run_latticefit_full(vals, n_null=5000, fixed_anchor=None):
    """Run latticefit across all standard bases, return dict of results."""
    rng = np.random.default_rng(42)
    results = {}

    for bname, base in BASES.items():
        for d in [2, 3, 4, 6, 8]:
            anchor = fixed_anchor if fixed_anchor is not None else vals.min()
            ks   = np.round(d * np.log(vals/anchor) / np.log(base))
            pred = anchor * base**(ks/d)
            res  = np.abs(np.log(vals/pred) / np.log(base))
            rms  = float(np.sqrt(np.mean(res**2)))
            max_rms = 0.5 / d

            # Null: draw log-uniform in same range AS SEEN BY THE LATTICE
            # i.e. uniform in log-space between min and max label
            k_min = np.floor(d * np.log(vals.min()/anchor) / np.log(base))
            k_max = np.ceil( d * np.log(vals.max()/anchor) / np.log(base))
            lo_k = k_min / d * np.log(base) + np.log(anchor)
            hi_k = k_max / d * np.log(base) + np.log(anchor)
            # Vectorised null using same anchor
            rand_mat = np.exp(rng.uniform(lo_k, hi_k, (n_null, len(vals))))
            ks_mat   = np.round(d * np.log(rand_mat/anchor) / np.log(base))
            pred_mat = anchor * base**(ks_mat/d)
            res_mat  = np.abs(np.log(rand_mat/pred_mat) / np.log(base))
            null_rms = np.sqrt(np.mean(res_mat**2, axis=1))

            p_val  = float(np.mean(null_rms <= rms))
            z_val  = float((null_rms.mean() - rms) / null_rms.std())
            effect = float((null_rms.mean() - rms) / null_rms.mean() * 100)

            key = (bname, d)
            results[key] = {
                "base": base, "base_name": bname, "denom": d,
                "rms": rms, "max_rms": max_rms,
                "null_mean": float(null_rms.mean()),
                "null_std":  float(null_rms.std()),
                "p_value": p_val, "z": z_val,
                "effect_pct": effect,
                "rms_pct": rms / max_rms * 100,
                "n": len(vals), "null_arr": null_rms,
            }

    # Find best by z-score
    best_key = max(results, key=lambda k: results[k]["z"])
    return results, best_key


# ── AI helpers ────────────────────────────────────────────────────
def call_claude(system_prompt, user_message, max_tokens=800):
    if not HAS_ANTHROPIC:
        return "Install anthropic package: pip install anthropic"
    import os
    api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    if not api_key:
        return ("AI assistant requires ANTHROPIC_API_KEY environment variable. "
                "Set it with: $env:ANTHROPIC_API_KEY = 'YOUR_API_KEY_HERE'")
    client = anthropic.Anthropic(api_key=api_key)
    msg = client.messages.create(
        model="claude-sonnet-4-5",
        max_tokens=max_tokens,
        system=system_prompt,
        messages=[{"role": "user", "content": user_message}],
    )
    return msg.content[0].text


def ai_parse_dataset(df_preview, filename):
    system = """You are a data scientist for LatticeFit.
Identify the best column for lattice fitting (positive, ≥1 order of magnitude range).
Flag binned data, currency symbols, European decimals, ceiling values, ID columns needing dedup.
Respond ONLY with JSON:
{"best_column":"name","preprocessing":[],"dedup_column":null,"dedup_reason":"",
 "reason":"one sentence","issues":"","alternative_columns":[]}"""
    col_info = []
    for col in df_preview.columns:
        try:
            vals = pd.to_numeric(
                df_preview[col].astype(str)
                .str.replace(',', '.').str.replace(r'[^0-9.\-]', '', regex=True),
                errors='coerce').dropna()
            pos = vals[vals > 0]
            if len(pos) > 3:
                decades = np.log10(pos.max()/pos.min()) if pos.min() > 0 else 0
                col_info.append(f"  {col}: n={len(pos)}, "
                                f"range={pos.min():.3g}–{pos.max():.3g}, "
                                f"decades={decades:.1f}")
            else:
                col_info.append(f"  {col}: non-numeric or <3 positive values")
        except:
            col_info.append(f"  {col}: parse error")
    user_msg = (f"File: {filename}\nColumns:\n{chr(10).join(col_info)}\n\n"
                f"First 3 rows:\n{df_preview.head(3).to_string()}")
    return call_claude(system, user_msg, max_tokens=300)


def ai_interpret(result, dataset_info, orders, validity_issues,
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
Be specific, accurate, and cautious. 3-5 sentences unless asked for more."""
    context = ""
    if chat_history:
        context = "\n\nPrevious conversation:\n"
        for msg in chat_history[-4:]:
            context += f"{msg['role'].title()}: {msg['content']}\n"
    base_names = {PHI: "φ (golden ratio)", 2**0.5: "√2", 2.0: "2",
                  np.e: "e", 10.0: "10"}
    b = result.get('base', 0)
    bn = min(base_names, key=lambda x: abs(x-b))
    bname = base_names[bn] if abs(bn-b) < 0.01 else f"{b:.5g}"
    issues_str = "; ".join(validity_issues) if validity_issues else "none"
    user_msg = (f"Dataset: {dataset_info.get('filename','?')}\n"
                f"Column: {dataset_info.get('column','?')}\n"
                f"n={result.get('n','?')}, {orders:.1f} decades\n"
                f"Best base: {bname} (d={result.get('denom','?')})\n"
                f"RMS={result.get('rms',0):.4f} ({result.get('rms_pct',0):.0f}% of max)\n"
                f"Effect size: {result.get('effect_pct',0):.1f}% RMS reduction\n"
                f"p={result.get('p_value',1):.4f}, z={result.get('z',0):+.2f}\n"
                f"Validity issues: {issues_str}\n{context}"
                f"\nInterpret these results.")
    return call_claude(system, user_msg, max_tokens=600)


def clean_column(series):
    s = series.astype(str)
    s = s.str.replace(r'[$₹£€¥]', '', regex=True)
    if s.str.contains(',').any() and not s.str.contains(r'\d\.\d').any():
        s = s.str.replace(',', '.')
    else:
        s = s.str.replace(r',(?=\d{3})', '', regex=True)
    return pd.to_numeric(s.str.strip(), errors='coerce')


# ── Session state ─────────────────────────────────────────────────
for key in ['fit_results', 'best_key', 'dataset_info', 'chat_history',
            'vals', 'validity', 'ai_suggestion']:
    if key not in st.session_state:
        st.session_state[key] = ([] if key == 'chat_history' else None)

# ── Sidebar ────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🔬 LatticeFit")
    st.caption("v0.3.0 · Discrete scaling analysis")
    st.markdown("---")

    n_null = st.slider("Null test trials", 1000, 20000, 5000, 1000,
                       help="More trials = more accurate p-values but slower")
    use_ai = st.toggle("AI assistant", value=True,
                       disabled=not HAS_ANTHROPIC,
                       help="Requires ANTHROPIC_API_KEY env variable")
    st.markdown("---")
    st.markdown(
        "📖 [GitHub](https://github.com/drmlgentry/latticefit)  \n"
        "📄 DOI: [10.5281/zenodo.19225731](https://doi.org/10.5281/zenodo.19225731)  \n"
        "🔒 Patent pending US 64/013,306"
    )

# ── Main tabs ─────────────────────────────────────────────────────
tab_analysis, tab_crossdomain, tab_validity_guide, tab_about = st.tabs([
    "📊 Analysis", "🌐 Cross-Domain Survey", "✅ Validity Guide", "ℹ️ About"
])

# ══════════════════════════════════════════════════════════════════
# TAB 1: ANALYSIS
# ══════════════════════════════════════════════════════════════════
with tab_analysis:
    st.header("Lattice Fitting Analysis")

    # ── Data input ────────────────────────────────────────────────
    c1, c2 = st.columns([3, 2])
    with c1:
        data_source = st.radio("Data source",
            ["Built-in demo", "Upload file", "Paste values"],
            horizontal=True)

    df = None
    filename = ""
    values, names = [], []

    # Built-in demos
    DEMOS = {
        "SM fermion masses (PDG 2024)": {
            "names": ["e","μ","τ","u","c","t","d","s","b"],
            "values": [5.11e-4,0.10566,1.77686,0.00216,1.275,172.76,
                       0.00467,0.0934,4.18],
            "unit": "GeV",
            "note": "9 fermion masses spanning 8 orders. φ-lattice, d=4.",
        },
        "Equal-tempered notes (acoustics)": {
            "names": [f"n{i}" for i in range(13)],
            "values": [440*2**(n/12) for n in range(13)],
            "unit": "Hz",
            "note": "Exact base-2^(1/12) lattice. Demonstrates perfect recovery.",
        },
        "Mammal body masses (ecology)": {
            "names": ["shrew","mouse","rat","rabbit","cat","dog",
                      "human","lion","horse","elephant","blue_whale"],
            "values": [0.003,0.02,0.2,2.0,4.0,30.0,70.0,
                       180.0,500.0,5000.0,150000.0],
            "unit": "kg",
            "note": "Body masses spanning 8 orders of magnitude.",
        },
        "Neutrino masses (predicted)": {
            "names": ["ν₁","ν₂","ν₃"],
            "values": [8.387e-3, 12.033e-3, 50.972e-3],
            "unit": "eV",
            "note": "φ-lattice prediction: q=(-149,-146,-134), sum=71.4 meV. "
                    "Testable by CMB-S4.",
        },
    }

    if data_source == "Built-in demo":
        demo_key = st.selectbox("Select dataset", list(DEMOS.keys()))
        demo = DEMOS[demo_key]
        names  = demo["names"]
        values = demo["values"]
        unit   = demo.get("unit", "")
        st.info(demo["note"])

    elif data_source == "Upload file":
        uploaded = st.file_uploader(
            "CSV, Excel, TSV",
            type=["csv","xlsx","xls","tsv"])
        unit = st.text_input("Unit (optional)", "")

        if uploaded:
            filename = uploaded.name
            try:
                if filename.endswith(('.xlsx','.xls')):
                    df = pd.read_excel(uploaded)
                elif filename.endswith('.tsv'):
                    df = pd.read_csv(uploaded, sep='\t', on_bad_lines='skip')
                else:
                    df = pd.read_csv(uploaded, on_bad_lines='skip')
                st.success(f"Loaded {df.shape[0]:,} rows × {df.shape[1]} columns")
                st.dataframe(df.head(4), use_container_width=True)
            except Exception as e:
                st.error(f"Could not read file: {e}")

            if df is not None:
                col_l, col_r = st.columns(2)
                with col_l:
                    if use_ai and HAS_ANTHROPIC:
                        if st.button("🤖 AI: identify best column"):
                            with st.spinner("Asking AI..."):
                                resp = ai_parse_dataset(df, filename)
                                try:
                                    st.session_state['ai_suggestion'] = json.loads(resp)
                                except:
                                    st.session_state['ai_suggestion'] = {
                                        "best_column": df.columns[0],
                                        "reason": resp,
                                        "preprocessing": [],
                                        "issues": "",
                                    }
                    if st.session_state['ai_suggestion']:
                        s = st.session_state['ai_suggestion']
                        st.info(f"**AI suggests:** `{s.get('best_column','?')}`  \n"
                                f"{s.get('reason','')}  \n"
                                f"Preprocessing: "
                                f"{', '.join(s.get('preprocessing',[])) or 'none'}")

                with col_r:
                    sugg_col = (st.session_state.get('ai_suggestion') or {}).get(
                        'best_column', df.columns[0])
                    idx = df.columns.tolist().index(sugg_col) \
                        if sugg_col in df.columns else 0
                    selected_col = st.selectbox("Column to analyze",
                                                df.columns.tolist(), index=idx)
                    dedup_col = st.selectbox("Deduplicate by (optional)",
                                             ["None"] + df.columns.tolist())

                if dedup_col != "None":
                    df = df.groupby(dedup_col)[selected_col].max().reset_index()
                raw = clean_column(df[selected_col]).dropna().values
                values = list(raw[raw > 0])
                names  = [str(i) for i in range(len(values))]

    else:  # Paste
        unit = st.text_input("Unit (optional)", "")
        pasted = st.text_area(
            "Paste data (name,value or value per line)",
            height=140,
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

    if not values:
        st.info("👈 Choose a demo dataset or upload/paste data to begin.")
        st.stop()

    x = np.array(values, dtype=float)
    x = x[x > 0]

    if len(x) < 3:
        st.error("Need at least 3 positive values.")
        st.stop()

    # ── Validity check ────────────────────────────────────────────
    passes, warnings_v, issues_v, orders = check_validity(x)

    with st.expander("✅ Validity check", expanded=bool(issues_v or warnings_v)):
        for p in passes:
            st.markdown(f'<div class="validity-ok">✓ {p}</div>',
                        unsafe_allow_html=True)
        for w in warnings_v:
            st.markdown(f'<div class="validity-warn">⚠️ {w}</div>',
                        unsafe_allow_html=True)
        for iss in issues_v:
            st.markdown(f'<div class="validity-fail">✗ {iss}</div>',
                        unsafe_allow_html=True)

    if issues_v:
        st.warning("This dataset has validity issues. Results may be unreliable — "
                   "see the Validity Guide tab for details.")

    # ── Fit parameters ────────────────────────────────────────────
    st.markdown("---")
    fc1, fc2, fc3 = st.columns(3)
    with fc1:
        auto_mode = st.checkbox("Auto-discover best base", value=False)
    with fc2:
        if not auto_mode:
            base_choice = st.selectbox("Base r", list(BASES.keys()))
            base_val = BASES[base_choice]
            denom = st.slider("Denominator d", 1, 8, 4)
        else:
            base_val, denom = PHI, 4
    with fc3:
        if not auto_mode:
            anchor_choice = st.selectbox(
                "Anchor A", ["Minimum", "First value",
                             "Geometric mean", "Custom"])
            if anchor_choice == "Custom":
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
                                             min_value=1e-30, format="%.6g")
            elif anchor_choice == "First value":
                anchor_val = float(x[0])
            elif anchor_choice == "Geometric mean":
                anchor_val = float(np.exp(np.mean(np.log(x))))
            else:
                anchor_val = float(x.min())
        else:
            anchor_val = float(x.min())

    run_btn = st.button("▶ Run LatticeFit", type="primary",
                        use_container_width=True)

    if not run_btn:
        st.stop()

    # ── Run fit ───────────────────────────────────────────────────
    with st.spinner(f"Fitting lattice to {len(x):,} values..."):
        fixed = anchor_val if not auto_mode and anchor_choice == 'Custom' else None
        all_results, best_key = run_latticefit_full(x, n_null=n_null, fixed_anchor=fixed)

    if auto_mode:
        r = all_results[best_key]
    else:
        # Use selected base + denom
        key = (base_choice, denom)
        r = all_results.get(key, all_results[best_key])
        if key not in all_results:
            st.warning(f"Base {base_choice} d={denom} not in results; "
                       f"showing best found.")

    st.session_state.update({
        'fit_results': all_results,
        'best_key': best_key,
        'vals': x,
        'dataset_info': {
            "filename": filename or data_source,
            "column": names[0] if len(names)==1 else "values",
            "n": len(x),
            "range": f"{x.min():.3g}–{x.max():.3g}",
            "decades": orders,
        },
        'chat_history': [],
        'validity': (passes, warnings_v, issues_v),
    })

    # ── Results ───────────────────────────────────────────────────
    st.header("Results")
    if not auto_mode and anchor_choice == "Custom":
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
    st.caption(anchor_note)

    # Metrics row
    p = r['p_value']
    z = r['z']
    eff = r['effect_pct']
    sig_label = ("✓✓ Strong (p<0.01)" if p < 0.01
                 else "✓ Marginal (p<0.05)" if p < 0.05
                 else "✗ Not significant")

    mc1, mc2, mc3, mc4, mc5 = st.columns(5)
    mc1.metric("Best base r",    r['base_name'],
               help="Base of geometric lattice")
    mc2.metric("Denominator d",  r['denom'],
               help="Lattice spaced at r^(k/d)")
    mc3.metric("RMS residual",   f"{r['rms']:.4f}",
               delta=f"{r['rms_pct']:.0f}% of max",
               delta_color="inverse")
    mc4.metric("Z-score",        f"{z:+.2f}σ",
               help="Sigma above null mean")
    mc5.metric("p-value",        f"{p:.4f}",
               delta=sig_label,
               delta_color="normal" if p < 0.05 else "off")

    # Effect size callout
    if eff > 5:
        st.success(f"**Effect size:** {eff:.1f}% RMS reduction below null mean. "
                   f"Z = {z:+.2f}σ")
    elif eff > 0:
        st.info(f"**Effect size:** {eff:.1f}% RMS reduction (marginal).")
    else:
        st.warning(f"**Effect size:** negative ({eff:.1f}%) — "
                   "data is *less* lattice-aligned than random.")

    # Bonferroni note if marginal
    if 0.01 < p < 0.05:
        st.caption(
            "⚠️ **Multiple testing note:** With multiple bases tested, "
            "a Bonferroni correction raises the significance threshold to p<0.017. "
            "This result does not survive that correction; interpret with caution."
        )

    result_tabs = st.tabs(["📈 Plot", "🔢 All bases", "📋 Export"])

    with result_tabs[0]:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        # Left: log histogram with lattice
        ax = axes[0]
        ax.hist(np.log10(x), bins=min(40, len(x)//2+1),
                color="#4f8ef7", edgecolor="white", linewidth=0.4, alpha=0.85)
        base_v = r['base']
        anchor_v = r.get('anchor', x.min())
        d_v = r['denom']
        lo10, hi10 = np.log10(x.min()), np.log10(x.max())
        for k in range(-200, 200):
            xl = np.log10(anchor_v) + k/d_v * np.log10(base_v)
            if lo10 <= xl <= hi10:
                ax.axvline(xl, color="orange", lw=0.8, alpha=0.6)
        ax.set_xlabel("log₁₀(value)", fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.set_title(f"Distribution with lattice (r={r['base_name']}, d={d_v})",
                     fontsize=10)
        ax.grid(True, alpha=0.2)

        # Right: null test
        ax = axes[1]
        null_arr = r.get('null_arr', np.array([r['null_mean']]))
        ax.hist(null_arr, bins=40, color="#aaaaaa", alpha=0.7,
                density=True, label="Null distribution")
        ax.axvline(r['rms'], color="#e74c3c", lw=2.5,
                   label=f"Observed RMS={r['rms']:.4f}")
        ax.axvline(r['null_mean'], color="#3498db", lw=2,
                   linestyle="--", label=f"Null mean={r['null_mean']:.4f}")
        ax.set_xlabel("RMS residual", fontsize=10)
        ax.set_title(f"Null test: p={p:.4f}, z={z:+.2f}σ", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

        fig.suptitle(f"LatticeFit — {filename or 'dataset'}  "
                     f"(n={len(x)}, {orders:.1f} decades)",
                     fontsize=11, fontweight="bold")
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        st.image(buf, use_container_width=True)
        st.download_button("⬇ Download plot", buf.getvalue(),
                           "latticefit_plot.png", "image/png")
        plt.close(fig)

    with result_tabs[1]:
        st.markdown("**All bases tested** — sorted by z-score:")
        rows = []
        for (bname, d), res in all_results.items():
            rows.append({
                "Base": bname, "d": d,
                "RMS": f"{res['rms']:.4f}",
                "% of max": f"{res['rms_pct']:.0f}%",
                "Effect": f"{res['effect_pct']:+.1f}%",
                "Z": f"{res['z']:+.2f}",
                "p": f"{res['p_value']:.4f}",
                "Sig": "✓" if res['p_value'] < 0.05 else "",
            })
        rows.sort(key=lambda r: float(r['Z']), reverse=True)
        st.dataframe(pd.DataFrame(rows), use_container_width=True,
                     hide_index=True)

    with result_tabs[2]:
        info = st.session_state['dataset_info']
        json_out = {
            "latticefit_version": "0.3.0",
            "dataset": info,
            "validity": {"passes": passes, "warnings": warnings_v,
                         "issues": issues_v},
            "best_result": {k: (float(v) if isinstance(v, (np.floating, float))
                                else (v.tolist() if hasattr(v, 'tolist') else v))
                            for k, v in r.items() if k != 'null_arr'},
        }
        st.download_button("⬇ Download JSON",
                           json.dumps(json_out, indent=2),
                           "latticefit_results.json", "application/json")
        st.markdown("**Cite LatticeFit:**")
        st.code(
            "Gentry, M. L. (2026). LatticeFit v0.3.0. "
            "https://github.com/drmlgentry/latticefit\n"
            "DOI: 10.5281/zenodo.19225731",
            language="text"
        )

    # ── AI interpretation ─────────────────────────────────────────
    if use_ai and HAS_ANTHROPIC:
        st.markdown("---")
        st.header("AI Interpretation")
        info = st.session_state['dataset_info']

        if not st.session_state['chat_history']:
            with st.spinner("Generating interpretation..."):
                interp = ai_interpret(r, info, orders, issues_v)
                st.session_state['chat_history'].append(
                    {"role": "assistant", "content": interp})

        for msg in st.session_state['chat_history']:
            with st.chat_message(msg['role']):
                st.write(msg['content'])

        if prompt := st.chat_input("Ask a follow-up question..."):
            st.session_state['chat_history'].append(
                {"role": "user", "content": prompt})
            with st.chat_message("user"):
                st.write(prompt)
            with st.chat_message("assistant"):
                with st.spinner("Thinking..."):
                    resp = ai_interpret(r, info, orders, issues_v,
                                        st.session_state['chat_history'])
                st.write(resp)
                st.session_state['chat_history'].append(
                    {"role": "assistant", "content": resp})

# ══════════════════════════════════════════════════════════════════
# TAB 2: CROSS-DOMAIN SURVEY
# ══════════════════════════════════════════════════════════════════
with tab_crossdomain:
    st.header("Cross-Domain Survey")
    st.markdown(
        "Results of applying LatticeFit (φ-lattice, d=4) across 14 datasets "
        "spanning physics, biology, chemistry, pharmacology, and finance. "
        "All datasets were tested with a log-uniform null (5,000 trials)."
    )

    # Summary counts
    n_genuine  = sum(1 for r in CROSS_DOMAIN if r[8]=="GENUINE")
    n_null_res = sum(1 for r in CROSS_DOMAIN if r[8]=="NULL")
    n_artifact = sum(1 for r in CROSS_DOMAIN if r[8]=="ARTIFACT")
    n_marginal = sum(1 for r in CROSS_DOMAIN if r[8]=="MARGINAL")

    sc1, sc2, sc3, sc4 = st.columns(4)
    sc1.metric("Genuine signals",  n_genuine,  help="p<0.05, valid criteria met")
    sc2.metric("Null results",     n_null_res, help="No lattice structure detected")
    sc3.metric("Marginal",         n_marginal, help="0.05<p<0.10")
    sc4.metric("Artifacts",        n_artifact, help="Narrow range or binned data")

    st.markdown("---")

    # Full table
    rows = []
    for (domain, n, rms, null_rms, z, p_str,
         orders_d, cat, status) in CROSS_DOMAIN:
        eff = (null_rms - rms) / null_rms * 100
        status_badge = {
            "GENUINE":    "✅ Genuine",
            "NULL":       "⬜ Null",
            "ARTIFACT":   "⚠️ Artifact",
            "MARGINAL":   "🟡 Marginal",
            "NEEDS NULL": "🔵 Needs physics null",
        }.get(status, status)
        rows.append({
            "Domain": domain,
            "n": f"{n:,}",
            "Orders": f"{orders_d:.1f}",
            "RMS": f"{rms:.4f}",
            "Effect": f"{eff:+.1f}%",
            "Z": f"{z:+.2f}σ",
            "p": p_str,
            "Category": cat,
            "Status": status_badge,
        })

    df_cross = pd.DataFrame(rows)
    st.dataframe(df_cross, use_container_width=True, hide_index=True)

    st.markdown("---")
    st.subheader("Pattern")
    st.markdown("""
**φ-lattice signals appear in:**
- Datasets governed by quantum arithmetic constraints (SM masses)
- Population-level biological frequency distributions (allele frequencies)
- Evolutionarily constrained size distributions (cetacean body masses)
- Single-series pharmacological data (EGFR SAR)

**φ-lattice signals are absent in:**
- Financial returns (governed by random walk / fat-tailed distributions)
- Crystal volumes (governed by packing geometry, not φ)
- Aggregated pharmacological data (heterogeneous assay conditions)
- Datasets governed by known non-φ physical laws (Gutenberg-Richter)

**Artifact classes identified:**
1. **Narrow range** (< 1 order): nuclear BE/A spans 0.97 orders → 97% of nuclei 
   fall in 2 adjacent lattice bins by necessity, not geometry.
2. **Binned data**: patent IC50 assays report 7 discrete values → 
   apparent structure from bin spacing, not biology.
    """)

    # Z-score bar chart
    st.subheader("Z-score comparison")
    fig, ax = plt.subplots(figsize=(10, 5))
    domains = [r[0][:30] for r in CROSS_DOMAIN]
    zscores = [r[4] for r in CROSS_DOMAIN]
    statuses = [r[8] for r in CROSS_DOMAIN]
    colors = {"GENUINE": "#2e7d32", "NULL": "#888888",
              "ARTIFACT": "#e65100", "MARGINAL": "#f57f17",
              "NEEDS NULL": "#1565c0"}
    bar_colors = [colors.get(s, "#888") for s in statuses]
    bars = ax.barh(domains, zscores, color=bar_colors, alpha=0.8)
    ax.axvline(0, color="black", lw=1)
    ax.axvline(1.96, color="green", lw=1.5, linestyle="--",
               label="p=0.05 (z=1.96)")
    ax.axvline(-1.96, color="green", lw=1.5, linestyle="--")
    ax.set_xlabel("Z-score (φ-lattice vs log-uniform null)", fontsize=10)
    ax.set_title("LatticeFit cross-domain survey", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2, axis="x")
    # Legend patches
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=v, label=k)
                       for k, v in colors.items()]
    ax.legend(handles=legend_elements, fontsize=8,
              loc="lower right")
    fig.tight_layout()
    buf2 = io.BytesIO()
    fig.savefig(buf2, format="png", dpi=150, bbox_inches="tight")
    buf2.seek(0)
    st.image(buf2, use_container_width=True)
    st.download_button("⬇ Download chart", buf2.getvalue(),
                       "cross_domain_zscores.png", "image/png")
    plt.close(fig)

# ══════════════════════════════════════════════════════════════════
# TAB 3: VALIDITY GUIDE
# ══════════════════════════════════════════════════════════════════
with tab_validity_guide:
    st.header("Validity Guide")
    st.markdown("""
Before trusting a LatticeFit result, verify all four criteria:

---

### 1. Minimum dynamic range ≥ 3 orders of magnitude

**Why:** When data spans less than ~2 lattice bins, any lattice will show
apparent clustering — the null distribution becomes non-uniform.

**Diagnostic:** `log10(max/min) ≥ 3`

**Counter-example:** Nuclear binding energies per nucleon (AME2020) span
only 0.97 orders → z=+7.18 but entirely explained by the fact that 97% of
nuclei fall in two adjacent bins near 8 MeV/nucleon.

**Rule:** Datasets with < 2 orders should be excluded from lattice analysis.

---

### 2. Data must not be binned or discretized

**Why:** Fixed bin boundaries align with lattice points by coincidence,
creating spurious signal.

**Diagnostic:** `len(unique values) / len(total values) > 0.5`

**Counter-example:** Patent IC50 data (CHEMBL3706050) has 7 discrete values
(1, 5.5, 20, 65, 200, 650, 1000 nM) → z=+9.81 from bin spacing alignment.

**Rule:** Screen for discretized data before analysis.

---

### 3. Use a physics-motivated null where available

**Why:** The log-uniform null assumes no prior structure. If a theoretical
model predicts the distribution (e.g., Gutenberg-Richter for earthquakes,
Bethe-Weizsäcker for nuclear masses), use residuals from that model.

**Counter-example:** Earthquake energies with GR-law null → z=-3.58 (null).
Without the GR correction, a naïve test might show spurious structure.

---

### 4. Effect size matters, not just z-score

**Why:** With large n, even tiny deviations from the null become statistically
significant. A 0.1% RMS reduction at n=100,000 may be statistically significant
but physically meaningless.

**Diagnostic:** Report `(null_RMS - obs_RMS) / null_RMS × 100%`

**Guideline:** Effect sizes below 2% should be interpreted with caution
regardless of p-value.

---

### 5. Multiple testing correction

When testing multiple bases (φ, √2, 2, e, 10) and denominators, apply
Bonferroni correction. With 5 bases × 5 denominators = 25 tests, the
corrected significance threshold is p < 0.002.

A result with p=0.04 for one base that does not survive Bonferroni correction
should be reported as exploratory, not confirmatory.
    """)

# ══════════════════════════════════════════════════════════════════
# TAB 4: ABOUT
# ══════════════════════════════════════════════════════════════════
with tab_about:
    st.header("About LatticeFit")
    st.markdown("""
**LatticeFit** tests whether a dataset of positive measurements clusters
near a discrete geometric lattice:

$$x_i \\approx A \\cdot r^{k_i/d}, \\quad k_i \\in \\mathbb{Z}$$

Given anchor $A$, base $r$, and denominator $d$, it:

1. Assigns the nearest integer label $k_i$ to each observation  
2. Computes logarithmic residuals $\\delta_i = |\\log_r(x_i/A) - k_i/d|$  
3. Tests whether the RMS residual is smaller than a log-uniform null  

---

### Neutrino mass prediction

The φ-lattice with anchor $m_e$ and denominator $d=4$ predicts:

| Neutrino | Mass | q |
|---|---|---|
| ν₁ | 8.387 meV | −149 |
| ν₂ | 12.033 meV | −146 |
| ν₃ | 50.972 meV | −134 |

**Sum: 71.4 meV** — testable by CMB-S4 (~2030).  
**Normal ordering only** — inverted ordering has no solution within 3σ.  
**Exact relation:** $m_3/m_2 = \\varphi^3$

---

### Algebraic bridge

The Alexander polynomial of the PMNS hyperbolic 3-manifold (m003) is:

$$\\Delta(t) = t^2 + 3t + 1$$

with roots $-\\varphi^2$ and $-\\varphi^{-2}$ (exact). Mahler measure:

$$\\log M(\\Delta) = 2\\log\\varphi = 2\\,\\mathrm{Reg}(\\mathbb{Q}(\\sqrt{5}))$$

This equals twice the SM mass lattice spacing — a purely topological
connection between the lepton sector and the golden ratio.

---

**Citation:**
```
Gentry, M. L. (2026). LatticeFit v0.3.0.
https://github.com/drmlgentry/latticefit
DOI: 10.5281/zenodo.19225731
```

**Patent pending:** US Provisional 64/013,306  
**ORCID:** 0009-0006-4550-2663
    """)

# ── Footer ────────────────────────────────────────────────────────
st.divider()
st.caption(
    "LatticeFit v0.3.0 · "
    "[GitHub](https://github.com/drmlgentry/latticefit) · "
    "DOI: [10.5281/zenodo.19225731](https://doi.org/10.5281/zenodo.19225731) · "
    "Patent pending US 64/013,306"
)
