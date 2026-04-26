"""
latticefit.cli
==============
Command-line interface for LatticeFit.

Usage examples
--------------
    # Basic fit with phi-lattice
    latticefit data.csv

    # Specify parameters explicitly
    latticefit data.csv --anchor 5.11e-4 --base phi --denom 4

    # Auto-discover best lattice
    latticefit data.csv --auto --top 5

    # With null tests
    latticefit data.csv --null 10000 --sector-null 0,1,2:3,4,5:6,7,8

    # Save plot and report
    latticefit data.csv --plot fit.png --report report.html

    # JSON output for pipeline integration
    latticefit data.csv --json

    # Run built-in Standard Model demo
    latticefit --demo

    # Pipe data directly
    echo "1.0,2.618,6.854,17.944" | latticefit -
"""

from __future__ import annotations
import argparse
import json
import sys
import os
import numpy as np
from datetime import datetime


# ── Colour helpers (works on most terminals) ───────────────────────
def _c(text, code): return f"\033[{code}m{text}\033[0m" if sys.stdout.isatty() else text
def green(t):  return _c(t, "32")
def yellow(t): return _c(t, "33")
def cyan(t):   return _c(t, "36")
def bold(t):   return _c(t, "1")
def dim(t):    return _c(t, "2")


# ── Built-in demo dataset ──────────────────────────────────────────
SM_DEMO = {
    "names":  ["e", "mu", "tau", "u", "c", "t", "d", "s", "b", "W", "Z", "H"],
    "values": [5.10999e-4, 0.105658, 1.77686,
               0.00216,    1.275,    172.76,
               0.00467,    0.0934,   4.18,
               80.379,     91.1876,  125.25],
    "anchor": 5.10999e-4,
    "base":   (1 + 5**0.5) / 2,
    "denom":  4,
    "sectors": [[0,1,2],[3,4,5],[6,7,8]],
}


def _load_csv(path: str):
    """Load name,value pairs from CSV or stdin (path='-')."""
    names, values = [], []
    src = sys.stdin if path == "-" else open(path)
    try:
        for line in src:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 2:
                try:
                    # Try name,value format
                    values.append(float(parts[-1]))
                    names.append(",".join(parts[:-1]))
                except ValueError:
                    # All numeric — treat first as value
                    values.append(float(parts[0]))
            else:
                values.append(float(parts[0]))
    finally:
        if path != "-":
            src.close()
    return names or None, np.array(values, dtype=float)


def _resolve_base(s: str) -> float:
    from latticefit import KNOWN, PHI
    BASES = {**KNOWN, "phi": PHI, "golden": PHI, "φ": PHI}
    b = BASES.get(s.lower())
    if b is not None:
        return b
    return float(s)


def _resolve_anchor(s: str, x: np.ndarray) -> float:
    if s == "first":      return float(x[0])
    if s == "last":       return float(x[-1])
    if s == "min":        return float(x.min())
    if s == "geomean":    return float(np.exp(np.mean(np.log(x))))
    if s == "median":     return float(np.median(x))
    return float(s)


def _parse_sectors(s: str) -> list[list[int]]:
    """Parse '0,1,2:3,4,5:6,7,8' into [[0,1,2],[3,4,5],[6,7,8]]."""
    return [[int(i) for i in grp.split(",")] for grp in s.split(":")]


def _progress(msg: str):
    print(dim(f"  {msg}…"), flush=True)


def _make_html_report(result, nulls: list, params: dict) -> str:
    """Generate a self-contained HTML report."""
    from latticefit import PHI
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")
    names = result.names or [str(i) for i in range(len(result.data))]

    rows = ""
    csv_lines = ["name,observed,k,predicted,residual,pct_error"]
    for n, xo, k, xp, res in zip(names, result.data, result.labels,
                                   result.predicted, result.residuals):
        pct = abs(xo - xp) / xo * 100
        bar_w = int(res / 0.125 * 120)
        rows += (
            f"<tr><td>{n}</td><td>{xo:.5g}</td><td>{k}</td>"
            f"<td>{xp:.5g}</td>"
            f"<td><div style='display:flex;align-items:center;gap:6px'>"
            f"<div style='width:{bar_w}px;height:10px;background:#4a90d9;"
            f"border-radius:3px'></div>{res:.4f}</div></td>"
            f"<td>{pct:.1f}%</td></tr>\n"
        )
        csv_lines.append(f"{n},{xo:.6g},{k},{xp:.6g},{res:.5f},{pct:.2f}")

    null_html = ""
    for nt in nulls:
        z = (nt.null_mean - nt.observed_rms) / max(nt.null_std, 1e-12)
        sig = "✓ significant" if nt.p_value < 0.05 else "~ marginal" if nt.p_value < 0.10 else "✗ not significant"
        null_html += f"""
        <div class="null-card">
          <h3>{nt.test_name} null test</h3>
          <table class="stats">
            <tr><td>Observed RMS</td><td><b>{nt.observed_rms:.5f}</b></td></tr>
            <tr><td>Null mean ± std</td><td>{nt.null_mean:.5f} ± {nt.null_std:.5f}</td></tr>
            <tr><td>Z-score</td><td>{z:.2f}</td></tr>
            <tr><td>p-value</td><td><b>{nt.p_value:.4f}</b> — {sig}</td></tr>
            <tr><td>Trials</td><td>{nt.n_trials:,}</td></tr>
          </table>
        </div>"""

    csv_content = "\n".join(csv_lines)
    base_name = {(1+5**0.5)/2: "φ (golden ratio)", 2.71828: "e",
                 2.0: "2", 10.0: "10"}.get(round(result.base, 4), f"{result.base:.6g}")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>LatticeFit Report</title>
<style>
  body {{ font-family: 'Helvetica Neue', sans-serif; max-width: 860px;
          margin: 40px auto; padding: 0 20px; color: #222; line-height:1.5 }}
  h1   {{ color: #1a3a5c; border-bottom: 2px solid #4a90d9; padding-bottom:8px }}
  h2   {{ color: #1a3a5c; margin-top: 32px }}
  h3   {{ color: #2a6496; margin: 0 0 8px 0 }}
  .meta {{ background:#f4f8ff; border-left:4px solid #4a90d9;
           padding:12px 16px; border-radius:4px; font-size:0.9em }}
  table {{ border-collapse:collapse; width:100%; margin:12px 0 }}
  th,td {{ padding:8px 12px; text-align:left; border-bottom:1px solid #e0e8f0; white-space:nowrap }}
  th    {{ background:#eaf2ff; font-weight:600 }}
  tr:hover {{ background:#f7faff }}
  .rms  {{ font-size:1.4em; color:#1a3a5c; font-weight:bold }}
  .null-card {{ background:#f9f9f9; border:1px solid #dde; border-radius:6px;
                padding:16px; margin:12px 0 }}
  .stats td:first-child {{ color:#555; width:160px }}
  .footer {{ margin-top:40px; font-size:0.8em; color:#888; border-top:1px solid #eee;
             padding-top:12px }}
</style>
</head>
<body>
<h1>LatticeFit Analysis Report</h1>
<div class="meta">
  Generated: {ts} &nbsp;|&nbsp;
  Base: {base_name} &nbsp;|&nbsp;
  Denominator: d = {result.denom} &nbsp;|&nbsp;
  Anchor: A = {result.anchor:.5g} &nbsp;|&nbsp;
  n = {len(result.data)} observations
</div>

<h2>Fit Quality</h2>
<p>RMS residual: <span class="rms">{result.rms:.5f}</span>
&nbsp; <span style="color:#666">(max possible = {0.5/result.denom:.4f},
fraction used = {result.rms/(0.5/result.denom)*100:.1f}%)</span></p>

<h2>Results Table</h2>
<table>
  <tr><th>Name</th><th>Observed</th><th>k</th><th>Predicted</th>
      <th>Residual δ</th><th>|Δ|/m</th></tr>
  {rows}
</table>
<details style="margin-top:8px">
  <summary style="cursor:pointer;color:#4a90d9;font-size:0.9em">
    ▼ Download results as CSV
  </summary>
  <textarea id="csvdata" style="width:100%;height:120px;font-family:monospace;
    font-size:0.8em;margin-top:6px;border:1px solid #dde;border-radius:4px;
    padding:8px" readonly>{csv_content}</textarea>
  <button onclick="
    var a=document.createElement('a');
    a.href='data:text/csv;charset=utf-8,'+
    encodeURIComponent(document.getElementById('csvdata').value);
    a.download='latticefit_results.csv';a.click();"
    style="margin-top:6px;padding:6px 14px;background:#4a90d9;color:white;
    border:none;border-radius:4px;cursor:pointer;font-size:0.9em">
    ⬇ Download CSV
  </button>
</details>

<h2>Statistical Validation</h2>
{null_html if null_html else '<p style="color:#888"><i>No null tests run. Use --null N to add.</i></p>'}

<div class="footer">
  LatticeFit 0.2.0 &nbsp;|&nbsp;
  M. L. Gentry, drmlgentry@protonmail.com &nbsp;|&nbsp;
  <a href="https://github.com/drmlgentry/latticefit">github.com/drmlgentry/latticefit</a>
</div>
</body>
</html>"""


def _json_output(result, nulls: list) -> dict:
    return {
        "latticefit_version": "0.2.0",
        "parameters": {
            "anchor": result.anchor,
            "base":   result.base,
            "denom":  result.denom,
        },
        "fit": {
            "rms":  result.rms,
            "n":    len(result.data),
            "max_possible_rms": 0.5 / result.denom,
        },
        "observations": [
            {
                "name":      (result.names[i] if result.names else str(i)),
                "observed":  float(result.data[i]),
                "label_k":   int(result.labels[i]),
                "predicted": float(result.predicted[i]),
                "residual":  float(result.residuals[i]),
                "pct_error": float(abs(result.data[i] - result.predicted[i])
                                   / result.data[i] * 100),
            }
            for i in range(len(result.data))
        ],
        "null_tests": [
            {
                "test":         nt.test_name,
                "p_value":      nt.p_value,
                "null_mean":    nt.null_mean,
                "null_std":     nt.null_std,
                "n_trials":     nt.n_trials,
            }
            for nt in nulls
        ],
    }



def _run_lucas(args, x, names):
    """Run Lucas-mode fit (base fixed to phi by theorem)."""
    from latticefit.lucas import fit_lucas, PRIME_LUCAS
    anchor = None if args.anchor == "first" else float(args.anchor)
    result = fit_lucas(x, anchor=anchor, q_grid=args.denom,
                       n_null=max(args.null, 1000))
    if not args.json:
        print(bold("\n── Lucas-mode fit (base = φ, theoretically derived) ──"))
        print(result.summary())
        print()
        if result.prime_lucas_hits:
            print(green("  Prime Lucas hits (prime dictionary of PMNS covering tower):"))
            for a in result.prime_lucas_hits:
                nm = names[x.tolist().index(a.value)] if names else f"{a.value:.5g}"
                print(f"    {nm}: q={a.q:+.2f}  L_{{k}}={a.lucas_integer}  "
                      f"delta={a.delta_pct:.3f}%")
        print(f"\n  Lucas-integer fraction: {result.lucas_integer_fraction:.3f}")
        print(f"  p-value: {result.p_value:.4f}  (null mean={result.null_mean:.5f})")
    else:
        import json
        print(json.dumps({
            "mode": "lucas",
            "anchor": result.anchor,
            "rms": result.rms,
            "p_value": result.p_value,
            "lucas_integer_fraction": result.lucas_integer_fraction,
            "observations": [
                {"value": a.value, "q": a.q, "lk": a.lk,
                 "is_integer_lucas": a.is_integer_lucas,
                 "lucas_integer": a.lucas_integer,
                 "is_prime_lucas": a.is_prime_lucas_,
                 "delta_pct": a.delta_pct}
                for a in result.assignments
            ]
        }, indent=2))


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="latticefit",
        description="Fit positive data to a geometric lattice  x ≈ A·r^(k/d).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  latticefit data.csv
  latticefit data.csv --auto --top 5
  latticefit data.csv --base phi --denom 4 --null 10000 --plot fit.png
  latticefit data.csv --sector-null 0,1,2:3,4,5 --null 50000
  latticefit data.csv --report report.html --json
  latticefit --demo
  echo "1.0,2.618,6.854" | latticefit -
  latticefit data.csv --lucas --anchor 5.11e-4
        """,
    )

    parser.add_argument("file", nargs="?", default=None,
                        help="CSV file (name,value or value per line); use - for stdin")
    parser.add_argument("--anchor", default="first",
                        help="Anchor A: first|last|min|geomean|median|<float> (default: first)")
    parser.add_argument("--base", default="phi",
                        help="Base r: phi|e|2|sqrt2|10|<float> (default: phi)")
    parser.add_argument("--denom", type=int, default=4,
                        help="Denominator d (default: 4)")
    parser.add_argument("--auto", action="store_true",
                        help="Auto-discover best (base, denom, anchor) combination")
    parser.add_argument("--top", type=int, default=3,
                        help="Number of top results to show with --auto (default: 3)")
    parser.add_argument("--null", type=int, default=0, metavar="N",
                        help="Run log-uniform null test with N trials")
    parser.add_argument("--sector-null", default=None, metavar="SECTORS",
                        help="Sector-anchor null: comma-separated indices per sector, "
                             "sectors separated by colon. E.g. 0,1,2:3,4,5:6,7,8")
    parser.add_argument("--plot", default=None, metavar="FILE",
                        help="Save fit plot to FILE (png/pdf/svg)")
    parser.add_argument("--null-plot", default=None, metavar="FILE",
                        help="Save null distribution plot to FILE")
    parser.add_argument("--report", default=None, metavar="FILE",
                        help="Save HTML report to FILE")
    parser.add_argument("--json", action="store_true",
                        help="Print JSON output (machine-readable)")
    parser.add_argument("--demo", action="store_true",
                        help="Run built-in Standard Model demo")
    parser.add_argument("--lucas", action="store_true",
                        help="Lucas mode: fix base=phi by theorem, report Lucas integer and prime Lucas hits")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Suppress progress messages")
    parser.add_argument("--version", action="version", version="latticefit 0.2.0")
    parser.add_argument("--cite", action="store_true",
                        help="Print citation information and exit")

    args = parser.parse_args(argv)

    # ── Cite mode ─────────────────────────────────────────────────
    if args.cite:
        print("""
To cite LatticeFit in academic work, please use:

  Gentry, M. L. (2026). LatticeFit: Discrete lattice fitting with
  statistical validation (v0.1.0). GitHub.
  https://github.com/drmlgentry/latticefit

BibTeX:
  @software{latticefit2026,
    author  = {Gentry, Marvin L.},
    title   = {{LatticeFit}: Discrete lattice fitting with statistical validation},
    year    = {2026},
    version = {0.1.0},
    url     = {https://github.com/drmlgentry/latticefit}
  }
""")
        return

    # ── Demo mode ─────────────────────────────────────────────────
    if args.demo:
        print(bold("\nLatticeFit — Standard Model Demo"))
        print(dim("  Fitting SM fermion masses to the φ-lattice\n"))
        d = SM_DEMO
        names  = d["names"]
        values = d["values"]
        anchor = d["anchor"]
        base   = d["base"]
        denom  = d["denom"]
        sectors = d["sectors"]
        null_n  = max(args.null, 10_000)
        sector_null_str = ":".join(",".join(str(i) for i in s) for s in sectors)
        # Recurse with demo params
        demo_argv = ["-", "--anchor", str(anchor), "--base", str(base),
                     "--denom", str(denom),
                     "--null", str(null_n),
                     "--sector-null", sector_null_str]
        if args.plot:       demo_argv += ["--plot", args.plot]
        if args.report:     demo_argv += ["--report", args.report]
        if args.json:       demo_argv += ["--json"]
        # Feed data via a temp file
        import tempfile, os
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv",
                                          delete=False) as f:
            for n, v in zip(names, values):
                f.write(f"{n},{v}\n")
            tmppath = f.name
        try:
            main([tmppath] + demo_argv[1:])
        finally:
            os.unlink(tmppath)
        return

    # ── Require file ──────────────────────────────────────────────
    if args.file is None:
        parser.print_help()
        sys.exit(0)

    # ── Load data ─────────────────────────────────────────────────
    if not args.quiet:
        src = "stdin" if args.file == "-" else args.file
        _progress(f"Loading {src}")
    try:
        names, x = _load_csv(args.file)
    except FileNotFoundError:
        print(f"Error: '{args.file}' not found.", file=sys.stderr); sys.exit(1)
    except ValueError as e:
        print(f"Error parsing data: {e}", file=sys.stderr); sys.exit(1)

    if len(x) == 0:
        print("Error: no data loaded.", file=sys.stderr); sys.exit(1)
    if np.any(x <= 0):
        print("Error: all values must be strictly positive.", file=sys.stderr); sys.exit(1)

    from latticefit import fit, discover, PHI, KNOWN
    from latticefit.stats import log_uniform_null, sector_anchor_null

    nulls = []

    # ── Fit ───────────────────────────────────────────────────────
    if args.auto:
        if not args.quiet:
            _progress(f"Auto-discovering lattice (top {args.top})")
        results = discover(x, names=names, top_n=args.top)
        if not args.json:
            print(bold("\n── Auto-discovery results ──"))
            for i, res in enumerate(results):
                medal = ["🥇","🥈","🥉","  ","  "][min(i,4)]
                print(f"\n{medal} Rank {i+1}  "
                      f"base={res.base:.5g}  d={res.denom}  "
                      f"A={res.anchor:.5g}  RMS={res.rms:.5f}")
                print(res.summary())
        result = results[0]
    else:
        try:
            base_val   = _resolve_base(args.base)
            anchor_val = _resolve_anchor(args.anchor, x)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr); sys.exit(1)

        if not args.quiet:
            _progress(f"Fitting  base={base_val:.5g}  d={args.denom}  "
                      f"A={anchor_val:.5g}")
        result = fit(x, anchor=anchor_val, base=base_val,
                     denom=args.denom, names=names)
        if not args.json:
            print(bold("\n── Fit results ──"))
            print(result.summary())
            max_rms = 0.5 / args.denom
            frac = result.rms / max_rms * 100
            print(f"\n{green('RMS = ' + f'{result.rms:.5f}')}"
                  f"  ({frac:.1f}% of max {max_rms:.4f})")

    # ── Lucas mode ────────────────────────────────────────────────
    if args.lucas:
        _run_lucas(args, x, names)
        return

    # ── Null tests ────────────────────────────────────────────────
    if args.null > 0:
        if not args.quiet:
            _progress(f"Log-uniform null ({args.null:,} trials)")
        nt = log_uniform_null(result, n_trials=args.null)
        nulls.append(nt)
        if not args.json:
            print(bold("\n── Log-uniform null ──"))
            print(nt.summary())

    if args.sector_null:
        try:
            sectors = _parse_sectors(args.sector_null)
        except ValueError:
            print("Error: --sector-null format is '0,1,2:3,4,5:6,7,8'",
                  file=sys.stderr); sys.exit(1)
        n_trials = max(args.null, 10_000)
        if not args.quiet:
            _progress(f"Sector-anchor null ({n_trials:,} trials, "
                      f"{len(sectors)} sectors)")
        nt2 = sector_anchor_null(result, sector_ids=sectors,
                                  n_trials=n_trials)
        nulls.append(nt2)
        if not args.json:
            print(bold("\n── Sector-anchor null ──"))
            print(nt2.summary())

    # ── JSON output ───────────────────────────────────────────────
    if args.json:
        print(json.dumps(_json_output(result, nulls), indent=2))

    # ── Plot ──────────────────────────────────────────────────────
    if args.plot:
        if not args.quiet:
            _progress(f"Saving plot → {args.plot}")
        from latticefit.plots import plot_fit
        plot_fit(result, outfile=args.plot)
        if not args.quiet:
            print(green(f"  ✓ Plot saved: {args.plot}"))

    if args.null_plot and nulls:
        if not args.quiet:
            _progress(f"Saving null plot → {args.null_plot}")
        from latticefit.plots import plot_null
        plot_null(nulls[0], outfile=args.null_plot)
        if not args.quiet:
            print(green(f"  ✓ Null plot saved: {args.null_plot}"))

    # ── HTML report ───────────────────────────────────────────────
    if args.report:
        if not args.quiet:
            _progress(f"Generating report → {args.report}")
        html = _make_html_report(result, nulls, vars(args))
        with open(args.report, "w", encoding="utf-8") as f:
            f.write(html)
        if not args.quiet:
            print(green(f"  ✓ Report saved: {args.report}"))


if __name__ == "__main__":
    main()
