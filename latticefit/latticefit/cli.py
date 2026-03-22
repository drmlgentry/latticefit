"""
latticefit.cli
==============
Command-line interface.

Usage
-----
    latticefit data.csv --anchor first --base phi --denom 4
    latticefit data.csv --auto
    latticefit data.csv --anchor 5.11e-4 --base 1.6180339887 --denom 4 --plot fit.png
"""

from __future__ import annotations
import argparse
import sys
import numpy as np


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="latticefit",
        description="Fit positive data to a geometric lattice x ≈ A·r^(k/d).",
    )
    parser.add_argument("file", help="CSV file (one positive number per line, "
                        "optionally two columns: name,value)")
    parser.add_argument("--anchor", default="first",
                        help="Anchor A: 'first', 'geomean', or a float value")
    parser.add_argument("--base", default="phi",
                        help="Base r: 'phi', 'e', '2', 'sqrt2', '10', "
                             "or a float value")
    parser.add_argument("--denom", type=int, default=4,
                        help="Denominator d (default: 4)")
    parser.add_argument("--auto", action="store_true",
                        help="Auto-discover best lattice parameters")
    parser.add_argument("--null", type=int, default=0,
                        help="Run log-uniform null test with N trials "
                             "(e.g. --null 10000)")
    parser.add_argument("--plot", default=None,
                        help="Save fit plot to this filename")
    parser.add_argument("--version", action="version", version="0.1.0")
    args = parser.parse_args(argv)

    # ── Load data ─────────────────────────────────────────────────
    names, values = [], []
    try:
        with open(args.file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(",")
                if len(parts) == 2:
                    names.append(parts[0].strip())
                    values.append(float(parts[1].strip()))
                else:
                    values.append(float(parts[0].strip()))
    except FileNotFoundError:
        print(f"Error: file '{args.file}' not found.", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error parsing data: {e}", file=sys.stderr)
        sys.exit(1)

    x = np.array(values)
    names = names if names else None

    from latticefit import fit, discover, PHI, KNOWN
    from latticefit.stats import log_uniform_null

    # ── Resolve parameters ────────────────────────────────────────
    BASES = {**KNOWN, "phi": PHI}
    base_val = BASES.get(args.base.lower())
    if base_val is None:
        try:
            base_val = float(args.base)
        except ValueError:
            print(f"Unknown base '{args.base}'. Use phi/e/2/sqrt2/10 "
                  "or a float.", file=sys.stderr)
            sys.exit(1)

    if args.anchor == "first":
        anchor_val = float(x[0])
    elif args.anchor == "geomean":
        anchor_val = float(np.exp(np.mean(np.log(x))))
    else:
        try:
            anchor_val = float(args.anchor)
        except ValueError:
            print(f"Unknown anchor '{args.anchor}'.", file=sys.stderr)
            sys.exit(1)

    # ── Fit ───────────────────────────────────────────────────────
    if args.auto:
        print("Auto-discovering lattice parameters…")
        results = discover(x, names=names, top_n=3)
        for i, res in enumerate(results):
            print(f"\n─── Rank {i+1} ───")
            print(res.summary())
        result = results[0]
    else:
        result = fit(x, anchor=anchor_val, base=base_val,
                     denom=args.denom, names=names)
        print(result.summary())

    # ── Null test ─────────────────────────────────────────────────
    if args.null > 0:
        print(f"\nRunning log-uniform null test ({args.null} trials)…")
        null = log_uniform_null(result, n_trials=args.null)
        print(null.summary())

    # ── Plot ─────────────────────────────────────────────────────
    if args.plot:
        from latticefit.plots import plot_fit
        plot_fit(result, outfile=args.plot)
        print(f"\nPlot saved to {args.plot}")


if __name__ == "__main__":
    main()
