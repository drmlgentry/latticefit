"""
Example 2: Musical Notes (Equal Temperament)
=============================================
The chromatic scale in equal temperament is an EXACT geometric lattice:

    f_n = 440 Hz × 2^(n/12)

This example demonstrates LatticeFit correctly recovering a known
lattice structure with near-zero residuals — a validation benchmark.
It also tests what happens when we use the WRONG base (phi instead of
2^(1/12)), showing the tool correctly identifies a poor fit.

This is the ideal demo for showing the tool works on known structure.
"""
import sys
sys.path.insert(0, "../")
import numpy as np
import latticefit

# Chromatic scale: A4 to A6 (25 notes, 2 octaves)
note_names = [
    "A4","A#4","B4","C5","C#5","D5","D#5","E5","F5","F#5","G5","G#5",
    "A5","A#5","B5","C6","C#6","D6","D#6","E6","F6","F#6","G6","G#6","A6"
]
freqs = np.array([440.0 * 2**(n/12) for n in range(25)])

print("=" * 60)
print("Musical Notes (Equal Temperament) — LatticeFit Analysis")
print("=" * 60)

# Fit with the CORRECT base: 2^(1/12) (semitone ratio)
semitone = 2 ** (1/12)
result_correct = latticefit.fit(freqs, anchor=440.0,
                                base=semitone, denom=1,
                                names=note_names)
print(f"\nCorrect base (2^1/12 = {semitone:.6f}):")
print(result_correct.summary())

# Fit with WRONG base: phi — should give large residuals
result_wrong = latticefit.fit(freqs, anchor=440.0,
                               base=latticefit.PHI, denom=4,
                               names=note_names)
print(f"\nWrong base (φ = {latticefit.PHI:.6f}):")
print(f"  RMS = {result_wrong.rms:.5f}  "
      f"(vs {result_correct.rms:.5f} for correct base)")
print(f"  → {result_wrong.rms/result_correct.rms:.0f}× worse fit")

# Auto-discover: should find 2^(1/12) or close
print("\nAuto-discovery:")
top = latticefit.discover(freqs, anchor_candidates=[440.0],
                           names=note_names, top_n=3) \
      if hasattr(latticefit.discover, '__wrapped__') else \
      latticefit.discover(freqs, names=note_names, top_n=3)
for i, r in enumerate(top):
    print(f"  #{i+1}: base={r.base:.6g}  d={r.denom}  "
          f"anchor={r.anchor:.4g}  RMS={r.rms:.6f}")

# Null test on correct fit
null = latticefit.log_uniform_null(result_correct, n_trials=10_000)
print()
print(null.summary())

# Plot
latticefit.plots.plot_fit(result_correct,
    outfile="musical_notes_fit.png",
    title="Musical notes (equal temperament)\n"
          r"$f_n = 440 \times 2^{n/12}$ — exact lattice recovery")
print("\nPlot saved: musical_notes_fit.png")
