#!/usr/bin/env pwsh
# step1_tag_latticefit.ps1
# Commit and tag LatticeFit v0.2.0

Set-Location "C:\dev\latticefit"

# Stage everything
git add -A

# Check status
Write-Host "Status:" -ForegroundColor Cyan
git status --short

# Commit
git commit -m "LatticeFit v0.2.0 — JOSS submission ready

New modules:
  latticefit/bootstrap.py  — bootstrap_ci(), propagate_uncertainty()
  latticefit/bundle.py     — generate_bundle(), bundle_from_csv()
  latticefit/models.py     — select_model() with AIC/BIC

Examples:
  examples/data/           — sm_masses.csv, musical_notes.csv
  examples/scripts/        — run_all_examples.py, prepare_benchmark_data.py
  validate_latticefit.py   — 18-dataset JOSS Table 1 reproduction

JOSS paper:
  paper.md                 — submission-ready, 18-dataset validation table
  paper.bib                — 12 bibliography entries

Patent: US provisional 64/013,306 (22 March 2026)
DOI pending JOSS review"

# Tag
git tag -a v0.2.0 -m "LatticeFit v0.2.0 — JOSS submission"

# Push with tags
git push origin main --tags

Write-Host ""
Write-Host "Done. v0.2.0 tagged and pushed." -ForegroundColor Green
Write-Host "Next: submit at https://joss.theoj.org/papers/new" -ForegroundColor Yellow
