# LatticeFit — JOSS Submission Checklist

## Pre-submission checklist (from joss.theoj.org)

### Software
- [x] Repository publicly available: https://github.com/drmlgentry/latticefit
- [x] OSI-approved license (MIT)
- [x] Version tagged: v0.2.0
- [x] Installable: `pip install latticefit`
- [x] Tests present (see tests/)
- [x] Documentation present (README.md, docstrings)
- [x] Examples/demos present (`latticefit --demo`)
- [x] CLI present (`latticefit`)
- [x] Python API present (`from latticefit import fit`)

### Paper (paper.md)
- [x] Title
- [x] Author list with ORCID
- [x] Affiliation
- [x] Summary (what the software does)
- [x] Statement of need (why it's needed, related software)
- [x] Algorithm description
- [x] Usage examples (Python API + CLI)
- [x] Validation / cross-domain results (Table 1, 18 datasets)
- [x] Reproducibility statement
- [x] Patent disclosure
- [x] References (paper.bib, 12 entries)

### Before submitting
- [ ] Run `pip install latticefit` from clean environment
- [ ] Run all demos: `latticefit --demo sm_masses` etc.
- [ ] Verify all tests pass: `pytest tests/`
- [ ] Tag release: `git tag v0.2.0 && git push --tags`
- [ ] Ensure paper.md + paper.bib are in repo root

## JOSS submission URL
https://joss.theoj.org/papers/new

## Submission fields
- Repository URL: https://github.com/drmlgentry/latticefit
- Branch: main
- Paper path: paper.md
- Language: Python
- Editor suggestion: leave blank (auto-assigned)
- Suggested reviewers: (optional — suggest 2-3 names)
  - Someone in computational statistics
  - Someone in econophysics or scaling laws
  - Someone in scientific Python tooling

## Suggested reviewers (search JOSS past papers)
Look for recent JOSS papers in:
  - scipy/statsmodels ecosystem
  - power law / heavy tail analysis
  - Any paper citing powerlaw package by Alstott et al.

## Expected timeline
- Initial screening: 1-2 weeks
- Review: 4-8 weeks
- Minor revisions: 1-2 weeks
- Total: ~3 months

## What reviewers will check
1. Software installs and runs
2. Tests pass
3. Paper accurately describes the software
4. Statement of need is convincing
5. Examples reproduce
6. License is present
