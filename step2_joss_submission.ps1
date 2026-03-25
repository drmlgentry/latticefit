#!/usr/bin/env pwsh
# step2_joss_submission.ps1
# Verify paper.md is in repo root and open JOSS submission page

Set-Location "C:\dev\latticefit"

# Verify required files exist
$required = @("paper.md", "paper.bib", "README.md", "LICENSE")
$missing = @()
foreach ($f in $required) {
    if (Test-Path $f) {
        Write-Host "  OK: $f" -ForegroundColor Green
    } else {
        Write-Host "  MISSING: $f" -ForegroundColor Red
        $missing += $f
    }
}

if ($missing.Count -gt 0) {
    Write-Host "Fix missing files before submitting." -ForegroundColor Red
    exit 1
}

# Verify paper.md has required JOSS frontmatter
$paper = Get-Content "paper.md" -Raw
$checks = @{
    "title:"       = $paper -match "^title:"
    "authors:"     = $paper -match "authors:"
    "orcid:"       = $paper -match "orcid:"
    "affiliations:"= $paper -match "affiliations:"
    "bibliography:"= $paper -match "bibliography:"
    "# Summary"    = $paper -match "# Summary"
    "# Statement"  = $paper -match "# Statement of Need"
}

Write-Host ""
Write-Host "paper.md checks:" -ForegroundColor Cyan
foreach ($k in $checks.Keys) {
    $v = $checks[$k]
    Write-Host "  $(if($v){'OK'}else{'MISSING'}): $k" -ForegroundColor $(if($v){'Green'}else{'Red'})
}

Write-Host ""
Write-Host "JOSS SUBMISSION STEPS:" -ForegroundColor Yellow
Write-Host "  1. Go to: https://joss.theoj.org/papers/new"
Write-Host "  2. Repository URL: https://github.com/drmlgentry/latticefit"
Write-Host "  3. Branch: main"
Write-Host "  4. Paper path: paper.md"
Write-Host "  5. Software language: Python"
Write-Host "  6. Click 'Submit'"
Write-Host ""
Write-Host "  Expected timeline: 4-8 weeks for review"
Write-Host "  Reviewers will run: pip install latticefit"
Write-Host "  and: python validate_latticefit.py"

# Open JOSS submission page
Start-Process "https://joss.theoj.org/papers/new"
