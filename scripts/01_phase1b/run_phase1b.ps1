param(
    [switch]$SkipRPackageInstall,
    [string]$PythonPath = "",
    [string]$RscriptPath = ""
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
if (-not $PythonPath) {
    $cmd = Get-Command python -ErrorAction SilentlyContinue
    if (-not $cmd) { throw "Python >=3.11 with numpy and pandas is required; pass -PythonPath" }
    $PythonPath = $cmd.Source
}
if (-not $RscriptPath) {
    $cmd = Get-Command Rscript -ErrorAction SilentlyContinue
    if (-not $cmd) { throw "R 4.5.x is required; pass -RscriptPath" }
    $RscriptPath = $cmd.Source
}
if (-not (Test-Path -LiteralPath $PythonPath)) { throw "Python executable not found: $PythonPath" }
if (-not (Test-Path -LiteralPath $RscriptPath)) { throw "Rscript executable not found: $RscriptPath" }

foreach ($required in @(
    "data\interim\phase1a\processed_matrices\GSE75010_BIOBANK_gene_level.tsv.gz",
    "data\interim\phase1a1\formal_matrices\GSE30186_formal_gene_level.tsv.gz",
    "results\01_phase1a1\formal_phase1b_matrix_registry.csv"
)) {
    if (-not (Test-Path -LiteralPath (Join-Path $repo $required))) { throw "Missing Phase 1A.1 prerequisite: $required" }
}

$env:R_LIBS_USER = Join-Path $repo "data\interim\phase1a1\Rlib"
if (-not $SkipRPackageInstall) { & $RscriptPath (Join-Path $PSScriptRoot "install_phase1b_r_packages.R") $repo }
& $RscriptPath (Join-Path $PSScriptRoot "run_phase1b_analysis.R") $repo
& $PythonPath (Join-Path $PSScriptRoot "build_phase1b_report.py")
& $PythonPath (Join-Path $PSScriptRoot "validate_phase1b_outputs.py")
Write-Host "Phase 1B PE disease-signature analysis complete: $repo"
