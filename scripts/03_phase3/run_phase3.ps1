param(
    [string]$PythonPath = "python",
    [string]$RscriptPath = "Rscript",
    [switch]$SkipDownload
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Set-Location $RepoRoot
$RequiredSoft = @("GSE182158", "GSE199071", "GSE117837")
foreach ($Accession in $RequiredSoft) {
    $SoftPath = Join-Path $RepoRoot "data\raw\${Accession}_family.soft.gz"
    if (-not (Test-Path -LiteralPath $SoftPath)) {
        throw "Missing Phase 0 metadata prerequisite: $SoftPath"
    }
}
if (-not $SkipDownload) {
    & $PythonPath scripts/03_phase3/download_phase3_inputs.py
    if ($LASTEXITCODE -ne 0) { throw "Phase 3 download failed" }
}
$RepoRLib = Join-Path $RepoRoot "data\interim\phase1a1\Rlib"
$UserRLib = Join-Path $env:LOCALAPPDATA "R\win-library\4.5"
$env:R_LIBS_USER = "${RepoRLib};${UserRLib}"
& $RscriptPath scripts/03_phase3/freeze_phase3_sender_inputs.R $RepoRoot
if ($LASTEXITCODE -ne 0) { throw "Phase 3 freeze failed" }
& $RscriptPath scripts/03_phase3/run_phase3_sender_analysis.R $RepoRoot
if ($LASTEXITCODE -ne 0) { throw "Phase 3 sender analysis failed" }
& $RscriptPath scripts/03_phase3/build_phase3_figures.R $RepoRoot
if ($LASTEXITCODE -ne 0) { throw "Phase 3 figure build failed" }
& $PythonPath scripts/03_phase3/validate_phase3_outputs.py
if ($LASTEXITCODE -ne 0) { throw "Phase 3 validation failed" }
Write-Host "Phase 3 complete. Phase 4 was not started."
