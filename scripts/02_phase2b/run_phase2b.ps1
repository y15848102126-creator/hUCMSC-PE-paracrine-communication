param(
    [string]$PythonPath = "python",
    [string]$RscriptPath = "Rscript"
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$priorRLibs = $env:R_LIBS_USER
$phase2bRLib = Join-Path $repo "data\interim\phase1a1\Rlib"
if (Test-Path -LiteralPath $phase2bRLib) { $env:R_LIBS_USER = $phase2bRLib }
Push-Location $repo
try {
    & $PythonPath "scripts/02_phase2b/freeze_phase2b_hypotheses.py"
    if ($LASTEXITCODE -ne 0) { throw "Phase 2B hypothesis freeze failed" }
    & $RscriptPath "scripts/02_phase2b/run_phase2b_analysis.R" $repo
    if ($LASTEXITCODE -ne 0) { throw "Phase 2B analysis failed" }
    & $RscriptPath "scripts/02_phase2b/build_phase2b_figures.R" $repo
    if ($LASTEXITCODE -ne 0) { throw "Phase 2B figure generation failed" }
    & $PythonPath "scripts/02_phase2b/validate_phase2b_outputs.py"
    if ($LASTEXITCODE -ne 0) { throw "Phase 2B structural validation failed" }
    & $RscriptPath "scripts/02_phase2b/validate_phase2b_meta.R" $repo
    if ($LASTEXITCODE -ne 0) { throw "Phase 2B REML validation failed" }
}
finally {
    Pop-Location
    $env:R_LIBS_USER = $priorRLibs
}
