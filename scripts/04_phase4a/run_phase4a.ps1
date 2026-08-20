param(
    [string]$PythonPath = "python",
    [string]$RscriptPath = "Rscript",
    [switch]$SkipDownload,
    [switch]$SkipDelivery
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..\..")).Path
Push-Location $repo
try {
    if (-not $SkipDownload) {
        & $PythonPath "scripts/04_phase4a/download_phase4a_resources.py"
        if ($LASTEXITCODE -ne 0) { throw "Phase 4A resource download failed" }
    }
    & $RscriptPath "scripts/04_phase4a/prepare_nichenet_resources.R" $repo
    if ($LASTEXITCODE -ne 0) { throw "NicheNet resource preparation failed" }
    & $PythonPath "scripts/04_phase4a/freeze_phase4a_inputs.py"
    if ($LASTEXITCODE -ne 0) { throw "Phase 4A freeze failed" }
    & $PythonPath "scripts/04_phase4a/run_phase4a_analysis.py"
    if ($LASTEXITCODE -ne 0) { throw "Phase 4A analysis failed" }
    & $RscriptPath "scripts/04_phase4a/make_phase4a_figures.R" $repo
    if ($LASTEXITCODE -ne 0) { throw "Phase 4A figure generation failed" }
    & $PythonPath "scripts/04_phase4a/validate_phase4a.py"
    if ($LASTEXITCODE -ne 0) { throw "Phase 4A validation failed" }
    if (-not $SkipDelivery) {
        & $PythonPath "scripts/04_phase4a/create_phase4a_delivery.py"
        if ($LASTEXITCODE -ne 0) { throw "Phase 4A delivery failed" }
        & $PythonPath "scripts/04_phase4a/validate_phase4a.py"
        if ($LASTEXITCODE -ne 0) { throw "Phase 4A delivery validation failed" }
    }
} finally {
    Pop-Location
}
