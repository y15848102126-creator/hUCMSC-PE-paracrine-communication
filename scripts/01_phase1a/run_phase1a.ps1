# HISTORICAL FREEZE RUNNER: Phase 1A.1 supersedes the withdrawn GSE30186 cohort-minimum shift-log route for formal analysis.
param(
    [switch]$SkipDownload,
    [switch]$Force,
    [string]$PythonPath = ""
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$python = $PythonPath
if (-not $python) {
    $pythonCommand = Get-Command python -ErrorAction SilentlyContinue
    if (-not $pythonCommand) {
        throw "Python >=3.11 with numpy and pandas is required; pass -PythonPath when python is not on PATH"
    }
    $python = $pythonCommand.Source
}
if (-not (Test-Path -LiteralPath $python)) {
    throw "Python executable not found: $python"
}

if (-not $SkipDownload) {
    $downloadArgs = @((Join-Path $PSScriptRoot "download_phase1a_data.py"))
    if ($Force) { $downloadArgs += "--force" }
    & $python @downloadArgs
}
& $python (Join-Path $PSScriptRoot "build_phase1a_freeze.py")
& $python (Join-Path $PSScriptRoot "validate_phase1a_outputs.py")

Write-Host "Phase 1A bulk data freeze complete: $repo"
