param(
    [switch]$SkipDownload,
    [switch]$SkipRPackageInstall,
    [switch]$Force,
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

if (-not $SkipDownload) {
    $args = @((Join-Path $PSScriptRoot "download_phase1a1_data.py"))
    if ($Force) { $args += "--force" }
    & $PythonPath @args
}
& $PythonPath (Join-Path $PSScriptRoot "prepare_phase1a1_inputs.py")

$rLib = Join-Path $repo "data\interim\phase1a1\Rlib"
$env:R_LIBS_USER = $rLib
if (-not $SkipRPackageInstall) {
    & $RscriptPath (Join-Path $PSScriptRoot "install_phase1a1_r_packages.R") $repo
}
& $RscriptPath (Join-Path $PSScriptRoot "preprocess_phase1a1.R") $repo
& $PythonPath (Join-Path $PSScriptRoot "build_phase1a1_amendment.py")
& $PythonPath (Join-Path $PSScriptRoot "validate_phase1a1_outputs.py")

Write-Host "Phase 1A.1 preprocessing amendment complete: $repo"
