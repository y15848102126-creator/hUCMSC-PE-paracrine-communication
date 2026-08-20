param(
    [switch]$SkipProcessed,
    [switch]$Force,
    [string]$PythonPath = ""
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$python = $PythonPath
if (-not $python) {
    $pythonCommand = Get-Command python -ErrorAction SilentlyContinue
    if (-not $pythonCommand) {
        throw "Python >=3.11 is required; pass -PythonPath when python is not on PATH"
    }
    $python = $pythonCommand.Source
}
if (-not (Test-Path -LiteralPath $python)) {
    throw "Python executable not found: $python"
}

$downloadArgs = @((Join-Path $PSScriptRoot "download_phase0b_metadata.py"))
if (-not $SkipProcessed) { $downloadArgs += "--include-processed" }
if ($Force) { $downloadArgs += "--force" }
& $python @downloadArgs
& $python (Join-Path $PSScriptRoot "build_phase0b_audit.py")
& $python (Join-Path $PSScriptRoot "validate_phase0b_outputs.py")

Write-Host "Phase 0B audit complete: $repo"
