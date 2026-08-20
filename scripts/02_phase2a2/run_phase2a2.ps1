param(
  [string]$PythonPath = "",
  [string]$RscriptPath = "C:\Program Files\R\R-4.5.3\bin\Rscript.exe",
  [switch]$SkipMetadataDownload
)
$ErrorActionPreference = "Stop"
$root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Set-Location $root
if (-not $PythonPath) {
  $cmd = Get-Command python -ErrorAction SilentlyContinue
  if (-not $cmd) { throw "Python 3.12 is required; pass -PythonPath" }
  $PythonPath = $cmd.Source
}
if (-not (Test-Path -LiteralPath $RscriptPath)) { throw "Rscript not found: $RscriptPath" }
if (-not $SkipMetadataDownload) { & $PythonPath scripts/02_phase2a2/download_phase2a2_metadata.py }
& $PythonPath scripts/02_phase2a2/build_phase2a2_audits.py
$env:R_LIBS_USER = (Resolve-Path "data/interim/phase1a1/Rlib").Path
& $RscriptPath scripts/02_phase2a2/run_phase2a2_corrected_analysis.R $root
& $PythonPath scripts/02_phase2a2/finalize_phase2a2.py
& $RscriptPath scripts/02_phase2a2/build_phase2a2_figures.R $root
& $PythonPath scripts/02_phase2a2/validate_phase2a2_outputs.py
