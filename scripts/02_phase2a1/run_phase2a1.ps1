param(
  [string]$PythonPath = "python",
  [string]$RscriptPath = "Rscript",
  [string]$CscPath = "C:\Windows\Microsoft.NET\Framework64\v4.0.30319\csc.exe",
  [switch]$SkipDownloads
)
$ErrorActionPreference = "Stop"
$root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Set-Location $root
$zip = "data\raw\phase0b\sc_PE_allcells_with_metadata_29-May-2023.txt.zip"
if (-not (Test-Path -LiteralPath $zip)) { throw "Missing audited Admati Figshare ZIP: $zip" }
if (-not (Test-Path -LiteralPath "data\interim\phase2a\admati_harmonized_pseudobulk_counts.csv")) {
  throw "Missing frozen Phase 2A pseudobulk cache; reproduce Phase 2A first"
}
if (-not $SkipDownloads) { & $PythonPath scripts/02_phase2a1/download_phase2a1_sources.py }
$interim = "data\interim\phase2a1"
New-Item -ItemType Directory -Force -Path $interim | Out-Null
& $CscPath /nologo /optimize+ /r:System.IO.Compression.dll /r:System.IO.Compression.FileSystem.dll /out:.\data\interim\phase2a1\AuditAdmatiCountLayer.exe .\scripts\02_phase2a1\AuditAdmatiCountLayer.cs
& .\data\interim\phase2a1\AuditAdmatiCountLayer.exe $zip data\interim\phase2a1\admati_count_layer_cell_audit.csv
$env:R_LIBS_USER = Join-Path $root "data\interim\phase1a1\Rlib"
& $RscriptPath scripts/02_phase2a1/run_phase2a1_analysis.R
& $PythonPath scripts/02_phase2a1/validate_phase2a1_outputs.py
