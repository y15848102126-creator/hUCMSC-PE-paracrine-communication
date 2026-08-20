# DEPRECATED / HISTORICAL ONLY: this runner executes the superseded Phase 2A count-likelihood workflow. It is excluded from the default public workflow.
param(
  [string]$PythonPath = "python",
  [string]$RscriptPath = "Rscript",
  [string]$CscPath = "C:\Windows\Microsoft.NET\Framework64\v4.0.30319\csc.exe",
  [switch]$SkipDownloads,
  [switch]$SkipRPackageInstall
)
$ErrorActionPreference = "Stop"
$root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Set-Location $root
if (-not (Test-Path -LiteralPath "data\raw\phase0b\sc_PE_allcells_with_metadata_29-May-2023.txt.zip")) { throw "Missing audited Admati Figshare ZIP" }
$env:R_LIBS_USER = Join-Path $root "data\interim\phase1a1\Rlib"
if (-not $SkipRPackageInstall) { & $RscriptPath scripts/02_phase2a/install_phase2a_r_packages.R $root }
& $PythonPath scripts/02_phase2a/freeze_admati_metadata.py
if (-not $SkipDownloads) { & $PythonPath scripts/02_phase2a/download_phase2a_resources.py }
& $PythonPath scripts/02_phase2a/prepare_phase2a_aggregation.py
& $CscPath /nologo /optimize+ /r:System.IO.Compression.dll /r:System.IO.Compression.FileSystem.dll /out:.\data\interim\phase2a\AggregatePseudobulk.exe .\scripts\02_phase2a\AggregatePseudobulk.cs
& .\data\interim\phase2a\AggregatePseudobulk.exe data\raw\phase0b\sc_PE_allcells_with_metadata_29-May-2023.txt.zip data\interim\phase2a\cell_group_ids.int32 data\interim\phase2a\pseudobulk_strata.tsv data\interim\phase2a\admati_harmonized_pseudobulk_counts.csv.gz data\interim\phase2a\pseudobulk_observed_totals.tsv
& $PythonPath scripts/02_phase2a/finalize_phase2a_aggregation.py
& $RscriptPath scripts/02_phase2a/run_phase2a_analysis.R
& $RscriptPath scripts/02_phase2a/build_phase2a_figures.R
& $PythonPath scripts/02_phase2a/add_phase2a_provenance.py
& $PythonPath scripts/02_phase2a/build_phase2a_report.py
& $PythonPath scripts/02_phase2a/validate_phase2a_outputs.py
