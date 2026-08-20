param(
    [switch]$Offline,
    [string]$PythonExecutable = $env:PHASE0_PYTHON
)

$ErrorActionPreference = 'Stop'
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
Set-Location $RepoRoot

function Resolve-WorkingPython {
    param([string]$Requested)
    $Candidates = @()
    if ($Requested) { $Candidates += $Requested }
    $Candidates += @('python', 'python3', 'py')
    foreach ($Candidate in $Candidates) {
        try {
            if ($Candidate -eq 'py') {
                & $Candidate -3 --version *> $null
                if ($LASTEXITCODE -eq 0) { return [pscustomobject]@{ Exe = $Candidate; Args = @('-3') } }
            } else {
                & $Candidate --version *> $null
                if ($LASTEXITCODE -eq 0) { return [pscustomobject]@{ Exe = $Candidate; Args = @() } }
            }
        } catch { }
    }
    throw 'No working Python >=3.11 found. Set PHASE0_PYTHON or -PythonExecutable.'
}

$PythonCommand = Resolve-WorkingPython -Requested $PythonExecutable
$Python = $PythonCommand.Exe
$PythonPrefix = $PythonCommand.Args

New-Item -ItemType Directory -Force -Path 'logs' | Out-Null
$LogPath = Join-Path $RepoRoot ('logs\phase0_' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '.log')
Start-Transcript -Path $LogPath | Out-Null
try {
    if (-not $Offline) {
        & $Python @PythonPrefix 'scripts\00_dataset_registry\download_metadata.py'
        if ($LASTEXITCODE -ne 0) { throw 'metadata download failed' }
    }
    & $Python @PythonPrefix 'scripts\00_dataset_registry\build_audit.py'
    if ($LASTEXITCODE -ne 0) { throw 'audit generation failed' }
    & $Python @PythonPrefix 'scripts\00_dataset_registry\validate_audit.py'
    if ($LASTEXITCODE -ne 0) { throw 'audit validation failed' }
    Write-Host 'Phase 0 audit completed successfully.'
} finally {
    Stop-Transcript | Out-Null
}
