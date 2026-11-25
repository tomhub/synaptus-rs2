param(
    [int]$Warmups = 2,
    [int]$Repeats = 9
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location $repoRoot

foreach ($tool in @('cargo', 'matlab')) {
    if (-not (Get-Command $tool -ErrorAction SilentlyContinue)) {
        throw "Required command '$tool' was not found on PATH. Install it and reopen PowerShell."
    }
}

$pythonLauncher = Get-Command 'py' -ErrorAction SilentlyContinue
if (-not $pythonLauncher) {
    $pythonLauncher = Get-Command 'python' -ErrorAction SilentlyContinue
}
if (-not $pythonLauncher) {
    throw 'Python 3 was not found. Install Python and make py.exe or python.exe available on PATH.'
}

$venvPython = Join-Path $repoRoot '.venv\Scripts\python.exe'
if (-not (Test-Path $venvPython)) {
    if ($pythonLauncher.Name -eq 'py.exe') {
        & $pythonLauncher.Source -3 -m venv (Join-Path $repoRoot '.venv')
    } else {
        & $pythonLauncher.Source -m venv (Join-Path $repoRoot '.venv')
    }
    if ($LASTEXITCODE -ne 0) { throw 'Creating the Python virtual environment failed.' }
}

$pythonVersion = & $venvPython -c "import sys; print('%d.%d' % sys.version_info[:2])"
if ($LASTEXITCODE -ne 0) { throw 'Could not determine the virtual environment Python version.' }
$versionParts = $pythonVersion -split '\.'
if ([int]$versionParts[0] -lt 3 -or ([int]$versionParts[0] -eq 3 -and [int]$versionParts[1] -lt 11)) {
    throw "Python 3.11 or newer is required for the pinned NumPy version; found $pythonVersion."
}

& $venvPython -m pip install -r (Join-Path $repoRoot 'requirements-benchmark.txt')
if ($LASTEXITCODE -ne 0) { throw 'Installing the benchmark Python dependency failed.' }

& $venvPython -I (Join-Path $repoRoot 'benchmarks\run_benchmark.py') `
    --engine matlab --warmups $Warmups --repeats $Repeats
if ($LASTEXITCODE -ne 0) { throw 'The benchmark did not complete successfully.' }

Write-Host ''
Write-Host "Report: $(Join-Path $repoRoot 'benchmarks\performance_report.html')"
Write-Host "Raw results: $(Join-Path $repoRoot 'benchmarks\benchmark_results.json')"