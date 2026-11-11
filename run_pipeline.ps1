param(
    [Parameter(ValueFromRemainingArguments=$true)]
    [string[]]$ArgsList
)
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$env:PYTHONPATH = $root
$py = Join-Path $root "environments\brd_venv\Scripts\python.exe"
if (-not $ArgsList) { $ArgsList = @("--run", "all") }
& $py (Join-Path $root "BRD_pipeline.py") @ArgsList
