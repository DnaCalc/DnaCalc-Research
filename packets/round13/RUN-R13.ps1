param([Parameter(Mandatory=$true)][string]$DnaCalcCheckout)
$ErrorActionPreference = 'Stop'
$packet = $PSScriptRoot
Push-Location $DnaCalcCheckout
try {
    $batches = Get-Content -LiteralPath (Join-Path $packet 'batches.json') -Raw | ConvertFrom-Json
    foreach ($batch in $batches) {
        $inputFile = Join-Path $packet $batch.path
        Write-Host ("R13: {0} / {1}: {2} calls" -f $batch.origin, $batch.function, $batch.calls)
        & cargo xtask witness read --func $batch.function --inputs $inputFile
        if ($LASTEXITCODE -ne 0) { throw "Reading failed: $($batch.path). Preserve all raw worker output before continuing." }
    }
} finally { Pop-Location }
