param([string]$OutputDirectory = (Join-Path (Split-Path -Parent $PSScriptRoot) 'excel-results-r15-native'))
$ErrorActionPreference = 'Stop'
$packet = $PSScriptRoot
$null = New-Item -ItemType Directory -Force -Path $OutputDirectory
$OutputDirectory = (Resolve-Path -LiteralPath $OutputDirectory).Path
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$exe = Join-Path $OutputDirectory ("raw-x87-probe-{0}.exe" -f $stamp)
& rustc --edition=2024 -O (Join-Path $packet 'raw-x87-probe.rs') -o $exe
if ($LASTEXITCODE -ne 0) { throw 'Native probe compilation failed' }
$result = Join-Path $OutputDirectory ("native-results-{0}.tsv" -f $stamp)
& $exe (Join-Path $packet 'native-leaf-requests.tsv') $result
if ($LASTEXITCODE -ne 0) { throw 'Native probe failed: preserve the partial result as a failure' }
Write-Host ("Native CPU results: {0}" -f $result)
Write-Host 'Return this file with the Excel readings. It contains no Excel answers.'
