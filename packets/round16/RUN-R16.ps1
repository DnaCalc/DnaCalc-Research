param(
    [Parameter(Mandatory=$true)][string]$DnaCalcCheckout,
    [string]$EvidenceCheckout = '',
    [string]$OutputDirectory = (Join-Path (Split-Path -Parent $PSScriptRoot) 'excel-results-r16')
)
$ErrorActionPreference = 'Stop'
function Write-Json([string]$Path, $Value) {
    $encoding = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($Path, (($Value | ConvertTo-Json -Depth 20) + "`n"), $encoding)
}
function Invoke-RecordedCargo([string[]]$Arguments, [string]$Log) {
    [void](Get-Command cargo -ErrorAction Stop)
    $preference = $ErrorActionPreference
    try {
        # Windows PowerShell 5.1 turns redirected native stderr into ErrorRecords.
        # Cargo's ordinary progress on stderr must not abort a successful command.
        $ErrorActionPreference = 'Continue'
        & cargo @Arguments 2>&1 | Tee-Object -FilePath $Log -ErrorAction Stop
        $code = $LASTEXITCODE
    } finally { $ErrorActionPreference = $preference }
    if ($code -ne 0) { throw "cargo $($Arguments -join ' ') exited $code; raw outputs are preserved." }
}
$packet = $PSScriptRoot
$checkout = (Resolve-Path -LiteralPath $DnaCalcCheckout).Path
if ($EvidenceCheckout -eq '') { $EvidenceCheckout = Join-Path (Split-Path -Parent $checkout) 'OxFunc-Evidence' }
$evidence = (Resolve-Path -LiteralPath $EvidenceCheckout).Path
if ([System.IO.Path]::GetFileName($evidence) -ine 'OxFunc-Evidence') {
    throw 'EvidenceCheckout must name the OxFunc-Evidence directory; xtask uses that repository name.'
}
$manifest = Get-Content -LiteralPath (Join-Path $packet 'manifest.json') -Raw | ConvertFrom-Json
foreach ($entry in $manifest.files_sha256.PSObject.Properties) {
    $got = (Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $packet $entry.Name)).Hash.ToLowerInvariant()
    if ($got -ne $entry.Value) { throw "Packet hash differs: $($entry.Name)" }
}
$section = Join-Path $checkout 'crates/calc-func/tests/fixtures/catalogue/sheet-dcj9k_r16.txt'
if ((Get-FileHash -Algorithm SHA256 -LiteralPath $section).Hash.ToLowerInvariant() -ne $manifest.catalogue_section_sha256) {
    throw 'The checkout does not contain this R16 catalogue section. Pull the published master first.'
}
$outputFull = [System.IO.Path]::GetFullPath($OutputDirectory)
foreach ($protected in @($packet, $checkout, $evidence)) {
    if ($outputFull -ieq $protected -or $outputFull.StartsWith($protected + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw 'Use an output directory outside the packet, application checkout and evidence checkout.'
    }
}
if (Test-Path -LiteralPath $outputFull) { throw 'Use a new output directory; preserve earlier attempts.' }
[void](New-Item -ItemType Directory -Path $outputFull)
$output = (Resolve-Path -LiteralPath $outputFull).Path
Copy-Item -LiteralPath $packet -Destination (Join-Path $output 'packet') -Recurse
$packetHashes = [ordered]@{}
Get-ChildItem -LiteralPath $packet -File -Recurse | ForEach-Object {
    $relative = $_.FullName.Substring($packet.Length + 1).Replace('\', '/')
    $packetHashes[$relative] = (Get-FileHash -Algorithm SHA256 -LiteralPath $_.FullName).Hash.ToLowerInvariant()
}
Write-Json (Join-Path $output 'packet-sha256.json') $packetHashes
$wanted = @{}
foreach ($batch in (Get-Content -LiteralPath (Join-Path $packet 'batches.json') -Raw | ConvertFrom-Json)) {
    $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $packet $batch.path)).Hash.ToLowerInvariant()
    $wanted[$hash] = $batch
}
$captureRoot = Join-Path $checkout 'target/c'
$captureBefore = @{}
if (Test-Path -LiteralPath $captureRoot) {
    Get-ChildItem -LiteralPath $captureRoot -Directory | ForEach-Object { $captureBefore[$_.FullName] = $true }
}
$readingRoot = Join-Path $evidence 'readings'
$readerRoot = Join-Path $checkout 'target/witness/read'
$frozenRoot = Join-Path $checkout 'target/witness/frozen'
$readingsBefore = @{}; $readerBefore = @{}; $frozenBefore = @{}
if (Test-Path -LiteralPath $readingRoot) {
    Get-ChildItem -LiteralPath $readingRoot -File -Recurse | ForEach-Object { $readingsBefore[$_.FullName] = $true }
}
if (Test-Path -LiteralPath $readerRoot) {
    Get-ChildItem -LiteralPath $readerRoot -Directory | ForEach-Object { $readerBefore[$_.FullName] = $true }
}
if (Test-Path -LiteralPath $frozenRoot) {
    Get-ChildItem -LiteralPath $frozenRoot -File | ForEach-Object { $frozenBefore[$_.FullName] = $true }
}
$oldLegacyRoot = $env:DNACALC_LEGACY_ROOT
$failure = $null; $found = @{}; $duplicates = @()
try {
    $env:DNACALC_LEGACY_ROOT = Split-Path -Parent $evidence
    Push-Location $checkout
    try {
        Invoke-RecordedCargo @('xtask', 'capture', '--pending', '--list') (Join-Path $output 'pending-before.txt')
        Invoke-RecordedCargo @('xtask', 'capture', '--pending', '--timeout', '900', '--attempts', '3') (Join-Path $output 'capture-console.txt')
        foreach ($batch in (Get-Content -LiteralPath (Join-Path $packet 'batches.json') -Raw | ConvertFrom-Json)) {
            Invoke-RecordedCargo @('xtask', 'witness', 'read', '--func', $batch.function, '--inputs', (Join-Path $packet $batch.path), '--timeout', '900') `
                (Join-Path $output ('read-' + $batch.function + '.txt'))
        }
    } finally { Pop-Location }
} catch { $failure = $_ }
finally {
    if ($null -eq $oldLegacyRoot) { Remove-Item Env:DNACALC_LEGACY_ROOT -ErrorAction SilentlyContinue }
    else { $env:DNACALC_LEGACY_ROOT = $oldLegacyRoot }
    foreach ($dir in @('catalogue-sections', 'worker/capture', 'worker/frozen', 'worker/read', 'OxFunc-Evidence/readings')) {
        [void](New-Item -ItemType Directory -Path (Join-Path $output $dir) -Force)
    }
    foreach ($name in @('dcj9k_r16', 'prefix_metadata', 'dcuar_cell')) {
        $source = Join-Path $checkout ('crates/calc-func/tests/captures/catalogue/sections/' + $name + '.capture.json')
        if (Test-Path -LiteralPath $source) { Copy-Item -LiteralPath $source -Destination (Join-Path $output 'catalogue-sections') }
    }
    if (Test-Path -LiteralPath $captureRoot) {
        Get-ChildItem -LiteralPath $captureRoot -Directory | Where-Object { -not $captureBefore.ContainsKey($_.FullName) } |
            ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $output 'worker/capture') -Recurse }
    }
    if (Test-Path -LiteralPath $readerRoot) {
        Get-ChildItem -LiteralPath $readerRoot -Directory | Where-Object { -not $readerBefore.ContainsKey($_.FullName) } |
            ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $output 'worker/read') -Recurse }
    }
    if (Test-Path -LiteralPath $frozenRoot) {
        Get-ChildItem -LiteralPath $frozenRoot -File | Where-Object { -not $frozenBefore.ContainsKey($_.FullName) } |
            ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $output 'worker/frozen') }
    }
    if (Test-Path -LiteralPath $readingRoot) {
        Get-ChildItem -LiteralPath $readingRoot -Filter '*.json' -File -Recurse |
            Where-Object { -not $readingsBefore.ContainsKey($_.FullName) } | ForEach-Object {
                $doc = Get-Content -LiteralPath $_.FullName -Raw | ConvertFrom-Json
                if ($null -eq $doc.capture_provenance -or $null -eq $doc.capture_provenance.frozen_model_predictions) { return }
                $hash = [string]$doc.capture_provenance.frozen_model_predictions.inputs_sha256
                if (-not $wanted.ContainsKey($hash)) { return }
                if ($found.ContainsKey($hash)) { $duplicates += $hash }
                $found[$hash] = $_.FullName
                $relative = $_.FullName.Substring($readingRoot.Length + 1)
                $destination = Join-Path (Join-Path $output 'OxFunc-Evidence/readings') $relative
                [void](New-Item -ItemType Directory -Path (Split-Path -Parent $destination) -Force)
                Copy-Item -LiteralPath $_.FullName -Destination $destination
            }
    }
    Write-Json (Join-Path $output 'run-status.json') ([ordered]@{
        failed = ($null -ne $failure -or $found.Count -ne $wanted.Count -or $duplicates.Count -gt 0)
        command_failure = $(if ($null -eq $failure) { $null } else { $failure.ToString() })
        requested_scalar_calls = $manifest.scalar_calls
        canonical_batches = $found.Count
        requested_batches = $wanted.Count
        duplicate_fresh_batches = @($duplicates)
        limits = 'Partial attempts and unchanged section copies are not automatically new semantic readings. Verify input hashes, section keys, raw results and settings before intake.'
    })
    $hashes = [ordered]@{}
    Get-ChildItem -LiteralPath $output -File -Recurse | ForEach-Object {
        $hashes[$_.FullName.Substring($output.Length + 1).Replace('\', '/')] =
            (Get-FileHash -Algorithm SHA256 -LiteralPath $_.FullName).Hash.ToLowerInvariant()
    }
    Write-Json (Join-Path $output 'results-sha256.json') $hashes
}
if ($null -ne $failure) { throw $failure }
if ($duplicates.Count -gt 0) { throw 'Duplicate fresh reader batches; both attempts are preserved.' }
if ($found.Count -ne $wanted.Count) { throw "Missing R16 reader batches: $($found.Count) of $($wanted.Count)." }
Write-Output "R16 preserved in $output. Predictions remain separate from Excel answers."
