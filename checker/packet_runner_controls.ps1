# Constructed control-flow tests for RUN-R16.ps1. This starts no Excel, creates
# no reading, and uses only an isolated synthetic application/evidence tree.
param(
    [Parameter(Mandatory=$true)][string]$Packet,
    [Parameter(Mandatory=$true)][string]$OutputDirectory
)
$ErrorActionPreference = 'Stop'
if (Test-Path -LiteralPath $OutputDirectory) { throw 'Use a fresh synthetic controls directory.' }
[void](New-Item -ItemType Directory -Path $OutputDirectory)
$root = (Resolve-Path -LiteralPath $OutputDirectory).Path
$packetPath = (Resolve-Path -LiteralPath $Packet).Path
$runner = Join-Path $packetPath 'RUN-R16.ps1'
$tokens = $null; $parseErrors = $null
[void][System.Management.Automation.Language.Parser]::ParseFile($runner, [ref]$tokens, [ref]$parseErrors)
if ($parseErrors.Count -gt 0) { throw ($parseErrors | Out-String) }
$global:R16MockCheckout = Join-Path $root 'application'
$global:R16MockEvidence = Join-Path $root 'override/OxFunc-Evidence'
[void](New-Item -ItemType Directory -Path (Join-Path $global:R16MockCheckout 'crates/calc-func/tests/fixtures/catalogue') -Force)
[void](New-Item -ItemType Directory -Path (Join-Path $global:R16MockEvidence 'readings/excel-20527') -Force)
Copy-Item -LiteralPath (Join-Path $packetPath 'sheet-dcj9k_r16.txt') `
    -Destination (Join-Path $global:R16MockCheckout 'crates/calc-func/tests/fixtures/catalogue/sheet-dcj9k_r16.txt')
$global:R16MockAttempt = 0
$global:R16MockFail = $false
$global:R16MockFailCapture = $false
$global:R16MockDuplicate = $false
$global:R16MockCommands = 0
function global:cargo {
    $arguments = @($args)
    if ($env:DNACALC_LEGACY_ROOT -ne (Split-Path -Parent $global:R16MockEvidence)) {
        throw 'Evidence override was not passed to xtask.'
    }
    $global:R16MockCommands += 1
    # This error-stream progress must not abort the runner; it models the
    # redirected ErrorRecord PowerShell 5.1 creates for Cargo stderr.
    Write-Error 'SYNTHETIC ordinary native stderr progress'
    $global:LASTEXITCODE = 0
    if ($arguments -contains '--list') { return }
    if ($arguments -contains 'capture') {
        $work = Join-Path $global:R16MockCheckout ('target/c/attempt-' + $global:R16MockAttempt)
        [void](New-Item -ItemType Directory -Path $work -Force)
        [System.IO.File]::WriteAllText((Join-Path $work 'result.json'), '{"synthetic":true}')
        $sections = Join-Path $global:R16MockCheckout 'crates/calc-func/tests/captures/catalogue/sections'
        [void](New-Item -ItemType Directory -Path $sections -Force)
        foreach ($name in @('dcj9k_r16','prefix_metadata','dcuar_cell')) {
            [System.IO.File]::WriteAllText((Join-Path $sections ($name + '.capture.json')), '{"synthetic":true}')
        }
        if ($global:R16MockFailCapture) { $global:LASTEXITCODE = 19 }
        return
    }
    $inputPath = $null
    for ($i = 0; $i -lt $arguments.Count - 1; $i++) {
        if ($arguments[$i] -eq '--inputs') { $inputPath = $arguments[$i + 1] }
    }
    if ($null -eq $inputPath) { throw 'Unexpected synthetic Cargo command.' }
    $doc = Get-Content -LiteralPath $inputPath -Raw | ConvertFrom-Json
    $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $inputPath).Hash.ToLowerInvariant()
    $stem = $doc.function + '-' + $global:R16MockAttempt
    $work = Join-Path $global:R16MockCheckout ('target/witness/read/' + $stem + '-1')
    [void](New-Item -ItemType Directory -Path $work -Force)
    [System.IO.File]::WriteAllText((Join-Path $work 'request.json'), '{"synthetic":true}')
    [System.IO.File]::WriteAllText((Join-Path $work 'result.json'), '{"synthetic":true}')
    $frozen = Join-Path $global:R16MockCheckout ('target/witness/frozen/' + $stem + '.json')
    [void](New-Item -ItemType Directory -Path (Split-Path -Parent $frozen) -Force)
    [System.IO.File]::WriteAllText($frozen, '{"synthetic":true}')
    if ($global:R16MockFail) { $global:LASTEXITCODE = 17; return }
    $reading = [ordered]@{
        synthetic = $true
        capture_provenance = [ordered]@{
            frozen_model_predictions = [ordered]@{
                inputs_sha256 = $hash
                frozen_path = $frozen
            }
        }
    }
    $encoding = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText((Join-Path $global:R16MockEvidence ('readings/excel-20527/' + $stem + '.json')),
        ($reading | ConvertTo-Json -Depth 10), $encoding)
    if ($global:R16MockDuplicate -and $doc.function -eq 'VAR.S') {
        [System.IO.File]::WriteAllText((Join-Path $global:R16MockEvidence ('readings/excel-20527/' + $stem + '-duplicate.json')),
            ($reading | ConvertTo-Json -Depth 10), $encoding)
    }
}
$oldLegacyRoot = $env:DNACALC_LEGACY_ROOT
$env:DNACALC_LEGACY_ROOT = 'SYNTHETIC original root'
$results = @()
try {
    foreach ($attempt in @(1,2)) {
        $global:R16MockAttempt = $attempt
        $output = Join-Path $root ('success-' + $attempt)
        & $runner -DnaCalcCheckout $global:R16MockCheckout -EvidenceCheckout $global:R16MockEvidence -OutputDirectory $output
        $status = Get-Content -LiteralPath (Join-Path $output 'run-status.json') -Raw | ConvertFrom-Json
        if ($status.failed -or $status.canonical_batches -ne 13) { throw 'Synthetic success/repeat lost a batch.' }
        if ($env:DNACALC_LEGACY_ROOT -ne 'SYNTHETIC original root') { throw 'Legacy-root environment was not restored.' }
        if (@(Get-ChildItem -LiteralPath (Join-Path $output 'worker/read') -Directory).Count -ne 13) { throw 'Fresh reader directories were not isolated.' }
        $results += [ordered]@{ synthetic=$true; control=('success/repeat-' + $attempt); passed=$true }
    }
    $global:R16MockAttempt = 3; $global:R16MockFail = $true
    $output = Join-Path $root 'failed-reader'
    $rejected = $false
    try { & $runner -DnaCalcCheckout $global:R16MockCheckout -EvidenceCheckout $global:R16MockEvidence -OutputDirectory $output }
    catch { if ($_.ToString() -notmatch 'exited 17') { throw }; $rejected=$true }
    if (-not $rejected) { throw 'Failed reader was accepted.' }
    $status = Get-Content -LiteralPath (Join-Path $output 'run-status.json') -Raw | ConvertFrom-Json
    if (-not $status.failed -or $status.canonical_batches -ne 0) { throw 'Partial reading was misreported.' }
    if (@(Get-ChildItem -LiteralPath (Join-Path $output 'worker/read') -Directory).Count -ne 1) { throw 'Failed raw worker was not retained.' }
    if (@(Get-ChildItem -LiteralPath (Join-Path $output 'worker/frozen') -File).Count -ne 1) { throw 'Failed frozen prediction was not retained.' }
    if ($env:DNACALC_LEGACY_ROOT -ne 'SYNTHETIC original root') { throw 'Failure did not restore the environment.' }
    $results += [ordered]@{ synthetic=$true; control='failed-reader-raw-preservation'; passed=$true }
    $global:R16MockAttempt = 4; $global:R16MockFail = $false; $global:R16MockFailCapture = $true
    $output = Join-Path $root 'failed-capture'; $rejected=$false
    try { & $runner -DnaCalcCheckout $global:R16MockCheckout -EvidenceCheckout $global:R16MockEvidence -OutputDirectory $output }
    catch { if ($_.ToString() -notmatch 'exited 19') { throw }; $rejected=$true }
    if (-not $rejected) { throw 'Failed capture was accepted.' }
    $status = Get-Content -LiteralPath (Join-Path $output 'run-status.json') -Raw | ConvertFrom-Json
    if (-not $status.failed -or $status.canonical_batches -ne 0) { throw 'Failed capture was misreported.' }
    if (@(Get-ChildItem -LiteralPath (Join-Path $output 'worker/capture') -Directory).Count -ne 1) { throw 'Failed capture raw output was not retained.' }
    if ($env:DNACALC_LEGACY_ROOT -ne 'SYNTHETIC original root') { throw 'Capture failure did not restore the environment.' }
    $results += [ordered]@{ synthetic=$true; control='failed-capture-raw-preservation'; passed=$true }
    $global:R16MockAttempt = 5; $global:R16MockFailCapture = $false; $global:R16MockDuplicate = $true
    $output = Join-Path $root 'duplicate-fresh'; $rejected=$false
    try { & $runner -DnaCalcCheckout $global:R16MockCheckout -EvidenceCheckout $global:R16MockEvidence -OutputDirectory $output }
    catch { if ($_.ToString() -notmatch 'Duplicate fresh reader batches') { throw }; $rejected=$true }
    if (-not $rejected) { throw 'Duplicate fresh batch was accepted.' }
    $status = Get-Content -LiteralPath (Join-Path $output 'run-status.json') -Raw | ConvertFrom-Json
    if (-not $status.failed -or $status.duplicate_fresh_batches.Count -ne 1) { throw 'Duplicate fresh batch was misreported.' }
    if (@(Get-ChildItem -LiteralPath (Join-Path $output 'OxFunc-Evidence/readings') -File -Recurse).Count -ne 14) { throw 'Both duplicate readings were not retained.' }
    $results += [ordered]@{ synthetic=$true; control='duplicate-fresh-both-preserved'; passed=$true }
    $forbidden = Join-Path $packetPath 'forbidden-synthetic-output'
    $rejected=$false
    try { & $runner -DnaCalcCheckout $global:R16MockCheckout -EvidenceCheckout $global:R16MockEvidence -OutputDirectory $forbidden }
    catch { if ($_.ToString() -notmatch 'outside the packet') { throw }; $rejected=$true }
    if (-not $rejected -or (Test-Path -LiteralPath $forbidden)) { throw 'Recursive output path was not rejected.' }
    $results += [ordered]@{ synthetic=$true; control='recursive-output-rejected'; passed=$true }
} finally {
    if ($null -eq $oldLegacyRoot) { Remove-Item Env:DNACALC_LEGACY_ROOT -ErrorAction SilentlyContinue }
    else { $env:DNACALC_LEGACY_ROOT = $oldLegacyRoot }
}
$report = [ordered]@{
    synthetic = $true
    complete = $true
    controls = @($results)
    mocked_cargo_commands = $global:R16MockCommands
    powershell = $PSVersionTable.PSVersion.ToString()
    runner_sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $runner).Hash.ToLowerInvariant()
    limits = 'Control-flow and packaging tests on constructed files only. No Excel started; no synthetic file is an answer or a WitnessSet. Linux PowerShell does not certify Windows 5.1 runtime equivalence.'
}
$encoding = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText((Join-Path $root 'controls.json'), (($report | ConvertTo-Json -Depth 20) + "`n"), $encoding)
Write-Output ($report | ConvertTo-Json -Depth 20)
exit 0
