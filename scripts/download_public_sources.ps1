param(
    [string]$Manifest = ".\pos_way_admet_benchmark\configs\download_manifest.json",
    [string]$OutputRoot = ".\pos_way_admet_benchmark\raw\public",
    [switch]$NoResume
)

$ErrorActionPreference = "Stop"

if (-not (Get-Command curl.exe -ErrorAction SilentlyContinue)) {
    throw "curl.exe is required but was not found."
}

$manifestPath = Resolve-Path -LiteralPath $Manifest
$manifestJson = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
$rootPath = Join-Path (Get-Location) $OutputRoot
New-Item -ItemType Directory -Force -Path $rootPath | Out-Null

$status = [System.Collections.Generic.List[object]]::new()

foreach ($file in $manifestJson.files) {
    $target = Join-Path $rootPath $file.path
    $targetDir = Split-Path -Parent $target
    New-Item -ItemType Directory -Force -Path $targetDir | Out-Null

    $alreadyExists = Test-Path -LiteralPath $target
    if ($alreadyExists -and ((Get-Item -LiteralPath $target).Length -gt 0)) {
        Write-Host "SKIP existing $($file.path)"
        $status.Add([pscustomobject]@{
            source = $file.source
            path = $file.path
            url = $file.url
            status = "skipped_existing"
            bytes = (Get-Item -LiteralPath $target).Length
        })
        continue
    }

    Write-Host "DOWNLOAD $($file.source) $($file.path)"
    $curlArgs = @(
        "--location",
        "--fail",
        "--retry", "5",
        "--retry-delay", "3",
        "--connect-timeout", "30",
        "--output", $target
    )
    if (-not $NoResume) {
        $curlArgs += @("--continue-at", "-")
    }
    $curlArgs += $file.url

    & curl.exe @curlArgs
    if ($LASTEXITCODE -ne 0) {
        Write-Warning "curl failed for $($file.url)"
        $status.Add([pscustomobject]@{
            source = $file.source
            path = $file.path
            url = $file.url
            status = "failed"
            bytes = $(if (Test-Path -LiteralPath $target) { (Get-Item -LiteralPath $target).Length } else { 0 })
        })
        continue
    }

    $status.Add([pscustomobject]@{
        source = $file.source
        path = $file.path
        url = $file.url
        status = "downloaded"
        bytes = (Get-Item -LiteralPath $target).Length
    })
}

$statusPath = Join-Path $rootPath "download_status.json"
$status | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding UTF8

$shaPath = Join-Path $rootPath "checksums_sha256.txt"
Get-ChildItem -LiteralPath $rootPath -Recurse -File |
    Where-Object { $_.Name -ne "checksums_sha256.txt" } |
    ForEach-Object {
        $hash = Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256
        $relative = Resolve-Path -LiteralPath $_.FullName -Relative
        "$($hash.Hash.ToLowerInvariant())  $relative"
    } |
    Set-Content -LiteralPath $shaPath -Encoding ASCII

Write-Host "DONE"
Write-Host "Status: $statusPath"
Write-Host "SHA256: $shaPath"
