param(
    [string]$IndexCsv = "pos_way_admet_benchmark\data\remaining_sources_3prop_2pos\pubchem_bioassay_csv_index.csv",
    [string]$OutDir = "pos_way_admet_benchmark\raw\public\pubchem\bioassay_csv",
    [double]$MaxShardMB = 1.0,
    [double]$MaxTotalMB = 50.0,
    [int]$MaxShards = 100,
    [string]$BaseUrl = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data",
    [int]$Retries = 5,
    [int]$TimeoutSec = 3600,
    [string]$Downloader = "curl"
)

$ErrorActionPreference = "Stop"

function Convert-SizeToMB {
    param([string]$SizeText)
    if ($null -eq $SizeText) {
        $text = ""
    }
    else {
        $text = $SizeText.Trim().ToUpperInvariant()
    }
    if ($text -match '^([0-9.]+)\s*([KMG])$') {
        $value = [double]$Matches[1]
        $unit = $Matches[2]
        if ($unit -eq "K") { return $value / 1024.0 }
        if ($unit -eq "M") { return $value }
        if ($unit -eq "G") { return $value * 1024.0 }
    }
    if ($text -match '^([0-9.]+)$') {
        return ([double]$Matches[1]) / 1024.0 / 1024.0
    }
    return 0.0
}

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$rows = Import-Csv -Path $IndexCsv
$selected = New-Object System.Collections.Generic.List[object]
$total = 0.0
foreach ($row in $rows) {
    $sizeMB = Convert-SizeToMB $row.size
    if ($sizeMB -le 0 -or $sizeMB -gt $MaxShardMB) {
        continue
    }
    $target = Join-Path $OutDir $row.name
    if (Test-Path $target) {
        continue
    }
    if ($selected.Count -ge $MaxShards) {
        break
    }
    if (($total + $sizeMB) -gt $MaxTotalMB) {
        break
    }
    $selected.Add([pscustomobject]@{
        name = $row.name
        size_mb = [math]::Round($sizeMB, 4)
        url = "$BaseUrl/$($row.name)"
        target = $target
    }) | Out-Null
    $total += $sizeMB
}

$downloaded = 0
$failed = 0
foreach ($item in $selected) {
    $ok = $false
    for ($attempt = 1; $attempt -le $Retries; $attempt++) {
        try {
            Write-Host "[download] $($item.name) $($item.size_mb) MB attempt=$attempt"
            $part = "$($item.target).part"
            if ($Downloader -eq "curl" -and (Get-Command curl.exe -ErrorAction SilentlyContinue)) {
                & curl.exe --fail --location --silent --show-error --retry 3 --retry-delay 3 --connect-timeout 60 --max-time $TimeoutSec --continue-at - --output $part $item.url
                if ($LASTEXITCODE -ne 0) {
                    throw "curl.exe failed with exit code $LASTEXITCODE"
                }
                Move-Item -Force -Path $part -Destination $item.target
            }
            else {
                Invoke-WebRequest -Uri $item.url -OutFile $part -TimeoutSec $TimeoutSec
                Move-Item -Force -Path $part -Destination $item.target
            }
            $ok = $true
            break
        }
        catch {
            Write-Warning "[retry] $($item.name): $($_.Exception.Message)"
            Start-Sleep -Seconds ([math]::Min(10, $attempt * 2))
        }
    }
    if ($ok) {
        $downloaded += 1
    }
    else {
        $failed += 1
    }
}

$manifest = [pscustomobject]@{
    index_csv = $IndexCsv
    out_dir = $OutDir
    max_shard_mb = $MaxShardMB
    max_total_mb = $MaxTotalMB
    max_shards = $MaxShards
    selected_shards = $selected.Count
    selected_total_mb = [math]::Round($total, 4)
    downloaded = $downloaded
    failed = $failed
    completed_at = (Get-Date).ToString("s")
}
$manifestPath = Join-Path $OutDir "download_pubchem_bioassay_shards_manifest.json"
$manifest | ConvertTo-Json -Depth 4 | Set-Content -Path $manifestPath -Encoding UTF8
$manifest | ConvertTo-Json -Depth 4
