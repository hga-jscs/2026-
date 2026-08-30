<#
逐项准备 AC 参考数据：复用已校验归档，必要时断点续传，只解出 file_list 及其中声明的
mappability、gene、centromere 文件；任何缺项或 tar 错误都会立即停止。
#>

param(
    [Parameter(Mandatory = $true)]
    [string]$DataRoot,
    [string]$DatasetManifest = ""
)

$ErrorActionPreference = "Stop"
$referenceRoot = Join-Path $DataRoot "reference_data"
$archiveRoot = Join-Path $DataRoot "reference_archives"
New-Item -ItemType Directory -Path $referenceRoot -Force | Out-Null
New-Item -ItemType Directory -Path $archiveRoot -Force | Out-Null

$references = @("GRCh37", "hg19", "GRCh38", "GRCh38_viral", "mm10")
foreach ($reference in $references) {
    $readyFile = Join-Path $referenceRoot "$reference\file_list.txt"
    if (Test-Path -LiteralPath $readyFile -PathType Leaf) {
        Write-Host "REFERENCE_REUSED ref=$reference"
        continue
    }

    $archive = Join-Path $archiveRoot "$reference.tar.gz"
    $url = "https://refs.ampliconrepository.org/data/module_support_files/AmpliconArchitect/$reference.tar.gz"
    $archiveValid = $false
    if (Test-Path -LiteralPath $archive -PathType Leaf) {
        & tar.exe -tzf $archive | Out-Null
        $archiveValid = ($LASTEXITCODE -eq 0)
    }
    if (-not $archiveValid) {
        Write-Host "REFERENCE_DOWNLOAD ref=$reference url=$url"
        & curl.exe -fL --retry 12 --retry-all-errors --retry-delay 5 -C - -o $archive $url
        if ($LASTEXITCODE -ne 0) {
            throw "Reference download failed for $reference with exit code $LASTEXITCODE"
        }
    }

    $fileListMember = "$reference/file_list.txt"
    $fileListText = (& tar.exe -xOzf $archive $fileListMember) -join "`n"
    if ($LASTEXITCODE -ne 0) {
        throw "Could not read file_list.txt for $reference"
    }
    $fileMap = @{}
    foreach ($line in ($fileListText -split "`n")) {
        if ($line -match '^\s*(\S+)\s+(.+?)\s*$') {
            $fileMap[$Matches[1]] = $Matches[2]
        }
    }
    foreach ($requiredKey in @("mapability_exclude_filename", "gene_filename", "centromere_filename")) {
        if (-not $fileMap.ContainsKey($requiredKey)) {
            throw "Missing $requiredKey in $reference/file_list.txt"
        }
    }
    $wantedMembers = @(
        $fileListMember,
        "$reference/$($fileMap['mapability_exclude_filename'])",
        "$reference/$($fileMap['gene_filename'])",
        "$reference/$($fileMap['centromere_filename'])"
    )
    & tar.exe -xzf $archive -C $referenceRoot @wantedMembers
    if ($LASTEXITCODE -ne 0) {
        throw "Required-file extraction failed for $reference"
    }
    if (-not (Test-Path -LiteralPath $readyFile -PathType Leaf)) {
        throw "Reference extraction did not create $readyFile"
    }
    $bytes = (Get-Item -LiteralPath $archive).Length
    Write-Host "REFERENCE_READY ref=$reference archive_bytes=$bytes"
}

if (-not $DatasetManifest) {
    $DatasetManifest = Join-Path (Split-Path -Parent $PSScriptRoot) "结果\数据集清单.json"
}
if (-not (Test-Path -LiteralPath $DatasetManifest -PathType Leaf)) {
    throw "Dataset manifest does not exist: $DatasetManifest"
}
$requiredReferences = Get-Content -LiteralPath $DatasetManifest -Raw -Encoding UTF8 |
    ConvertFrom-Json |
    Where-Object { $_.status -eq "READY" } |
    ForEach-Object {
        if ($_.reference_genome -eq "hg38") { "GRCh38" }
        elseif ($_.reference_genome -eq "GRCm38") { "mm10" }
        else { [string]$_.reference_genome }
    } |
    Sort-Object -Unique
$missingReferences = @(
    $requiredReferences | Where-Object {
        -not (Test-Path -LiteralPath (Join-Path $referenceRoot "$_\file_list.txt") -PathType Leaf)
    }
)
if ($missingReferences.Count -gt 0) {
    throw "Reference set is incomplete for the manifest: $($missingReferences -join ', ')"
}

Write-Host "REFERENCES_COMPLETE count=$($references.Count) manifest_required=$($requiredReferences.Count)"
