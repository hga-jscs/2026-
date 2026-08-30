<#
统一入口：在已有公开数据上执行 AC、最大环状占比合并与自动一致性检查。
代码只使用 task2 内冻结的 AC 和算法源码；检查点只有在输入、代码、配置和运行参数指纹
全部一致时才复用。-ForceAC 与 -ForceLP 可分别强制重算。
#>

param(
    [string]$PythonExe = "D:\anaconda3\python.exe",
    [string]$WslPython = "/usr/bin/python3",
    [Parameter(Mandatory = $true)]
    [string]$AcRoot,
    [string]$DataRoot = "",
    [string]$ResultDir = "",
    [switch]$ForceAC,
    [switch]$ForceLP
)

$ErrorActionPreference = "Stop"
$codeRoot = $PSScriptRoot
$repositoryRoot = Split-Path -Parent $codeRoot

if (-not $DataRoot) {
    $DataRoot = Join-Path $repositoryRoot "data"
}
if (-not $ResultDir) {
    $ResultDir = Join-Path $repositoryRoot "results"
}

$datasetManifest = Join-Path $ResultDir "数据集清单.json"
$acSource = $AcRoot
$algorithmSource = Join-Path $repositoryRoot "source"

foreach ($required in @($PythonExe, $DataRoot, $datasetManifest, $acSource, $algorithmSource)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Required path does not exist: $required"
    }
}

New-Item -ItemType Directory -Force -Path $ResultDir | Out-Null

Write-Host "STEP 1/3  AC execution or fingerprint-verified reuse"
$acArguments = @(
    (Join-Path $codeRoot "run_ac_all.py"),
    "--dataset-manifest", $datasetManifest,
    "--data-root", $DataRoot,
    "--output-dir", $ResultDir,
    "--ac-source", $acSource,
    "--ac-python", $WslPython,
    "--jobs", "8",
    "--no-bfbarchitect"
)
if ($ForceAC) { $acArguments += "--force" }
& $PythonExe @acArguments
if ($LASTEXITCODE -ne 0) { throw "AC stage failed with exit code $LASTEXITCODE" }

Write-Host "STEP 2/3  Maximum cyclic ratio and four-column merge"
$lpArguments = @(
    (Join-Path $codeRoot "run_lwcn_and_merge.py"),
    "--dataset-manifest", $datasetManifest,
    "--data-root", $DataRoot,
    "--output-dir", $ResultDir,
    "--jobs", "8"
)
if ($ForceLP) { $lpArguments += "--force" }
& $PythonExe @lpArguments
if ($LASTEXITCODE -ne 0) { throw "Maximum cyclic ratio stage failed with exit code $LASTEXITCODE" }

Write-Host "STEP 3/3  Automated consistency check"
& $PythonExe (Join-Path $codeRoot "verify_all_results.py") `
    --data-root $DataRoot `
    --output-dir $ResultDir `
    --ac-source $acSource
if ($LASTEXITCODE -ne 0) { throw "Verification failed with exit code $LASTEXITCODE" }

Write-Host "TASK2_FULL_PIPELINE_PASSED"
