<#
统一入口：在已有公开数据和检查点上重新执行 AC 调度、最大环状占比合并与独立验收。
默认不重复下载约 8 GB 归档，也不重算已完整的 28,142 个检查点；传入 -ForceLP 时只强制
重算 LP。AC 的缺失/不完整项目会由 run_ac_all.py 自动继续。
#>

param(
    [string]$PythonExe = "D:\anaconda3\python.exe",
    [string]$DataRoot = "",
    [string]$ResultDir = "",
    [switch]$ForceLP
)

$ErrorActionPreference = "Stop"
$codeRoot = $PSScriptRoot
$task2Root = Split-Path -Parent (Split-Path -Parent $codeRoot)
$exchangeRoot = Split-Path -Parent $task2Root
$researchRoot = Split-Path -Parent $exchangeRoot

if (-not $DataRoot) {
    $DataRoot = Join-Path $exchangeRoot ".task1_data"
}
if (-not $ResultDir) {
    $ResultDir = Join-Path (Split-Path -Parent $codeRoot) "结果"
}

$datasetManifest = Join-Path $ResultDir "数据集清单.json"
$acSource = Join-Path $task2Root "02_AC判定论文与流程图\AmpliconClassifier_注释版"
$algorithmRoot = Join-Path $researchRoot "8月17日交流之后的任务\algorithm_revised"
$baselineRoot = Join-Path $researchRoot "8月17日交流之后的任务\原始算法的代码实现"

foreach ($required in @($PythonExe, $DataRoot, $datasetManifest, $acSource, $algorithmRoot, $baselineRoot)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Required path does not exist: $required"
    }
}

New-Item -ItemType Directory -Force -Path $ResultDir | Out-Null
$env:PYTHONPATH = "$(Join-Path $algorithmRoot 'src');$(Join-Path $baselineRoot 'src')"

Write-Host "STEP 1/3  AC project dispatch and checkpoint validation"
& $PythonExe (Join-Path $codeRoot "run_ac_all.py") `
    --dataset-manifest $datasetManifest `
    --data-root $DataRoot `
    --output-dir $ResultDir `
    --ac-source $acSource `
    --ac-python $PythonExe `
    --jobs 8 `
    --no-bfbarchitect
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

Write-Host "STEP 3/3  Independent acceptance"
& $PythonExe (Join-Path $codeRoot "verify_all_results.py") `
    --data-root $DataRoot `
    --output-dir $ResultDir
if ($LASTEXITCODE -ne 0) { throw "Verification failed with exit code $LASTEXITCODE" }

Write-Host "TASK2_FULL_PIPELINE_PASSED"
