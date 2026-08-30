# AmpliconClassifier 与原图最大环状 LWCN 全数据流程

本目录是 task2 的可复核代码交付。它下载并整理 AmpliconRepository 公开 AA 结果，批量调用 AmpliconClassifier（AC），在原始 breakpoint graph 上求最大可行环状长度加权拷贝数（LWCN），合并四列结果，并执行独立验收。

## 目录

| 路径 | 内容 |
|---|---|
| `pipeline/` | 数据下载、解包、AC 批处理、LP 批处理、合并与验收入口 |
| `src/original_graph_lwcn/` | 当前原图 LP 实现 |
| `src/cyclic_lwcn/` | CoRAL graph/cycles 解析与对照实现 |
| `tests/` | 原图 LP 单元测试 |
| `results/` | 2026-08-30 全量运行清单、最终 CSV 与验收 JSON |
| `server/` | 实验室服务器的隔离部署、受限运行与验收脚本 |

`pipeline/逐行注释索引.md` 逐段对应所有可执行代码；源文件本身保留了模块、函数、关键分支与安全边界注释。

## 已复核结果

- 公开项目：32。
- 配对 graph/cycles：28,142。
- AC 完成：28,142；LP 完成：28,142。
- LP 状态：28,131 个 `OPTIMAL`，11 个 `TRIVIAL_OPTIMAL_ZERO`。
- 最大流平衡残差：`9.594493466380527e-08`；最大下界/上界违反分别为 `3.597120556975142e-08`、`9.977164339147748e-08`，均不超过 `1e-7` 验收容差。
- 最终文件：`results/all_ac_lwcn_results.csv`；机器验收：`results/verification.json`。

28,030 行带解析警告，主要来自上游 AA 文本中的额外行、长度约定差异或无法匹配的端点；这些行仍通过 LP 状态、有限数值和 `1e-7` 可行性验收。警告不是静默丢弃，完整内容保留在每个项目的 LP checkpoint 中（原始数据与 checkpoint 未提交到 GitHub）。

## 环境

1. Python 3.10 或更新版本。
2. 安装 `requirements.txt`。
3. 单独安装 AmpliconClassifier，并准备其 `GRCh38` 参考数据；本仓库没有复制上游 AC 仓库或基因组参考文件。
4. 将仓库 `src` 加入 `PYTHONPATH`。

PowerShell 示例：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
$env:PYTHONPATH = (Resolve-Path .\src).Path
.\.venv\Scripts\python.exe -m pytest .\tests -q
```

上述命令执行 4 个无需外部数据的单元测试，并在未配置数据时跳过 112 对 CoRAL 集成测试。若已有 CoRAL 输入，先设置 `$env:CORAL_DATA_ROOT` 为包含 112 对 graph/cycles 的目录，再运行 pytest。

## 全流程入口

```powershell
.\pipeline\run_all.ps1 `
  -PythonExe .\.venv\Scripts\python.exe `
  -AcRoot D:\path\to\AmpliconClassifier `
  -DataRoot D:\path\to\data_root `
  -OutputDir D:\path\to\task2_results
```

脚本按“公开项目清单 → 下载/安全解包 → AC → 原图 LP → 四列 CSV → 独立验收”执行，并复用已完成的合法 checkpoint。下载数据可能较大，先检查磁盘空间。服务器脚本固定到 task2 私有目录，运行前应按本机路径修改或复核，不要在未知目录直接执行。

## 解释边界

`lwcn` 是给定 breakpoint graph 容量和流平衡约束下的最大环状 LWCN；`classification` 是 AC 对 graph 与 cycles 等证据的分类。仅凭 graph 得到正的 `lwcn` 不能替代 AC 的 `Cyclic/ecDNA` 判定，也不能恢复原始 cycles 分解。
