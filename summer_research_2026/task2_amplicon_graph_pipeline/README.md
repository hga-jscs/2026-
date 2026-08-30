# Amplicon graph maximum cyclic ratio pipeline

This directory contains the reproducible code and compact results for task2. It pairs AmpliconRepository graph/cycles files, runs AmpliconClassifier, solves the original-graph linear program from task1 task 3, writes the required four-column table, and performs automated provenance and numerical consistency checks.

## Output semantics

- `sample` is `project_id::original_sample`, so `(sample, amplicon)` is unique across projects.
- `lwcn` is the maximum cyclic length-weighted copy-number **ratio** in `[0,1]`, as required by the task; it is not an edge copy number.
- `classification` is AmpliconClassifier's `amplicon_decomposition_class`, not its `ecDNA+` field.
- The extended provenance table separately retains the absolute maximum cyclic LWCN, ratio, `ecDNA+`, `BFB+`, `FAN+`, graph/cycles hashes, parser messages, solver status and residuals.

## Contents

| Path | Purpose |
|---|---|
| `pipeline/` | Download, extraction, AC, LP, merge and automated-check scripts |
| `source/` | `original_graph_lwcn` and required `cyclic_lwcn` source |
| `tests/` | Six tests and the 112-pair CoRAL regression data |
| `results/` | Four-column results, provenance, manifests and automated check |
| `server/` | Limited Linux deployment, run and verification scripts |

## Verified snapshot

- 32 public projects and 28,142 project-scoped amplicon records.
- AC: 28,142 classified records, generated with `--no_bfbarchitect`.
- LP: 28,131 `OPTIMAL` and 11 `TRIVIAL_OPTIMAL_ZERO` records.
- Required four-column key duplicates: 0.
- Unsupported graph formats: 0. Legal AA source records, length conventions and header variants are INFO; only 2 records retain a structured WARNING.
- The automated check re-hashes all 32 archives and all 56,284 graph/cycles source files, checks current code/config/input fingerprints, recomputes ratios, and checks LP feasibility and duality-gap scale.
- The Linux evidence is a current-code compatibility smoke test on 112 CoRAL amplicons, not a full 28,142-record server run.

These checks establish file provenance, identity consistency and numerical feasibility. They are not an independent biological validation.

## Run

Python 3.10 or newer is required. Install `requirements.txt`, install AmpliconClassifier separately, prepare its GRCh38 reference data, and then run:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
$env:PYTHONPATH = (Resolve-Path .\source).Path
.\.venv\Scripts\python.exe -m pytest .\tests -q

.\pipeline\run_all.ps1 `
  -PythonExe .\.venv\Scripts\python.exe `
  -AcRoot D:\path\to\AmpliconClassifier `
  -DataRoot D:\path\to\data `
  -ResultDir D:\path\to\results
```

The input archive and reference data are intentionally not committed. See `pipeline/逐行注释索引.md` and the local server guide for code locations and deployment details.
