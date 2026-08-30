"""按公开项目运行当前 AmpliconClassifier，并持续写入可恢复运行清单。

每个项目独立执行；若 profile 行数已与数据清单一致，则安全复用该完整检查点。
不完整项目在 WSL 中使用隔离 scratch 运行，完成后核对行数再标记 COMPLETE。
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def to_wsl(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"expected a Windows drive path: {resolved}")
    tail = resolved.as_posix().split(":", 1)[1].lstrip("/")
    return f"/mnt/{drive}/{tail}"


def row_count(tsv_path: Path) -> int:
    if not tsv_path.exists():
        return -1
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def write_run_manifest(rows: list[dict[str, object]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "ac_run_manifest.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    columns = [
        "project_id",
        "project_name",
        "reference_genome",
        "bfbarchitect_mode",
        "expected_amplicons",
        "classified_amplicons",
        "status",
        "started_at_utc",
        "finished_at_utc",
        "elapsed_seconds",
        "exit_code",
        "error",
    ]
    with (output_dir / "ac_run_manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the current AmpliconClassifier once per public project.")
    parser.add_argument("--dataset-manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ac-source", required=True, type=Path)
    parser.add_argument("--ac-python", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--bfb-solver", choices=("auto", "gurobi", "mosek", "cbc"), default="auto")
    parser.add_argument("--no-bfbarchitect", action="store_true")
    parser.add_argument("--wsl-scratch-root", default="/tmp/task2_ac")
    args = parser.parse_args()

    projects = json.loads(args.dataset_manifest.read_text(encoding="utf-8"))
    projects = [row for row in projects if row.get("status") == "READY"]
    ac_script = args.ac_source / "amplicon_classifier.py"
    reference_root = args.data_root / "reference_data"
    ac_runs = args.data_root / "ac_runs"
    rows: list[dict[str, object]] = []
    print(f"AC_PROJECTS count={len(projects)}", flush=True)

    for index, project in enumerate(projects, start=1):
        project_id = str(project["project_id"])
        project_name = str(project["project_name"])
        reference = str(project["reference_genome"])
        if reference == "hg38":
            reference = "GRCh38"
        elif reference == "GRCm38":
            reference = "mm10"
        expected = int(project["paired_amplicon_count"])
        project_input = args.data_root / "extracted" / project_id
        project_output = ac_runs / project_id
        project_output.mkdir(parents=True, exist_ok=True)
        profile = project_output / f"{project_id}_amplicon_classification_profiles.tsv"
        existing = row_count(profile)
        record: dict[str, object] = {
            "project_id": project_id,
            "project_name": project_name,
            "reference_genome": reference,
            "bfbarchitect_mode": "disabled" if args.no_bfbarchitect else f"enabled:{args.bfb_solver}",
            "expected_amplicons": expected,
            "classified_amplicons": max(existing, 0),
            "status": "PENDING",
            "started_at_utc": "",
            "finished_at_utc": "",
            "elapsed_seconds": "",
            "exit_code": "",
            "error": "",
        }
        print(
            f"AC_PROJECT {index}/{len(projects)} id={project_id} ref={reference} expected={expected}",
            flush=True,
        )
        if existing == expected:
            record["status"] = "REUSED_COMPLETE"
            rows.append(record)
            write_run_manifest(rows, args.output_dir)
            print(f"AC_REUSED id={project_id} rows={existing}", flush=True)
            continue

        if not (reference_root / reference / "file_list.txt").exists():
            record["status"] = "ERROR"
            record["error"] = f"missing reference data: {reference}"
            rows.append(record)
            write_run_manifest(rows, args.output_dir)
            print(f"AC_ERROR id={project_id} error={record['error']}", file=sys.stderr, flush=True)
            continue

        if not project_id.isalnum():
            raise ValueError(f"unsafe project id for scratch path: {project_id}")
        scratch_dir = f"{args.wsl_scratch_root.rstrip('/')}/{project_id}"
        scratch_prefix = f"{scratch_dir}/{project_id}"
        if subprocess.call(["wsl.exe", "--", "mkdir", "-p", scratch_dir]) != 0:
            record["status"] = "ERROR"
            record["error"] = f"could not create WSL scratch directory: {scratch_dir}"
            rows.append(record)
            write_run_manifest(rows, args.output_dir)
            print(f"AC_ERROR id={project_id} error={record['error']}", file=sys.stderr, flush=True)
            continue

        command = [
            "wsl.exe",
            "--",
            "/usr/bin/env",
            f"AA_DATA_REPO={to_wsl(reference_root)}",
            f"PYTHONPATH={to_wsl(args.ac_source)}",
            to_wsl(args.ac_python),
            to_wsl(ac_script),
            "--ref",
            reference,
            "--AA_results",
            to_wsl(project_input),
            "-o",
            scratch_prefix,
            "--jobs",
            str(max(1, min(args.jobs, expected))),
            "--bfb_solver",
            args.bfb_solver,
            "--bfb_threads",
            "1",
            "--no_results_table",
        ]
        if args.no_bfbarchitect:
            command.append("--no_bfbarchitect")
        record["started_at_utc"] = utc_now()
        start = time.perf_counter()
        # AmpliconClassifier writes its detailed log to the output prefix.  Keep
        # the batch console compact so a full public-repository run does not
        # accumulate millions of progress lines in the orchestration session.
        exit_code = subprocess.call(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        elapsed = time.perf_counter() - start
        if exit_code == 0:
            copy_specs = [
                (f"{scratch_prefix}_amplicon_classification_profiles.tsv", profile, True),
                (f"{scratch_prefix}.log", project_output / f"{project_id}.log", True),
                (
                    f"{scratch_prefix}_foldback_qc_filtered.txt",
                    project_output / f"{project_id}_foldback_qc_filtered.txt",
                    False,
                ),
            ]
            for source, destination, required in copy_specs:
                exists = subprocess.call(["wsl.exe", "--", "test", "-f", source]) == 0
                if not exists:
                    if required:
                        exit_code = 3
                    continue
                if subprocess.call(["wsl.exe", "--", "cp", source, to_wsl(destination)]) != 0:
                    exit_code = 3
            if exit_code == 0:
                subprocess.call(["wsl.exe", "--", "rm", "-rf", "--", scratch_dir])
        classified = row_count(profile)
        record.update(
            {
                "classified_amplicons": max(classified, 0),
                "finished_at_utc": utc_now(),
                "elapsed_seconds": round(elapsed, 3),
                "exit_code": exit_code,
            }
        )
        if exit_code == 0 and classified == expected:
            record["status"] = "COMPLETE"
            print(f"AC_DONE id={project_id} rows={classified} seconds={elapsed:.3f}", flush=True)
        else:
            record["status"] = "ERROR"
            record["error"] = f"exit={exit_code}; expected={expected}; classified={classified}"
            print(f"AC_ERROR id={project_id} error={record['error']}", file=sys.stderr, flush=True)
        rows.append(record)
        write_run_manifest(rows, args.output_dir)

    errors = sum(row["status"] == "ERROR" for row in rows)
    classified_total = sum(int(row["classified_amplicons"]) for row in rows)
    print(f"AC_COMPLETE projects={len(rows)} errors={errors} rows={classified_total}", flush=True)
    return 0 if errors == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
