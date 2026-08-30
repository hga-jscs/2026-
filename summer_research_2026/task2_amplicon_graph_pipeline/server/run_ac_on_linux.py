#!/usr/bin/env python3
"""在 Linux 上运行单个项目的 AmpliconClassifier，并写入可核对元数据。"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from pipeline_provenance import (
    ac_runtime_tree_sha256,
    canonical_sha256,
    paired_input_snapshot,
    sha256,
    tree_sha256,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def row_count(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def reference_sha256(reference_dir: Path) -> str:
    file_list = reference_dir / "file_list.txt"
    paths = [file_list]
    for line in file_list.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 2 and not fields[1].startswith("#"):
            candidate = reference_dir / fields[1]
            if candidate.is_file():
                paths.append(candidate)
    return canonical_sha256(
        [
            {
                "path": path.relative_to(reference_dir).as_posix(),
                "sha256": sha256(path),
            }
            for path in sorted(set(paths))
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--deploy-root", required=True, type=Path)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--reference", default="GRCh38")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()

    root = args.deploy_root.resolve()
    project_id = args.project_id
    ac_source = root / "AmpliconClassifier"
    input_dir = root / "data_root" / "extracted" / project_id
    reference_dir = root / "reference_data" / args.reference
    output_dir = root / "data_root" / "ac_runs" / project_id
    scratch = root / "scratch" / project_id
    profile = output_dir / f"{project_id}_amplicon_classification_profiles.tsv"
    log_path = output_dir / "batch_console.log"
    metadata_path = output_dir / "run_metadata.json"
    manifest = json.loads((root / "runner" / "coral_manifest.json").read_text(encoding="utf-8"))
    project = next(row for row in manifest if row["project_id"] == project_id)
    expected = int(project["paired_amplicon_count"])

    output_dir.mkdir(parents=True, exist_ok=True)
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    prefix = scratch / project_id

    version = json.loads((ac_source / "版本溯源.json").read_text(encoding="utf-8"))
    input_snapshot = paired_input_snapshot(input_dir)
    source_hash = ac_runtime_tree_sha256(ac_source)
    config_hash = sha256(ac_source / "ampclasslib" / "default_config.json")
    payload = {
        "project_id": project_id,
        "reference_genome": args.reference,
        "bfbarchitect_mode": "disabled",
        "bfb_solver": "auto",
        "jobs": min(args.jobs, expected),
        "input_snapshot_sha256": input_snapshot["snapshot_sha256"],
        "reference_snapshot_sha256": reference_sha256(reference_dir),
        "ac_source_sha256": source_hash,
        "default_config_sha256": config_hash,
        "upstream_commit": version["upstream_commit"],
        "pipeline_sha256": tree_sha256(root / "runner", ["run_ac_on_linux.py", "pipeline_provenance.py"]),
        "python_executable": sys.executable,
    }
    fingerprint = canonical_sha256(payload)
    command = [
        sys.executable,
        str(ac_source / "amplicon_classifier.py"),
        "--ref", args.reference,
        "--AA_results", str(input_dir),
        "-o", str(prefix),
        "--jobs", str(min(args.jobs, expected)),
        "--bfb_solver", "auto",
        "--bfb_threads", "1",
        "--no_results_table",
        "--no_bfbarchitect",
    ]
    environment = os.environ.copy()
    environment.update(
        {
            "AA_DATA_REPO": str(root / "reference_data"),
            "PYTHONPATH": str(ac_source),
            "MPLBACKEND": "Agg",
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
    )
    started = utc_now()
    start = time.perf_counter()
    with log_path.open("wb") as log_handle:
        exit_code = subprocess.call(command, env=environment, stdout=log_handle, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - start
    source_profile = Path(f"{prefix}_amplicon_classification_profiles.tsv")
    source_log = Path(f"{prefix}.log")
    if exit_code == 0 and source_profile.is_file():
        shutil.copy2(source_profile, profile)
        if source_log.is_file():
            shutil.copy2(source_log, output_dir / f"{project_id}.log")
    classified = row_count(profile) if profile.is_file() else -1
    status = "COMPLETE" if exit_code == 0 and classified == expected else "ERROR"
    metadata = {
        "status": status,
        "project_id": project_id,
        "hostname": socket.gethostname(),
        "started_at_utc": started,
        "finished_at_utc": utc_now(),
        "elapsed_seconds": round(elapsed, 3),
        "exit_code": exit_code,
        "execution_status": "RECOMPUTED_ON_LINUX_SERVER",
        "run_fingerprint": fingerprint,
        **payload,
        "upstream_repository": version["upstream_repository"],
        "profile_sha256": sha256(profile) if profile.is_file() else "",
        "console_log_sha256": sha256(log_path),
        "classified_amplicons": max(classified, 0),
    }
    metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"SERVER_AC_{status} rows={classified} exit_code={exit_code}")
    return 0 if status == "COMPLETE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
