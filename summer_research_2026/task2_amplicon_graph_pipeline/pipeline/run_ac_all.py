"""Run AmpliconClassifier per public project with fingerprint-checked reuse."""

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

from pipeline_provenance import (
    ac_runtime_tree_sha256,
    canonical_sha256,
    paired_input_snapshot,
    sha256,
    tree_sha256,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def to_wsl(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"expected a Windows drive path: {resolved}")
    tail = resolved.as_posix().split(":", 1)[1].lstrip("/")
    return f"/mnt/{drive}/{tail}"


def wsl_python(value: str) -> str:
    """Accept either a native Linux interpreter or a Windows interpreter path."""

    if value.startswith("/"):
        return value
    return to_wsl(Path(value))


def row_count(tsv_path: Path) -> int:
    if not tsv_path.exists():
        return -1
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def read_version_record(ac_source: Path) -> dict[str, str]:
    path = ac_source / "版本溯源.json"
    if not path.is_file():
        raise FileNotFoundError(f"missing AC version record: {path}")
    record = json.loads(path.read_text(encoding="utf-8"))
    required = {"upstream_repository", "upstream_commit"}
    if not required.issubset(record):
        raise ValueError(f"incomplete AC version record: {path}")
    return {key: str(value) for key, value in record.items()}


def reference_snapshot_sha256(reference_dir: Path) -> str:
    """Hash file_list and every existing file it names."""

    file_list = reference_dir / "file_list.txt"
    if not file_list.is_file():
        raise FileNotFoundError(f"missing reference file list: {file_list}")
    paths = [file_list]
    for raw_line in file_list.read_text(encoding="utf-8").splitlines():
        fields = raw_line.split()
        if len(fields) < 2 or fields[1].startswith("#"):
            continue
        candidate = reference_dir / fields[1]
        if candidate.is_file():
            paths.append(candidate)
    records = [
        {"path": path.relative_to(reference_dir).as_posix(), "sha256": sha256(path)}
        for path in sorted(set(paths), key=lambda item: item.relative_to(reference_dir).as_posix())
    ]
    return canonical_sha256(records)


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
        "execution_status",
        "started_at_utc",
        "finished_at_utc",
        "validated_at_utc",
        "elapsed_seconds",
        "exit_code",
        "input_snapshot_sha256",
        "reference_snapshot_sha256",
        "ac_source_sha256",
        "default_config_sha256",
        "upstream_commit",
        "run_fingerprint",
        "profile_sha256",
        "error",
    ]
    with (output_dir / "ac_run_manifest.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the frozen AmpliconClassifier once per public project."
    )
    parser.add_argument("--dataset-manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ac-source", required=True, type=Path)
    parser.add_argument(
        "--ac-python",
        required=True,
        help="Python interpreter used inside WSL; accepts /usr/bin/python3 or a Windows path.",
    )
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--bfb-solver", choices=("auto", "gurobi", "mosek", "cbc"), default="auto")
    parser.add_argument("--no-bfbarchitect", action="store_true")
    parser.add_argument("--force", action="store_true", help="Run AC even when fingerprints match.")
    parser.add_argument("--wsl-scratch-root", default="/tmp/task2_ac")
    args = parser.parse_args()

    projects = json.loads(args.dataset_manifest.read_text(encoding="utf-8"))
    projects = [row for row in projects if row.get("status") == "READY"]
    ac_script = args.ac_source / "amplicon_classifier.py"
    default_config = args.ac_source / "ampclasslib" / "default_config.json"
    reference_root = args.data_root / "reference_data"
    ac_runs = args.data_root / "ac_runs"
    version_record = read_version_record(args.ac_source)
    ac_python = wsl_python(args.ac_python)
    ac_source_sha256 = ac_runtime_tree_sha256(args.ac_source)
    default_config_sha256 = sha256(default_config)
    pipeline_sha256 = tree_sha256(
        Path(__file__).resolve().parent,
        ["run_ac_all.py", "pipeline_provenance.py"],
    )
    reference_hashes: dict[str, str] = {}
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
        metadata_path = project_output / "run_metadata.json"
        console_log = project_output / "batch_console.log"
        bfbarchitect_mode = "disabled" if args.no_bfbarchitect else f"enabled:{args.bfb_solver}"
        validated_at = utc_now()
        record: dict[str, object] = {
            "project_id": project_id,
            "project_name": project_name,
            "reference_genome": reference,
            "bfbarchitect_mode": bfbarchitect_mode,
            "expected_amplicons": expected,
            "classified_amplicons": max(row_count(profile), 0),
            "status": "PENDING",
            "execution_status": "",
            "started_at_utc": "",
            "finished_at_utc": "",
            "validated_at_utc": validated_at,
            "elapsed_seconds": "",
            "exit_code": "",
            "input_snapshot_sha256": "",
            "reference_snapshot_sha256": "",
            "ac_source_sha256": ac_source_sha256,
            "default_config_sha256": default_config_sha256,
            "upstream_commit": version_record["upstream_commit"],
            "run_fingerprint": "",
            "profile_sha256": sha256(profile) if profile.is_file() else "",
            "error": "",
        }
        print(
            f"AC_PROJECT {index}/{len(projects)} id={project_id} ref={reference} expected={expected}",
            flush=True,
        )
        try:
            reference_file_list = reference_root / reference / "file_list.txt"
            if not reference_file_list.is_file():
                raise FileNotFoundError(f"missing reference data: {reference}")
            if reference not in reference_hashes:
                reference_hashes[reference] = reference_snapshot_sha256(reference_root / reference)
            input_snapshot = paired_input_snapshot(project_input)
            if int(input_snapshot["pair_count"]) != expected:
                raise ValueError(
                    f"input snapshot pair count {input_snapshot['pair_count']} != expected {expected}"
                )
            fingerprint_payload = {
                "project_id": project_id,
                "reference_genome": reference,
                "bfbarchitect_mode": bfbarchitect_mode,
                "bfb_solver": args.bfb_solver,
                "jobs": max(1, min(args.jobs, expected)),
                "input_snapshot_sha256": input_snapshot["snapshot_sha256"],
                "reference_snapshot_sha256": reference_hashes[reference],
                "ac_source_sha256": ac_source_sha256,
                "default_config_sha256": default_config_sha256,
                "upstream_commit": version_record["upstream_commit"],
                "pipeline_sha256": pipeline_sha256,
                "python_executable": ac_python,
            }
            run_fingerprint = canonical_sha256(fingerprint_payload)
            record.update(
                {
                    "input_snapshot_sha256": input_snapshot["snapshot_sha256"],
                    "reference_snapshot_sha256": reference_hashes[reference],
                    "run_fingerprint": run_fingerprint,
                }
            )

            existing = row_count(profile)
            reuse = False
            if metadata_path.is_file() and not args.force:
                previous = json.loads(metadata_path.read_text(encoding="utf-8"))
                reuse = (
                    previous.get("status") == "COMPLETE"
                    and previous.get("run_fingerprint") == run_fingerprint
                    and existing == expected
                    and previous.get("profile_sha256") == sha256(profile)
                )
                if reuse:
                    record.update(
                        {
                            "classified_amplicons": existing,
                            "status": "COMPLETE",
                            "execution_status": "REUSED_FINGERPRINT_MATCHED",
                            "started_at_utc": previous.get("started_at_utc", ""),
                            "finished_at_utc": previous.get("finished_at_utc", ""),
                            "elapsed_seconds": previous.get("elapsed_seconds", 0),
                            "exit_code": 0,
                            "profile_sha256": previous["profile_sha256"],
                        }
                    )
                    rows.append(record)
                    write_run_manifest(rows, args.output_dir)
                    print(f"AC_REUSED id={project_id} rows={existing}", flush=True)
                    continue

            if not project_id.isalnum():
                raise ValueError(f"unsafe project id for scratch path: {project_id}")
            if not args.wsl_scratch_root.startswith("/tmp/"):
                raise ValueError("WSL scratch root must stay below /tmp")
            scratch_dir = f"{args.wsl_scratch_root.rstrip('/')}/{project_id}"
            scratch_prefix = f"{scratch_dir}/{project_id}"
            subprocess.run(["wsl.exe", "--", "rm", "-rf", "--", scratch_dir], check=True)
            subprocess.run(["wsl.exe", "--", "mkdir", "-p", scratch_dir], check=True)
            command = [
                "wsl.exe",
                "--",
                "/usr/bin/env",
                f"AA_DATA_REPO={to_wsl(reference_root)}",
                f"PYTHONPATH={to_wsl(args.ac_source)}",
                ac_python,
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
            started_at = utc_now()
            start = time.perf_counter()
            with console_log.open("wb") as log_handle:
                exit_code = subprocess.call(command, stdout=log_handle, stderr=subprocess.STDOUT)
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
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    if subprocess.call(["wsl.exe", "--", "cp", source, to_wsl(destination)]) != 0:
                        exit_code = 3
            classified = row_count(profile)
            finished_at = utc_now()
            status = "COMPLETE" if exit_code == 0 and classified == expected else "ERROR"
            profile_hash = sha256(profile) if profile.is_file() else ""
            metadata = {
                "status": status,
                "project_id": project_id,
                "started_at_utc": started_at,
                "finished_at_utc": finished_at,
                "elapsed_seconds": round(elapsed, 3),
                "exit_code": exit_code,
                "execution_status": "RECOMPUTED",
                "run_fingerprint": run_fingerprint,
                **fingerprint_payload,
                "upstream_repository": version_record["upstream_repository"],
                "profile_sha256": profile_hash,
                "console_log_sha256": sha256(console_log),
                "command": [
                    "amplicon_classifier.py",
                    "--ref",
                    reference,
                    "--AA_results",
                    f"<data_root>/extracted/{project_id}",
                    "--jobs",
                    str(max(1, min(args.jobs, expected))),
                    "--no_results_table",
                    *( ["--no_bfbarchitect"] if args.no_bfbarchitect else [] ),
                ],
            }
            metadata_path.write_text(
                json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
            )
            record.update(
                {
                    "classified_amplicons": max(classified, 0),
                    "status": status,
                    "execution_status": "RECOMPUTED",
                    "started_at_utc": started_at,
                    "finished_at_utc": finished_at,
                    "elapsed_seconds": round(elapsed, 3),
                    "exit_code": exit_code,
                    "profile_sha256": profile_hash,
                }
            )
            if status == "COMPLETE":
                subprocess.run(["wsl.exe", "--", "rm", "-rf", "--", scratch_dir], check=False)
                print(f"AC_DONE id={project_id} rows={classified} seconds={elapsed:.3f}", flush=True)
            else:
                record["error"] = f"exit={exit_code}; expected={expected}; classified={classified}"
                print(f"AC_ERROR id={project_id} error={record['error']}", file=sys.stderr, flush=True)
        except Exception as exc:
            record["status"] = "ERROR"
            record["execution_status"] = "FAILED_BEFORE_OR_DURING_RUN"
            record["exit_code"] = 2
            record["error"] = f"{type(exc).__name__}: {exc}"
            print(f"AC_ERROR id={project_id} error={record['error']}", file=sys.stderr, flush=True)
        rows.append(record)
        write_run_manifest(rows, args.output_dir)

    errors = sum(row["status"] == "ERROR" for row in rows)
    classified_total = sum(int(row["classified_amplicons"]) for row in rows)
    print(f"AC_COMPLETE projects={len(rows)} errors={errors} rows={classified_total}", flush=True)
    return 0 if errors == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
