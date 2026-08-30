"""Automated consistency checks for the bounded CoRAL Linux server rerun."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ID = "6a581fe4b219fa804843a951"
EXPECTED_ROWS = 112
TOLERANCE = 1e-7
ACCEPTED_STATUSES = {"OPTIMAL", "TRIVIAL_OPTIMAL_ZERO"}


def read_rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--deploy-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.deploy_root.resolve()
    profile_path = (
        root
        / "data_root"
        / "ac_runs"
        / PROJECT_ID
        / f"{PROJECT_ID}_amplicon_classification_profiles.tsv"
    )
    checkpoint_path = root / "data_root" / "lp_runs" / f"{PROJECT_ID}.csv"
    final_path = root / "results" / "server_coral_112_ac_lwcn_ratio.csv"
    ac_metadata_path = root / "data_root" / "ac_runs" / PROJECT_ID / "run_metadata.json"
    lp_metadata_path = root / "data_root" / "lp_runs" / f"{PROJECT_ID}.metadata.json"
    errors: list[str] = []

    for path in (profile_path, checkpoint_path, final_path, ac_metadata_path, lp_metadata_path):
        if not path.is_file():
            errors.append(f"missing file: {path}")
    if errors:
        print("TASK2 SERVER VERIFICATION FAILED")
        return 2

    profile = read_rows(profile_path, delimiter="\t")
    checkpoint = read_rows(checkpoint_path)
    final = read_rows(final_path)
    if len(profile) != EXPECTED_ROWS:
        errors.append(f"AC rows={len(profile)}")
    if len(checkpoint) != EXPECTED_ROWS:
        errors.append(f"LP rows={len(checkpoint)}")
    if len(final) != EXPECTED_ROWS:
        errors.append(f"final rows={len(final)}")
    if final and list(final[0]) != ["sample", "amplicon", "lwcn", "classification"]:
        errors.append("unexpected final CSV header")

    profile_map = {
        (row["sample_name"], row["amplicon_number"]): row["amplicon_decomposition_class"]
        for row in profile
    }
    checkpoint_map = {(row["sample"], row["amplicon"]): row for row in checkpoint}
    final_map = {
        (row["sample"].split("::", 1)[-1], row["amplicon"]): row for row in final
    }
    if set(profile_map) != set(checkpoint_map) or set(profile_map) != set(final_map):
        errors.append("sample/amplicon keys differ across AC, LP and final CSV")

    maxima = {
        "max_balance_residual": 0.0,
        "max_lower_bound_violation": 0.0,
        "max_upper_bound_violation": 0.0,
    }
    for key, row in checkpoint_map.items():
        if row["lp_status"] not in ACCEPTED_STATUSES:
            errors.append(f"unaccepted LP status for {key}: {row['lp_status']}")
        for field in maxima:
            value = float(row[field])
            if not math.isfinite(value):
                errors.append(f"non-finite {field} for {key}")
            maxima[field] = max(maxima[field], value)
            if value > TOLERANCE:
                errors.append(f"{field} exceeds tolerance for {key}: {value}")
        if key in final_map:
            if final_map[key]["classification"] != profile_map[key]:
                errors.append(f"classification mismatch for {key}")
            lwcn = float(final_map[key]["lwcn"])
            if not math.isfinite(lwcn) or lwcn < -TOLERANCE or lwcn > 1.0 + TOLERANCE:
                errors.append(f"invalid maximum cyclic LWCN ratio for {key}: {lwcn}")

    ac_metadata = json.loads(ac_metadata_path.read_text(encoding="utf-8"))
    lp_metadata = json.loads(lp_metadata_path.read_text(encoding="utf-8"))
    if ac_metadata.get("status") != "COMPLETE" or int(ac_metadata.get("exit_code", -1)) != 0:
        errors.append("AC metadata is not a completed zero-exit run")
    if lp_metadata.get("status") != "COMPLETE" or int(lp_metadata.get("exit_code", -1)) != 0:
        errors.append("LP metadata is not a completed zero-exit run")
    if ac_metadata.get("profile_sha256") != sha256(profile_path):
        errors.append("AC profile fingerprint mismatch")
    if lp_metadata.get("checkpoint_sha256") != sha256(checkpoint_path):
        errors.append("LP checkpoint fingerprint mismatch")

    report = {
        "check_type": "automated_consistency_check",
        "status": "PASSED" if not errors else "FAILED",
        "verified_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "project_id": PROJECT_ID,
        "hostname": ac_metadata.get("hostname", ""),
        "ac_started_at_utc": ac_metadata.get("started_at_utc", ""),
        "ac_finished_at_utc": ac_metadata.get("finished_at_utc", ""),
        "ac_run_fingerprint": ac_metadata.get("run_fingerprint", ""),
        "lp_run_fingerprint": lp_metadata.get("lp_run_fingerprint", ""),
        "ac_rows": len(profile),
        "lp_rows": len(checkpoint),
        "final_rows": len(final),
        "tolerance": TOLERANCE,
        **maxima,
        "ac_source_sha256": sha256(root / "AmpliconClassifier" / "amplicon_classifier.py"),
        "lp_source_sha256": sha256(
            root
            / "runner"
            / "算法源码"
            / "original_graph_lwcn"
            / "original_graph_linear_program.py"
        ),
        "final_csv_sha256": sha256(final_path),
        "errors": errors,
    }
    result_dir = root / "results"
    (result_dir / "server_automated_check.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    summary = (
        f"TASK2 SERVER AUTOMATED CONSISTENCY CHECK {report['status']}\n"
        f"AC={len(profile)}; LP={len(checkpoint)}; CSV={len(final)}\n"
        f"max_balance_residual={maxima['max_balance_residual']:.6g}; "
        f"max_lower_bound_violation={maxima['max_lower_bound_violation']:.6g}; "
        f"max_upper_bound_violation={maxima['max_upper_bound_violation']:.6g}\n"
    )
    (result_dir / "server_automated_check.txt").write_text(summary, encoding="utf-8")
    print(summary, end="")
    if errors:
        for error in errors[:30]:
            print(f"- {error}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
