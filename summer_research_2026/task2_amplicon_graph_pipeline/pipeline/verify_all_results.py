"""Automated consistency checks for the complete AC + cyclic-ratio snapshot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from pipeline_provenance import (
    ac_runtime_tree_sha256,
    canonical_sha256,
    paired_input_snapshot,
    sha256,
    tree_sha256,
)


FEASIBILITY_TOLERANCE = 1e-7
DUALITY_GAP_RATIO_TOLERANCE = 1e-12
ACCEPTED_LP_STATUSES = {"OPTIMAL", "TRIVIAL_OPTIMAL_ZERO"}
ALLOWED_CLASSIFICATIONS = {
    "Cyclic",
    "Complex-non-cyclic",
    "Linear",
    "No-FSCNA",
    "Invalid",
    "Virus",
}


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def fail_if(condition: bool, message: str, failures: list[str]) -> None:
    if condition:
        failures.append(message)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Automated provenance, identity, hash, range and residual consistency checks."
    )
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ac-source", type=Path)
    args = parser.parse_args()

    code_root = Path(__file__).resolve().parent
    repository_root = code_root.parent
    ac_source = args.ac_source or repository_root / "AmpliconClassifier"
    failures: list[str] = []
    dataset_rows = read_csv(args.output_dir / "数据集清单.csv")
    ac_rows = read_csv(args.output_dir / "ac_run_manifest.csv")
    lp_rows = read_csv(args.output_dir / "lp_run_manifest.csv")
    final_path = args.output_dir / "全量AC与环状LWCN结果.csv"
    extended_path = args.output_dir / "全量扩展溯源结果.csv"
    with final_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        final_header = reader.fieldnames
        final_rows = list(reader)
    with extended_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        extended_header = reader.fieldnames or []
        extended_rows = list(reader)

    fail_if(final_header != ["sample", "amplicon", "lwcn", "classification"],
            f"unexpected four-column header: {final_header}", failures)
    fail_if(not dataset_rows, "dataset manifest is empty", failures)
    fail_if(any(row["status"] != "READY" for row in dataset_rows),
            "at least one dataset is not READY", failures)
    expected_total = sum(int(row["paired_amplicon_count"]) for row in dataset_rows)
    expected_projects = len(dataset_rows)

    archive_total = 0
    archive_hashes_verified = 0
    snapshot_by_project: dict[str, dict[str, object]] = {}
    for row in dataset_rows:
        project_id = row["project_id"]
        archive = args.data_root / "archives" / f"{project_id}.tar.gz"
        if not archive.is_file():
            failures.append(f"archive missing: {project_id}")
            continue
        actual_size = archive.stat().st_size
        archive_total += actual_size
        fail_if(actual_size != int(row["archive_bytes"]),
                f"archive size mismatch: {project_id}", failures)
        actual_hash = sha256(archive)
        fail_if(actual_hash != row["sha256"], f"archive SHA-256 mismatch: {project_id}", failures)
        if actual_hash == row["sha256"]:
            archive_hashes_verified += 1
        snapshot = paired_input_snapshot(args.data_root / "extracted" / project_id)
        snapshot_by_project[project_id] = snapshot
        fail_if(int(snapshot["pair_count"]) != int(row["paired_amplicon_count"]),
                f"input pair count mismatch: {project_id}", failures)

    fail_if(len(ac_rows) != expected_projects,
            f"AC project count {len(ac_rows)} != {expected_projects}", failures)
    fail_if(any(row["status"] != "COMPLETE" for row in ac_rows),
            "at least one AC project is incomplete", failures)
    ac_expected = sum(int(row["expected_amplicons"]) for row in ac_rows)
    ac_actual = sum(int(row["classified_amplicons"]) for row in ac_rows)
    fail_if(ac_expected != expected_total or ac_actual != expected_total,
            f"AC count mismatch: expected={ac_expected}, actual={ac_actual}, input={expected_total}", failures)

    current_ac_source_hash = ac_runtime_tree_sha256(ac_source)
    current_config_hash = sha256(ac_source / "ampclasslib" / "default_config.json")
    ac_metadata_verified = 0
    for manifest in ac_rows:
        project_id = manifest["project_id"]
        failure_count_before_project = len(failures)
        metadata_path = args.data_root / "ac_runs" / project_id / "run_metadata.json"
        profile = args.data_root / "ac_runs" / project_id / f"{project_id}_amplicon_classification_profiles.tsv"
        if not metadata_path.is_file():
            failures.append(f"AC metadata missing: {project_id}")
            continue
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        fail_if(metadata.get("status") != "COMPLETE", f"AC metadata incomplete: {project_id}", failures)
        fail_if(int(metadata.get("exit_code", -1)) != 0, f"AC exit code is not zero: {project_id}", failures)
        fail_if(not metadata.get("started_at_utc") or not metadata.get("finished_at_utc"),
                f"AC timestamps missing: {project_id}", failures)
        fail_if(metadata.get("ac_source_sha256") != current_ac_source_hash,
                f"AC source fingerprint mismatch: {project_id}", failures)
        fail_if(metadata.get("default_config_sha256") != current_config_hash,
                f"AC config fingerprint mismatch: {project_id}", failures)
        fail_if(metadata.get("profile_sha256") != sha256(profile),
                f"AC profile hash mismatch: {project_id}", failures)
        snapshot = snapshot_by_project.get(project_id, {})
        fail_if(metadata.get("input_snapshot_sha256") != snapshot.get("snapshot_sha256"),
                f"AC input fingerprint mismatch: {project_id}", failures)
        fingerprint_payload = {
            "project_id": project_id,
            "reference_genome": metadata["reference_genome"],
            "bfbarchitect_mode": metadata["bfbarchitect_mode"],
            "bfb_solver": metadata["bfb_solver"],
            "jobs": metadata["jobs"],
            "input_snapshot_sha256": metadata["input_snapshot_sha256"],
            "reference_snapshot_sha256": metadata["reference_snapshot_sha256"],
            "ac_source_sha256": metadata["ac_source_sha256"],
            "default_config_sha256": metadata["default_config_sha256"],
            "upstream_commit": metadata["upstream_commit"],
            "pipeline_sha256": metadata["pipeline_sha256"],
            "python_executable": metadata["python_executable"],
        }
        fail_if(canonical_sha256(fingerprint_payload) != metadata.get("run_fingerprint"),
                f"AC run fingerprint is internally inconsistent: {project_id}", failures)
        if len(failures) == failure_count_before_project:
            ac_metadata_verified += 1

    fail_if(len(lp_rows) != expected_projects,
            f"LP project count {len(lp_rows)} != {expected_projects}", failures)
    fail_if(any(row["status"] != "COMPLETE" for row in lp_rows),
            "at least one LP project is incomplete", failures)
    lp_expected = sum(int(row["expected_amplicons"]) for row in lp_rows)
    lp_actual = sum(int(row["solved_amplicons"]) for row in lp_rows)
    fail_if(lp_expected != expected_total or lp_actual != expected_total,
            f"LP count mismatch: expected={lp_expected}, actual={lp_actual}, input={expected_total}", failures)

    current_lp_source_hash = tree_sha256(repository_root / "source")
    current_lp_pipeline_hash = tree_sha256(
        code_root, ["run_lwcn_and_merge.py", "pipeline_provenance.py"]
    )
    checkpoint_rows: list[dict[str, str]] = []
    lp_status_counts: Counter[str] = Counter()
    max_balance = max_lower = max_upper = max_gap = max_gap_ratio = 0.0
    parse_info_rows = parse_warning_rows = 0
    source_files_verified = 0
    for project in dataset_rows:
        project_id = project["project_id"]
        checkpoint = args.data_root / "lp_runs" / f"{project_id}.csv"
        metadata_path = args.data_root / "lp_runs" / f"{project_id}.metadata.json"
        if not checkpoint.is_file() or not metadata_path.is_file():
            failures.append(f"LP checkpoint or metadata missing: {project_id}")
            continue
        rows = read_csv(checkpoint)
        checkpoint_rows.extend(rows)
        fail_if(len(rows) != int(project["paired_amplicon_count"]),
                f"LP checkpoint count mismatch: {project_id}", failures)
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        fail_if(metadata.get("status") != "COMPLETE", f"LP metadata incomplete: {project_id}", failures)
        fail_if(int(metadata.get("exit_code", -1)) != 0, f"LP exit code is not zero: {project_id}", failures)
        fail_if(metadata.get("lp_source_sha256") != current_lp_source_hash,
                f"LP source fingerprint mismatch: {project_id}", failures)
        fail_if(metadata.get("pipeline_sha256") != current_lp_pipeline_hash,
                f"LP pipeline fingerprint mismatch: {project_id}", failures)
        fail_if(metadata.get("checkpoint_sha256") != sha256(checkpoint),
                f"LP checkpoint hash mismatch: {project_id}", failures)
        fail_if(metadata.get("input_snapshot_sha256") != snapshot_by_project[project_id]["snapshot_sha256"],
                f"LP input fingerprint mismatch: {project_id}", failures)
        lp_payload = {
            "project_id": project_id,
            "input_snapshot_sha256": metadata["input_snapshot_sha256"],
            "ac_profile_sha256": metadata["ac_profile_sha256"],
            "ac_run_fingerprint": metadata["ac_run_fingerprint"],
            "lp_source_sha256": metadata["lp_source_sha256"],
            "pipeline_sha256": metadata["pipeline_sha256"],
            "feasibility_tolerance": metadata["feasibility_tolerance"],
        }
        fail_if(canonical_sha256(lp_payload) != metadata.get("lp_run_fingerprint"),
                f"LP run fingerprint is internally inconsistent: {project_id}", failures)
        pair_map = {
            (row["sample"], row["amplicon"]): row
            for row in snapshot_by_project[project_id]["pairs"]
        }
        for row in rows:
            key = (row["sample"], row["amplicon"])
            pair = pair_map.get(key)
            if pair is None:
                failures.append(f"checkpoint identity missing from input snapshot: {project_id}/{key}")
                continue
            # ``paired_input_snapshot`` above has already re-read and hashed both
            # source files in this same verification run.  Reuse those freshly
            # computed values here instead of reading all 56,284 files twice.
            graph_hash = pair["graph_sha256"]
            cycles_hash = pair["cycles_sha256"]
            fail_if(row["graph_sha256"] != graph_hash or row["cycles_sha256"] != cycles_hash,
                    f"source file hash mismatch: {project_id}/{key}", failures)
            if row["graph_sha256"] == graph_hash and row["cycles_sha256"] == cycles_hash:
                source_files_verified += 2
            lp_status_counts[row["lp_status"]] += 1
            fail_if(row["lp_status"] not in ACCEPTED_LP_STATUSES,
                    f"non-optimal LP status: {project_id}/{key}", failures)
            parse_info_rows += int(row.get("parse_info_count") or 0) > 0
            parse_warning_rows += int(row.get("parse_warning_count") or 0) > 0
            balance = float(row["max_balance_residual"])
            lower = float(row["max_lower_bound_violation"])
            upper = float(row["max_upper_bound_violation"])
            gap = float(row["certified_gap"])
            max_balance, max_lower, max_upper, max_gap = (
                max(max_balance, balance), max(max_lower, lower),
                max(max_upper, upper), max(max_gap, abs(gap)),
            )
            for field, value in (
                ("balance", balance), ("lower", lower), ("upper", upper)
            ):
                fail_if(not math.isfinite(value) or value > FEASIBILITY_TOLERANCE,
                        f"{field} residual exceeds tolerance: {project_id}/{key}={value}", failures)
            total = float(row["total_length_weighted_copy_number"])
            maximum = float(row["maximum_cyclic_length_weighted_copy_number"])
            ratio = float(row["lwcn"])
            gap_ratio = abs(gap) / max(1.0, abs(total), abs(maximum))
            max_gap_ratio = max(max_gap_ratio, gap_ratio)
            fail_if(not math.isfinite(gap) or gap < 0.0 or gap_ratio > DUALITY_GAP_RATIO_TOLERANCE,
                    f"relative duality gap exceeds tolerance: {project_id}/{key}={gap_ratio}", failures)
            expected_ratio = maximum / total if total > 0 else 0.0
            fail_if(not (math.isfinite(total) and math.isfinite(maximum) and math.isfinite(ratio)),
                    f"non-finite LP objective: {project_id}/{key}", failures)
            fail_if(maximum < -FEASIBILITY_TOLERANCE or maximum > total + FEASIBILITY_TOLERANCE,
                    f"absolute cyclic LWCN outside capacity: {project_id}/{key}", failures)
            fail_if(ratio < -FEASIBILITY_TOLERANCE or ratio > 1 + FEASIBILITY_TOLERANCE,
                    f"cyclic ratio outside [0,1]: {project_id}/{key}", failures)
            fail_if(abs(ratio - expected_ratio) > 2e-9,
                    f"cyclic ratio does not match numerator/denominator: {project_id}/{key}", failures)

    fail_if(len(checkpoint_rows) != expected_total,
            f"checkpoint row count {len(checkpoint_rows)} != {expected_total}", failures)
    fail_if(len(extended_rows) != expected_total,
            f"extended row count {len(extended_rows)} != {expected_total}", failures)
    fail_if(len(final_rows) != expected_total,
            f"four-column row count {len(final_rows)} != {expected_total}", failures)
    required_extended = {
        "project_id", "sample", "project_scoped_sample", "amplicon", "graph_sha256",
        "cycles_sha256", "maximum_cyclic_length_weighted_copy_number", "lwcn",
        "classification", "ecDNA+", "BFB+", "FAN+", "lp_status", "parse_warning_types",
    }
    fail_if(not required_extended.issubset(extended_header),
            f"extended table missing columns: {sorted(required_extended - set(extended_header))}", failures)

    checkpoint_map = {
        (row["project_id"], row["sample"], row["amplicon"]): row for row in checkpoint_rows
    }
    extended_map = {
        (row["project_id"], row["sample"], row["amplicon"]): row for row in extended_rows
    }
    fail_if(len(checkpoint_map) != expected_total, "checkpoint project-scoped identities are not unique", failures)
    fail_if(len(extended_map) != expected_total, "extended project-scoped identities are not unique", failures)
    fail_if(set(checkpoint_map) != set(extended_map), "checkpoint and extended identities differ", failures)
    common_columns = set(extended_header) & set(checkpoint_rows[0] if checkpoint_rows else {})
    for key in set(checkpoint_map) & set(extended_map):
        if any(checkpoint_map[key].get(column, "") != extended_map[key].get(column, "") for column in common_columns):
            failures.append(f"extended table differs from checkpoint: {key}")
            break

    final_map = {(row["sample"], row["amplicon"]): row for row in final_rows}
    fail_if(len(final_map) != expected_total, "four-column project-scoped identities are not unique", failures)
    for row in extended_rows:
        final = final_map.get((row["project_scoped_sample"], row["amplicon"]))
        if final is None:
            failures.append(f"four-column row missing: {row['project_scoped_sample']}/{row['amplicon']}")
            break
        if final["lwcn"] != row["lwcn"] or final["classification"] != row["classification"]:
            failures.append(f"four-column value mismatch: {row['project_scoped_sample']}/{row['amplicon']}")
            break
        fail_if(row["classification"] not in ALLOWED_CLASSIFICATIONS,
                f"unexpected decomposition class: {row['classification']}", failures)

    parser_audit = json.loads((args.output_dir / "解析信息审计.json").read_text(encoding="utf-8"))
    fail_if(int(parser_audit["record_count"]) != expected_total,
            "parser audit row count mismatch", failures)
    fail_if(int(parser_audit["unsupported_format_records"]) != 0,
            "parser audit contains unsupported formats", failures)
    audit_warning_rows = int(parser_audit["records_with_warning"])
    fail_if(audit_warning_rows != parse_warning_rows,
            "parser warning count differs between audit and checkpoints", failures)

    classification_counts = Counter(row["classification"] for row in final_rows)
    summary = {
        "check_type": "automated_consistency_check",
        "verified_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "PASSED" if not failures else "FAILED",
        "public_project_count": expected_projects,
        "archive_total_bytes": archive_total,
        "archive_sha256_verified_count": archive_hashes_verified,
        "paired_amplicon_record_count": expected_total,
        "project_scoped_identity_unique": len(final_map) == expected_total,
        "ac_classified_count": ac_actual,
        "ac_metadata_verified_count": ac_metadata_verified,
        "lp_optimal_count": len(checkpoint_rows),
        "lp_status_counts": dict(sorted(lp_status_counts.items())),
        "source_graph_and_cycles_sha256_verified_count": source_files_verified,
        "parse_info_rows": parse_info_rows,
        "parse_warning_rows": parse_warning_rows,
        "final_csv_row_count": len(final_rows),
        "extended_csv_row_count": len(extended_rows),
        "max_balance_residual": max_balance,
        "max_lower_bound_violation": max_lower,
        "max_upper_bound_violation": max_upper,
        "max_certified_gap": max_gap,
        "max_certified_gap_ratio": max_gap_ratio,
        "feasibility_tolerance": FEASIBILITY_TOLERANCE,
        "duality_gap_ratio_tolerance": DUALITY_GAP_RATIO_TOLERANCE,
        "classification_semantics": "AmpliconClassifier amplicon_decomposition_class",
        "lwcn_semantics": "maximum cyclic length-weighted copy-number ratio",
        "classification_counts": dict(sorted(classification_counts.items())),
        "failures": failures,
    }
    json_path = args.output_dir / "全量自动一致性检查.json"
    text_path = args.output_dir / "全量自动一致性检查.txt"
    json_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    message = (
        f"TASK2 FULL AUTOMATED CONSISTENCY CHECK {summary['status']}\n"
        f"projects={expected_projects}; project_scoped_records={expected_total}; "
        f"AC={ac_actual}; LP={len(checkpoint_rows)}; CSV={len(final_rows)}\n"
        f"archive_hashes={archive_hashes_verified}/{expected_projects}; "
        f"source_file_hashes={source_files_verified}/{expected_total * 2}; "
        f"parser_warning_rows={parse_warning_rows}\n"
        f"max_balance_residual={max_balance:.6g}; max_lower_bound_violation={max_lower:.6g}; "
        f"max_upper_bound_violation={max_upper:.6g}; max_certified_gap_ratio={max_gap_ratio:.6g}\n"
    )
    if failures:
        message += "failures=" + " | ".join(failures[:50]) + "\n"
    text_path.write_text(message, encoding="utf-8")
    print(message, end="")
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
