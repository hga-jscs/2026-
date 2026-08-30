"""独立验收全量 AC + 最大环状占比结果；任一计数、状态或残差不符即失败。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

FEASIBILITY_TOLERANCE = 1e-7
ACCEPTED_LP_STATUSES = {"OPTIMAL", "TRIVIAL_OPTIMAL_ZERO"}


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Independent acceptance checks for AC and maximum cyclic-ratio results."
    )
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    failures: list[str] = []
    dataset_rows = read_csv(args.output_dir / "数据集清单.csv")
    ac_rows = read_csv(args.output_dir / "ac_run_manifest.csv")
    lp_rows = read_csv(args.output_dir / "lp_run_manifest.csv")
    result_path = args.output_dir / "全量AC与环状LWCN结果.csv"
    with result_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames
        result_rows = list(reader)

    if header != ["sample", "amplicon", "lwcn", "classification"]:
        failures.append(f"unexpected result header: {header}")
    if len(dataset_rows) != 32:
        failures.append(f"dataset count is {len(dataset_rows)}, expected 32")
    if any(row["status"] != "READY" for row in dataset_rows):
        failures.append("at least one dataset is not READY")
    expected_total = sum(int(row["paired_amplicon_count"]) for row in dataset_rows)
    archive_total = 0
    for row in dataset_rows:
        archive = args.data_root / "archives" / f"{row['project_id']}.tar.gz"
        if not archive.is_file():
            failures.append(f"archive missing: {row['project_id']}")
            continue
        actual_size = archive.stat().st_size
        declared_size = int(row["archive_bytes"])
        archive_total += actual_size
        if actual_size != declared_size:
            failures.append(f"archive size mismatch: {row['project_id']}")
        if len(row["sha256"]) != 64:
            failures.append(f"bad SHA-256 field: {row['project_id']}")

    if len(ac_rows) != 32:
        failures.append(f"AC project count is {len(ac_rows)}, expected 32")
    if any(row["status"] not in {"COMPLETE", "REUSED_COMPLETE"} for row in ac_rows):
        failures.append("at least one AC project is incomplete")
    ac_expected = sum(int(row["expected_amplicons"]) for row in ac_rows)
    ac_actual = sum(int(row["classified_amplicons"]) for row in ac_rows)
    if ac_expected != expected_total or ac_actual != expected_total:
        failures.append(f"AC count mismatch: expected={ac_expected}, actual={ac_actual}, input={expected_total}")

    if len(lp_rows) != 32:
        failures.append(f"LP project count is {len(lp_rows)}, expected 32")
    if any(row["status"] != "COMPLETE" for row in lp_rows):
        failures.append("at least one LP project is incomplete")
    lp_expected = sum(int(row["expected_amplicons"]) for row in lp_rows)
    lp_actual = sum(int(row["solved_amplicons"]) for row in lp_rows)
    if lp_expected != expected_total or lp_actual != expected_total:
        failures.append(f"LP count mismatch: expected={lp_expected}, actual={lp_actual}, input={expected_total}")

    max_balance = 0.0
    max_lower = 0.0
    max_upper = 0.0
    solver_rows = 0
    lp_status_counts: Counter[str] = Counter()
    parse_warning_rows = 0
    for project in dataset_rows:
        checkpoint = args.data_root / "lp_runs" / f"{project['project_id']}.csv"
        if not checkpoint.is_file():
            failures.append(f"LP checkpoint missing: {project['project_id']}")
            continue
        rows = read_csv(checkpoint)
        if len(rows) != int(project["paired_amplicon_count"]):
            failures.append(f"LP checkpoint count mismatch: {project['project_id']}")
        for row in rows:
            solver_rows += 1
            lp_status_counts[row["lp_status"]] += 1
            if row["lp_status"] not in ACCEPTED_LP_STATUSES:
                failures.append(f"non-optimal LP status: {project['project_id']}")
            if int(row.get("parse_warning_count") or 0) > 0:
                parse_warning_rows += 1
            max_balance = max(max_balance, float(row["max_balance_residual"]))
            max_lower = max(max_lower, float(row["max_lower_bound_violation"]))
            max_upper = max(max_upper, float(row["max_upper_bound_violation"]))

    if max_balance > FEASIBILITY_TOLERANCE:
        failures.append(f"balance residual too large: {max_balance}")
    if max_lower > FEASIBILITY_TOLERANCE:
        failures.append(f"lower-bound violation too large: {max_lower}")
    if max_upper > FEASIBILITY_TOLERANCE:
        failures.append(f"upper-bound violation too large: {max_upper}")
    if solver_rows != expected_total:
        failures.append(f"solver row count mismatch: {solver_rows}/{expected_total}")

    for index, row in enumerate(result_rows, start=2):
        if not row["sample"] or not row["amplicon"] or not row["classification"]:
            failures.append(f"blank required field at result row {index}")
            break
        try:
            lwcn = float(row["lwcn"])
            if not math.isfinite(lwcn) or lwcn < -1e-10 or lwcn > 1.0 + 1e-10:
                raise ValueError
        except ValueError:
            failures.append(f"invalid lwcn at result row {index}: {row['lwcn']}")
            break
    if len(result_rows) != expected_total:
        failures.append(f"final CSV row count mismatch: {len(result_rows)}/{expected_total}")

    classification_counts = Counter(row["classification"] for row in result_rows)
    summary = {
        "verified_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "PASSED" if not failures else "FAILED",
        "public_project_count": len(dataset_rows),
        "archive_total_bytes": archive_total,
        "paired_amplicon_count": expected_total,
        "ac_classified_count": ac_actual,
        "lp_optimal_count": solver_rows,
        "lp_status_counts": dict(sorted(lp_status_counts.items())),
        "parse_warning_rows": parse_warning_rows,
        "final_csv_row_count": len(result_rows),
        "max_balance_residual": max_balance,
        "max_lower_bound_violation": max_lower,
        "max_upper_bound_violation": max_upper,
        "feasibility_tolerance": FEASIBILITY_TOLERANCE,
        "classification_counts": dict(sorted(classification_counts.items())),
        "failures": failures,
    }
    (args.output_dir / "全量验收.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    message = (
        f"TASK2 FULL RESULTS VERIFICATION {summary['status']}\n"
        f"projects={len(dataset_rows)}; paired_amplicons={expected_total}; "
        f"AC={ac_actual}; LP_OPTIMAL={solver_rows}; CSV={len(result_rows)}\n"
        f"max_balance_residual={max_balance:.6g}; "
        f"max_lower_bound_violation={max_lower:.6g}; "
        f"max_upper_bound_violation={max_upper:.6g}\n"
    )
    if failures:
        message += "failures=" + " | ".join(failures) + "\n"
    (args.output_dir / "全量验收.txt").write_text(message, encoding="utf-8")
    print(message, end="")
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
