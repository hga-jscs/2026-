"""Solve the maximum cyclic LWCN ratio and build traceable task2 result tables.

The required four-column CSV keeps exactly ``sample,amplicon,lwcn,classification``.
``lwcn`` is the maximum cyclic length-weighted copy-number ratio requested by the
user. To prevent cross-project collisions, its ``sample`` value is scoped as
``project_id::original_sample``. A separate provenance table retains original
sample names, source hashes, AC mechanism fields, the absolute LP numerator and
all parser/solver evidence.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ALGORITHM_SOURCE = Path(__file__).resolve().parents[1] / "source"
if str(ALGORITHM_SOURCE) not in sys.path:
    sys.path.insert(0, str(ALGORITHM_SOURCE))

from cyclic_lwcn.model import BreakpointEdge, BreakpointGraph, SequenceEdge
from cyclic_lwcn.parser import CoRALParseError, _as_nonnegative_float, _infer_endnodes, parse_node
from original_graph_lwcn.original_graph_linear_program import solve_original_graph_linear_program

from pipeline_provenance import (
    canonical_sha256,
    filesystem_path,
    paired_input_snapshot,
    sha256,
    tree_sha256,
)


FEASIBILITY_TOLERANCE = 1e-7
ACCEPTED_LP_STATUSES = {"OPTIMAL", "TRIVIAL_OPTIMAL_ZERO"}
ALLOWED_CLASSIFICATIONS = {
    "Cyclic",
    "Complex-non-cyclic",
    "Linear",
    "No-FSCNA",
    "Invalid",
    "Virus",
}
MESSAGE_PATTERN = re.compile(r"^(INFO|WARNING)\|([A-Z0-9_]+)\|(.*)$")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_ac_profile(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    """Read all AC fields and reject duplicate project-local identities."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "sample_name",
        "amplicon_number",
        "amplicon_decomposition_class",
        "ecDNA+",
        "BFB+",
        "FAN+",
        "ecDNA_amplicons",
    }
    fields = set(rows[0]) if rows else set()
    if not required.issubset(fields):
        raise ValueError(f"missing required AC columns in {path}: {sorted(required - fields)}")
    result: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["sample_name"], row["amplicon_number"])
        if key in result:
            raise ValueError(f"duplicate AC identity in {path}: {key}")
        if row["amplicon_decomposition_class"] not in ALLOWED_CLASSIFICATIONS:
            raise ValueError(
                f"unexpected AC decomposition class in {path}: "
                f"{row['amplicon_decomposition_class']!r}"
            )
        result[key] = row
    return result


def as_nonnegative_compatible(
    token: str,
    label: str,
    line_number: int,
    messages: list[str],
) -> float:
    """Clamp only solver-scale numerical noise and record it as a warning."""

    try:
        value = float(token)
    except ValueError:
        return _as_nonnegative_float(token, label, line_number)
    if -FEASIBILITY_TOLERANCE <= value < 0.0:
        messages.append(
            "WARNING|NUMERICAL_NOISE_CLAMP|"
            f"line={line_number}; field={label}; original={value:.12g}; replacement=0"
        )
        return 0.0
    return _as_nonnegative_float(token, label, line_number)


def parse_graph_file_compatible(path: Path) -> BreakpointGraph:
    """Strictly map supported CoRAL/AA records to the LP graph.

    Supported AA ``source`` records are validated and then omitted because cycle
    load on an edge incident to CoRAL's source nodes is fixed to zero. Known
    headers and path constraints are explicitly categorized. Any other record
    type raises ``UNSUPPORTED_FORMAT`` instead of being silently ignored.
    """

    source_path = path.resolve()
    sequence_edges: list[SequenceEdge] = []
    breakpoint_edges: list[BreakpointEdge] = []
    nodes = {}
    sequence_by_node: dict[str, SequenceEdge] = {}
    type_counts = {"concordant": 0, "discordant": 0}
    messages: list[str] = []
    aa_length_rows = 0
    source_edge_rows = 0
    path_constraint_rows = 0
    legacy_header_rows = 0
    comma_prefixed_header_rows = 0
    with source_path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(("SequenceEdge:", "BreakpointEdge:", "PathConstraint:")):
                continue
            if line.startswith((",SequenceEdge:", ",BreakpointEdge:", ",PathConstraint:")):
                comma_prefixed_header_rows += 1
                continue
            if line.startswith(("Sequence edge:", "Breakpoint edge:")):
                legacy_header_rows += 1
                continue
            fields = line.split()
            record_type = fields[0]
            if record_type == "sequence":
                if len(fields) < 7:
                    raise CoRALParseError(
                        f"UNSUPPORTED_FORMAT line {line_number}: sequence record needs at least 7 fields"
                    )
                start_node = parse_node(fields[1])
                end_node = parse_node(fields[2])
                copy_number = as_nonnegative_compatible(
                    fields[3], "sequence CN", line_number, messages
                )
                try:
                    declared_length = int(fields[5])
                except ValueError as exc:
                    raise CoRALParseError(
                        f"UNSUPPORTED_FORMAT line {line_number}: sequence length is not an integer"
                    ) from exc
                coordinate_length = end_node.position - start_node.position + 1
                if start_node.chromosome != end_node.chromosome or coordinate_length <= 0:
                    raise CoRALParseError(
                        f"UNSUPPORTED_FORMAT line {line_number}: invalid sequence interval "
                        f"{fields[1]}..{fields[2]}"
                    )
                if declared_length == coordinate_length - 1:
                    aa_length_rows += 1
                elif declared_length != coordinate_length:
                    raise CoRALParseError(
                        f"UNSUPPORTED_FORMAT line {line_number}: declared length {declared_length} "
                        f"matches neither inclusive ({coordinate_length}) nor AA "
                        f"end-start ({coordinate_length - 1}) convention"
                    )
                edge = SequenceEdge(
                    edge_id=f"e{len(sequence_edges) + 1}",
                    start_node=start_node.raw,
                    end_node=end_node.raw,
                    chromosome=start_node.chromosome,
                    start=start_node.position,
                    end=end_node.position,
                    copy_number=copy_number,
                    length=declared_length,
                )
                for node in (start_node, end_node):
                    if node.raw in sequence_by_node:
                        raise CoRALParseError(
                            f"UNSUPPORTED_FORMAT line {line_number}: node {node.raw} "
                            "belongs to multiple sequence edges"
                        )
                    nodes[node.raw] = node
                    sequence_by_node[node.raw] = edge
                sequence_edges.append(edge)
            elif record_type in type_counts:
                if len(fields) < 3 or "->" not in fields[1]:
                    raise CoRALParseError(
                        f"UNSUPPORTED_FORMAT line {line_number}: malformed {record_type} record"
                    )
                left, right = fields[1].split("->", maxsplit=1)
                node1, node2 = parse_node(left), parse_node(right)
                copy_number = as_nonnegative_compatible(
                    fields[2], "breakpoint CN", line_number, messages
                )
                type_counts[record_type] += 1
                prefix = "c" if record_type == "concordant" else "d"
                breakpoint_edges.append(
                    BreakpointEdge(
                        edge_id=f"{prefix}{type_counts[record_type]}",
                        edge_type=record_type,
                        node1=node1.raw,
                        node2=node2.raw,
                        copy_number=copy_number,
                    )
                )
                nodes.setdefault(node1.raw, node1)
                nodes.setdefault(node2.raw, node2)
            elif record_type == "source":
                if len(fields) < 3 or "->" not in fields[1]:
                    raise CoRALParseError(
                        f"UNSUPPORTED_FORMAT line {line_number}: malformed source record"
                    )
                left, right = fields[1].split("->", maxsplit=1)
                if re.fullmatch(r".+:-1[+-]", left) is None:
                    raise CoRALParseError(
                        f"UNSUPPORTED_FORMAT line {line_number}: source endpoint is not coordinate -1"
                    )
                parse_node(right)
                as_nonnegative_compatible(fields[2], "source-edge CN", line_number, messages)
                source_edge_rows += 1
            elif record_type == "path_constraint":
                path_constraint_rows += 1
            else:
                raise CoRALParseError(
                    f"UNSUPPORTED_FORMAT line {line_number}: unrecognized record {record_type!r}"
                )

    if source_edge_rows:
        messages.append(
            "INFO|SOURCE_EDGE_RECORDS|"
            f"count={source_edge_rows}; validated AA source connections; cyclic load fixed to zero"
        )
    if aa_length_rows:
        messages.append(
            "INFO|AA_LENGTH_CONVENTION|"
            f"count={aa_length_rows}; used declared AA end-start genomic lengths"
        )
    if path_constraint_rows:
        messages.append(
            "INFO|PATH_CONSTRAINT_RECORDS|"
            f"count={path_constraint_rows}; intentionally excluded from the graph-capacity LP"
        )
    if legacy_header_rows:
        messages.append(
            "INFO|LEGACY_AA_HEADERS|"
            f"count={legacy_header_rows}; recognized legacy human-readable headers"
        )
    if comma_prefixed_header_rows:
        messages.append(
            "INFO|AA_COMMA_PREFIXED_HEADERS|"
            f"count={comma_prefixed_header_rows}; recognized comma-prefixed AA headers"
        )
    if not sequence_edges:
        if not breakpoint_edges and not nodes:
            messages.append(
                "INFO|EMPTY_GRAPH|count=1; no edge records; exact maximum cyclic ratio is zero"
            )
            return BreakpointGraph(
                source_path=source_path,
                nodes={},
                sequence_edges=(),
                breakpoint_edges=(),
                sequence_by_node={},
                endnodes=frozenset(),
                interval_count=0,
                warnings=messages,
            )
        raise CoRALParseError(f"UNSUPPORTED_FORMAT: no sequence edges found in {source_path}")
    unknown_nodes = sorted(set(nodes) - set(sequence_by_node))
    if unknown_nodes:
        messages.append(
            "WARNING|UNMATCHED_BREAKPOINT_ENDPOINT|"
            f"count={len(unknown_nodes)}; circular load forced to zero at endpoints="
            + ",".join(unknown_nodes[:5])
        )
    endnodes, interval_count = _infer_endnodes(sequence_edges)
    return BreakpointGraph(
        source_path=source_path,
        nodes=nodes,
        sequence_edges=tuple(sequence_edges),
        breakpoint_edges=tuple(breakpoint_edges),
        sequence_by_node=sequence_by_node,
        endnodes=endnodes,
        interval_count=interval_count,
        warnings=messages,
    )


def message_fields(messages: list[str]) -> dict[str, object]:
    """Summarize parser information without conflating INFO and WARNING."""

    info_types: set[str] = set()
    warning_types: set[str] = set()
    for message in messages:
        match = MESSAGE_PATTERN.fullmatch(message)
        if match is None:
            raise ValueError(f"unstructured parser message: {message}")
        severity, code, _ = match.groups()
        (warning_types if severity == "WARNING" else info_types).add(code)
    return {
        "parse_info_count": sum(message.startswith("INFO|") for message in messages),
        "parse_info_types": ";".join(sorted(info_types)),
        "parse_warning_count": sum(message.startswith("WARNING|") for message in messages),
        "parse_warning_types": ";".join(sorted(warning_types)),
        "highest_parse_severity": "WARNING" if warning_types else "INFO" if info_types else "NONE",
        "parse_messages": " | ".join(messages),
    }


def valid_violation(value: float | None) -> bool:
    return value is not None and math.isfinite(value) and value <= FEASIBILITY_TOLERANCE


def solve_one_graph(task: dict[str, object]) -> tuple[dict[str, object] | None, str | None]:
    """Solve one graph and return both result semantics and provenance."""

    graph_path = Path(str(task["graph_path"]))
    try:
        graph = parse_graph_file_compatible(graph_path)
        total_lwcn = sum(edge.copy_number * edge.length for edge in graph.sequence_edges)
        if not graph.sequence_edges and not graph.breakpoint_edges:
            lp_status = "TRIVIAL_OPTIMAL_ZERO"
            maximum_lwcn = ratio = balance = lower = upper = gap = 0.0
        else:
            result = solve_original_graph_linear_program(graph)
            if result.status != "OPTIMAL":
                raise RuntimeError(f"LP status={result.status}: {result.solver_message}")
            if result.maximum_cyclic_lwcn is None or result.maximum_cyclic_ratio is None:
                raise ValueError("LP did not return both absolute and ratio objectives")
            maximum_lwcn = float(result.maximum_cyclic_lwcn)
            ratio = float(result.maximum_cyclic_ratio)
            if not math.isfinite(maximum_lwcn) or maximum_lwcn < -FEASIBILITY_TOLERANCE:
                raise ValueError(f"invalid maximum cyclic LWCN: {maximum_lwcn}")
            if not math.isfinite(ratio) or ratio < -FEASIBILITY_TOLERANCE or ratio > 1 + FEASIBILITY_TOLERANCE:
                raise ValueError(f"maximum cyclic ratio is outside [0,1]: {ratio}")
            if not valid_violation(result.max_balance_residual):
                raise ValueError(f"balance residual={result.max_balance_residual}")
            if not valid_violation(result.max_lower_bound_violation):
                raise ValueError(f"lower-bound violation={result.max_lower_bound_violation}")
            if not valid_violation(result.max_upper_bound_violation):
                raise ValueError(f"upper-bound violation={result.max_upper_bound_violation}")
            lp_status = result.status
            balance = float(result.max_balance_residual)
            lower = float(result.max_lower_bound_violation)
            upper = float(result.max_upper_bound_violation)
            gap = float(result.certified_gap or 0.0)
        ratio = min(1.0, max(0.0, ratio))
        maximum_lwcn = max(0.0, maximum_lwcn)
        ac_row = dict(task["ac_row"])
        sample = str(task["sample"])
        project_id = str(task["project_id"])
        if "::" in sample:
            raise ValueError("original sample name contains reserved project-scope delimiter '::'")
        row: dict[str, object] = {
            "project_id": project_id,
            "project_name": task["project_name"],
            "sample": sample,
            "project_scoped_sample": f"{project_id}::{sample}",
            "amplicon": task["amplicon"],
            "graph_relative_path": task["graph_relative_path"],
            "graph_sha256": task["graph_sha256"],
            "cycles_relative_path": task["cycles_relative_path"],
            "cycles_sha256": task["cycles_sha256"],
            "total_length_weighted_copy_number": format(total_lwcn, ".12g"),
            "maximum_cyclic_length_weighted_copy_number": format(maximum_lwcn, ".12g"),
            "lwcn": format(ratio, ".12g"),
            "classification": ac_row["amplicon_decomposition_class"],
            "ecDNA+": ac_row["ecDNA+"],
            "BFB+": ac_row["BFB+"],
            "FAN+": ac_row["FAN+"],
            "ecDNA_amplicons": ac_row["ecDNA_amplicons"],
            "lp_status": lp_status,
            "max_balance_residual": balance,
            "max_lower_bound_violation": lower,
            "max_upper_bound_violation": upper,
            "certified_gap": gap,
            "input_snapshot_sha256": task["input_snapshot_sha256"],
            "ac_run_fingerprint": task["ac_run_fingerprint"],
            "ac_source_sha256": task["ac_source_sha256"],
            "ac_default_config_sha256": task["ac_default_config_sha256"],
            "ac_upstream_commit": task["ac_upstream_commit"],
            "bfbarchitect_mode": task["bfbarchitect_mode"],
            "lp_source_sha256": task["lp_source_sha256"],
        }
        row.update(message_fields(graph.warnings))
        return row, None
    except Exception as exc:
        return None, f"{graph_path}: {type(exc).__name__}: {exc}"


CHECKPOINT_FIELDS = [
    "project_id", "project_name", "sample", "project_scoped_sample", "amplicon",
    "graph_relative_path", "graph_sha256", "cycles_relative_path", "cycles_sha256",
    "total_length_weighted_copy_number", "maximum_cyclic_length_weighted_copy_number",
    "lwcn", "classification", "ecDNA+", "BFB+", "FAN+", "ecDNA_amplicons",
    "lp_status", "max_balance_residual", "max_lower_bound_violation",
    "max_upper_bound_violation", "certified_gap", "parse_info_count", "parse_info_types",
    "parse_warning_count", "parse_warning_types", "highest_parse_severity", "parse_messages",
    "input_snapshot_sha256", "ac_run_fingerprint", "ac_source_sha256",
    "ac_default_config_sha256", "ac_upstream_commit", "bfbarchitect_mode", "lp_source_sha256",
]


def write_checkpoint(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CHECKPOINT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_checkpoint(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def solve_project(
    project: dict[str, object],
    project_input: Path,
    profile: Path,
    snapshot: dict[str, object],
    ac_metadata: dict[str, object],
    lp_source_sha256: str,
    jobs: int,
) -> tuple[list[dict[str, object]], list[str], float]:
    """Build one-to-one tasks from the frozen input snapshot and solve them."""

    classifications = read_ac_profile(profile)
    tasks: list[dict[str, object]] = []
    errors: list[str] = []
    seen: set[tuple[str, str]] = set()
    for pair in list(snapshot["pairs"]):
        key = (str(pair["sample"]), str(pair["amplicon"]))
        if key in seen:
            errors.append(f"duplicate project-local graph identity: {key}")
            continue
        seen.add(key)
        ac_row = classifications.get(key)
        if ac_row is None:
            errors.append(f"AC classification missing for {key}")
            continue
        tasks.append(
            {
                **pair,
                "project_id": str(project["project_id"]),
                "project_name": str(project["project_name"]),
                "graph_path": str(filesystem_path(project_input) / pair["graph_relative_path"]),
                "ac_row": ac_row,
                "input_snapshot_sha256": snapshot["snapshot_sha256"],
                "ac_run_fingerprint": ac_metadata["run_fingerprint"],
                "ac_source_sha256": ac_metadata["ac_source_sha256"],
                "ac_default_config_sha256": ac_metadata["default_config_sha256"],
                "ac_upstream_commit": ac_metadata["upstream_commit"],
                "bfbarchitect_mode": ac_metadata["bfbarchitect_mode"],
                "lp_source_sha256": lp_source_sha256,
            }
        )
    if set(classifications) != seen:
        extra = sorted(set(classifications) - seen)[:10]
        errors.append(f"AC profile contains identities without selected graph pairs: {extra}")

    rows: list[dict[str, object]] = []
    start = time.perf_counter()
    worker_count = max(1, min(jobs, len(tasks)))
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_count) as pool:
        for index, (row, error) in enumerate(
            pool.map(solve_one_graph, tasks, chunksize=10), start=1
        ):
            if row is not None:
                rows.append(row)
            if error is not None:
                errors.append(error)
            if index % 500 == 0:
                print(
                    f"LP_PROGRESS id={project['project_id']} solved={len(rows)} "
                    f"attempted={index} errors={len(errors)}",
                    flush=True,
                )
    return rows, errors, time.perf_counter() - start


def write_manifest(rows: list[dict[str, object]], output_dir: Path) -> None:
    (output_dir / "lp_run_manifest.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    columns = [
        "project_id", "project_name", "expected_amplicons", "solved_amplicons", "status",
        "execution_status", "started_at_utc", "finished_at_utc", "elapsed_seconds", "exit_code",
        "input_snapshot_sha256", "ac_run_fingerprint", "lp_source_sha256", "lp_run_fingerprint",
        "error_count", "errors",
    ]
    with (output_dir / "lp_run_manifest.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_parser_audit(rows: list[dict[str, object]], output_dir: Path) -> None:
    category_rows: dict[tuple[str, str], set[str]] = defaultdict(set)
    category_occurrences: Counter[tuple[str, str]] = Counter()
    category_examples: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in rows:
        identity = f"{row['project_scoped_sample']}/{row['amplicon']}"
        messages = str(row.get("parse_messages") or "")
        if not messages:
            continue
        for message in messages.split(" | "):
            match = MESSAGE_PATTERN.fullmatch(message)
            if match is None:
                raise ValueError(f"unstructured parser message in checkpoint: {message}")
            severity, code, detail = match.groups()
            key = (severity, code)
            category_rows[key].add(identity)
            count_match = re.search(r"(?:^|;)count=(\d+)", detail)
            category_occurrences[key] += int(count_match.group(1)) if count_match else 1
            if len(category_examples[key]) < 3:
                category_examples[key].append(f"{identity}: {detail}")
    handling = {
        "SOURCE_EDGE_RECORDS": "validated; excluded because source-edge cyclic load is zero",
        "AA_LENGTH_CONVENTION": "accepted declared AA end-start length",
        "PATH_CONSTRAINT_RECORDS": "documented model scope; not used by this LP",
        "LEGACY_AA_HEADERS": "recognized header; no data loss",
        "AA_COMMA_PREFIXED_HEADERS": "recognized AA header variant; no data loss",
        "EMPTY_GRAPH": "exact zero result",
        "NUMERICAL_NOISE_CLAMP": "kept only within 1e-7 and disclosed",
        "UNMATCHED_BREAKPOINT_ENDPOINT": "retained with WARNING; cyclic load forced to zero",
    }
    categories = [
        {
            "severity": key[0],
            "code": key[1],
            "affected_record_count": len(category_rows[key]),
            "record_or_edge_occurrence_count": category_occurrences[key],
            "examples": category_examples[key],
            "handling": handling.get(key[1], "unsupported category"),
        }
        for key in sorted(category_rows)
    ]
    report = {
        "record_count": len(rows),
        "records_with_info": sum(int(row.get("parse_info_count") or 0) > 0 for row in rows),
        "records_with_warning": sum(int(row.get("parse_warning_count") or 0) > 0 for row in rows),
        "unsupported_format_records": 0,
        "categories": categories,
    }
    (output_dir / "解析信息审计.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Solve maximum cyclic LWCN ratios and write scoped/provenance results."
    )
    parser.add_argument("--dataset-manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--force", action="store_true", help="Recompute every LP checkpoint.")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    projects = json.loads(args.dataset_manifest.read_text(encoding="utf-8"))
    projects = [row for row in projects if row.get("status") == "READY"]
    code_root = Path(__file__).resolve().parent
    lp_source_sha256 = tree_sha256(code_root.parent / "source")
    pipeline_sha256 = tree_sha256(code_root, ["run_lwcn_and_merge.py", "pipeline_provenance.py"])
    all_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    total_errors = 0
    print(f"LP_PROJECTS count={len(projects)}", flush=True)

    for index, project in enumerate(projects, start=1):
        project_id = str(project["project_id"])
        expected = int(project["paired_amplicon_count"])
        project_input = args.data_root / "extracted" / project_id
        profile = args.data_root / "ac_runs" / project_id / f"{project_id}_amplicon_classification_profiles.tsv"
        ac_metadata_path = args.data_root / "ac_runs" / project_id / "run_metadata.json"
        checkpoint = args.data_root / "lp_runs" / f"{project_id}.csv"
        metadata_path = args.data_root / "lp_runs" / f"{project_id}.metadata.json"
        print(f"LP_PROJECT {index}/{len(projects)} id={project_id} expected={expected}", flush=True)
        started_at = utc_now()
        errors: list[str] = []
        rows: list[dict[str, object]] = []
        elapsed = 0.0
        execution_status = "RECOMPUTED"
        try:
            if not ac_metadata_path.is_file():
                raise FileNotFoundError(f"AC fingerprint metadata missing: {ac_metadata_path}")
            ac_metadata = json.loads(ac_metadata_path.read_text(encoding="utf-8"))
            if ac_metadata.get("status") != "COMPLETE":
                raise ValueError(f"AC metadata is not COMPLETE: {ac_metadata_path}")
            snapshot = paired_input_snapshot(project_input)
            if int(snapshot["pair_count"]) != expected:
                raise ValueError(f"input snapshot pair count {snapshot['pair_count']} != expected {expected}")
            fingerprint_payload = {
                "project_id": project_id,
                "input_snapshot_sha256": snapshot["snapshot_sha256"],
                "ac_profile_sha256": sha256(profile),
                "ac_run_fingerprint": ac_metadata["run_fingerprint"],
                "lp_source_sha256": lp_source_sha256,
                "pipeline_sha256": pipeline_sha256,
                "feasibility_tolerance": FEASIBILITY_TOLERANCE,
            }
            lp_run_fingerprint = canonical_sha256(fingerprint_payload)
            reuse = False
            if checkpoint.is_file() and metadata_path.is_file() and not args.force:
                previous = json.loads(metadata_path.read_text(encoding="utf-8"))
                candidate_rows = read_checkpoint(checkpoint)
                reuse = (
                    previous.get("lp_run_fingerprint") == lp_run_fingerprint
                    and previous.get("status") == "COMPLETE"
                    and len(candidate_rows) == expected
                    and all(row.get("lp_status") in ACCEPTED_LP_STATUSES for row in candidate_rows)
                )
                if reuse:
                    rows = candidate_rows
                    execution_status = "REUSED_FINGERPRINT_MATCHED"
                    print(f"LP_REUSED id={project_id} rows={len(rows)}", flush=True)
            if not reuse:
                rows, errors, elapsed = solve_project(
                    project, project_input, profile, snapshot, ac_metadata, lp_source_sha256, args.jobs
                )
                write_checkpoint(checkpoint, rows)
            status = "COMPLETE" if len(rows) == expected and not errors else "ERROR"
            metadata = {
                "status": status,
                "project_id": project_id,
                "started_at_utc": started_at,
                "finished_at_utc": utc_now(),
                "elapsed_seconds": round(elapsed, 3),
                "exit_code": 0 if status == "COMPLETE" else 2,
                "execution_status": execution_status,
                "lp_run_fingerprint": lp_run_fingerprint,
                **fingerprint_payload,
                "checkpoint_sha256": sha256(checkpoint) if checkpoint.is_file() else "",
            }
            metadata_path.write_text(
                json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
            )
        except Exception as exc:
            errors.append(f"{type(exc).__name__}: {exc}")
            status = "ERROR"
            lp_run_fingerprint = ""
            snapshot = {"snapshot_sha256": ""}
            ac_metadata = {"run_fingerprint": ""}
        total_errors += len(errors) + (0 if len(rows) == expected else 1)
        manifest_rows.append(
            {
                "project_id": project_id,
                "project_name": str(project["project_name"]),
                "expected_amplicons": expected,
                "solved_amplicons": len(rows),
                "status": status,
                "execution_status": execution_status,
                "started_at_utc": started_at,
                "finished_at_utc": utc_now(),
                "elapsed_seconds": round(elapsed, 3),
                "exit_code": 0 if status == "COMPLETE" else 2,
                "input_snapshot_sha256": snapshot.get("snapshot_sha256", ""),
                "ac_run_fingerprint": ac_metadata.get("run_fingerprint", ""),
                "lp_source_sha256": lp_source_sha256,
                "lp_run_fingerprint": lp_run_fingerprint,
                "error_count": len(errors),
                "errors": " | ".join(errors[:20]),
            }
        )
        all_rows.extend(rows)
        write_manifest(manifest_rows, args.output_dir)
        print(f"LP_DONE id={project_id} rows={len(rows)} errors={len(errors)} seconds={elapsed:.3f}", flush=True)

    all_rows.sort(
        key=lambda row: (
            str(row["project_id"]), str(row["sample"]), str(row["amplicon"]),
            str(row["graph_relative_path"]),
        )
    )
    scoped_keys = [(row["project_scoped_sample"], row["amplicon"]) for row in all_rows]
    duplicate_scoped_keys = [key for key, count in Counter(scoped_keys).items() if count > 1]
    if duplicate_scoped_keys:
        raise RuntimeError(f"duplicate scoped result identities: {duplicate_scoped_keys[:10]}")

    final_path = args.output_dir / "全量AC与环状LWCN结果.csv"
    with final_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", "amplicon", "lwcn", "classification"],
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in all_rows:
            writer.writerow(
                {
                    "sample": row["project_scoped_sample"],
                    "amplicon": row["amplicon"],
                    "lwcn": row["lwcn"],
                    "classification": row["classification"],
                }
            )
    with (args.output_dir / "全量扩展溯源结果.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=CHECKPOINT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)
    write_parser_audit(all_rows, args.output_dir)

    print(
        f"LP_COMPLETE projects={len(manifest_rows)} rows={len(all_rows)} errors={total_errors}",
        flush=True,
    )
    return 0 if total_errors == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
