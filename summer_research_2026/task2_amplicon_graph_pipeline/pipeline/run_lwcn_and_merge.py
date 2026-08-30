"""对每张公开 graph 求最大环状 LWCN，并与 AC 分类一对一合并。

解析器同时接受 CoRAL inclusive length 与 AA ``end-start`` 长度口径；数值噪声只在
``1e-7`` 容差内夹到零。每个项目写独立检查点，最终 CSV 只包含用户要求的四列。
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

from cyclic_lwcn.model import BreakpointEdge, BreakpointGraph, SequenceEdge
from cyclic_lwcn.parser import (
    CoRALParseError,
    _as_nonnegative_float,
    _infer_endnodes,
    parse_node,
)
from original_graph_lwcn.original_graph_linear_program import (
    solve_original_graph_linear_program,
)

FEASIBILITY_TOLERANCE = 1e-7
ACCEPTED_LP_STATUSES = {"OPTIMAL", "TRIVIAL_OPTIMAL_ZERO"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_ac_profile(path: Path) -> dict[tuple[str, str], str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample_name", "amplicon_number", "amplicon_decomposition_class"}
    if rows and not required.issubset(rows[0]):
        raise ValueError(f"missing required AC columns in {path}")
    return {
        (row["sample_name"], row["amplicon_number"]): row["amplicon_decomposition_class"]
        for row in rows
    }


def graph_identity(graph_path: Path) -> tuple[str, str]:
    suffix = "_graph.txt"
    prefix = graph_path.name[: -len(suffix)]
    marker = prefix.rfind("_amplicon")
    if marker < 1 or not prefix[marker + 1 :].removeprefix("amplicon").isdigit():
        raise ValueError(f"cannot parse sample/amplicon from {graph_path.name}")
    return prefix[:marker], prefix[marker + 1 :]


def paired_graphs(project_input: Path) -> list[Path]:
    graphs: list[Path] = []
    for graph in project_input.rglob("*_graph.txt"):
        lowered_parts = tuple(part.lower() for part in graph.parts)
        if any(
            part == "bfbarchitect_outputs" or part.endswith("_bfbarchitect_outputs")
            for part in lowered_parts
        ):
            continue
        lowered = graph.as_posix().lower()
        if "bpg_converted" in lowered or "/_classification/files/" in lowered:
            continue
        if graph.name.endswith("_features_to_graph.txt") or graph.name.endswith("_feature_to_graph.txt"):
            continue
        cycles = graph.with_name(graph.name[: -len("_graph.txt")] + "_cycles.txt")
        if cycles.is_file():
            graphs.append(graph)
    return sorted(graphs, key=lambda path: path.as_posix())


def valid_number(
    value: float | None,
    upper: float = FEASIBILITY_TOLERANCE,
) -> bool:
    return value is not None and math.isfinite(value) and value <= upper


def as_nonnegative_compatible(
    token: str,
    label: str,
    line_number: int,
    warnings: list[str],
) -> float:
    try:
        value = float(token)
    except ValueError:
        return _as_nonnegative_float(token, label, line_number)
    if -FEASIBILITY_TOLERANCE <= value < 0.0:
        warnings.append(
            f"Clamped numerical-noise {label} {value:.12g} to zero on line {line_number}"
        )
        return 0.0
    return _as_nonnegative_float(token, label, line_number)


def parse_graph_file_compatible(path: Path) -> BreakpointGraph:
    """Parse CoRAL or AA graph length conventions without changing LP semantics.

    CoRAL exports the inclusive coordinate length, while AmpliconArchitect archives
    commonly store ``end - start``.  The LP uses the declared source length in both
    cases; every other parser validation remains strict.
    """

    source_path = path.resolve()
    sequence_edges: list[SequenceEdge] = []
    breakpoint_edges: list[BreakpointEdge] = []
    nodes = {}
    sequence_by_node: dict[str, SequenceEdge] = {}
    type_counts = {"concordant": 0, "discordant": 0}
    warnings: list[str] = []
    aa_length_rows = 0
    with source_path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith(("SequenceEdge:", "BreakpointEdge:", "PathConstraint:")):
                continue
            fields = line.split()
            record_type = fields[0]
            if record_type == "sequence":
                if len(fields) < 7:
                    raise CoRALParseError(
                        f"Line {line_number}: sequence record needs at least 7 fields"
                    )
                start_node = parse_node(fields[1])
                end_node = parse_node(fields[2])
                copy_number = as_nonnegative_compatible(
                    fields[3], "sequence CN", line_number, warnings
                )
                try:
                    declared_length = int(fields[5])
                except ValueError as exc:
                    raise CoRALParseError(
                        f"Line {line_number}: sequence length is not an integer"
                    ) from exc
                coordinate_length = end_node.position - start_node.position + 1
                if start_node.chromosome != end_node.chromosome or coordinate_length <= 0:
                    raise CoRALParseError(
                        f"Line {line_number}: invalid sequence interval {fields[1]}..{fields[2]}"
                    )
                if declared_length == coordinate_length - 1:
                    aa_length_rows += 1
                elif declared_length != coordinate_length:
                    raise CoRALParseError(
                        f"Line {line_number}: declared length {declared_length} matches neither "
                        f"inclusive ({coordinate_length}) nor AA ({coordinate_length - 1}) convention"
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
                            f"Line {line_number}: node {node.raw} belongs to multiple sequence edges"
                        )
                    nodes[node.raw] = node
                    sequence_by_node[node.raw] = edge
                sequence_edges.append(edge)
            elif record_type in type_counts:
                if len(fields) < 3 or "->" not in fields[1]:
                    raise CoRALParseError(f"Line {line_number}: malformed {record_type} record")
                left, right = fields[1].split("->", maxsplit=1)
                node1, node2 = parse_node(left), parse_node(right)
                copy_number = as_nonnegative_compatible(
                    fields[2], "breakpoint CN", line_number, warnings
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
            elif record_type == "path_constraint":
                continue
            else:
                warnings.append(f"Ignored line {line_number}: {line[:100]}")
    if not sequence_edges:
        if not breakpoint_edges and not nodes:
            warnings.append("Empty graph: exact trivial cyclic LWCN is zero")
            return BreakpointGraph(
                source_path=source_path,
                nodes={},
                sequence_edges=(),
                breakpoint_edges=(),
                sequence_by_node={},
                endnodes=frozenset(),
                interval_count=0,
                warnings=warnings,
            )
        raise CoRALParseError(f"No sequence edges found in {source_path}")
    unknown_nodes = sorted(set(nodes) - set(sequence_by_node))
    if unknown_nodes:
        warnings.append(
            f"Kept {len(unknown_nodes)} breakpoint endpoint(s) without sequence incidence; "
            "LP balance forces their circular load to zero: "
            + ", ".join(unknown_nodes[:5])
        )
    if aa_length_rows:
        warnings.append(f"Accepted AA end-start length convention on {aa_length_rows} sequence rows")
    endnodes, interval_count = _infer_endnodes(sequence_edges)
    return BreakpointGraph(
        source_path=source_path,
        nodes=nodes,
        sequence_edges=tuple(sequence_edges),
        breakpoint_edges=tuple(breakpoint_edges),
        sequence_by_node=sequence_by_node,
        endnodes=endnodes,
        interval_count=interval_count,
        warnings=warnings,
    )


def solve_one_graph(
    task: tuple[str, str, str, str, str]
) -> tuple[dict[str, object] | None, str | None]:
    project_id, graph_path_text, sample, amplicon, classification = task
    graph_path = Path(graph_path_text)
    try:
        graph = parse_graph_file_compatible(graph_path)
        if not graph.sequence_edges and not graph.breakpoint_edges:
            return (
                {
                    "sample": sample,
                    "amplicon": amplicon,
                    "lwcn": "0",
                    "classification": classification,
                    "project_id": project_id,
                    "lp_status": "TRIVIAL_OPTIMAL_ZERO",
                    "max_balance_residual": 0.0,
                    "max_lower_bound_violation": 0.0,
                    "max_upper_bound_violation": 0.0,
                    "certified_gap": 0.0,
                    "parse_warning_count": len(graph.warnings),
                    "parse_warnings": " | ".join(graph.warnings),
                },
                None,
            )
        result = solve_original_graph_linear_program(graph)
        if result.status != "OPTIMAL":
            raise RuntimeError(f"LP status={result.status}: {result.solver_message}")
        if result.maximum_cyclic_lwcn is None or not math.isfinite(result.maximum_cyclic_lwcn):
            raise ValueError("LP returned a non-finite cyclic LWCN")
        if not valid_number(result.max_balance_residual):
            raise ValueError(f"balance residual={result.max_balance_residual}")
        if not valid_number(result.max_lower_bound_violation):
            raise ValueError(f"lower-bound violation={result.max_lower_bound_violation}")
        if not valid_number(result.max_upper_bound_violation):
            raise ValueError(f"upper-bound violation={result.max_upper_bound_violation}")
        return (
            {
                "sample": sample,
                "amplicon": amplicon,
                "lwcn": format(result.maximum_cyclic_lwcn, ".12g"),
                "classification": classification,
                "project_id": project_id,
                "lp_status": result.status,
                "max_balance_residual": result.max_balance_residual,
                "max_lower_bound_violation": result.max_lower_bound_violation,
                "max_upper_bound_violation": result.max_upper_bound_violation,
                "certified_gap": result.certified_gap,
                "parse_warning_count": len(graph.warnings),
                "parse_warnings": " | ".join(graph.warnings),
            },
            None,
        )
    except Exception as exc:
        return None, f"{graph_path}: {type(exc).__name__}: {exc}"


def solve_project(
    project_id: str,
    project_input: Path,
    ac_profile: Path,
    checkpoint: Path,
    jobs: int,
) -> tuple[list[dict[str, object]], list[str], float]:
    classifications = read_ac_profile(ac_profile)
    graphs = paired_graphs(project_input)
    rows: list[dict[str, object]] = []
    errors: list[str] = []
    tasks: list[tuple[str, str, str, str, str]] = []
    for graph_path in graphs:
        try:
            sample, amplicon = graph_identity(graph_path)
            classification = classifications.get((sample, amplicon))
            if classification is None:
                raise KeyError(f"AC classification missing for {sample}/{amplicon}")
            tasks.append((project_id, str(graph_path), sample, amplicon, classification))
        except Exception as exc:
            errors.append(f"{graph_path}: {type(exc).__name__}: {exc}")

    start = time.perf_counter()
    worker_count = max(1, min(jobs, len(tasks)))
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_count) as pool:
        outcomes = pool.map(solve_one_graph, tasks, chunksize=10)
        for index, (row, error) in enumerate(outcomes, start=1):
            if row is not None:
                rows.append(row)
            if error is not None:
                errors.append(error)
            if index % 500 == 0:
                print(
                    f"LP_PROGRESS id={project_id} solved={len(rows)} attempted={index} errors={len(errors)}",
                    flush=True,
                )
    elapsed = time.perf_counter() - start
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample",
        "amplicon",
        "lwcn",
        "classification",
        "project_id",
        "lp_status",
        "max_balance_residual",
        "max_lower_bound_violation",
        "max_upper_bound_violation",
        "certified_gap",
        "parse_warning_count",
        "parse_warnings",
    ]
    with checkpoint.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows, errors, elapsed


def read_checkpoint(path: Path) -> list[dict[str, object]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_summary(rows: list[dict[str, object]], output_dir: Path) -> None:
    (output_dir / "lp_run_manifest.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    columns = [
        "project_id",
        "project_name",
        "expected_amplicons",
        "solved_amplicons",
        "status",
        "elapsed_seconds",
        "finished_at_utc",
        "error_count",
        "errors",
    ]
    with (output_dir / "lp_run_manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Solve cyclic-LWCN LPs and merge with AC classifications.")
    parser.add_argument("--dataset-manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--force", action="store_true", help="Recompute all LP checkpoints.")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    projects = json.loads(args.dataset_manifest.read_text(encoding="utf-8"))
    projects = [row for row in projects if row.get("status") == "READY"]
    all_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    total_errors = 0
    print(f"LP_PROJECTS count={len(projects)}", flush=True)

    for index, project in enumerate(projects, start=1):
        project_id = str(project["project_id"])
        project_name = str(project["project_name"])
        expected = int(project["paired_amplicon_count"])
        checkpoint = args.data_root / "lp_runs" / f"{project_id}.csv"
        profile = (
            args.data_root
            / "ac_runs"
            / project_id
            / f"{project_id}_amplicon_classification_profiles.tsv"
        )
        print(f"LP_PROJECT {index}/{len(projects)} id={project_id} expected={expected}", flush=True)
        errors: list[str] = []
        elapsed = 0.0
        if checkpoint.exists() and not args.force:
            rows = read_checkpoint(checkpoint)
            if len(rows) == expected and all(
                row.get("lp_status") in ACCEPTED_LP_STATUSES for row in rows
            ):
                print(f"LP_REUSED id={project_id} rows={len(rows)}", flush=True)
            else:
                rows, errors, elapsed = solve_project(
                    project_id,
                    args.data_root / "extracted" / project_id,
                    profile,
                    checkpoint,
                    args.jobs,
                )
        else:
            rows, errors, elapsed = solve_project(
                project_id,
                args.data_root / "extracted" / project_id,
                profile,
                checkpoint,
                args.jobs,
            )
        status = "COMPLETE" if len(rows) == expected and not errors else "ERROR"
        total_errors += len(errors) + (0 if len(rows) == expected else 1)
        run_rows.append(
            {
                "project_id": project_id,
                "project_name": project_name,
                "expected_amplicons": expected,
                "solved_amplicons": len(rows),
                "status": status,
                "elapsed_seconds": round(elapsed, 3),
                "finished_at_utc": utc_now(),
                "error_count": len(errors),
                "errors": " | ".join(errors[:20]),
            }
        )
        all_rows.extend(rows)
        write_summary(run_rows, args.output_dir)
        print(
            f"LP_DONE id={project_id} rows={len(rows)} errors={len(errors)} seconds={elapsed:.3f}",
            flush=True,
        )

    final_columns = ["sample", "amplicon", "lwcn", "classification"]
    with (args.output_dir / "全量AC与环状LWCN结果.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=final_columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    print(
        f"LP_COMPLETE projects={len(run_rows)} rows={len(all_rows)} errors={total_errors}",
        flush=True,
    )
    return 0 if total_errors == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
