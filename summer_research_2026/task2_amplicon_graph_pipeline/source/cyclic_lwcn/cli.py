"""Command-line interface for single-file and batch workflows."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .analysis import run_batch_analysis
from .lwcns import calculate_structure_lwcn, check_decomposition_against_graph
from .parser import parse_cycles_file, parse_graph_file
from .state_lp import solve_cyclic_lwcn_upper_bound


def _emit(payload: dict[str, object], output: str | None) -> None:
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    if output:
        Path(output).resolve().write_text(text + "\n", encoding="utf-8")
    print(text)


def _cycle_command(args: argparse.Namespace) -> int:
    decomposition = parse_cycles_file(args.cycles)
    metrics = calculate_structure_lwcn(decomposition)
    payload: dict[str, object] = {
        "cycles_file": str(decomposition.source_path),
        "cycle_lwcn": metrics.cycle_lwcn,
        "path_lwcn": metrics.path_lwcn,
        "structure_lwcn": metrics.structure_lwcn,
        "cycle_count": metrics.cycle_count,
        "path_count": metrics.path_count,
    }
    if args.graph:
        graph = parse_graph_file(args.graph)
        check = check_decomposition_against_graph(graph, decomposition, metrics)
        payload.update(
            {
                "graph_lwcn": check.graph_lwcn,
                "actual_cycle_ratio": (
                    metrics.cycle_lwcn / check.graph_lwcn
                    if check.graph_lwcn > 0
                    else 0.0
                ),
                "decomposition_check": {
                    **check.__dict__,
                    "unmatched_graph_segments": list(check.unmatched_graph_segments),
                    "unmatched_cycle_segments": list(check.unmatched_cycle_segments),
                },
            }
        )
    _emit(payload, args.output)
    return 0


def _upper_command(args: argparse.Namespace) -> int:
    graph = parse_graph_file(args.graph)
    result = solve_cyclic_lwcn_upper_bound(graph)
    payload = {"graph_file": str(graph.source_path), **result.as_dict(args.include_flow)}
    _emit(payload, args.output)
    return 0 if result.status == "OPTIMAL" else 2


def _analyze_command(args: argparse.Namespace) -> int:
    summary = run_batch_analysis(args.data_root, args.output_dir)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["algorithm_error_count"] == 0 else 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cyclic-lwcn",
        description="Calculate CoRAL cycle LWCN and its state-graph LP upper bound.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    cycle = subparsers.add_parser("cycle", help="calculate LWCN from a cycles file")
    cycle.add_argument("--cycles", required=True, help="CoRAL *_cycles.txt")
    cycle.add_argument("--graph", help="optional matching *_graph.txt for validation")
    cycle.add_argument("--output", help="optional JSON output path")
    cycle.set_defaults(handler=_cycle_command)

    upper = subparsers.add_parser(
        "upper-bound", help="solve the graph-only maximum cyclic LWCN"
    )
    upper.add_argument("--graph", required=True, help="CoRAL *_graph.txt")
    upper.add_argument("--output", help="optional JSON output path")
    upper.add_argument(
        "--include-flow", action="store_true", help="include positive state-arc flows"
    )
    upper.set_defaults(handler=_upper_command)

    analyze = subparsers.add_parser(
        "analyze", help="analyze all original CoRAL graph/cycle pairs"
    )
    analyze.add_argument("--data-root", required=True, help="extracted project data root")
    analyze.add_argument("--output-dir", required=True, help="analysis output directory")
    analyze.set_defaults(handler=_analyze_command)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except Exception as exc:  # Deliberate fail-fast CLI boundary; never emits fake values.
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
