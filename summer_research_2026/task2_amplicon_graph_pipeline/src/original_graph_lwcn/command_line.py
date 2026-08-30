"""原图 LP 的单图求解与全量评测命令行入口。"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from cyclic_lwcn.parser import parse_graph_file

from .benchmarking import run_coral_benchmark
from .original_graph_linear_program import solve_original_graph_linear_program


def _write_json(payload: dict[str, object], output: str | None) -> None:
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    if output:
        output_path = Path(output).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(text + "\n", encoding="utf-8")
    print(text)


def _solve_command(args: argparse.Namespace) -> int:
    graph = parse_graph_file(args.graph)
    result = solve_original_graph_linear_program(graph)
    _write_json(result.as_dict(args.include_edge_load), args.output)
    return 0 if result.status == "OPTIMAL" else 2


def _benchmark_command(args: argparse.Namespace) -> int:
    summary = run_coral_benchmark(
        args.data_root,
        args.output_dir,
        solve_original_graph_linear_program,
        state_match_tolerance_ratio=1e-10,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["algorithm_error_count"] == 0 else 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="original-graph-lwcn",
        description="直接在原始 CoRAL 断点图上求最大环状 LWCN。",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    solve = subparsers.add_parser("solve", help="求解单个 *_graph.txt")
    solve.add_argument("--graph", required=True, help="CoRAL *_graph.txt")
    solve.add_argument("--output", help="可选 JSON 输出文件")
    solve.add_argument(
        "--include-edge-load",
        action="store_true",
        help="在 JSON 中包含每条正环流物理边的 CN",
    )
    solve.set_defaults(handler=_solve_command)

    benchmark = subparsers.add_parser("benchmark", help="运行 112 张官方图全量对照")
    benchmark.add_argument("--data-root", required=True, help="CoRAL_cell_lines 根目录")
    benchmark.add_argument("--output-dir", required=True, help="CSV/JSON/PNG 输出目录")
    benchmark.set_defaults(handler=_benchmark_command)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except Exception as exc:
        # CLI 边界只负责把真实异常转成非零退出码，不返回替代值或空结果。
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
