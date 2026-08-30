"""在 AmpliconRepository CoRAL 数据上统一评测 task2 的三种求解器。"""

from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter
from typing import Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cyclic_lwcn.analysis import coral_reconstruction_status, discover_official_pairs
from cyclic_lwcn.lwcns import calculate_structure_lwcn
from cyclic_lwcn.parser import parse_cycles_file, parse_graph_file
from cyclic_lwcn.state_lp import solve_cyclic_lwcn_upper_bound

from .solver_result import CyclicLWCNSolverResult


Solver = Callable[[object], CyclicLWCNSolverResult]


def _sample_name(graph_path: Path) -> str:
    """从 CoRAL 的 ``*_reconstruction_results`` 目录稳定提取样本名。"""

    for parent in graph_path.parents:
        if parent.name.endswith("_reconstruction_results"):
            return parent.name.removesuffix("_reconstruction_results")
    return graph_path.parent.name


def _analyze_one_pair(
    graph_path: Path,
    cycles_path: Path,
    solver: Solver,
) -> dict[str, object]:
    """求解一张图，并用同一输入重新计算状态图基准与 CoRAL 已提取值。"""

    started = perf_counter()
    graph = parse_graph_file(graph_path)
    decomposition = parse_cycles_file(cycles_path)
    extracted = calculate_structure_lwcn(decomposition)
    reference = solve_cyclic_lwcn_upper_bound(graph)
    if reference.status != "OPTIMAL":
        raise RuntimeError(
            f"状态图基准未达到 OPTIMAL：{graph_path.name}: {reference.solver_message}"
        )

    result = solver(graph)
    official_status = coral_reconstruction_status(graph_path)
    total = float(graph.total_lwcn)
    actual_ratio = extracted.cycle_lwcn / total if total > 0 else 0.0
    method_ratio = result.maximum_cyclic_ratio
    reference_ratio = reference.maximum_cyclic_ratio
    ratio_difference = (
        None
        if method_ratio is None or reference_ratio is None
        else float(method_ratio - reference_ratio)
    )
    actual_excess = (
        None
        if result.maximum_cyclic_lwcn is None
        else float(extracted.cycle_lwcn - result.maximum_cyclic_lwcn)
    )
    actual_minus_upper = (
        None
        if result.cyclic_upper_certificate is None
        else float(extracted.cycle_lwcn - result.cyclic_upper_certificate)
    )
    comparison_eligible = official_status == "SUCCESS" and result.succeeded

    return {
        "sample": _sample_name(graph_path),
        "amplicon": graph_path.name.removesuffix("_graph.txt"),
        "graph_path": str(graph_path.resolve()),
        "cycles_path": str(cycles_path.resolve()),
        "coral_reconstruction_status": official_status,
        "comparison_eligible": comparison_eligible,
        "algorithm": result.algorithm,
        "status": result.status,
        "solver_status": result.solver_status,
        "solver_message": result.solver_message,
        "sequence_edges": len(graph.sequence_edges),
        "breakpoint_edges": len(graph.breakpoint_edges),
        "variables": result.variable_count,
        "equalities": result.equality_count,
        "graph_lwcn": total,
        "actual_cycle_lwcn": float(extracted.cycle_lwcn),
        "actual_cycle_ratio": actual_ratio,
        "method_cycle_lwcn": result.maximum_cyclic_lwcn,
        "method_cycle_ratio": method_ratio,
        "state_lp_cycle_lwcn": reference.maximum_cyclic_lwcn,
        "state_lp_cycle_ratio": reference_ratio,
        "method_minus_state_lwcn": (
            None
            if result.maximum_cyclic_lwcn is None
            or reference.maximum_cyclic_lwcn is None
            else float(result.maximum_cyclic_lwcn - reference.maximum_cyclic_lwcn)
        ),
        "method_minus_state_ratio": ratio_difference,
        "actual_minus_method_lwcn": actual_excess,
        "actual_minus_method_ratio": (
            None if actual_excess is None or total <= 0 else actual_excess / total
        ),
        "actual_minus_upper_certificate_lwcn": actual_minus_upper,
        "actual_minus_upper_certificate_ratio": (
            None
            if actual_minus_upper is None or total <= 0
            else actual_minus_upper / total
        ),
        "solve_seconds": result.solve_seconds,
        "analysis_seconds": perf_counter() - started,
        "max_balance_residual": result.max_balance_residual,
        "max_lower_bound_violation": result.max_lower_bound_violation,
        "max_upper_bound_violation": result.max_upper_bound_violation,
        "cyclic_lower_certificate": result.cyclic_lower_certificate,
        "cyclic_upper_certificate": result.cyclic_upper_certificate,
        "certified_gap": result.certified_gap,
        "certified_gap_ratio": (
            None
            if result.certified_gap is None or total <= 0
            else result.certified_gap / total
        ),
        "diagnostics_json": json.dumps(
            result.diagnostics, ensure_ascii=False, sort_keys=True
        ),
    }


def summarize_benchmark(
    frame: pd.DataFrame,
    *,
    state_match_tolerance_ratio: float,
) -> dict[str, object]:
    """把逐图结果汇总成可机器验收的指标。"""

    succeeded = frame[frame["status"].isin(["OPTIMAL", "CONVERGED", "CONVERGED_CERTIFIED"])].copy()
    comparable = succeeded[succeeded["comparison_eligible"] == True].copy()  # noqa: E712
    ratio_error = pd.to_numeric(
        succeeded["method_minus_state_ratio"], errors="coerce"
    ).abs()
    actual_excess_ratio = pd.to_numeric(
        comparable["actual_minus_method_ratio"], errors="coerce"
    )
    actual_upper_excess_ratio = pd.to_numeric(
        comparable["actual_minus_upper_certificate_ratio"], errors="coerce"
    )
    gap_ratio = pd.to_numeric(succeeded["certified_gap_ratio"], errors="coerce")
    balance = pd.to_numeric(succeeded["max_balance_residual"], errors="coerce")
    lower_violation = pd.to_numeric(
        succeeded["max_lower_bound_violation"], errors="coerce"
    )
    upper_violation = pd.to_numeric(
        succeeded["max_upper_bound_violation"], errors="coerce"
    )
    return {
        "algorithm": str(frame["algorithm"].dropna().iloc[0]) if len(frame) else None,
        "graph_count": int(len(frame)),
        "solved_count": int(len(succeeded)),
        "algorithm_error_count": int(len(frame) - len(succeeded)),
        "coral_reconstruction_success_count": int(
            (frame["coral_reconstruction_status"] == "SUCCESS").sum()
        ),
        "coral_reconstruction_failure_count": int(
            (frame["coral_reconstruction_status"] == "FAILURE").sum()
        ),
        "comparison_count": int(len(comparable)),
        "state_match_tolerance_ratio": state_match_tolerance_ratio,
        "state_lp_match_count": int((ratio_error <= state_match_tolerance_ratio).sum()),
        "max_abs_state_difference_ratio": (
            float(ratio_error.max()) if len(ratio_error.dropna()) else None
        ),
        "median_abs_state_difference_ratio": (
            float(ratio_error.median()) if len(ratio_error.dropna()) else None
        ),
        "actual_exceeds_method_count": int((actual_excess_ratio > 1e-8).sum()),
        "max_actual_minus_method_ratio": (
            float(actual_excess_ratio.max())
            if len(actual_excess_ratio.dropna())
            else None
        ),
        "actual_exceeds_upper_certificate_count": int(
            (actual_upper_excess_ratio > 1e-8).sum()
        ),
        "max_actual_minus_upper_certificate_ratio": (
            float(actual_upper_excess_ratio.max())
            if len(actual_upper_excess_ratio.dropna())
            else None
        ),
        "max_balance_residual": float(balance.max()) if len(balance.dropna()) else None,
        "max_lower_bound_violation": (
            float(lower_violation.max()) if len(lower_violation.dropna()) else None
        ),
        "max_upper_bound_violation": (
            float(upper_violation.max()) if len(upper_violation.dropna()) else None
        ),
        "max_certified_gap_ratio": (
            float(gap_ratio.max()) if len(gap_ratio.dropna()) else None
        ),
        "median_certified_gap_ratio": (
            float(gap_ratio.median()) if len(gap_ratio.dropna()) else None
        ),
        "total_solve_seconds": float(succeeded["solve_seconds"].sum()),
        "median_solve_seconds": float(succeeded["solve_seconds"].median()),
        "max_solve_seconds": float(succeeded["solve_seconds"].max()),
        "status_counts": {
            str(status): int(count)
            for status, count in frame["status"].value_counts().items()
        },
    }


def _write_comparison_figure(frame: pd.DataFrame, output_path: Path) -> None:
    """生成不依赖中文字体的三联图，供三份中文报告共同引用。"""

    ok = frame[
        frame["status"].isin(["OPTIMAL", "CONVERGED", "CONVERGED_CERTIFIED"])
    ].copy()
    error = pd.to_numeric(ok["method_minus_state_ratio"], errors="coerce").dropna()
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})
    figure, axes = plt.subplots(1, 3, figsize=(12, 3.6), constrained_layout=True)
    axes[0].scatter(
        ok["state_lp_cycle_ratio"], ok["method_cycle_ratio"], s=15, alpha=0.7
    )
    axes[0].plot([0, 1], [0, 1], "--", color="black", linewidth=1)
    axes[0].set(
        xlabel="State-graph LP ratio",
        ylabel="Current method ratio",
        title="Agreement",
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
    )
    axes[1].hist(error, bins=20, color="#2f6f9f", edgecolor="white")
    axes[1].set(
        xlabel="Method - state LP ratio",
        ylabel="Amplicons",
        title="Numerical difference",
    )
    axes[2].scatter(
        ok["variables"], ok["solve_seconds"], s=15, alpha=0.7, color="#b45f06"
    )
    axes[2].set(
        xlabel="Original-edge variables",
        ylabel="Solve seconds",
        title="Runtime scaling",
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def run_coral_benchmark(
    data_root: str | Path,
    output_dir: str | Path,
    solver: Solver,
    *,
    state_match_tolerance_ratio: float,
) -> dict[str, object]:
    """对全部官方 graph/cycles 配对运行指定求解器并写出 CSV/JSON/PNG。"""

    pairs = discover_official_pairs(data_root)
    if not pairs:
        raise FileNotFoundError(f"未在 {Path(data_root).resolve()} 下找到 CoRAL 配对")
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for graph_path, cycles_path in pairs:
        try:
            rows.append(_analyze_one_pair(graph_path, cycles_path, solver))
        except Exception as exc:
            # 批量层保留错误记录，但绝不填入虚假数值；最终验收会因 error_count 非零而失败。
            rows.append(
                {
                    "sample": _sample_name(graph_path),
                    "amplicon": graph_path.name.removesuffix("_graph.txt"),
                    "graph_path": str(graph_path.resolve()),
                    "cycles_path": str(cycles_path.resolve()),
                    "coral_reconstruction_status": coral_reconstruction_status(graph_path),
                    "comparison_eligible": False,
                    "algorithm": getattr(solver, "__name__", "unknown_solver"),
                    "status": type(exc).__name__,
                    "solver_status": "EXCEPTION",
                    "solver_message": str(exc),
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(
        output / "per_amplicon_results.csv", index=False, encoding="utf-8-sig"
    )
    issues = frame[~frame["status"].isin(["OPTIMAL", "CONVERGED", "CONVERGED_CERTIFIED"])].copy()
    issues.to_csv(output / "issues.csv", index=False, encoding="utf-8-sig")
    summary = summarize_benchmark(
        frame, state_match_tolerance_ratio=state_match_tolerance_ratio
    )
    (output / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    _write_comparison_figure(frame, output / "method_comparison.png")
    return summary
