"""Batch analysis of official AmpliconRepository CoRAL reconstructions."""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .balance import UnbalancedGraphError
from .lwcns import calculate_structure_lwcn, check_decomposition_against_graph
from .parser import CoRALParseError, parse_cycles_file, parse_graph_file
from .state_lp import solve_cyclic_lwcn_upper_bound


def discover_official_pairs(data_root: str | Path) -> list[tuple[Path, Path]]:
    """Find original CoRAL graph/cycle pairs, excluding annotations and BFB output."""

    root = Path(data_root).resolve()
    pairs: list[tuple[Path, Path]] = []
    for graph_path in root.rglob("*_graph.txt"):
        normalized = str(graph_path).replace("/", "\\").lower()
        if "\\results\\samples\\" not in normalized:
            continue
        if "_reconstruction_results" not in normalized:
            continue
        cycle_path = graph_path.with_name(
            graph_path.name.removesuffix("_graph.txt") + "_cycles.txt"
        )
        if cycle_path.is_file():
            pairs.append((graph_path, cycle_path))
    return sorted(pairs)


def _sample_name(graph_path: Path) -> str:
    for parent in graph_path.parents:
        if parent.name.endswith("_reconstruction_results"):
            return parent.name.removesuffix("_reconstruction_results")
    return graph_path.parent.name


@lru_cache(maxsize=None)
def _summary_status_map(summary_path: Path) -> dict[int, str]:
    statuses: dict[int, str] = {}
    current_amplicon: int | None = None
    for raw_line in summary_path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        match = re.fullmatch(r"AmpliconID\s*=\s*(\d+)", line)
        if match:
            current_amplicon = int(match.group(1))
            continue
        match = re.fullmatch(r"Cycle Decomposition Status:\s*(\S+)", line)
        if match and current_amplicon is not None:
            statuses[current_amplicon] = match.group(1).upper()
    return statuses


def coral_reconstruction_status(graph_path: Path) -> str:
    """Read the authoritative per-amplicon SUCCESS/FAILURE status."""

    match = re.search(r"_amplicon(\d+)_graph\.txt$", graph_path.name)
    if not match:
        return "UNKNOWN"
    sample = _sample_name(graph_path)
    summary_path = graph_path.parent / f"{sample}_summary.txt"
    if not summary_path.is_file():
        return "UNKNOWN"
    return _summary_status_map(summary_path).get(int(match.group(1)), "UNKNOWN")


def analyze_pair(graph_path: Path, cycle_path: Path) -> dict[str, object]:
    started = perf_counter()
    sample = _sample_name(graph_path)
    amplicon = graph_path.name.removesuffix("_graph.txt")
    base: dict[str, object] = {
        "sample": sample,
        "amplicon": amplicon,
        "graph_path": str(graph_path),
        "cycles_path": str(cycle_path),
        "coral_reconstruction_status": coral_reconstruction_status(graph_path),
    }
    try:
        graph = parse_graph_file(graph_path)
        decomposition = parse_cycles_file(cycle_path)
        structure = calculate_structure_lwcn(decomposition)
        decomposition_check = check_decomposition_against_graph(
            graph, decomposition, structure
        )
        upper = solve_cyclic_lwcn_upper_bound(graph)
        if upper.status != "OPTIMAL":
            return {
                **base,
                "status": upper.status,
                "error": upper.solver_message,
                "analysis_seconds": perf_counter() - started,
            }
        assert upper.maximum_cyclic_lwcn is not None
        assert upper.maximum_cyclic_ratio is not None
        total = graph.total_lwcn
        actual_ratio = structure.cycle_lwcn / total if total > 0 else 0.0
        excess = structure.cycle_lwcn - upper.maximum_cyclic_lwcn
        excess_ratio = excess / total if total > 0 else 0.0
        comparison_eligible = base["coral_reconstruction_status"] == "SUCCESS"
        violation = (
            comparison_eligible and excess_ratio > 1e-8 and excess > 1e-3
        )
        return {
            **base,
            "status": "ACTUAL_EXCEEDS_UPPER_BOUND" if violation else "OK",
            "error": "" if not violation else "Extracted cycle LWCN exceeds LP upper bound",
            "comparison_eligible": comparison_eligible,
            "sequence_edges": len(graph.sequence_edges),
            "breakpoint_edges": len(graph.breakpoint_edges),
            "intervals": graph.interval_count,
            "cycles": structure.cycle_count,
            "paths": structure.path_count,
            "graph_lwcn": total,
            "actual_cycle_lwcn": structure.cycle_lwcn,
            "actual_path_lwcn": structure.path_lwcn,
            "actual_cycle_ratio": actual_ratio,
            "upper_cycle_lwcn": upper.maximum_cyclic_lwcn,
            "upper_cycle_ratio": upper.maximum_cyclic_ratio,
            "headroom_lwcn": upper.maximum_cyclic_lwcn - structure.cycle_lwcn,
            "headroom_ratio": upper.maximum_cyclic_ratio - actual_ratio,
            "actual_minus_upper_lwcn": excess,
            "actual_minus_upper_ratio": excess_ratio,
            "decomposition_lwcn_residual": decomposition_check.lwcn_residual,
            "max_segment_copy_residual": decomposition_check.max_segment_copy_residual,
            "unmatched_graph_segments": len(decomposition_check.unmatched_graph_segments),
            "unmatched_cycle_segments": len(decomposition_check.unmatched_cycle_segments),
            "input_max_internal_balance_residual": upper.max_input_internal_balance_residual,
            "lp_minimum_walk_lwcn": upper.minimum_walk_lwcn,
            "lp_primal_dual_gap": upper.primal_dual_gap,
            "lp_max_state_residual": upper.max_state_residual,
            "lp_max_terminal_residual": upper.max_terminal_residual,
            "lp_max_capacity_violation": upper.max_capacity_violation,
            "lp_max_remaining_balance_residual": upper.max_remaining_balance_residual,
            "lp_arcs": upper.arc_count,
            "lp_solve_seconds": upper.solve_seconds,
            "analysis_seconds": perf_counter() - started,
        }
    except (CoRALParseError, UnbalancedGraphError, OSError, ValueError) as exc:
        return {
            **base,
            "status": type(exc).__name__,
            "error": str(exc),
            "analysis_seconds": perf_counter() - started,
        }


def _finite_series(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame[column], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()


def summarize_results(frame: pd.DataFrame) -> dict[str, object]:
    ok = frame[frame["status"] == "OK"].copy()
    comparable = ok[ok["comparison_eligible"] == True].copy()  # noqa: E712
    total_lwcn = float(comparable["graph_lwcn"].sum()) if not comparable.empty else 0.0
    actual_total = float(comparable["actual_cycle_lwcn"].sum()) if not comparable.empty else 0.0
    upper_total = float(comparable["upper_cycle_lwcn"].sum()) if not comparable.empty else 0.0
    headroom = _finite_series(comparable, "headroom_ratio") if not comparable.empty else pd.Series(dtype=float)
    actual = _finite_series(comparable, "actual_cycle_ratio") if not comparable.empty else pd.Series(dtype=float)
    upper = _finite_series(comparable, "upper_cycle_ratio") if not comparable.empty else pd.Series(dtype=float)
    correlation = float(actual.corr(upper)) if len(actual) >= 2 else None
    return {
        "pair_count": int(len(frame)),
        "lp_solved_ok_count": int(len(ok)),
        "algorithm_error_count": int((frame["status"] != "OK").sum()),
        "coral_reconstruction_success_count": int(
            (frame["coral_reconstruction_status"] == "SUCCESS").sum()
        ),
        "coral_reconstruction_failure_count": int(
            (frame["coral_reconstruction_status"] == "FAILURE").sum()
        ),
        "comparison_count": int(len(comparable)),
        "excluded_from_comparison_count": int(len(ok) - len(comparable)),
        "actual_exceeds_upper_bound_count": int(
            (frame["status"] == "ACTUAL_EXCEEDS_UPPER_BOUND").sum()
        ),
        "sample_count": int(comparable["sample"].nunique()) if not comparable.empty else 0,
        "all_lp_graph_lwcn": float(ok["graph_lwcn"].sum()) if not ok.empty else 0.0,
        "comparison_graph_lwcn": total_lwcn,
        "comparison_actual_cycle_lwcn": actual_total,
        "comparison_upper_cycle_lwcn": upper_total,
        "weighted_actual_cycle_ratio": actual_total / total_lwcn if total_lwcn else None,
        "weighted_upper_cycle_ratio": upper_total / total_lwcn if total_lwcn else None,
        "actual_as_fraction_of_upper_total": actual_total / upper_total if upper_total else None,
        "mean_headroom_ratio": float(headroom.mean()) if len(headroom) else None,
        "median_headroom_ratio": float(headroom.median()) if len(headroom) else None,
        "p95_headroom_ratio": float(headroom.quantile(0.95)) if len(headroom) else None,
        "actual_upper_ratio_pearson": correlation,
        "tight_within_1e_6_count": int((headroom.abs() <= 1e-6).sum()),
        "max_actual_minus_upper_ratio": float(
            _finite_series(comparable, "actual_minus_upper_ratio").max()
        ) if "actual_minus_upper_ratio" in frame and _finite_series(frame, "actual_minus_upper_ratio").size else None,
        "max_input_internal_balance_residual": float(
            _finite_series(ok, "input_max_internal_balance_residual").max()
        ) if not ok.empty else None,
        "max_decomposition_lwcn_abs_residual": float(
            _finite_series(comparable, "decomposition_lwcn_residual").abs().max()
        ) if not comparable.empty else None,
        "weighted_decomposition_completeness": (
            1.0
            - float(comparable["decomposition_lwcn_residual"].sum()) / total_lwcn
            if total_lwcn
            else None
        ),
        "median_decomposition_completeness": float(
            (1.0 - comparable["decomposition_lwcn_residual"] / comparable["graph_lwcn"]).median()
        ) if not comparable.empty else None,
        "decomposition_incomplete_over_1e_6_count": int(
            (
                comparable["decomposition_lwcn_residual"].abs()
                / comparable["graph_lwcn"]
                > 1e-6
            ).sum()
        ) if not comparable.empty else 0,
        "max_segment_copy_residual": float(
            _finite_series(comparable, "max_segment_copy_residual").max()
        ) if not comparable.empty else None,
        "max_lp_primal_dual_gap": float(
            _finite_series(ok, "lp_primal_dual_gap").max()
        ) if not ok.empty else None,
        "max_lp_state_residual": float(
            _finite_series(ok, "lp_max_state_residual").max()
        ) if not ok.empty else None,
        "max_lp_terminal_residual": float(
            _finite_series(ok, "lp_max_terminal_residual").max()
        ) if not ok.empty else None,
        "max_lp_capacity_violation": float(
            _finite_series(ok, "lp_max_capacity_violation").max()
        ) if not ok.empty else None,
        "max_lp_remaining_balance_residual": float(
            _finite_series(ok, "lp_max_remaining_balance_residual").max()
        ) if not ok.empty else None,
        "total_lp_solve_seconds": float(ok["lp_solve_seconds"].sum()) if not ok.empty else 0.0,
        "median_lp_solve_seconds": float(ok["lp_solve_seconds"].median()) if not ok.empty else None,
        "max_lp_solve_seconds": float(ok["lp_solve_seconds"].max()) if not ok.empty else None,
        "status_counts": {str(k): int(v) for k, v in frame["status"].value_counts().items()},
    }


def _write_figure(frame: pd.DataFrame, output_path: Path) -> None:
    ok = frame[
        (frame["status"] == "OK") & (frame["comparison_eligible"] == True)  # noqa: E712
    ].copy()
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})
    figure, axes = plt.subplots(1, 3, figsize=(12, 3.6), constrained_layout=True)
    axes[0].scatter(
        ok["upper_cycle_ratio"], ok["actual_cycle_ratio"], s=18, alpha=0.7
    )
    axes[0].plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=1)
    axes[0].set(xlabel="LP upper ratio", ylabel="Extracted cycle ratio", title="Bound check")
    axes[0].set_xlim(-0.02, 1.02)
    axes[0].set_ylim(-0.02, 1.02)
    axes[1].hist(ok["headroom_ratio"], bins=20, color="#2f6f9f", edgecolor="white")
    axes[1].set(xlabel="Upper ratio - extracted ratio", ylabel="Amplicons", title="Upper-bound headroom")
    axes[2].scatter(ok["sequence_edges"], ok["lp_solve_seconds"], s=18, alpha=0.7, color="#b45f06")
    axes[2].set(xlabel="Sequence edges", ylabel="LP solve seconds", title="Runtime scaling")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def run_batch_analysis(data_root: str | Path, output_dir: str | Path) -> dict[str, object]:
    pairs = discover_official_pairs(data_root)
    if not pairs:
        raise FileNotFoundError(
            f"No original CoRAL graph/cycle pairs found below {Path(data_root).resolve()}"
        )
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows = [analyze_pair(graph, cycles) for graph, cycles in pairs]
    frame = pd.DataFrame(rows)
    frame.to_csv(output / "per_amplicon_results.csv", index=False, encoding="utf-8-sig")
    issues = frame[
        (frame["status"] != "OK")
        | (frame["coral_reconstruction_status"] != "SUCCESS")
    ].copy()
    if not issues.empty:
        issues.insert(
            0,
            "issue_type",
            np.where(
                issues["status"] != "OK",
                issues["status"],
                "CORAL_RECONSTRUCTION_" + issues["coral_reconstruction_status"].astype(str),
            ),
        )
    issues.to_csv(output / "issues.csv", index=False, encoding="utf-8-sig")
    numeric = [
        "graph_lwcn",
        "actual_cycle_lwcn",
        "upper_cycle_lwcn",
        "lp_solve_seconds",
    ]
    comparable = frame[
        (frame["status"] == "OK") & (frame["comparison_eligible"] == True)  # noqa: E712
    ]
    sample_summary = comparable.groupby("sample", as_index=False)[numeric].sum()
    if not sample_summary.empty:
        sample_summary["actual_cycle_ratio"] = (
            sample_summary["actual_cycle_lwcn"] / sample_summary["graph_lwcn"]
        )
        sample_summary["upper_cycle_ratio"] = (
            sample_summary["upper_cycle_lwcn"] / sample_summary["graph_lwcn"]
        )
    sample_summary.to_csv(output / "per_sample_summary.csv", index=False, encoding="utf-8-sig")
    summary = summarize_results(frame)
    (output / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    _write_figure(frame, output / "performance_overview.png")
    return summary
