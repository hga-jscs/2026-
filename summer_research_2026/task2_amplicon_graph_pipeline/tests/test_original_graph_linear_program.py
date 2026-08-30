"""原图 LP 的合成算例与 112 张真实图回归测试。"""

import os
from pathlib import Path

import pytest

from cyclic_lwcn.analysis import discover_official_pairs
from cyclic_lwcn.parser import parse_graph_file
from cyclic_lwcn.state_lp import solve_cyclic_lwcn_upper_bound
from original_graph_lwcn import solve_original_graph_linear_program


def _write_graph(
    tmp_path: Path,
    sequences: list[tuple[int, int, float]],
    breakpoints: list[tuple[str, str, str, float]],
) -> Path:
    """写出最小 CoRAL graph 测试文件；测试数据本身也遵守真实文本格式。"""

    graph = tmp_path / "synthetic_graph.txt"
    lines = [
        "SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberReadsMapped"
    ]
    for start, end, copy_number in sequences:
        length = end - start + 1
        lines.append(
            f"sequence\tchr1:{start}-\tchr1:{end}+\t{copy_number}\t0\t{length}\t0"
        )
    lines.append(
        "BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfReadPairs"
    )
    for edge_type, node1, node2, copy_number in breakpoints:
        lines.append(f"{edge_type}\t{node1}->{node2}\t{copy_number}\t0")
    graph.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return graph


def test_pure_linear_graph_has_no_circular_load(tmp_path: Path) -> None:
    graph = parse_graph_file(_write_graph(tmp_path, [(1, 100, 1)], []))
    result = solve_original_graph_linear_program(graph)
    assert result.status == "OPTIMAL"
    assert result.maximum_cyclic_lwcn == pytest.approx(0.0)


def test_single_segment_cycle_uses_the_whole_sequence(tmp_path: Path) -> None:
    path = _write_graph(
        tmp_path,
        [(1, 100, 1)],
        [("discordant", "chr1:1-", "chr1:100+", 1)],
    )
    result = solve_original_graph_linear_program(parse_graph_file(path))
    assert result.maximum_cyclic_lwcn == pytest.approx(100.0)
    assert result.maximum_cyclic_ratio == pytest.approx(1.0)


def test_mixed_path_and_cycle_matches_hand_solution(tmp_path: Path) -> None:
    path = _write_graph(
        tmp_path,
        [(1, 3, 2), (4, 8, 2)],
        [
            ("concordant", "chr1:3+", "chr1:4-", 2),
            ("discordant", "chr1:8+", "chr1:1-", 1),
        ],
    )
    result = solve_original_graph_linear_program(parse_graph_file(path))
    assert result.total_lwcn == pytest.approx(16.0)
    assert result.maximum_cyclic_lwcn == pytest.approx(8.0)
    assert result.maximum_cyclic_ratio == pytest.approx(0.5)
    assert result.max_balance_residual <= 1e-9


def test_foldback_self_loops_have_incidence_two(tmp_path: Path) -> None:
    path = _write_graph(
        tmp_path,
        [(1, 100, 2)],
        [
            ("discordant", "chr1:1-", "chr1:1-", 1),
            ("discordant", "chr1:100+", "chr1:100+", 1),
        ],
    )
    result = solve_original_graph_linear_program(parse_graph_file(path))
    assert result.maximum_cyclic_lwcn == pytest.approx(200.0)
    assert result.maximum_cyclic_ratio == pytest.approx(1.0)


def test_all_112_official_graphs_match_the_state_graph_lp() -> None:
    configured_root = os.environ.get("CORAL_DATA_ROOT")
    if not configured_root:
        pytest.skip("set CORAL_DATA_ROOT to run the 112-pair CoRAL integration test")
    data_root = Path(configured_root).expanduser().resolve()
    pairs = discover_official_pairs(data_root)
    assert len(pairs) == 112
    maximum_ratio_difference = 0.0
    for graph_path, _ in pairs:
        graph = parse_graph_file(graph_path)
        direct = solve_original_graph_linear_program(graph)
        state = solve_cyclic_lwcn_upper_bound(graph)
        assert direct.status == "OPTIMAL"
        assert state.status == "OPTIMAL"
        assert direct.maximum_cyclic_ratio is not None
        assert state.maximum_cyclic_ratio is not None
        maximum_ratio_difference = max(
            maximum_ratio_difference,
            abs(direct.maximum_cyclic_ratio - state.maximum_cyclic_ratio),
        )
    assert maximum_ratio_difference <= 1e-10
