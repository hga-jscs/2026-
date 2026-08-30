"""Length-weighted copy-number calculations for extracted structures."""

from __future__ import annotations

from dataclasses import dataclass

from .model import BreakpointGraph, CycleDecomposition


@dataclass(frozen=True)
class StructureLWCNResult:
    cycle_lwcn: float
    path_lwcn: float
    structure_lwcn: float
    cycle_count: int
    path_count: int
    segment_total_copy: dict[int, float]
    segment_cycle_copy: dict[int, float]


@dataclass(frozen=True)
class DecompositionCheck:
    graph_lwcn: float
    decomposition_lwcn: float
    lwcn_residual: float
    max_segment_copy_residual: float
    unmatched_graph_segments: tuple[str, ...]
    unmatched_cycle_segments: tuple[int, ...]


def calculate_structure_lwcn(
    decomposition: CycleDecomposition,
) -> StructureLWCNResult:
    """Calculate cycle/path LWCN, counting every repeated segment traversal."""

    cycle_lwcn = 0.0
    path_lwcn = 0.0
    cycle_count = 0
    path_count = 0
    total_copy = {segment_id: 0.0 for segment_id in decomposition.segments}
    cycle_copy = {segment_id: 0.0 for segment_id in decomposition.segments}
    for structure in decomposition.structures:
        length = sum(
            decomposition.segments[traversal.segment_id].length
            for traversal in structure.traversals
        )
        lwcn = structure.copy_number * length
        if structure.kind == "cycle":
            cycle_lwcn += lwcn
            cycle_count += 1
        else:
            path_lwcn += lwcn
            path_count += 1
        for traversal in structure.traversals:
            total_copy[traversal.segment_id] += structure.copy_number
            if structure.kind == "cycle":
                cycle_copy[traversal.segment_id] += structure.copy_number
    return StructureLWCNResult(
        cycle_lwcn=cycle_lwcn,
        path_lwcn=path_lwcn,
        structure_lwcn=cycle_lwcn + path_lwcn,
        cycle_count=cycle_count,
        path_count=path_count,
        segment_total_copy=total_copy,
        segment_cycle_copy=cycle_copy,
    )


def check_decomposition_against_graph(
    graph: BreakpointGraph,
    decomposition: CycleDecomposition,
    metrics: StructureLWCNResult | None = None,
) -> DecompositionCheck:
    """Compare aggregate structure segment CN with graph sequence-edge CN."""

    metrics = metrics or calculate_structure_lwcn(decomposition)
    graph_by_coordinate = {
        (e.chromosome, e.start, e.end): e for e in graph.sequence_edges
    }
    cycles_by_coordinate = {
        (s.chromosome, s.start, s.end): s for s in decomposition.segments.values()
    }
    unmatched_graph = sorted(
        edge.edge_id
        for coordinate, edge in graph_by_coordinate.items()
        if coordinate not in cycles_by_coordinate
    )
    unmatched_cycles = sorted(
        segment.segment_id
        for coordinate, segment in cycles_by_coordinate.items()
        if coordinate not in graph_by_coordinate
    )
    max_residual = 0.0
    for coordinate in graph_by_coordinate.keys() & cycles_by_coordinate.keys():
        edge = graph_by_coordinate[coordinate]
        segment = cycles_by_coordinate[coordinate]
        max_residual = max(
            max_residual,
            abs(edge.copy_number - metrics.segment_total_copy[segment.segment_id]),
        )
    return DecompositionCheck(
        graph_lwcn=graph.total_lwcn,
        decomposition_lwcn=metrics.structure_lwcn,
        lwcn_residual=graph.total_lwcn - metrics.structure_lwcn,
        max_segment_copy_residual=max_residual,
        unmatched_graph_segments=tuple(unmatched_graph),
        unmatched_cycle_segments=tuple(unmatched_cycles),
    )

