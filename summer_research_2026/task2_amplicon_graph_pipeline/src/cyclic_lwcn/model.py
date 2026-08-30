"""Typed in-memory representation of CoRAL graphs and cycle files."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class GenomicNode:
    raw: str
    chromosome: str
    position: int
    strand: str


@dataclass(frozen=True)
class SequenceEdge:
    edge_id: str
    start_node: str
    end_node: str
    chromosome: str
    start: int
    end: int
    copy_number: float
    length: int

    def other(self, node: str) -> str:
        if node == self.start_node:
            return self.end_node
        if node == self.end_node:
            return self.start_node
        raise KeyError(f"Node {node!r} is not incident to {self.edge_id}")


@dataclass(frozen=True)
class BreakpointEdge:
    edge_id: str
    edge_type: str
    node1: str
    node2: str
    copy_number: float


@dataclass
class BreakpointGraph:
    source_path: Path
    nodes: dict[str, GenomicNode]
    sequence_edges: tuple[SequenceEdge, ...]
    breakpoint_edges: tuple[BreakpointEdge, ...]
    sequence_by_node: dict[str, SequenceEdge]
    endnodes: frozenset[str]
    interval_count: int
    warnings: list[str] = field(default_factory=list)
    @property
    def total_lwcn(self) -> float:
        return sum(edge.copy_number * edge.length for edge in self.sequence_edges)


@dataclass(frozen=True)
class Segment:
    segment_id: int
    chromosome: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class SegmentTraversal:
    segment_id: int
    strand: str


@dataclass(frozen=True)
class Structure:
    structure_id: str
    kind: str
    copy_number: float
    traversals: tuple[SegmentTraversal, ...]
    raw_line: str


@dataclass
class CycleDecomposition:
    source_path: Path
    segments: dict[int, Segment]
    structures: tuple[Structure, ...]
    warnings: list[str] = field(default_factory=list)
