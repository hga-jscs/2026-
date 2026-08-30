"""Strict parsers for CoRAL ``*_graph.txt`` and ``*_cycles.txt`` files."""

from __future__ import annotations

import re
from pathlib import Path

from .model import (
    BreakpointEdge,
    BreakpointGraph,
    CycleDecomposition,
    GenomicNode,
    Segment,
    SegmentTraversal,
    SequenceEdge,
    Structure,
)


_NODE_RE = re.compile(r"^(?P<chrom>.+):(?P<position>\d+)(?P<strand>[+-])$")
_TRAVERSAL_RE = re.compile(r"^(?P<segment>\d+)(?P<strand>[+-])$")


class CoRALParseError(ValueError):
    """Raised when a CoRAL input is malformed or internally inconsistent."""


def parse_node(raw: str) -> GenomicNode:
    match = _NODE_RE.fullmatch(raw.strip())
    if not match:
        raise CoRALParseError(f"Invalid genomic endpoint: {raw!r}")
    return GenomicNode(
        raw=raw.strip(),
        chromosome=match.group("chrom"),
        position=int(match.group("position")),
        strand=match.group("strand"),
    )


def _as_nonnegative_float(raw: str, label: str, line_number: int) -> float:
    try:
        value = float(raw)
    except ValueError as exc:
        raise CoRALParseError(
            f"Line {line_number}: {label} is not numeric: {raw!r}"
        ) from exc
    if value < 0:
        raise CoRALParseError(f"Line {line_number}: {label} must be nonnegative")
    return value


def _infer_endnodes(
    sequence_edges: list[SequenceEdge],
) -> tuple[frozenset[str], int]:
    if not sequence_edges:
        return frozenset(), 0
    ordered = sorted(sequence_edges, key=lambda e: (e.chromosome, e.start, e.end))
    endnodes: set[str] = set()
    interval_count = 0
    first = ordered[0]
    interval_start = first.start_node
    interval_end = first.end_node
    previous_chromosome = first.chromosome
    previous_end = first.end
    for edge in ordered[1:]:
        contiguous = (
            edge.chromosome == previous_chromosome
            and edge.start == previous_end + 1
        )
        if not contiguous:
            endnodes.update((interval_start, interval_end))
            interval_count += 1
            interval_start = edge.start_node
        interval_end = edge.end_node
        previous_chromosome = edge.chromosome
        previous_end = edge.end
    endnodes.update((interval_start, interval_end))
    return frozenset(endnodes), interval_count + 1


def parse_graph_file(path: str | Path) -> BreakpointGraph:
    """Parse a CoRAL breakpoint graph and infer interval boundary endnodes.

    The inference follows CoRAL's own rule: adjacent sequence edges belong to the
    same interval only when they share a chromosome and the next start equals
    the previous end plus one.
    """

    source_path = Path(path).resolve()
    sequence_edges: list[SequenceEdge] = []
    breakpoint_edges: list[BreakpointEdge] = []
    nodes: dict[str, GenomicNode] = {}
    sequence_by_node: dict[str, SequenceEdge] = {}
    type_counts = {"concordant": 0, "discordant": 0}
    warnings: list[str] = []

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
                copy_number = _as_nonnegative_float(fields[3], "sequence CN", line_number)
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
                if declared_length != coordinate_length:
                    raise CoRALParseError(
                        f"Line {line_number}: declared length {declared_length} does not match "
                        f"inclusive coordinate length {coordinate_length}"
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
                    raise CoRALParseError(
                        f"Line {line_number}: malformed {record_type} record"
                    )
                left, right = fields[1].split("->", maxsplit=1)
                node1, node2 = parse_node(left), parse_node(right)
                copy_number = _as_nonnegative_float(fields[2], "breakpoint CN", line_number)
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
        raise CoRALParseError(f"No sequence edges found in {source_path}")
    unknown_nodes = sorted(set(nodes) - set(sequence_by_node))
    if unknown_nodes:
        preview = ", ".join(unknown_nodes[:5])
        raise CoRALParseError(
            f"Breakpoint endpoint(s) do not belong to a sequence edge: {preview}"
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
        warnings=warnings,
    )


def _parse_semicolon_fields(line: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in line.split(";"):
        if "=" in field:
            key, value = field.split("=", maxsplit=1)
            result[key.strip()] = value.strip()
    return result


def parse_cycles_file(path: str | Path) -> CycleDecomposition:
    """Parse CoRAL cycle/path structures, retaining repeated segment traversals."""

    source_path = Path(path).resolve()
    segments: dict[int, Segment] = {}
    structures: list[Structure] = []
    warnings: list[str] = []
    with source_path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if fields[0] == "Segment":
                if len(fields) < 5:
                    raise CoRALParseError(
                        f"Line {line_number}: malformed segment definition"
                    )
                try:
                    segment_id = int(fields[1])
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as exc:
                    raise CoRALParseError(
                        f"Line {line_number}: invalid segment number or coordinate"
                    ) from exc
                if segment_id == 0 or end < start:
                    raise CoRALParseError(
                        f"Line {line_number}: invalid segment {segment_id}"
                    )
                segment = Segment(segment_id, fields[2], start, end)
                if segment_id in segments and segments[segment_id] != segment:
                    raise CoRALParseError(
                        f"Line {line_number}: conflicting definition of segment {segment_id}"
                    )
                segments[segment_id] = segment
                continue
            if not (line.startswith("Cycle=") or line.startswith("Path=")):
                continue
            parsed = _parse_semicolon_fields(line)
            kind = "cycle" if line.startswith("Cycle=") else "path"
            id_key = "Cycle" if kind == "cycle" else "Path"
            required = (id_key, "Copy_count", "Segments")
            missing = [key for key in required if key not in parsed]
            if missing:
                raise CoRALParseError(
                    f"Line {line_number}: missing structure field(s): {', '.join(missing)}"
                )
            copy_number = _as_nonnegative_float(
                parsed["Copy_count"], "structure copy count", line_number
            )
            traversals: list[SegmentTraversal] = []
            for token in filter(None, parsed["Segments"].split(",")):
                match = _TRAVERSAL_RE.fullmatch(token.strip())
                if not match:
                    raise CoRALParseError(
                        f"Line {line_number}: invalid segment traversal {token!r}"
                    )
                segment_id = int(match.group("segment"))
                if segment_id == 0:
                    if kind != "path":
                        raise CoRALParseError(
                            f"Line {line_number}: segment 0 is only valid in paths"
                        )
                    continue
                traversals.append(
                    SegmentTraversal(segment_id, match.group("strand"))
                )
            structures.append(
                Structure(
                    structure_id=parsed[id_key],
                    kind=kind,
                    copy_number=copy_number,
                    traversals=tuple(traversals),
                    raw_line=line,
                )
            )

    referenced = {t.segment_id for s in structures for t in s.traversals}
    missing_segments = sorted(referenced - set(segments))
    if missing_segments:
        raise CoRALParseError(
            f"Undefined segment ID(s) used by structures: {missing_segments[:10]}"
        )
    if not structures:
        warnings.append("No extracted cycles or paths were present")
    return CycleDecomposition(source_path, segments, tuple(structures), warnings)
