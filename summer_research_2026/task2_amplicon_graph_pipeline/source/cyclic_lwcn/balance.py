"""Validate the exact CN balance required by the state-graph LP."""

from __future__ import annotations

from dataclasses import dataclass

from .model import BreakpointGraph


@dataclass(frozen=True)
class BalanceReport:
    terminal_capacities: dict[str, float]
    node_residuals: dict[str, float]
    max_internal_abs_residual: float
    most_negative_endnode_residual: float
    tolerance: float


class UnbalancedGraphError(ValueError):
    """Raised when input capacities do not admit an exact base-model decomposition."""


def _node_tolerance(sequence_cn: float, absolute: float, relative: float) -> float:
    return absolute + relative * max(1.0, abs(sequence_cn))


def validate_and_infer_terminals(
    graph: BreakpointGraph,
    *,
    absolute_tolerance: float = 1e-7,
    relative_tolerance: float = 1e-10,
) -> BalanceReport:
    """Validate local balance and infer terminal half-edge CN at interval ends.

    At an internal endpoint, sequence CN must equal incident breakpoint CN. At
    a CoRAL interval end, a nonnegative residual is the terminal capacity used
    by open paths. Values within tolerance are clipped only for arithmetic; no
    optimization-based data repair or fabricated capacity is performed.
    """

    breakpoint_incidence = {node: 0.0 for node in graph.nodes}
    for edge in graph.breakpoint_edges:
        breakpoint_incidence[edge.node1] += edge.copy_number
        breakpoint_incidence[edge.node2] += edge.copy_number

    residuals: dict[str, float] = {}
    terminals: dict[str, float] = {}
    invalid_internal: list[tuple[str, float, float]] = []
    invalid_endnodes: list[tuple[str, float, float]] = []
    max_internal = 0.0
    most_negative_endnode = 0.0
    for node, sequence in graph.sequence_by_node.items():
        residual = sequence.copy_number - breakpoint_incidence[node]
        residuals[node] = residual
        tolerance = _node_tolerance(
            sequence.copy_number, absolute_tolerance, relative_tolerance
        )
        if node in graph.endnodes:
            most_negative_endnode = min(most_negative_endnode, residual)
            if residual < -tolerance:
                invalid_endnodes.append((node, residual, tolerance))
            terminals[node] = max(0.0, residual) if abs(residual) > tolerance else 0.0
        else:
            max_internal = max(max_internal, abs(residual))
            if abs(residual) > tolerance:
                invalid_internal.append((node, residual, tolerance))

    if invalid_internal or invalid_endnodes:
        details: list[str] = []
        if invalid_internal:
            node, residual, tolerance = max(
                invalid_internal, key=lambda item: abs(item[1])
            )
            details.append(
                f"internal node {node}: residual={residual:.12g}, tolerance={tolerance:.3g}"
            )
        if invalid_endnodes:
            node, residual, tolerance = min(
                invalid_endnodes, key=lambda item: item[1]
            )
            details.append(
                f"endnode {node}: negative residual={residual:.12g}, tolerance={tolerance:.3g}"
            )
        raise UnbalancedGraphError("; ".join(details))

    return BalanceReport(
        terminal_capacities=terminals,
        node_residuals=residuals,
        max_internal_abs_residual=max_internal,
        most_negative_endnode_residual=most_negative_endnode,
        tolerance=absolute_tolerance,
    )

