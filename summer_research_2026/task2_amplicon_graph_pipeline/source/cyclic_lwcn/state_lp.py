"""Exact state-graph LP for the base-model maximum cyclic LWCN."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix

from .balance import BalanceReport, validate_and_infer_terminals
from .model import BreakpointGraph


@dataclass(frozen=True)
class StateArc:
    arc_id: str
    tail: str | None
    head: str | None
    sequence_edge_id: str | None
    breakpoint_edge_id: str | None
    terminal_node: str | None
    cost: float


@dataclass
class UpperBoundResult:
    status: str
    solver_status: int
    solver_message: str
    total_lwcn: float
    minimum_walk_lwcn: float | None
    maximum_cyclic_lwcn: float | None
    maximum_cyclic_ratio: float | None
    dual_walk_bound: float | None
    cyclic_lower_certificate: float | None
    cyclic_upper_certificate: float | None
    primal_dual_gap: float | None
    solve_seconds: float
    state_count: int
    arc_count: int
    equality_count: int
    capacity_count: int
    max_state_residual: float | None
    max_terminal_residual: float | None
    max_capacity_violation: float | None
    min_remaining_capacity: float | None
    max_remaining_balance_residual: float | None
    max_input_internal_balance_residual: float
    terminal_capacity_total: float
    balance_report: BalanceReport
    arcs: tuple[StateArc, ...]
    flow: np.ndarray | None = None

    def as_dict(self, include_flow: bool = False) -> dict[str, object]:
        record: dict[str, object] = {
            "status": self.status,
            "solver_status": self.solver_status,
            "solver_message": self.solver_message,
            "total_lwcn": self.total_lwcn,
            "minimum_walk_lwcn": self.minimum_walk_lwcn,
            "maximum_cyclic_lwcn": self.maximum_cyclic_lwcn,
            "maximum_cyclic_ratio": self.maximum_cyclic_ratio,
            "dual_walk_bound": self.dual_walk_bound,
            "cyclic_lower_certificate": self.cyclic_lower_certificate,
            "cyclic_upper_certificate": self.cyclic_upper_certificate,
            "primal_dual_gap": self.primal_dual_gap,
            "solve_seconds": self.solve_seconds,
            "state_count": self.state_count,
            "arc_count": self.arc_count,
            "equality_count": self.equality_count,
            "capacity_count": self.capacity_count,
            "max_state_residual": self.max_state_residual,
            "max_terminal_residual": self.max_terminal_residual,
            "max_capacity_violation": self.max_capacity_violation,
            "min_remaining_capacity": self.min_remaining_capacity,
            "max_remaining_balance_residual": self.max_remaining_balance_residual,
            "max_input_internal_balance_residual": self.max_input_internal_balance_residual,
            "terminal_capacity_total": self.terminal_capacity_total,
        }
        if include_flow and self.flow is not None:
            record["positive_arcs"] = [
                {
                    "arc_id": arc.arc_id,
                    "tail": arc.tail or "S",
                    "head": arc.head or "T",
                    "sequence_edge_id": arc.sequence_edge_id,
                    "breakpoint_edge_id": arc.breakpoint_edge_id,
                    "terminal_node": arc.terminal_node,
                    "cost": arc.cost,
                    "flow": float(value),
                }
                for arc, value in zip(self.arcs, self.flow, strict=True)
                if value > 1e-10
            ]
        return record


def _build_arcs(
    graph: BreakpointGraph, balance: BalanceReport
) -> tuple[StateArc, ...]:
    arcs: list[StateArc] = []
    for edge in graph.breakpoint_edges:
        orientations = [(edge.node1, edge.node2)]
        if edge.node1 != edge.node2:
            orientations.append((edge.node2, edge.node1))
        for direction, (tail, breakpoint_arrival) in enumerate(orientations, start=1):
            sequence = graph.sequence_by_node[breakpoint_arrival]
            arcs.append(
                StateArc(
                    arc_id=f"{edge.edge_id}:dir{direction}",
                    tail=tail,
                    head=sequence.other(breakpoint_arrival),
                    sequence_edge_id=sequence.edge_id,
                    breakpoint_edge_id=edge.edge_id,
                    terminal_node=None,
                    cost=float(sequence.length),
                )
            )
    for node, capacity in sorted(balance.terminal_capacities.items()):
        if capacity <= 0:
            continue
        sequence = graph.sequence_by_node[node]
        arcs.append(
            StateArc(
                arc_id=f"terminal:{node}:source",
                tail=None,
                head=sequence.other(node),
                sequence_edge_id=sequence.edge_id,
                breakpoint_edge_id=None,
                terminal_node=node,
                cost=float(sequence.length),
            )
        )
        arcs.append(
            StateArc(
                arc_id=f"terminal:{node}:sink",
                tail=node,
                head=None,
                sequence_edge_id=None,
                breakpoint_edge_id=None,
                terminal_node=node,
                cost=0.0,
            )
        )
    return tuple(arcs)


def _sparse_matrix(
    row_count: int, column_count: int, entries: list[tuple[int, int, float]]
) -> csr_matrix:
    if not entries:
        return csr_matrix((row_count, column_count), dtype=float)
    rows, columns, values = zip(*entries, strict=True)
    return coo_matrix(
        (values, (rows, columns)), shape=(row_count, column_count), dtype=float
    ).tocsr()


def solve_cyclic_lwcn_upper_bound(
    graph: BreakpointGraph,
    *,
    absolute_balance_tolerance: float = 1e-7,
    relative_balance_tolerance: float = 1e-10,
) -> UpperBoundResult:
    """Solve the minimum-walk LP and return ``D - W_min``.

    This is exact for the aggregate base model: all CN is explained, repeated
    nodes/edges are permitted, and no per-structure long-read or traversal-count
    constraints are imposed. It is a safe upper bound when those full-CoRAL
    structure-level constraints are added.
    """

    balance = validate_and_infer_terminals(
        graph,
        absolute_tolerance=absolute_balance_tolerance,
        relative_tolerance=relative_balance_tolerance,
    )
    arcs = _build_arcs(graph, balance)
    nodes = tuple(sorted(graph.nodes))
    terminals = tuple(
        node for node, capacity in sorted(balance.terminal_capacities.items()) if capacity > 0
    )
    sequence_edges = tuple(graph.sequence_edges)
    breakpoint_edges = tuple(graph.breakpoint_edges)
    node_row = {node: index for index, node in enumerate(nodes)}
    terminal_row = {
        node: len(nodes) + index for index, node in enumerate(terminals)
    }
    seq_row = {edge.edge_id: index for index, edge in enumerate(sequence_edges)}
    bp_row = {
        edge.edge_id: len(sequence_edges) + index
        for index, edge in enumerate(breakpoint_edges)
    }

    equality_entries: list[tuple[int, int, float]] = []
    capacity_entries: list[tuple[int, int, float]] = []
    for column, arc in enumerate(arcs):
        if arc.tail is not None:
            equality_entries.append((node_row[arc.tail], column, -1.0))
        if arc.head is not None:
            equality_entries.append((node_row[arc.head], column, 1.0))
        if arc.terminal_node is not None:
            equality_entries.append((terminal_row[arc.terminal_node], column, 1.0))
        if arc.sequence_edge_id is not None:
            capacity_entries.append((seq_row[arc.sequence_edge_id], column, 1.0))
        if arc.breakpoint_edge_id is not None:
            capacity_entries.append((bp_row[arc.breakpoint_edge_id], column, 1.0))

    equality_count = len(nodes) + len(terminals)
    capacity_count = len(sequence_edges) + len(breakpoint_edges)
    a_eq = _sparse_matrix(equality_count, len(arcs), equality_entries)
    a_ub = _sparse_matrix(capacity_count, len(arcs), capacity_entries)
    b_eq = np.zeros(equality_count, dtype=float)
    for node, row in terminal_row.items():
        b_eq[row] = balance.terminal_capacities[node]
    b_ub = np.array(
        [edge.copy_number for edge in sequence_edges]
        + [edge.copy_number for edge in breakpoint_edges],
        dtype=float,
    )
    objective = np.array([arc.cost for arc in arcs], dtype=float)

    start = perf_counter()
    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=(0, None),
        method="highs",
    )
    elapsed = perf_counter() - start
    common = dict(
        total_lwcn=graph.total_lwcn,
        solve_seconds=elapsed,
        state_count=len(nodes),
        arc_count=len(arcs),
        equality_count=equality_count,
        capacity_count=capacity_count,
        max_input_internal_balance_residual=balance.max_internal_abs_residual,
        terminal_capacity_total=sum(balance.terminal_capacities.values()),
        balance_report=balance,
        arcs=arcs,
    )
    if not result.success:
        return UpperBoundResult(
            status="NO_EXACT_OPEN_WALK_DECOMPOSITION",
            solver_status=int(result.status),
            solver_message=str(result.message),
            minimum_walk_lwcn=None,
            maximum_cyclic_lwcn=None,
            maximum_cyclic_ratio=None,
            dual_walk_bound=None,
            cyclic_lower_certificate=None,
            cyclic_upper_certificate=None,
            primal_dual_gap=None,
            max_state_residual=None,
            max_terminal_residual=None,
            max_capacity_violation=None,
            min_remaining_capacity=None,
            max_remaining_balance_residual=None,
            flow=None,
            **common,
        )

    flow = np.asarray(result.x, dtype=float)
    equality_activity = np.asarray(a_eq @ flow).ravel()
    capacity_activity = np.asarray(a_ub @ flow).ravel()
    state_residual = equality_activity[: len(nodes)]
    terminal_residual = equality_activity[len(nodes) :] - b_eq[len(nodes) :]
    remaining = b_ub - capacity_activity
    max_capacity_violation = max(0.0, float(np.max(-remaining, initial=0.0)))

    remaining_sequence = {
        edge.edge_id: remaining[index] for index, edge in enumerate(sequence_edges)
    }
    remaining_breakpoint = {
        edge.edge_id: remaining[len(sequence_edges) + index]
        for index, edge in enumerate(breakpoint_edges)
    }
    remaining_bp_incidence = {node: 0.0 for node in nodes}
    for edge in breakpoint_edges:
        value = remaining_breakpoint[edge.edge_id]
        remaining_bp_incidence[edge.node1] += value
        remaining_bp_incidence[edge.node2] += value
    remaining_balance_residuals = []
    for node in nodes:
        seq = graph.sequence_by_node[node]
        remaining_balance_residuals.append(
            remaining_sequence[seq.edge_id] - remaining_bp_incidence[node]
        )

    primal = float(result.fun)
    dual = None
    if hasattr(result.eqlin, "marginals") and hasattr(result.ineqlin, "marginals"):
        dual = float(
            b_eq @ np.asarray(result.eqlin.marginals)
            + b_ub @ np.asarray(result.ineqlin.marginals)
        )
    gap = None if dual is None else max(0.0, primal - dual)
    total = graph.total_lwcn
    maximum = total - primal
    scale_tolerance = 1e-8 * max(1.0, total)
    if -scale_tolerance <= maximum < 0:
        maximum = 0.0
    if total <= 0:
        ratio = 0.0
    else:
        ratio = maximum / total
    return UpperBoundResult(
        status="OPTIMAL",
        solver_status=int(result.status),
        solver_message=str(result.message),
        minimum_walk_lwcn=primal,
        maximum_cyclic_lwcn=maximum,
        maximum_cyclic_ratio=ratio,
        dual_walk_bound=dual,
        cyclic_lower_certificate=maximum,
        cyclic_upper_certificate=None if dual is None else total - dual,
        primal_dual_gap=gap,
        max_state_residual=float(np.max(np.abs(state_residual), initial=0.0)),
        max_terminal_residual=float(np.max(np.abs(terminal_residual), initial=0.0)),
        max_capacity_violation=max_capacity_violation,
        min_remaining_capacity=float(np.min(remaining, initial=0.0)),
        max_remaining_balance_residual=float(
            np.max(np.abs(remaining_balance_residuals), initial=0.0)
        ),
        flow=flow,
        **common,
    )

