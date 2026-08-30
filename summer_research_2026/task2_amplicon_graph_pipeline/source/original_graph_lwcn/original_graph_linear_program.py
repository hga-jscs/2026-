"""在原断点图物理边上直接求最大环状 LWCN 的稀疏 LP。"""

from __future__ import annotations

from time import perf_counter

import numpy as np
from scipy.optimize import linprog

from cyclic_lwcn.model import BreakpointGraph

from .graph_cycle_model import (
    build_original_graph_cycle_model,
    validate_circular_edge_load,
)
from .solver_result import CyclicLWCNSolverResult


def solve_original_graph_linear_program(
    graph: BreakpointGraph,
) -> CyclicLWCNSolverResult:
    """求解原图最大平衡环流。

    数学模型为：

    ``max w^T z``

    ``s.t. A z = 0, 0 <= z <= c``

    其中 ``z`` 是每条原物理边分配给闭合交替游走的 CN，``A`` 是端点上的两色平衡矩阵，
    ``w`` 只在序列边上等于片段长度。SciPy 的 ``linprog`` 只接受最小化形式，所以代码传入
    ``-w``；返回后再次用原目标和原约束独立复算，不直接相信求解器摘要。
    """

    model = build_original_graph_cycle_model(graph)
    started = perf_counter()
    result = linprog(
        -model.objective_coefficients,
        A_eq=model.balance_matrix,
        b_eq=np.zeros(model.balance_matrix.shape[0], dtype=float),
        bounds=list(zip(np.zeros_like(model.capacities), model.capacities, strict=True)),
        method="highs",
    )
    elapsed = perf_counter() - started

    if not result.success:
        return CyclicLWCNSolverResult(
            algorithm="original_graph_lp",
            status="SOLVER_FAILURE",
            solver_status=int(result.status),
            solver_message=str(result.message),
            total_lwcn=model.total_lwcn,
            maximum_cyclic_lwcn=None,
            maximum_cyclic_ratio=None,
            solve_seconds=elapsed,
            variable_count=len(model.variables),
            equality_count=model.balance_matrix.shape[0],
            max_balance_residual=None,
            max_lower_bound_violation=None,
            max_upper_bound_violation=None,
            cyclic_lower_certificate=None,
            cyclic_upper_certificate=None,
            certified_gap=None,
            variable_ids=model.variable_ids,
            circular_edge_load=None,
        )

    load = np.asarray(result.x, dtype=float)
    validation = validate_circular_edge_load(model, load)

    # HiGHS 对上界约束给出的 marginal 属于最小化问题。等式右端为 0、下界也为 0，
    # 因此“上界 × 上界边际值”就是最小化对偶目标；取负后得到最大化问题的对偶上界。
    dual_upper = None
    if hasattr(result, "upper") and hasattr(result.upper, "marginals"):
        dual_minimum = float(
            model.capacities @ np.asarray(result.upper.marginals, dtype=float)
        )
        dual_upper = -dual_minimum
    maximum = validation.maximum_cyclic_lwcn
    gap = None if dual_upper is None else max(0.0, dual_upper - maximum)

    return CyclicLWCNSolverResult(
        algorithm="original_graph_lp",
        status="OPTIMAL",
        solver_status=int(result.status),
        solver_message=str(result.message),
        total_lwcn=model.total_lwcn,
        maximum_cyclic_lwcn=maximum,
        maximum_cyclic_ratio=validation.maximum_cyclic_ratio,
        solve_seconds=elapsed,
        variable_count=len(model.variables),
        equality_count=model.balance_matrix.shape[0],
        max_balance_residual=validation.max_balance_residual,
        max_lower_bound_violation=validation.max_lower_bound_violation,
        max_upper_bound_violation=validation.max_upper_bound_violation,
        cyclic_lower_certificate=maximum,
        cyclic_upper_certificate=dual_upper,
        certified_gap=gap,
        variable_ids=model.variable_ids,
        circular_edge_load=load,
        diagnostics={
            "simplex_or_ipm_iterations": int(getattr(result, "nit", 0)),
            "sequence_variable_count": len(graph.sequence_edges),
            "breakpoint_variable_count": len(graph.breakpoint_edges),
        },
    )

