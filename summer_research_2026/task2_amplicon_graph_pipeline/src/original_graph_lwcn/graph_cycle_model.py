"""把 CoRAL 原断点图直接写成“最大平衡环流”模型。

核心思想
--------
原图的边天然分成两种颜色：序列边和断点边。一个非负边负载能够分解为闭合交替游走，
当且仅当每个端点上两种颜色的关联负载相等。因此不需要把“断点边 + 下一条序列边”
压缩成状态弧，只需为每条原物理边设置一个“分配给环的 CN”变量并写局部平衡式。

对于断点自环，边在同一个端点出现两次，所以平衡矩阵中的关联系数必须是 2；这是实现中
最容易被漏掉、也会直接破坏 foldback 结果的细节。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags

from cyclic_lwcn.model import BreakpointGraph


@dataclass(frozen=True)
class OriginalEdgeVariable:
    """原断点图中一个环状 CN 变量的物理含义。"""

    index: int
    edge_id: str
    edge_kind: str
    node1: str
    node2: str
    capacity: float
    lwcn_coefficient: float


@dataclass(frozen=True)
class OriginalGraphCycleModel:
    """未做状态展开的原图线性模型。"""

    graph: BreakpointGraph
    variables: tuple[OriginalEdgeVariable, ...]
    capacities: np.ndarray
    objective_coefficients: np.ndarray
    balance_matrix: csr_matrix
    balance_nodes: tuple[str, ...]
    total_lwcn: float

    @property
    def variable_ids(self) -> tuple[str, ...]:
        return tuple(variable.edge_id for variable in self.variables)


@dataclass(frozen=True)
class ScaledCycleModel:
    """供非线性方法使用的无量纲模型。

    令 ``z_j = capacity_j * x_j``，则 ``0 <= x_j <= 1``。目标再除以图总 LWCN，
    直接得到环状比例。行归一化只改善数值条件，不改变平衡等式的解集。
    """

    original: OriginalGraphCycleModel
    active_indices: np.ndarray
    normalized_objective: np.ndarray
    normalized_balance_matrix: csr_matrix
    kept_balance_rows: np.ndarray
    row_scales: np.ndarray

    def to_original_load(self, scaled_load: np.ndarray) -> np.ndarray:
        """把正容量变量的 [0,1] 比例还原成全部原边的实际 CN。"""

        load = np.zeros(len(self.original.variables), dtype=float)
        capacities = self.original.capacities[self.active_indices]
        load[self.active_indices] = capacities * np.asarray(scaled_load, dtype=float)
        return load


@dataclass(frozen=True)
class CircularLoadValidation:
    """对一个候选环流重新计算的独立可行性指标。"""

    maximum_cyclic_lwcn: float
    maximum_cyclic_ratio: float
    max_balance_residual: float
    max_lower_bound_violation: float
    max_upper_bound_violation: float


def build_original_graph_cycle_model(graph: BreakpointGraph) -> OriginalGraphCycleModel:
    """在原物理边上建立最大环状 LWCN 模型。

    变量顺序固定为“全部序列边在前、全部断点边在后”，使 CSV、JSON 和不同算法间的
    边级结果可以逐项对齐。每个端点 ``v`` 的等式为：

    ``环状序列 CN(v) - 断点边环状关联 CN 总和(v) = 0``。
    """

    nodes = tuple(sorted(graph.nodes))
    node_row = {node: row for row, node in enumerate(nodes)}
    variables: list[OriginalEdgeVariable] = []
    entries: list[tuple[int, int, float]] = []

    # 序列边在其两个端点各贡献一次“红色”关联量，目标系数就是闭区间长度。
    for edge in graph.sequence_edges:
        column = len(variables)
        variables.append(
            OriginalEdgeVariable(
                index=column,
                edge_id=edge.edge_id,
                edge_kind="sequence",
                node1=edge.start_node,
                node2=edge.end_node,
                capacity=float(edge.copy_number),
                lwcn_coefficient=float(edge.length),
            )
        )
        entries.append((node_row[edge.start_node], column, 1.0))
        entries.append((node_row[edge.end_node], column, 1.0))

    # 断点边作为另一种颜色，以负号进入平衡式。自环的两次 append 会由稀疏矩阵自动求和为 -2。
    for edge in graph.breakpoint_edges:
        column = len(variables)
        variables.append(
            OriginalEdgeVariable(
                index=column,
                edge_id=edge.edge_id,
                edge_kind=edge.edge_type,
                node1=edge.node1,
                node2=edge.node2,
                capacity=float(edge.copy_number),
                lwcn_coefficient=0.0,
            )
        )
        entries.append((node_row[edge.node1], column, -1.0))
        entries.append((node_row[edge.node2], column, -1.0))

    if entries:
        rows, columns, values = zip(*entries, strict=True)
        balance = coo_matrix(
            (values, (rows, columns)),
            shape=(len(nodes), len(variables)),
            dtype=float,
        ).tocsr()
    else:
        balance = csr_matrix((len(nodes), len(variables)), dtype=float)

    capacities = np.array([variable.capacity for variable in variables], dtype=float)
    objective = np.array(
        [variable.lwcn_coefficient for variable in variables], dtype=float
    )
    return OriginalGraphCycleModel(
        graph=graph,
        variables=tuple(variables),
        capacities=capacities,
        objective_coefficients=objective,
        balance_matrix=balance,
        balance_nodes=nodes,
        total_lwcn=float(graph.total_lwcn),
    )


def _independent_balance_rows(model: OriginalGraphCycleModel) -> np.ndarray:
    """按连通分量删除恰好一条冗余的二部图平衡式。

    这里不用每张图都做昂贵的稠密 QR。对无向关联矩阵，行相关向量必须沿每条边交替变号；
    因而一个连通分量仅在底层图为二部图时产生一维行冗余。若分量含奇环或自环，则所有行
    独立。这个判定只影响非线性求解器的数值坐标，不改变原始平衡方程。
    """

    adjacency: dict[str, list[str]] = {node: [] for node in model.balance_nodes}
    has_self_loop: set[str] = set()
    for variable in model.variables:
        adjacency[variable.node1].append(variable.node2)
        adjacency[variable.node2].append(variable.node1)
        if variable.node1 == variable.node2:
            has_self_loop.add(variable.node1)

    row_by_node = {node: row for row, node in enumerate(model.balance_nodes)}
    keep = np.ones(len(model.balance_nodes), dtype=bool)
    color: dict[str, int] = {}
    for start in model.balance_nodes:
        if start in color:
            continue
        color[start] = 0
        stack = [start]
        component = []
        bipartite = True
        while stack:
            node = stack.pop()
            component.append(node)
            if node in has_self_loop:
                bipartite = False
            for neighbor in adjacency[node]:
                if neighbor not in color:
                    color[neighbor] = 1 - color[node]
                    stack.append(neighbor)
                elif color[neighbor] == color[node]:
                    bipartite = False
        if bipartite and component:
            # 删除哪个端点都不改变方程组；固定删排序最小者可保证完全可复现。
            dropped = min(component)
            keep[row_by_node[dropped]] = False
    return np.flatnonzero(keep)


def build_scaled_cycle_model(
    model: OriginalGraphCycleModel,
    *,
    positive_capacity_tolerance: float = 0.0,
) -> ScaledCycleModel:
    """构造非线性求解使用的 [0,1] 变量、独立平衡行和归一化目标。"""

    active = np.flatnonzero(model.capacities > positive_capacity_tolerance)
    kept_rows = _independent_balance_rows(model)
    if active.size:
        scaled_balance = model.balance_matrix[kept_rows][:, active] @ diags(
            model.capacities[active]
        )
        scaled_balance = scaled_balance.tocsr()
    else:
        scaled_balance = csr_matrix((len(kept_rows), 0), dtype=float)

    # 每行除以二范数，使“大 CN 边”和“小 CN 边”共同出现时不会让停止条件失真。
    row_norms = np.sqrt(np.asarray(scaled_balance.power(2).sum(axis=1)).ravel())
    row_scales = np.ones_like(row_norms)
    nonzero = row_norms > 0
    row_scales[nonzero] = 1.0 / row_norms[nonzero]
    normalized_balance = diags(row_scales) @ scaled_balance

    if model.total_lwcn > 0:
        normalized_objective = (
            model.objective_coefficients[active] * model.capacities[active]
        ) / model.total_lwcn
    else:
        normalized_objective = np.zeros(len(active), dtype=float)
    return ScaledCycleModel(
        original=model,
        active_indices=active,
        normalized_objective=np.asarray(normalized_objective, dtype=float),
        normalized_balance_matrix=normalized_balance.tocsr(),
        kept_balance_rows=kept_rows,
        row_scales=row_scales,
    )


def validate_circular_edge_load(
    model: OriginalGraphCycleModel,
    circular_edge_load: np.ndarray,
) -> CircularLoadValidation:
    """不依赖求解器状态，直接从原矩阵复核边负载。"""

    load = np.asarray(circular_edge_load, dtype=float)
    if load.shape != model.capacities.shape:
        raise ValueError(
            f"边负载维度 {load.shape} 与模型维度 {model.capacities.shape} 不一致"
        )
    balance_residual = np.asarray(model.balance_matrix @ load).ravel()
    maximum = float(model.objective_coefficients @ load)
    ratio = maximum / model.total_lwcn if model.total_lwcn > 0 else 0.0
    return CircularLoadValidation(
        maximum_cyclic_lwcn=maximum,
        maximum_cyclic_ratio=ratio,
        max_balance_residual=float(
            np.max(np.abs(balance_residual), initial=0.0)
        ),
        max_lower_bound_violation=max(0.0, float(np.max(-load, initial=0.0))),
        max_upper_bound_violation=max(
            0.0, float(np.max(load - model.capacities, initial=0.0))
        ),
    )
