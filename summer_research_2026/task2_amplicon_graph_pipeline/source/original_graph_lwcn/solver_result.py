"""三种 task2 求解器共同使用的结果协议。

统一结果结构的目的不是强迫三种算法拥有相同内部过程，而是保证批量评测时比较的是
同一组量：环状 LWCN、可行性残差、上下界证书和真实耗时。这样可以防止某个算法只
返回一个看似合理的数字，却没有说明该数字是否来自可行解、是否已经收敛。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class CyclicLWCNSolverResult:
    """一个求解器对单张断点图给出的完整、可核验结果。"""

    algorithm: str
    status: str
    solver_status: int | str
    solver_message: str
    total_lwcn: float
    maximum_cyclic_lwcn: float | None
    maximum_cyclic_ratio: float | None
    solve_seconds: float
    variable_count: int
    equality_count: int
    max_balance_residual: float | None
    max_lower_bound_violation: float | None
    max_upper_bound_violation: float | None
    cyclic_lower_certificate: float | None
    cyclic_upper_certificate: float | None
    certified_gap: float | None
    variable_ids: tuple[str, ...] = ()
    circular_edge_load: np.ndarray | None = None
    diagnostics: dict[str, object] = field(default_factory=dict)

    @property
    def succeeded(self) -> bool:
        """只有明确的最优/收敛状态才被批量评测视为成功。"""

        return self.status in {"OPTIMAL", "CONVERGED", "CONVERGED_CERTIFIED"}

    def as_dict(self, include_edge_load: bool = False) -> dict[str, object]:
        """转换成可直接写入 JSON 的字典。"""

        payload: dict[str, object] = {
            "algorithm": self.algorithm,
            "status": self.status,
            "solver_status": self.solver_status,
            "solver_message": self.solver_message,
            "total_lwcn": self.total_lwcn,
            "maximum_cyclic_lwcn": self.maximum_cyclic_lwcn,
            "maximum_cyclic_ratio": self.maximum_cyclic_ratio,
            "solve_seconds": self.solve_seconds,
            "variable_count": self.variable_count,
            "equality_count": self.equality_count,
            "max_balance_residual": self.max_balance_residual,
            "max_lower_bound_violation": self.max_lower_bound_violation,
            "max_upper_bound_violation": self.max_upper_bound_violation,
            "cyclic_lower_certificate": self.cyclic_lower_certificate,
            "cyclic_upper_certificate": self.cyclic_upper_certificate,
            "certified_gap": self.certified_gap,
            "diagnostics": self.diagnostics,
        }
        if include_edge_load and self.circular_edge_load is not None:
            payload["circular_edge_load"] = {
                edge_id: float(value)
                for edge_id, value in zip(
                    self.variable_ids, self.circular_edge_load, strict=True
                )
                if value > 1e-10
            }
        return payload

