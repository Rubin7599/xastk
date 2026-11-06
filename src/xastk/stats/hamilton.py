"""Hamilton (F) test for nested EXAFS models.

References
----------
W. C. Hamilton, Acta Cryst. 18, 502 (1965).
L. Downward, C. H. Booth, W. W. Lukens, F. Bridges, AIP Conf. Proc. 882, 129 (2007).
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

from scipy import stats


@dataclass(frozen=True)
class HamiltonResult:
    """Outcome of comparing a simple model A with a nested, richer model B."""

    F: float
    dof_num: int          # b = p_B - p_A, number of added variables
    dof_den: float        # N_idp - p_B
    confidence: float     # F-distribution CDF: confidence that B is a real improvement
    p_value: float        # 1 - confidence
    ratio: float          # Hamilton R ratio = sqrt(chi2_A / chi2_B)
    level: float
    significant: bool

    def as_dict(self) -> dict:
        return asdict(self)

    def __str__(self) -> str:
        verdict = "significant" if self.significant else "not significant"
        return (
            f"Hamilton test: F = {self.F:.3f} (b = {self.dof_num}, N_idp - p_B = {self.dof_den:.2f}), "
            f"confidence = {self.confidence:.3f} -> {verdict} at {self.level:.0%}"
        )


def hamilton_test(
    chi2_a: float,
    nvar_a: int,
    chi2_b: float,
    nvar_b: int,
    n_idp: float,
    level: float = 0.95,
) -> HamiltonResult:
    """F-test of whether nested model B (more variables) improves on model A.

    ``chi2_a`` and ``chi2_b`` must be computed on the same data, windows and
    noise estimate (Larch ``chi_square`` satisfies this when both fits share
    one transform). The test is only valid for nested models: B must reduce
    to A for some setting of its extra variables.
    """
    b = int(nvar_b) - int(nvar_a)
    if b <= 0:
        raise ValueError("model B must have more variables than model A (nested test)")
    dof = float(n_idp) - float(nvar_b)
    if dof <= 0:
        raise ValueError(f"no degrees of freedom left: N_idp = {n_idp:.2f}, p_B = {nvar_b}")
    if chi2_a <= 0 or chi2_b <= 0:
        raise ValueError("chi-square values must be positive")
    F = ((chi2_a - chi2_b) / b) / (chi2_b / dof)
    conf = float(stats.f.cdf(F, b, dof)) if F > 0 else 0.0
    return HamiltonResult(
        F=float(F),
        dof_num=b,
        dof_den=dof,
        confidence=conf,
        p_value=1.0 - conf,
        ratio=math.sqrt(chi2_a / chi2_b),
        level=level,
        significant=conf >= level,
    )


def hamilton_from_results(result_a, result_b, level: float = 0.95) -> HamiltonResult:
    """Hamilton test from two Larch ``feffit`` results (A nested in B).

    Uses ``chi_square``, ``nvarys`` and ``n_independent`` as computed by Larch.
    """
    nidp_a = float(result_a.n_independent)
    nidp_b = float(result_b.n_independent)
    if abs(nidp_a - nidp_b) > 1e-6:
        raise ValueError("the two fits use different fit ranges (N_idp differs)")
    return hamilton_test(
        result_a.chi_square, result_a.nvarys, result_b.chi_square, result_b.nvarys, nidp_b, level
    )
