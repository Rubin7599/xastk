"""Multivariate curve resolution by alternating least squares (MCR-ALS).

Model:  D (n_t x n_E) = C (n_t x n) S (n x n_E) + E

Constraints available per iteration:
  - non-negativity of C and/or S (clipping after an unconstrained solve)
  - closure: each row of C sums to one (mass balance of the absorbing element)
  - kinetic hard model: C is replaced by the best-fitting kinetic model
    (hybrid hard-soft MCR, de Juan et al.), which strongly reduces
    rotational ambiguity

Rotational ambiguity warning: a good lack of fit does NOT show that C and S
are correct. See docs/technical-report.md, section on operando MCR.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .kinetics import KineticModel


@dataclass
class MCRResult:
    C: np.ndarray
    S: np.ndarray
    lack_of_fit: float          # percent, ||D - CS|| / ||D|| * 100
    n_iter: int
    converged: bool
    kinetic_params: np.ndarray | None = None

    def residual(self, D: np.ndarray) -> np.ndarray:
        return D - self.C @ self.S


def _lof(D, C, S):
    return 100.0 * np.linalg.norm(D - C @ S) / np.linalg.norm(D)


def mcr_als(
    D: np.ndarray,
    S0: np.ndarray | None = None,
    C0: np.ndarray | None = None,
    *,
    nonneg_C: bool = True,
    nonneg_S: bool = True,
    closure: bool = True,
    kinetic: KineticModel | None = None,
    t: np.ndarray | None = None,
    max_iter: int = 500,
    tol: float = 1e-9,
) -> MCRResult:
    """Run MCR-ALS starting from initial spectra ``S0`` or concentrations ``C0``.

    Parameters
    ----------
    D : data matrix, rows = spectra in time order, columns = energy points
    S0 : initial component spectra (n x n_E), e.g. selected rows of D
    C0 : initial concentrations (n_t x n); used if S0 is None
    kinetic : optional :class:`KineticModel`; requires ``t``
    tol : stop when the relative change of lack of fit is below tol
    """
    D = np.asarray(D, float)
    if (S0 is None) == (C0 is None):
        raise ValueError("give exactly one of S0 or C0")
    if kinetic is not None and t is None:
        raise ValueError("a kinetic constraint needs the time axis t")
    if S0 is not None:
        S = np.asarray(S0, float).copy()
        C = None
    else:
        C = np.asarray(C0, float).copy()
        S = np.linalg.lstsq(C, D, rcond=None)[0]
    if kinetic is not None and S.shape[0] != kinetic.n_components:
        raise ValueError(f"kinetic model has {kinetic.n_components} components, S has {S.shape[0]}")

    kp = None
    prev = np.inf
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        C = np.linalg.lstsq(S.T, D.T, rcond=None)[0].T
        if nonneg_C:
            C = np.clip(C, 0.0, None)
        if closure:
            sums = C.sum(axis=1, keepdims=True)
            sums[sums == 0] = 1.0
            C = C / sums
        if kinetic is not None:
            kp, C = kinetic.fit(t, C, kp)
        S = np.linalg.lstsq(C, D, rcond=None)[0]
        if nonneg_S:
            S = np.clip(S, 0.0, None)
        lof = _lof(D, C, S)
        if abs(prev - lof) <= tol * max(lof, 1e-30):
            converged = True
            break
        prev = lof
    return MCRResult(C=C, S=S, lack_of_fit=_lof(D, C, S), n_iter=it, converged=converged, kinetic_params=kp)
