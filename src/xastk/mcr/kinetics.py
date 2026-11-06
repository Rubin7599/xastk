"""Kinetic hard models used as constraints inside MCR-ALS.

Each model maps a time (or potential / charge) axis and rate parameters to a
concentration matrix C (n_t x n_components) whose rows sum to one.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


def first_order(t: np.ndarray, k: float) -> np.ndarray:
    """A -> B with rate k."""
    a = np.exp(-k * t)
    return np.column_stack([a, 1.0 - a])


def consecutive(t: np.ndarray, k1: float, k2: float) -> np.ndarray:
    """A -> B -> C (first-order, consecutive), handles k1 == k2."""
    a = np.exp(-k1 * t)
    if abs(k2 - k1) < 1e-9 * max(k1, k2, 1e-30):
        b = k1 * t * np.exp(-k1 * t)
    else:
        b = k1 / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t))
    return np.column_stack([a, b, 1.0 - a - b])


@dataclass
class KineticModel:
    """A kinetic hard model with positive rate parameters.

    Parameters are fitted in log space, which keeps them positive and makes
    rates spanning orders of magnitude well conditioned.
    """

    func: Callable[..., np.ndarray]
    param_names: Sequence[str]
    p0: Sequence[float]
    n_components: int

    def __call__(self, t: np.ndarray, params: Sequence[float]) -> np.ndarray:
        return self.func(t, *params)

    def fit(self, t: np.ndarray, C: np.ndarray, p_start: Sequence[float] | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Least-squares fit of the model to concentration profiles C.

        Returns (params, C_model).
        """
        start = np.log(np.asarray(p_start if p_start is not None else self.p0, float))

        def res(logp):
            return (self.func(t, *np.exp(logp)) - C).ravel()

        sol = least_squares(res, start)
        p = np.exp(sol.x)
        return p, self.func(t, *p)


FIRST_ORDER = KineticModel(first_order, ["k"], [0.1], 2)
CONSECUTIVE = KineticModel(consecutive, ["k1", "k2"], [0.2, 0.02], 3)
