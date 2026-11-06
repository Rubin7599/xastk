"""Operando deconvolution without pure references: PCA rank + constrained MCR-ALS."""

from .als import MCRResult, mcr_als
from .kinetics import CONSECUTIVE, FIRST_ORDER, KineticModel, consecutive, first_order
from .rank import RankEstimate, pca_rank

__all__ = [
    "MCRResult",
    "mcr_als",
    "KineticModel",
    "CONSECUTIVE",
    "FIRST_ORDER",
    "consecutive",
    "first_order",
    "RankEstimate",
    "pca_rank",
]
